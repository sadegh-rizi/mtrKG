import os
import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON


def _resolve_query_input(query_input):
    query_path = None

    if isinstance(query_input, str) and os.path.isfile(query_input):
        query_path = query_input
        print(f"Reading SPARQL query from file: {query_path}")
        with open(query_path, "r", encoding="utf-8-sig") as file:
            query_string = file.read()
    else:
        query_string = query_input

    if not isinstance(query_string, str):
        raise TypeError(
            "query_input must be a SPARQL string or a valid path to a query file."
        )

    clean_query = (
        query_string.strip()
        .replace("\xa0", " ")
        .encode("ascii", "ignore")
        .decode("utf-8")
    )
    return clean_query, query_path


def _csv_path_from_query_path(query_path):
    base_path, _ = os.path.splitext(query_path)
    return f"{base_path}.csv"


def query_graphdb(
    query_input,
    repo_name="mtrKG",
    host="localhost",
    port=7200,
    show_query=True,
    save_csv=True,
):
    """
    Send a SPARQL query to a GraphDB repository and return results as a Pandas DataFrame.

    Args:
        query_input (str): Raw SPARQL string or path to a query file (e.g., 'query.rq').
        repo_name (str): Target GraphDB repository name.
        host (str): GraphDB host address.
        port (int): GraphDB port.
        show_query (bool): Print the final SPARQL query text before execution.
        save_csv (bool): Save SELECT results to CSV when query_input is a file path.
    """
    query_string, query_path = _resolve_query_input(query_input)

    if show_query:
        print("SPARQL query being executed:")
        print(query_string)

    # 1. Define the GraphDB endpoint URL
    # GraphDB always exposes repositories at /repositories/{name}
    endpoint_url = f"http://{host}:{port}/repositories/{repo_name}"

    # 2. Initialize the wrapper
    sparql = SPARQLWrapper(endpoint_url)
    sparql.setQuery(query_string)

    # 3. Request JSON format (Best for converting to Pandas)
    sparql.setReturnFormat(JSON)

    print(f"Sending query to GraphDB ({repo_name})...")

    try:
        # 4. Execute the query
        results = sparql.query().convert()

        if "boolean" in results:
            print(f"ASK query result: {results['boolean']}")
            return results["boolean"]

        # 5. Parse the standard W3C SPARQL JSON format
        # Extract column names (the variables in your SELECT statement)
        variables = results.get("head", {}).get("vars", [])

        # Extract the rows
        bindings = results.get("results", {}).get("bindings", [])

        # Build the data table
        data = []
        for row in bindings:
            parsed_row = []
            for var in variables:
                if var in row:
                    # Extract just the string/number value, ignoring the datatype metadata
                    parsed_row.append(row[var]["value"])
                else:
                    parsed_row.append(None)
            data.append(parsed_row)

        # 6. Convert to Pandas DataFrame
        df = pd.DataFrame(data, columns=variables)
        print(f"Success! Retrieved {len(df)} rows.")

        if save_csv and query_path:
            csv_path = _csv_path_from_query_path(query_path)
            df.to_csv(csv_path, index=False)
            print(f"Saved results to CSV: {csv_path}")

        return df

    except Exception as e:
        print(f"GraphDB Connection/Query Error: {e}")
        print("-> Make sure GraphDB is running and your repository name is correct!")
        return None
