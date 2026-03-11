import os
import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON


def query_graphdb(query_input, repo_name="mtrKG", host="localhost", port=7200):
    """
    Sends a SPARQL query to a local Ontotext GraphDB repository and
    returns the results as a clean Pandas DataFrame.

    Args:
        query_input (str): Either a raw SPARQL string OR a path to a file (e.g., 'query.rq').
        repo_name (str): Target GraphDB repository name.
        host (str): GraphDB host address.
        port (int): GraphDB port.
    """

    # ==========================================
    # NEW: CHECK IF INPUT IS A FILE OR A STRING
    # ==========================================
    if isinstance(query_input, str) and os.path.isfile(query_input):
        print(f"Reading SPARQL query from file: {query_input}")
        with open(query_input, 'r', encoding='utf-8-sig') as file:
            query_string = file.read()
    else:
        # If it's not a valid file path, assume it is a raw SPARQL string
        query_string = query_input
    query_string = query_string.strip().replace('\xa0', ' ').encode('ascii', 'ignore').decode('utf-8')
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

        # 5. Parse the standard W3C SPARQL JSON format
        # Extract column names (the variables in your SELECT statement)
        variables = results["head"]["vars"]

        # Extract the rows
        bindings = results["results"]["bindings"]

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
        return df

    except Exception as e:
        print(f"GraphDB Connection/Query Error: {e}")
        print("-> Make sure GraphDB is running and your repository name is correct!")
        return None