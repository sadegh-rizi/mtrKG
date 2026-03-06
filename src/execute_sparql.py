import os
import pandas as pd
from rdflib import Graph

def execute_sparql(g, query_input):
    """
    Executes an arbitrary SPARQL query against an rdflib Graph.

    Args:
        g (rdflib.Graph): Your populated Knowledge Graph.
        query_input (str): Either a raw SPARQL string, OR a path to a file
                           (e.g., 'query.txt', 'my_query.rq').

    Returns:
        pd.DataFrame: A Pandas DataFrame containing the query results (for SELECT queries),
                      or the raw boolean/graph output for ASK/CONSTRUCT queries.
    """
    # 1. Determine if input is a file path or a raw query string
    # We check if it looks like a file extension to avoid OS errors on massive strings
    if isinstance(query_input, str) and query_input.endswith(('.txt', '.rq', '.sparql')):
        if os.path.exists(query_input):
            with open(query_input, 'r') as file:
                query_string = file.read()
                print(f"Loaded query from file: {query_input}")
        else:
            print(f"Error: File '{query_input}' not found.")
            return None
    else:
        query_string = query_input

    print("Executing SPARQL query...")

    try:
        # 2. Run the query
        results = g.query(query_string)

        # 3. Handle SELECT queries (The most common type, returns tabular data)
        if results.type == 'SELECT':
            # Extract variable names (the column headers from your SELECT ?var1 ?var2 line)
            columns = [str(var) for var in results.vars]

            data = []
            for row in results:
                # Convert rdflib URIRefs and Literals into standard Python strings/numbers
                clean_row = []
                for item in row:
                    if item is None:
                        clean_row.append(None)
                    else:
                        # .toPython() automatically converts XSD.float to Python floats, XSD.integer to ints, etc.
                        clean_row.append(item.toPython())
                data.append(clean_row)

            # Create the DataFrame
            df = pd.DataFrame(data, columns=columns)
            print(f"Success! Query returned {len(df)} rows.")
            return df

        # 4. Handle ASK queries (Returns True/False)
        elif results.type == 'ASK':
            print(f"ASK Query Result: {results.askAnswer}")
            return results.askAnswer

        # 5. Handle CONSTRUCT/DESCRIBE queries (Returns a sub-graph)
        elif results.type in ['CONSTRUCT', 'DESCRIBE']:
            print(f"Success! Query returned a sub-graph with {len(results.graph)} triples.")
            return results.graph

    except Exception as e:
        print(f"SPARQL Execution Error: {e}")
        return None