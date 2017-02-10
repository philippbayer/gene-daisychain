# This class provides functionality to start, stop or restart a project-specific neo4j graph database instance
import os
import subprocess


class DBRUnner:

    def __init__(self, send_data):
        # Define path to neo4j binary
        self.send_data = send_data

    # Reply to request send from a user app
    # User_request is a list produced by the "_" split command
    # e.g. [ProjectID, START]
    def evaluate_user_request(self, user_request):
        # Check for correct syntax: Proj_ID + command
        if len(user_request) != 2 or not user_request[0].isdigit():
            self.send_data("-5")
        # Define path to neo4j binary
        neo4j_binary = os.path.join("Projects", str(user_request[0]), "proj_graph_db", "bin", "neo4j")
        # Start db?
        if user_request[1] == "START":
            self.start(neo4j_binary)
        elif user_request[1] == "STOP":
            self.stop(neo4j_binary)
        elif user_request[1] == "RESTART":
            self.restart(neo4j_binary)
        elif user_request[1] == "STATUS":
            self.get_status(neo4j_binary)
        else:
            self.send_data("-5")

    def get_status(self,neo4j_binary):
        # Retrieve status (running/not running) from neo4j instance
        status = subprocess.run([neo4j_binary], "status", stdout=subprocess.PIPE)
        self.send_data(status)

    def start(self,neo4j_binary):
        # Start up the project neo4j database
        self.send_data("Starting Project DB")
        subprocess.run([neo4j_binary], "start")

    def stop(self,neo4j_binary):
        # Start up the project neo4j database
        self.send_data("Stopping Project DB")
        subprocess.run([neo4j_binary], "stop")

    def restart(self,neo4j_binary):
        self.send_data("Restarting Project DB")
        subprocess.run([neo4j_binary], "stop")
        subprocess.run([neo4j_binary], "start")