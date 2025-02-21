from flask import Flask, render_template

app = Flask(__name__)

@app.route('/')
def main():
    return render_template('main.html')

@app.route('/pagina1')
def pagina1():
    return render_template('pagina1.html')

@app.route('/pagina2')
def pagina2():
    return render_template('pagina2.html')

if __name__ == '__main__':
    app.run(debug=True)
