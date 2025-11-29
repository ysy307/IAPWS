# 実行ファイルの定義
add_executable(main main.F90)

# ライブラリをリンク
target_link_libraries(main PRIVATE iapws)

# テスト設定
add_test(NAME run_main COMMAND main)