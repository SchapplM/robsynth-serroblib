% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t55 = sin(pkin(9));
t56 = cos(pkin(9));
t64 = sin(qJ(1));
t65 = cos(qJ(1));
t24 = -t55 * t64 - t56 * t65;
t25 = t55 * t65 - t56 * t64;
t50 = g(1) * t24 + g(2) * t25;
t43 = g(3) * t40 - t38 * t50;
t71 = t38 * pkin(8);
t67 = g(3) * t38;
t63 = t24 * t38;
t62 = t24 * t40;
t61 = t25 * t38;
t60 = t25 * t40;
t37 = sin(qJ(5));
t59 = t37 * t40;
t39 = cos(qJ(5));
t58 = t39 * t40;
t57 = pkin(1) * t65 + qJ(2) * t64;
t54 = pkin(2) * t65 + t57;
t10 = -t24 * t59 - t25 * t39;
t6 = -t24 * t39 + t25 * t59;
t53 = g(1) * t6 + g(2) * t10;
t52 = -pkin(1) * t64 + qJ(2) * t65;
t51 = g(1) * t25 - g(2) * t24;
t7 = t24 * t37 + t25 * t58;
t48 = -pkin(3) * t24 + pkin(7) * t25 + t54;
t47 = -g(1) * t10 + g(2) * t6 + t37 * t67;
t11 = -t24 * t58 + t25 * t37;
t46 = -g(1) * t11 + g(2) * t7 + t39 * t67;
t45 = -pkin(4) * t62 - pkin(8) * t63 + t48;
t44 = -pkin(2) * t64 + t52;
t42 = pkin(3) * t25 + pkin(7) * t24 + t44;
t41 = pkin(4) * t60 + pkin(8) * t61 + t42;
t27 = g(1) * t65 + g(2) * t64;
t26 = g(1) * t64 - g(2) * t65;
t17 = pkin(8) * t62;
t15 = pkin(8) * t60;
t12 = t51 * t38;
t5 = -t40 * t50 - t67;
t4 = t43 * t39;
t3 = t43 * t37;
t2 = -g(1) * t7 - g(2) * t11;
t1 = [0, 0, 0, 0, 0, 0, t26, t27, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t27, -g(1) * t52 - g(2) * t57, 0, 0, 0, 0, 0, 0, -t51, t50, 0, -g(1) * t44 - g(2) * t54, 0, 0, 0, 0, 0, 0, -t51 * t40, t12, -t50, -g(1) * t42 - g(2) * t48, 0, 0, 0, 0, 0, 0, t2, t53, -t12, -g(1) * t41 - g(2) * t45, 0, 0, 0, 0, 0, 0, t2, -t12, -t53, -g(1) * (pkin(5) * t7 + qJ(6) * t6 + t41) - g(2) * (pkin(5) * t11 + qJ(6) * t10 + t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t5, -g(1) * (pkin(4) * t63 - t17) - g(2) * (pkin(4) * t61 - t15) - g(3) * (-pkin(4) * t40 - t71) 0, 0, 0, 0, 0, 0, t4, -t5, t3, g(3) * t71 + g(1) * t17 + g(2) * t15 + t43 * (pkin(5) * t39 + qJ(6) * t37 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t46, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, t46, -g(1) * (-pkin(5) * t10 + qJ(6) * t11) - g(2) * (pkin(5) * t6 - qJ(6) * t7) - (pkin(5) * t37 - qJ(6) * t39) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47;];
taug_reg  = t1;
