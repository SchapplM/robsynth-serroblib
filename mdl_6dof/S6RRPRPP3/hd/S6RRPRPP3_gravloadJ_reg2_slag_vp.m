% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t37 = pkin(9) + qJ(4);
t31 = sin(t37);
t32 = cos(t37);
t44 = cos(qJ(1));
t70 = t44 * t32;
t42 = sin(qJ(1));
t43 = cos(qJ(2));
t72 = t42 * t43;
t12 = t31 * t72 + t70;
t71 = t44 * t31;
t13 = t32 * t72 - t71;
t83 = -t13 * pkin(4) - t12 * qJ(5);
t79 = g(2) * t42;
t23 = g(1) * t44 + t79;
t41 = sin(qJ(2));
t16 = -g(3) * t43 + t23 * t41;
t82 = pkin(4) * t32;
t81 = g(1) * t42;
t40 = -pkin(8) - qJ(3);
t77 = pkin(5) - t40;
t76 = t31 * t41;
t75 = t32 * t41;
t74 = t41 * t40;
t38 = sin(pkin(9));
t73 = t42 * t38;
t39 = cos(pkin(9));
t30 = t39 * pkin(3) + pkin(2);
t22 = t43 * t30;
t69 = t44 * t38;
t68 = t44 * t39;
t67 = -pkin(4) - qJ(6);
t66 = t44 * pkin(1) + t42 * pkin(7);
t65 = qJ(5) * t31;
t63 = t44 * t74;
t34 = t44 * pkin(7);
t62 = pkin(3) * t69 + t42 * t74 + t34;
t61 = t77 * t44;
t60 = -pkin(1) - t22;
t59 = -t12 * pkin(4) + t13 * qJ(5);
t14 = -t42 * t32 + t43 * t71;
t15 = t42 * t31 + t43 * t70;
t58 = -t14 * pkin(4) + t15 * qJ(5);
t57 = -t30 - t65;
t56 = pkin(3) * t73 + t44 * t22 + t66;
t55 = g(3) * (t22 + (t65 + t82) * t43);
t4 = g(1) * t12 - g(2) * t14;
t5 = g(1) * t13 - g(2) * t15;
t54 = -g(2) * t44 + t81;
t53 = t43 * pkin(2) + t41 * qJ(3);
t49 = t23 * t43;
t48 = t15 * pkin(4) + t14 * qJ(5) + t56;
t2 = g(1) * t14 + g(2) * t12 + g(3) * t76;
t46 = g(1) * t15 + g(2) * t13 + g(3) * t75;
t45 = t60 * t42 + t62;
t20 = qJ(5) * t75;
t18 = t54 * t41;
t17 = g(3) * t41 + t49;
t7 = t16 * t32;
t6 = t16 * t31;
t1 = [0, 0, 0, 0, 0, 0, t54, t23, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t43, -t18, -t23, -g(1) * (-t42 * pkin(1) + t34) - g(2) * t66, 0, 0, 0, 0, 0, 0, -g(1) * (-t39 * t72 + t69) - g(2) * (t43 * t68 + t73) -g(1) * (t38 * t72 + t68) - g(2) * (t42 * t39 - t43 * t69) t18, -g(1) * t34 - g(2) * (t53 * t44 + t66) - (-pkin(1) - t53) * t81, 0, 0, 0, 0, 0, 0, t5, -t4, t18, -g(1) * t45 - g(2) * (t56 - t63) 0, 0, 0, 0, 0, 0, t18, -t5, t4, -g(1) * (t45 + t83) - g(2) * (t48 - t63) 0, 0, 0, 0, 0, 0, t18, t4, t5, -g(1) * (-t13 * qJ(6) + t62 + t83) - g(2) * (t15 * qJ(6) + t41 * t61 + t48) - (-t41 * pkin(5) + t60) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t39, -t16 * t38, -t17, -g(3) * t53 + t23 * (pkin(2) * t41 - qJ(3) * t43) 0, 0, 0, 0, 0, 0, t7, -t6, -t17, -g(3) * (t22 - t74) + t23 * (t30 * t41 + t40 * t43) 0, 0, 0, 0, 0, 0, -t17, -t7, t6, -t55 + t40 * t49 + (g(3) * t40 + t23 * (-t57 + t82)) * t41, 0, 0, 0, 0, 0, 0, -t17, t6, t7, -t55 + (-g(3) * qJ(6) * t32 - g(1) * t61 - t77 * t79) * t43 + (-g(3) * t77 + t23 * (-t67 * t32 - t57)) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t46, -g(1) * t58 - g(2) * t59 - g(3) * (-pkin(4) * t76 + t20) 0, 0, 0, 0, 0, 0, 0, -t46, t2, -g(1) * (-t14 * qJ(6) + t58) - g(2) * (-t12 * qJ(6) + t59) - g(3) * (t67 * t76 + t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46;];
taug_reg  = t1;
