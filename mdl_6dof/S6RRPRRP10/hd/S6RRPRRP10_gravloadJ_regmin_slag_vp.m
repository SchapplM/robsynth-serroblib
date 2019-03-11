% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t58 = sin(qJ(2));
t59 = sin(qJ(1));
t61 = cos(qJ(2));
t80 = cos(pkin(6));
t90 = cos(qJ(1));
t69 = t80 * t90;
t36 = t58 * t69 + t59 * t61;
t52 = pkin(11) + qJ(4);
t49 = sin(t52);
t50 = cos(t52);
t54 = sin(pkin(6));
t77 = t54 * t90;
t21 = t36 * t50 - t49 * t77;
t35 = t59 * t58 - t61 * t69;
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t7 = t21 * t57 - t35 * t60;
t8 = t21 * t60 + t35 * t57;
t74 = t59 * t80;
t37 = t90 * t58 + t61 * t74;
t97 = -g(1) * t37 - g(2) * t35;
t38 = -t58 * t74 + t90 * t61;
t84 = t54 * t59;
t24 = t38 * t49 - t50 * t84;
t75 = -t36 * t49 - t50 * t77;
t85 = t54 * t58;
t64 = -g(3) * (-t49 * t85 + t80 * t50) - g(2) * t75 + g(1) * t24;
t91 = g(3) * t54;
t87 = t50 * t57;
t86 = t50 * t60;
t83 = t54 * t61;
t82 = t60 * t61;
t81 = t90 * pkin(1) + pkin(8) * t84;
t79 = t57 * t83;
t53 = sin(pkin(11));
t78 = t53 * t84;
t76 = -t59 * pkin(1) + pkin(8) * t77;
t73 = t53 * t77;
t25 = t38 * t50 + t49 * t84;
t11 = t25 * t57 - t37 * t60;
t72 = -g(1) * t7 + g(2) * t11;
t71 = g(1) * t75 + g(2) * t24;
t70 = g(1) * t35 - g(2) * t37;
t29 = t80 * t49 + t50 * t85;
t18 = t29 * t57 + t54 * t82;
t1 = g(1) * t11 + g(2) * t7 + g(3) * t18;
t12 = t25 * t60 + t37 * t57;
t19 = t29 * t60 - t79;
t66 = g(1) * t12 + g(2) * t8 + g(3) * t19;
t13 = -t35 * t87 - t36 * t60;
t15 = -t37 * t87 - t38 * t60;
t26 = t50 * t79 - t60 * t85;
t65 = g(1) * t15 + g(2) * t13 + g(3) * t26;
t63 = g(1) * t25 + g(2) * t21 + g(3) * t29;
t17 = g(3) * t83 + t97;
t62 = g(1) * t38 + g(2) * t36 + g(3) * t85;
t56 = -pkin(9) - qJ(3);
t55 = cos(pkin(11));
t48 = t55 * pkin(3) + pkin(2);
t27 = (t50 * t82 + t57 * t58) * t54;
t16 = -t37 * t86 + t38 * t57;
t14 = -t35 * t86 + t36 * t57;
t6 = t17 * t49;
t5 = t64 * t60;
t4 = t64 * t57;
t3 = -g(1) * t16 - g(2) * t14 - g(3) * t27;
t2 = g(1) * t8 - g(2) * t12;
t9 = [0, g(1) * t59 - g(2) * t90, g(1) * t90 + g(2) * t59, 0, 0, 0, 0, 0, g(1) * t36 - g(2) * t38, -t70, -g(1) * (-t36 * t55 + t73) - g(2) * (t38 * t55 + t78) -g(1) * (t36 * t53 + t55 * t77) - g(2) * (-t38 * t53 + t55 * t84) t70, -g(1) * (-t36 * pkin(2) - t35 * qJ(3) + t76) - g(2) * (t38 * pkin(2) + t37 * qJ(3) + t81) 0, 0, 0, 0, 0, g(1) * t21 - g(2) * t25, t71, 0, 0, 0, 0, 0, t2, t72, t2, -t71, -t72, -g(1) * (pkin(3) * t73 - pkin(4) * t21 - pkin(5) * t8 + pkin(10) * t75 - qJ(6) * t7 + t35 * t56 - t36 * t48 + t76) - g(2) * (pkin(3) * t78 + t25 * pkin(4) + t12 * pkin(5) + t24 * pkin(10) + t11 * qJ(6) - t37 * t56 + t38 * t48 + t81); 0, 0, 0, 0, 0, 0, 0, 0, -t17, t62, -t17 * t55, t17 * t53, -t62, -g(1) * (-t37 * pkin(2) + t38 * qJ(3)) - g(2) * (-t35 * pkin(2) + t36 * qJ(3)) - (pkin(2) * t61 + qJ(3) * t58) * t91, 0, 0, 0, 0, 0, -t17 * t50, t6, 0, 0, 0, 0, 0, t3, t65, t3, -t6, -t65, -g(1) * (t16 * pkin(5) + t15 * qJ(6) - t38 * t56) - g(2) * (t14 * pkin(5) + t13 * qJ(6) - t36 * t56) - g(3) * (t27 * pkin(5) + t26 * qJ(6)) + t56 * t58 * t91 + (-t61 * t91 - t97) * (pkin(4) * t50 + pkin(10) * t49 + t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t63, 0, 0, 0, 0, 0, t5, -t4, t5, -t63, t4, -pkin(10) * t63 + t64 * (pkin(5) * t60 + qJ(6) * t57 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t66, t1, 0, -t66, -g(1) * (-t11 * pkin(5) + t12 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - g(3) * (-t18 * pkin(5) + t19 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t9;
