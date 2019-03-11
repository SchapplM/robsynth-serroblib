% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t39 = sin(qJ(2));
t40 = sin(qJ(1));
t42 = cos(qJ(2));
t43 = cos(qJ(1));
t60 = cos(pkin(6));
t57 = t43 * t60;
t21 = t39 * t57 + t40 * t42;
t58 = t40 * t60;
t23 = -t39 * t58 + t43 * t42;
t77 = -g(1) * t23 - g(2) * t21;
t34 = pkin(11) + qJ(6);
t31 = sin(t34);
t32 = cos(t34);
t20 = t40 * t39 - t42 * t57;
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t36 = sin(pkin(6));
t63 = t36 * t43;
t50 = -t20 * t38 + t41 * t63;
t76 = t21 * t32 + t31 * t50;
t75 = -t21 * t31 + t32 * t50;
t72 = g(3) * t36;
t69 = t31 * t38;
t68 = t32 * t38;
t35 = sin(pkin(11));
t67 = t35 * t38;
t66 = t36 * t39;
t65 = t36 * t40;
t64 = t36 * t42;
t37 = cos(pkin(11));
t62 = t37 * t38;
t61 = t38 * t39;
t59 = g(3) * (pkin(2) * t64 + qJ(3) * t66);
t49 = t20 * t41 + t38 * t63;
t22 = t43 * t39 + t42 * t58;
t5 = -t22 * t41 + t38 * t65;
t56 = g(1) * t49 + g(2) * t5;
t55 = g(1) * t20 - g(2) * t22;
t54 = g(1) * t21 - g(2) * t23;
t53 = g(1) * t43 + g(2) * t40;
t52 = pkin(4) * t38 - qJ(5) * t41;
t51 = t43 * pkin(1) + t23 * pkin(2) + pkin(8) * t65 + t22 * qJ(3);
t18 = t60 * t38 + t41 * t64;
t47 = g(1) * t5 - g(2) * t49 + g(3) * t18;
t19 = -t38 * t64 + t60 * t41;
t6 = t22 * t38 + t41 * t65;
t46 = g(1) * t6 - g(2) * t50 + g(3) * t19;
t45 = -t40 * pkin(1) - t21 * pkin(2) + pkin(8) * t63 - t20 * qJ(3);
t4 = -g(1) * t22 - g(2) * t20 + g(3) * t64;
t44 = g(3) * t66 - t77;
t16 = t22 * pkin(2);
t14 = t20 * pkin(2);
t3 = t44 * t41;
t2 = t23 * t31 + t6 * t32;
t1 = t23 * t32 - t6 * t31;
t7 = [0, g(1) * t40 - g(2) * t43, t53, 0, 0, 0, 0, 0, t54, -t55, -t53 * t36, -t54, t55, -g(1) * t45 - g(2) * t51, 0, 0, 0, 0, 0, -g(1) * t50 - g(2) * t6, t56, -g(1) * (-t21 * t35 + t37 * t50) - g(2) * (t23 * t35 + t6 * t37) -g(1) * (-t21 * t37 - t35 * t50) - g(2) * (t23 * t37 - t6 * t35) -t56, -g(1) * (pkin(3) * t63 + pkin(4) * t50 - t21 * pkin(9) + qJ(5) * t49 + t45) - g(2) * (pkin(3) * t65 + t6 * pkin(4) + t23 * pkin(9) + t5 * qJ(5) + t51) 0, 0, 0, 0, 0, -g(1) * t75 - g(2) * t2, g(1) * t76 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t4, t44, 0, t4, -t44, -g(1) * (t23 * qJ(3) - t16) - g(2) * (t21 * qJ(3) - t14) - t59, 0, 0, 0, 0, 0, -t44 * t38, -t3, -g(1) * (-t22 * t35 + t23 * t62) - g(2) * (-t20 * t35 + t21 * t62) - (t35 * t42 + t37 * t61) * t72, -g(1) * (-t22 * t37 - t23 * t67) - g(2) * (-t20 * t37 - t21 * t67) - (-t35 * t61 + t37 * t42) * t72, t3, -g(1) * (-t22 * pkin(9) - t16) - g(2) * (-t20 * pkin(9) - t14) - t59 - (pkin(9) * t42 + t52 * t39) * t72 + t77 * (qJ(3) + t52) 0, 0, 0, 0, 0, -g(1) * (-t22 * t31 + t23 * t68) - g(2) * (-t20 * t31 + t21 * t68) - (t31 * t42 + t32 * t61) * t72, -g(1) * (-t22 * t32 - t23 * t69) - g(2) * (-t20 * t32 - t21 * t69) - (-t31 * t61 + t32 * t42) * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t46, t47 * t37, -t47 * t35, -t46, -g(1) * (-t5 * pkin(4) + t6 * qJ(5)) - g(2) * (pkin(4) * t49 - qJ(5) * t50) - g(3) * (-t18 * pkin(4) + t19 * qJ(5)) 0, 0, 0, 0, 0, t47 * t32, -t47 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t76 - g(3) * (-t19 * t31 + t32 * t66) g(1) * t2 - g(2) * t75 - g(3) * (-t19 * t32 - t31 * t66);];
taug_reg  = t7;
