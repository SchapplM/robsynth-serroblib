% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t41 = sin(pkin(11));
t46 = sin(qJ(2));
t50 = cos(qJ(2));
t68 = cos(pkin(11));
t32 = -t41 * t50 - t46 * t68;
t47 = sin(qJ(1));
t51 = cos(qJ(1));
t43 = cos(pkin(6));
t56 = -t41 * t46 + t50 * t68;
t52 = t56 * t43;
t13 = t47 * t32 + t51 * t52;
t44 = sin(qJ(6));
t48 = cos(qJ(6));
t25 = t32 * t43;
t14 = -t25 * t51 + t47 * t56;
t45 = sin(qJ(4));
t49 = cos(qJ(4));
t42 = sin(pkin(6));
t76 = t42 * t51;
t6 = t14 * t45 + t49 * t76;
t88 = -t13 * t48 + t44 * t6;
t87 = t13 * t44 + t48 * t6;
t70 = t51 * t46;
t71 = t47 * t50;
t29 = -t43 * t71 - t70;
t77 = t42 * t50;
t86 = -g(1) * t29 - g(3) * t77;
t16 = t32 * t51 - t47 * t52;
t23 = t56 * t42;
t85 = -g(1) * t16 - g(2) * t13 - g(3) * t23;
t78 = t42 * t47;
t75 = t44 * t45;
t74 = t45 * t48;
t72 = t47 * t46;
t69 = t51 * t50;
t65 = t43 * t69;
t7 = t14 * t49 - t45 * t76;
t26 = t43 * t46 * pkin(2) + (-pkin(8) - qJ(3)) * t42;
t40 = pkin(2) * t50 + pkin(1);
t64 = -t26 * t47 + t40 * t51;
t15 = -t25 * t47 - t51 * t56;
t10 = -t15 * t45 - t49 * t78;
t62 = -g(1) * t6 + g(2) * t10;
t11 = -t15 * t49 + t45 * t78;
t61 = -g(1) * t7 + g(2) * t11;
t60 = g(1) * t51 + g(2) * t47;
t59 = g(1) * t47 - g(2) * t51;
t58 = -t26 * t51 - t40 * t47;
t24 = t32 * t42;
t18 = -t24 * t45 - t43 * t49;
t55 = g(1) * t10 + g(2) * t6 + g(3) * t18;
t19 = -t24 * t49 + t43 * t45;
t54 = g(1) * t11 + g(2) * t7 + g(3) * t19;
t33 = pkin(2) * t65;
t30 = -t43 * t72 + t69;
t28 = -t43 * t70 - t71;
t27 = -t65 + t72;
t22 = -g(3) * t43 - t42 * t59;
t5 = t10 * t44 - t16 * t48;
t4 = t10 * t48 + t16 * t44;
t3 = t85 * t49;
t2 = t85 * t45;
t1 = [0, t59, t60, 0, 0, 0, 0, 0, -g(1) * t28 - g(2) * t30, -g(1) * t27 - g(2) * t29, -t60 * t42, -g(1) * t58 - g(2) * t64, 0, 0, 0, 0, 0, -t61, t62, -g(1) * t13 + g(2) * t16, t61, -t62, -g(1) * (-pkin(3) * t14 - pkin(4) * t7 + pkin(9) * t13 - qJ(5) * t6 + t58) - g(2) * (-pkin(3) * t15 + pkin(4) * t11 - pkin(9) * t16 + qJ(5) * t10 + t64) 0, 0, 0, 0, 0, g(1) * t88 - g(2) * t5, g(1) * t87 - g(2) * t4; 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t27 + t86, g(3) * t42 * t46 + g(1) * t30 - g(2) * t28, 0, -g(2) * t33 + (g(2) * t72 + t86) * pkin(2), 0, 0, 0, 0, 0, t3, -t2, g(1) * t15 - g(2) * t14 + g(3) * t24, -t3, t2, -g(1) * (pkin(2) * t29 - t15 * pkin(9)) - g(2) * (-pkin(2) * t72 + pkin(9) * t14 + t33) - g(3) * (pkin(2) * t77 - t24 * pkin(9)) + t85 * (pkin(4) * t49 + qJ(5) * t45 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (-t15 * t48 + t16 * t75) - g(2) * (t13 * t75 + t14 * t48) - g(3) * (t23 * t75 - t24 * t48) -g(1) * (t15 * t44 + t16 * t74) - g(2) * (t13 * t74 - t14 * t44) - g(3) * (t23 * t74 + t24 * t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t54, 0, -t55, -t54, -g(1) * (-pkin(4) * t10 + qJ(5) * t11) - g(2) * (-pkin(4) * t6 + qJ(5) * t7) - g(3) * (-pkin(4) * t18 + qJ(5) * t19) 0, 0, 0, 0, 0, -t54 * t44, -t54 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t87 - g(3) * (t18 * t48 + t23 * t44) g(1) * t5 + g(2) * t88 - g(3) * (-t18 * t44 + t23 * t48);];
taug_reg  = t1;
