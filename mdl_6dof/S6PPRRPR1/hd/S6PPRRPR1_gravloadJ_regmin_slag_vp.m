% Calculate minimal parameter regressor of gravitation load for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% taug_reg [6x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t58 = sin(pkin(12));
t59 = sin(pkin(11));
t49 = t59 * t58;
t62 = cos(pkin(12));
t63 = cos(pkin(11));
t56 = t63 * t62;
t65 = cos(pkin(6));
t39 = -t65 * t56 + t49;
t60 = sin(pkin(7));
t61 = sin(pkin(6));
t53 = t61 * t60;
t64 = cos(pkin(7));
t76 = t39 * t64 + t63 * t53;
t51 = t59 * t62;
t54 = t63 * t58;
t40 = t65 * t51 + t54;
t50 = t59 * t61;
t75 = t40 * t64 - t60 * t50;
t74 = t62 * t64 * t61 + t65 * t60;
t32 = sin(qJ(3));
t52 = t61 * t58;
t70 = cos(qJ(3));
t14 = t32 * t52 - t74 * t70;
t20 = t65 * t54 + t51;
t6 = t20 * t32 + t76 * t70;
t21 = -t65 * t49 + t56;
t8 = t21 * t32 + t75 * t70;
t45 = g(1) * t8 + g(2) * t6 + g(3) * t14;
t28 = pkin(13) + qJ(6);
t26 = sin(t28);
t33 = cos(qJ(4));
t69 = t26 * t33;
t27 = cos(t28);
t68 = t27 * t33;
t29 = sin(pkin(13));
t67 = t29 * t33;
t30 = cos(pkin(13));
t66 = t30 * t33;
t55 = t63 * t61;
t15 = t74 * t32 + t70 * t52;
t31 = sin(qJ(4));
t38 = -t62 * t53 + t65 * t64;
t10 = t15 * t31 - t38 * t33;
t34 = t39 * t60 - t64 * t55;
t7 = t20 * t70 - t76 * t32;
t2 = t7 * t31 - t34 * t33;
t35 = t40 * t60 + t64 * t50;
t9 = t21 * t70 - t75 * t32;
t4 = t9 * t31 - t35 * t33;
t47 = g(1) * t4 + g(2) * t2 + g(3) * t10;
t11 = t15 * t33 + t38 * t31;
t3 = t34 * t31 + t7 * t33;
t5 = t35 * t31 + t9 * t33;
t46 = g(1) * t5 + g(2) * t3 + g(3) * t11;
t44 = g(1) * t9 + g(2) * t7 + g(3) * t15;
t18 = -g(1) * t50 + g(2) * t55 - g(3) * t65;
t1 = t45 * t31;
t12 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t45, t44, 0, 0, 0, 0, 0, t45 * t33, -t1, -g(1) * (t9 * t29 - t8 * t66) - g(2) * (t7 * t29 - t6 * t66) - g(3) * (-t14 * t66 + t15 * t29) -g(1) * (t9 * t30 + t8 * t67) - g(2) * (t7 * t30 + t6 * t67) - g(3) * (t14 * t67 + t15 * t30) t1, -t44 * pkin(9) + t45 * (pkin(4) * t33 + qJ(5) * t31 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (t9 * t26 - t8 * t68) - g(2) * (t7 * t26 - t6 * t68) - g(3) * (-t14 * t68 + t15 * t26) -g(1) * (t9 * t27 + t8 * t69) - g(2) * (t7 * t27 + t6 * t69) - g(3) * (t14 * t69 + t15 * t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t46, t47 * t30, -t47 * t29, -t46, -g(1) * (-t4 * pkin(4) + t5 * qJ(5)) - g(2) * (-t2 * pkin(4) + t3 * qJ(5)) - g(3) * (-t10 * pkin(4) + t11 * qJ(5)) 0, 0, 0, 0, 0, t47 * t27, -t47 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t5 * t26 + t8 * t27) - g(2) * (-t3 * t26 + t6 * t27) - g(3) * (-t11 * t26 + t14 * t27) -g(1) * (-t8 * t26 - t5 * t27) - g(2) * (-t6 * t26 - t3 * t27) - g(3) * (-t11 * t27 - t14 * t26);];
taug_reg  = t12;
