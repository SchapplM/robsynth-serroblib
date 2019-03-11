% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t41 = sin(qJ(2));
t31 = t41 * qJ(3);
t78 = pkin(1) + t31;
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t22 = g(1) * t45 + g(2) * t42;
t44 = cos(qJ(2));
t9 = g(3) * t41 + t22 * t44;
t77 = pkin(2) + pkin(8);
t76 = pkin(2) * t41;
t75 = g(1) * t42;
t71 = g(3) * t44;
t35 = t44 * pkin(2);
t40 = sin(qJ(4));
t70 = t40 * t44;
t69 = t42 * t40;
t43 = cos(qJ(4));
t68 = t42 * t43;
t67 = t44 * t45;
t66 = t45 * t40;
t65 = t45 * t43;
t61 = qJ(3) * t44;
t24 = t42 * t61;
t58 = pkin(4) * t70;
t64 = t42 * t58 + t24;
t26 = t45 * t61;
t63 = t45 * t58 + t26;
t62 = t35 + t31;
t60 = qJ(5) * t40;
t59 = qJ(5) * t43;
t57 = -qJ(6) + t77;
t56 = t44 * t59;
t14 = -t41 * t65 + t69;
t15 = t41 * t66 + t68;
t54 = -t14 * pkin(4) + t15 * qJ(5);
t16 = t41 * t68 + t66;
t17 = -t41 * t69 + t65;
t53 = t16 * pkin(4) - t17 * qJ(5);
t52 = pkin(2) * t67 + t42 * pkin(7) + t78 * t45;
t51 = g(1) * t16 + g(2) * t14;
t50 = -g(2) * t45 + t75;
t49 = pkin(5) * t40 - t59;
t48 = g(3) * (t41 * t40 * pkin(4) + t44 * pkin(8) + t62);
t36 = t45 * pkin(7);
t47 = t45 * pkin(3) + t17 * pkin(4) + t16 * qJ(5) + t36;
t2 = g(1) * t14 - g(2) * t16 + t43 * t71;
t3 = -g(1) * t15 + g(2) * t17 + g(3) * t70;
t46 = t42 * pkin(3) + t15 * pkin(4) + pkin(8) * t67 + t14 * qJ(5) + t52;
t19 = -g(2) * t67 + t44 * t75;
t18 = t50 * t41;
t8 = t22 * t41 - t71;
t7 = t9 * t43;
t6 = t9 * t40;
t5 = -g(1) * t17 - g(2) * t15;
t1 = [0, t50, t22, 0, 0, 0, 0, 0, t19, -t18, -t22, -t19, t18, -g(1) * t36 - g(2) * t52 - (-t78 - t35) * t75, 0, 0, 0, 0, 0, t5, t51, t5, t19, -t51, -g(1) * t47 - g(2) * t46 - (-t77 * t44 - t78) * t75, t5, -t51, -t19, -g(1) * (t17 * pkin(5) + t47) - g(2) * (t15 * pkin(5) - qJ(6) * t67 + t46) - (-t57 * t44 - t78) * t75; 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, -t8, -t9, -g(1) * (-t45 * t76 + t26) - g(2) * (-t42 * t76 + t24) - g(3) * t62, 0, 0, 0, 0, 0, -t6, -t7, -t6, t8, t7, -g(1) * (-t45 * t56 + t63) - g(2) * (-t42 * t56 + t64) - t48 + (g(3) * t59 + t22 * t77) * t41, -t6, t7, -t8, -g(1) * t63 - g(2) * t64 - t48 + (g(3) * qJ(6) - t22 * t49) * t44 + (-g(3) * t49 + t22 * t57) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t3, t2, 0, t3, -g(1) * t54 - g(2) * t53 - (-pkin(4) * t43 - t60) * t71, t2, t3, 0, -g(1) * (-t14 * pkin(5) + t54) - g(2) * (t16 * pkin(5) + t53) - (-t60 + (-pkin(4) - pkin(5)) * t43) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9;];
taug_reg  = t1;
