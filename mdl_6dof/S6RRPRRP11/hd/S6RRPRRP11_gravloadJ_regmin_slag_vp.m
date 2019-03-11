% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t36 = sin(qJ(2));
t26 = t36 * qJ(3);
t39 = cos(qJ(2));
t47 = t39 * pkin(2) + t26;
t37 = sin(qJ(1));
t40 = cos(qJ(1));
t18 = g(1) * t40 + g(2) * t37;
t34 = qJ(4) + qJ(5);
t25 = cos(t34);
t50 = t40 * t25;
t24 = sin(t34);
t58 = t37 * t24;
t3 = t36 * t50 - t58;
t51 = t40 * t24;
t57 = t37 * t25;
t5 = t36 * t57 + t51;
t63 = g(3) * t39;
t1 = -g(1) * t3 - g(2) * t5 + t25 * t63;
t8 = g(3) * t36 + t18 * t39;
t67 = g(1) * t37;
t61 = t36 * t37;
t60 = t36 * t40;
t35 = sin(qJ(4));
t16 = t35 * pkin(4) + pkin(5) * t24;
t59 = t37 * t16;
t56 = t37 * t35;
t38 = cos(qJ(4));
t55 = t37 * t38;
t33 = -qJ(6) - pkin(9) - pkin(8);
t54 = t39 * t33;
t53 = t39 * t40;
t52 = t40 * t16;
t49 = t40 * t35;
t48 = t40 * t38;
t17 = t38 * pkin(4) + pkin(5) * t25;
t46 = qJ(3) * t39;
t44 = pkin(2) * t53 + t37 * pkin(7) + (pkin(1) + t26) * t40;
t43 = -g(2) * t40 + t67;
t42 = -pkin(1) - t47;
t30 = t40 * pkin(7);
t21 = t40 * t46;
t19 = t37 * t46;
t15 = pkin(3) + t17;
t14 = t43 * t39;
t13 = t43 * t36;
t12 = -t36 * t56 + t48;
t11 = t36 * t55 + t49;
t10 = t36 * t49 + t55;
t9 = t36 * t48 - t56;
t7 = t18 * t36 - t63;
t6 = -t36 * t58 + t50;
t4 = t36 * t51 + t57;
t2 = g(1) * t4 - g(2) * t6 - t24 * t63;
t20 = [0, t43, t18, 0, 0, 0, 0, 0, t14, -t13, -t18, -t14, t13, -g(1) * t30 - g(2) * t44 - t42 * t67, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t14, -g(1) * (t40 * t15 + t30) - g(2) * (-t33 * t53 + t36 * t52 + t44) + (-g(1) * (-t36 * t16 + t42 + t54) - g(2) * t15) * t37; 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, -t7, -t8, -g(1) * (-pkin(2) * t60 + t21) - g(2) * (-pkin(2) * t61 + t19) - g(3) * t47, 0, 0, 0, 0, 0, -t8 * t35, -t8 * t38, 0, 0, 0, 0, 0, -t8 * t24, -t8 * t25, t7, -g(1) * (t39 * t52 + t21) - g(2) * (t39 * t59 + t19) - g(3) * (t47 - t54) + (-g(3) * t16 + t18 * (pkin(2) - t33)) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11 + t38 * t63, g(1) * t10 - g(2) * t12 - t35 * t63, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t17 * t60 - t59) - g(2) * (t17 * t61 + t52) + t17 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg  = t20;
