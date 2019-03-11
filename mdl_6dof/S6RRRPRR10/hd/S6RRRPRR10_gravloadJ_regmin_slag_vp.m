% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t31 = sin(qJ(1));
t33 = cos(qJ(3));
t34 = cos(qJ(2));
t29 = sin(qJ(3));
t35 = cos(qJ(1));
t52 = t35 * t29;
t20 = -t31 * t33 + t34 * t52;
t51 = t35 * t33;
t21 = t31 * t29 + t34 * t51;
t27 = qJ(5) + qJ(6);
t25 = sin(t27);
t26 = cos(t27);
t3 = t20 * t26 - t21 * t25;
t43 = t25 * t33 - t26 * t29;
t30 = sin(qJ(2));
t57 = g(3) * t30;
t53 = t31 * t34;
t18 = t29 * t53 + t51;
t19 = t33 * t53 - t52;
t73 = t18 * t26 - t19 * t25;
t77 = g(1) * t3 + g(2) * t73 - t43 * t57;
t28 = sin(qJ(5));
t32 = cos(qJ(5));
t41 = t28 * t33 - t29 * t32;
t6 = t20 * t32 - t21 * t28;
t74 = t18 * t32 - t19 * t28;
t76 = g(1) * t6 + g(2) * t74 - t41 * t57;
t47 = g(1) * t35 + g(2) * t31;
t68 = g(3) * t34 - t47 * t30;
t40 = t28 * t29 + t32 * t33;
t44 = t18 * t28 + t19 * t32;
t7 = t20 * t28 + t21 * t32;
t70 = g(1) * t7 + g(2) * t44 + t40 * t57;
t4 = t20 * t25 + t21 * t26;
t42 = t25 * t29 + t26 * t33;
t45 = t18 * t25 + t19 * t26;
t2 = g(1) * t4 + g(2) * t45 + t42 * t57;
t48 = g(1) * t18 - g(2) * t20;
t46 = g(1) * t31 - g(2) * t35;
t39 = t34 * pkin(2) + t30 * pkin(8) + pkin(1);
t38 = pkin(3) * t33 + qJ(4) * t29 + pkin(2);
t5 = g(1) * t20 + g(2) * t18 + t29 * t57;
t37 = g(1) * t21 + g(2) * t19 + t33 * t57;
t22 = t46 * t30;
t14 = t47 * t34 + t57;
t10 = t68 * t33;
t9 = t68 * t29;
t8 = g(1) * t19 - g(2) * t21;
t1 = [0, t46, t47, 0, 0, 0, 0, 0, t46 * t34, -t22, 0, 0, 0, 0, 0, t8, -t48, t8, t22, t48, -g(1) * (-t19 * pkin(3) - t18 * qJ(4)) - g(2) * (t21 * pkin(3) + t20 * qJ(4)) + (-g(1) * pkin(7) - g(2) * t39) * t35 + (-g(2) * pkin(7) + g(1) * t39) * t31, 0, 0, 0, 0, 0, g(1) * t44 - g(2) * t7, g(1) * t74 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t45 - g(2) * t4, g(1) * t73 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, -t68, t14, 0, 0, 0, 0, 0, -t10, t9, -t10, -t14, -t9 (-t47 * pkin(8) - g(3) * t38) * t34 + (-g(3) * pkin(8) + t47 * t38) * t30, 0, 0, 0, 0, 0, -t68 * t40, t68 * t41, 0, 0, 0, 0, 0, -t68 * t42, t68 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t37, t5, 0, -t37, -g(1) * (-t20 * pkin(3) + t21 * qJ(4)) - g(2) * (-t18 * pkin(3) + t19 * qJ(4)) - (-pkin(3) * t29 + qJ(4) * t33) * t57, 0, 0, 0, 0, 0, t76, -t70, 0, 0, 0, 0, 0, t77, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t70, 0, 0, 0, 0, 0, -t77, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, t2;];
taug_reg  = t1;
