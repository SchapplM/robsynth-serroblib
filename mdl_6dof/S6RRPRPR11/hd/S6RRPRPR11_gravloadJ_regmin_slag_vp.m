% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t32 = sin(qJ(2));
t24 = t32 * qJ(3);
t35 = cos(qJ(2));
t45 = t35 * pkin(2) + t24;
t33 = sin(qJ(1));
t36 = cos(qJ(1));
t15 = g(1) * t36 + g(2) * t33;
t31 = sin(qJ(4));
t34 = cos(qJ(4));
t48 = t33 * t34;
t11 = t31 * t36 + t32 * t48;
t56 = g(3) * t35;
t47 = t34 * t36;
t49 = t33 * t31;
t9 = t32 * t47 - t49;
t64 = -g(1) * t9 - g(2) * t11 + t34 * t56;
t8 = g(3) * t32 + t15 * t35;
t62 = pkin(4) * t31;
t61 = g(1) * t33;
t30 = -qJ(5) - pkin(8);
t54 = t30 * t35;
t53 = t31 * t35;
t52 = t32 * t36;
t23 = qJ(4) + pkin(10) + qJ(6);
t20 = sin(t23);
t51 = t33 * t20;
t21 = cos(t23);
t50 = t33 * t21;
t46 = t35 * t36;
t44 = qJ(3) * t35;
t43 = pkin(4) * t53;
t41 = t31 * t52;
t40 = pkin(2) * t46 + t33 * pkin(7) + (pkin(1) + t24) * t36;
t39 = -g(2) * t36 + t61;
t38 = -pkin(1) - t45;
t27 = t36 * pkin(7);
t22 = pkin(4) * t34 + pkin(3);
t18 = t36 * t44;
t16 = t33 * t44;
t14 = t39 * t35;
t13 = t39 * t32;
t12 = -t32 * t49 + t47;
t10 = t41 + t48;
t7 = t15 * t32 - t56;
t6 = t21 * t36 - t32 * t51;
t5 = t20 * t36 + t32 * t50;
t4 = t20 * t52 + t50;
t3 = t21 * t52 - t51;
t2 = g(1) * t4 - g(2) * t6 - t20 * t56;
t1 = -g(1) * t3 - g(2) * t5 + t21 * t56;
t17 = [0, t39, t15, 0, 0, 0, 0, 0, t14, -t13, -t15, -t14, t13, -g(1) * t27 - g(2) * t40 - t38 * t61, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t14, -g(1) * (t22 * t36 + t27) - g(2) * (pkin(4) * t41 - t30 * t46 + t40) + (-g(1) * (-t32 * t62 + t38 + t54) - g(2) * t22) * t33, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, -t7, -t8, -g(1) * (-pkin(2) * t52 + t18) - g(2) * (-pkin(2) * t32 * t33 + t16) - g(3) * t45, 0, 0, 0, 0, 0, -t8 * t31, -t8 * t34, t7, -g(1) * (t36 * t43 + t18) - g(2) * (t33 * t43 + t16) - g(3) * (t45 - t54) + (-g(3) * t62 + t15 * (pkin(2) - t30)) * t32, 0, 0, 0, 0, 0, -t8 * t20, -t8 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, g(1) * t10 - g(2) * t12 - g(3) * t53, 0, t64 * pkin(4), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t17;
