% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR5
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:43:18
% EndTime: 2019-05-07 10:43:19
% DurationCPUTime: 0.25s
% Computational Cost: add. (248->61), mult. (288->80), div. (0->0), fcn. (300->10), ass. (0->51)
t37 = cos(qJ(2));
t32 = qJ(2) + qJ(3);
t27 = sin(t32);
t29 = cos(t32);
t45 = t29 * pkin(3) + t27 * qJ(4);
t57 = t37 * pkin(2) + t45;
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t19 = g(1) * t38 + g(2) * t35;
t8 = g(3) * t27 + t19 * t29;
t56 = pkin(3) * t27;
t54 = g(3) * t29;
t31 = qJ(5) + qJ(6);
t26 = sin(t31);
t53 = t35 * t26;
t28 = cos(t31);
t52 = t35 * t28;
t33 = sin(qJ(5));
t51 = t35 * t33;
t36 = cos(qJ(5));
t50 = t35 * t36;
t49 = t38 * t26;
t48 = t38 * t28;
t47 = t38 * t33;
t46 = t38 * t36;
t44 = qJ(4) * t29;
t34 = sin(qJ(2));
t43 = -pkin(2) * t34 - t56;
t42 = g(1) * t35 - g(2) * t38;
t41 = pkin(1) + t57;
t39 = -pkin(8) - pkin(7);
t21 = t38 * t44;
t20 = t35 * t44;
t18 = -t27 * t51 + t46;
t17 = t27 * t50 + t47;
t16 = t27 * t47 + t50;
t15 = t27 * t46 - t51;
t14 = t42 * t29;
t13 = t42 * t27;
t12 = -t27 * t53 + t48;
t11 = t27 * t52 + t49;
t10 = t27 * t49 + t52;
t9 = t27 * t48 - t53;
t7 = t19 * t27 - t54;
t6 = t8 * t36;
t5 = t8 * t33;
t4 = t8 * t28;
t3 = t8 * t26;
t2 = g(1) * t10 - g(2) * t12 - t26 * t54;
t1 = -g(1) * t9 - g(2) * t11 + t28 * t54;
t22 = [0, t42, t19, 0, 0, 0, 0, 0, t42 * t37, -t42 * t34, 0, 0, 0, 0, 0, t14, -t13, -t19, -t14, t13 (g(1) * t39 - g(2) * t41) * t38 + (g(1) * t41 + g(2) * t39) * t35, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t16, g(1) * t17 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t37 + t19 * t34, g(3) * t34 + t19 * t37, 0, 0, 0, 0, 0, t7, t8, 0, -t7, -t8, -g(1) * (t43 * t38 + t21) - g(2) * (t43 * t35 + t20) - g(3) * t57, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, -t7, -t8, -g(1) * (-t38 * t56 + t21) - g(2) * (-t35 * t56 + t20) - g(3) * t45, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17 + t36 * t54, g(1) * t16 - g(2) * t18 - t33 * t54, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t22;
