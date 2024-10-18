% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR14_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:37
% EndTime: 2024-09-27 18:43:37
% DurationCPUTime: 0.25s
% Computational Cost: add. (532->49), mult. (274->77), div. (0->0), fcn. (276->20), ass. (0->58)
t62 = g(3) * sin(pkin(5));
t45 = cos(pkin(5));
t46 = sin(qJ(3));
t61 = t45 * t46;
t48 = cos(qJ(3));
t60 = t45 * t48;
t42 = qJ(3) + qJ(4);
t36 = pkin(5) + t42;
t33 = qJ(5) + t36;
t59 = sin(t33) / 0.2e1;
t58 = sin(t36) / 0.2e1;
t57 = pkin(5) - t42;
t56 = -qJ(5) + t57;
t55 = sin(t57);
t50 = sin(t56);
t19 = t59 - t50 / 0.2e1;
t41 = qJ(5) + t42;
t35 = cos(t41);
t43 = qJ(1) + qJ(2);
t38 = sin(t43);
t40 = cos(t43);
t54 = t40 * t19 + t38 * t35;
t53 = t38 * t19 - t40 * t35;
t21 = t58 - t55 / 0.2e1;
t39 = cos(t42);
t52 = t40 * t21 + t38 * t39;
t51 = t38 * t21 - t40 * t39;
t49 = cos(qJ(1));
t47 = sin(qJ(1));
t37 = sin(t42);
t34 = sin(t41);
t32 = cos(t57);
t31 = cos(t56);
t30 = cos(t36) / 0.2e1;
t27 = cos(t33) / 0.2e1;
t24 = g(1) * t40 + g(2) * t38;
t23 = g(1) * t38 - g(2) * t40;
t22 = t32 / 0.2e1 + t30;
t20 = t31 / 0.2e1 + t27;
t18 = -t38 * t61 + t40 * t48;
t17 = -t38 * t60 - t40 * t46;
t16 = -t38 * t48 - t40 * t61;
t15 = t38 * t46 - t40 * t60;
t14 = -t38 * t22 - t40 * t37;
t13 = -t40 * t22 + t38 * t37;
t12 = -t38 * t20 - t40 * t34;
t11 = -t40 * t20 + t38 * t34;
t10 = -g(1) * t16 - g(2) * t18;
t9 = -g(1) * t15 - g(2) * t17;
t8 = g(1) * t52 + g(2) * t51;
t7 = -g(1) * t13 - g(2) * t14;
t6 = g(1) * t54 + g(2) * t53;
t5 = -g(1) * t11 - g(2) * t12;
t4 = -g(1) * t51 + g(2) * t52 - g(3) * (t30 - t32 / 0.2e1);
t3 = -g(1) * t14 + g(2) * t13 - g(3) * (t58 + t55 / 0.2e1);
t2 = -g(1) * t53 + g(2) * t54 - g(3) * (t27 - t31 / 0.2e1);
t1 = -g(1) * t12 + g(2) * t11 - g(3) * (t59 + t50 / 0.2e1);
t25 = [0, g(1) * t47 - g(2) * t49, g(1) * t49 + g(2) * t47, 0, t23, t24, 0, 0, 0, 0, 0, t10, t9, 0, 0, 0, 0, 0, t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, t23, t24, 0, 0, 0, 0, 0, t10, t9, 0, 0, 0, 0, 0, t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t17 + g(2) * t15 - t48 * t62, g(1) * t18 - g(2) * t16 + t46 * t62, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t25;
