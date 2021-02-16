% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:46:37
% EndTime: 2021-01-16 03:46:39
% DurationCPUTime: 0.39s
% Computational Cost: add. (352->86), mult. (588->156), div. (0->0), fcn. (732->14), ass. (0->47)
t22 = sin(pkin(11));
t27 = sin(qJ(2));
t43 = cos(pkin(11));
t44 = cos(pkin(6));
t36 = t44 * t43;
t52 = cos(qJ(2));
t11 = t22 * t52 + t27 * t36;
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t23 = sin(pkin(6));
t37 = t23 * t43;
t45 = t23 * t29;
t46 = t23 * t27;
t38 = t22 * t44;
t9 = t27 * t38 - t43 * t52;
t54 = -g(1) * (t22 * t45 + t9 * t26) - g(2) * (-t11 * t26 - t29 * t37) - g(3) * (-t26 * t46 + t44 * t29);
t53 = g(3) * t23;
t20 = qJ(3) + pkin(12);
t17 = cos(t20);
t21 = qJ(5) + qJ(6);
t18 = sin(t21);
t51 = t17 * t18;
t19 = cos(t21);
t50 = t17 * t19;
t25 = sin(qJ(5));
t49 = t17 * t25;
t28 = cos(qJ(5));
t48 = t17 * t28;
t47 = t22 * t23;
t42 = t17 * t52;
t41 = t23 * t52;
t40 = t25 * t52;
t39 = t28 * t52;
t16 = sin(t20);
t6 = t16 * t47 - t9 * t17;
t35 = g(1) * (t9 * t16 + t17 * t47) + g(2) * (-t11 * t16 - t17 * t37) + g(3) * (-t16 * t46 + t44 * t17);
t33 = -g(1) * t9 + g(2) * t11 + g(3) * t46;
t10 = t22 * t27 - t52 * t36;
t12 = t43 * t27 + t52 * t38;
t30 = g(1) * t12 + g(2) * t10 - g(3) * t41;
t24 = qJ(4) + pkin(8);
t15 = t29 * pkin(3) + pkin(2);
t8 = t44 * t16 + t17 * t46;
t4 = t11 * t17 - t16 * t37;
t2 = -g(1) * (-t12 * t18 - t6 * t19) - g(2) * (-t10 * t18 - t4 * t19) - g(3) * (t18 * t41 - t8 * t19);
t1 = -g(1) * (t12 * t19 - t6 * t18) - g(2) * (t10 * t19 - t4 * t18) - g(3) * (-t8 * t18 - t19 * t41);
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t30, t33, 0, 0, 0, 0, 0, t30 * t29, -t30 * t26, t30 * t17, -t30 * t16, -t33, -g(1) * (-t12 * t15 - t9 * t24) - g(2) * (-t10 * t15 + t11 * t24) - (t52 * t15 + t24 * t27) * t53, 0, 0, 0, 0, 0, -g(1) * (-t12 * t48 - t9 * t25) - g(2) * (-t10 * t48 + t11 * t25) - (t17 * t39 + t25 * t27) * t53, -g(1) * (t12 * t49 - t9 * t28) - g(2) * (t10 * t49 + t11 * t28) - (-t17 * t40 + t27 * t28) * t53, 0, 0, 0, 0, 0, -g(1) * (-t12 * t50 - t9 * t18) - g(2) * (-t10 * t50 + t11 * t18) - (t18 * t27 + t19 * t42) * t53, -g(1) * (t12 * t51 - t9 * t19) - g(2) * (t10 * t51 + t11 * t19) - (-t18 * t42 + t19 * t27) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -g(1) * (-t26 * t47 + t9 * t29) - g(2) * (-t11 * t29 + t26 * t37) - g(3) * (-t44 * t26 - t27 * t45), -t35, g(1) * t6 + g(2) * t4 + g(3) * t8, 0, t54 * pkin(3), 0, 0, 0, 0, 0, -t35 * t28, t35 * t25, 0, 0, 0, 0, 0, -t35 * t19, t35 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t28 - t6 * t25) - g(2) * (t10 * t28 - t4 * t25) - g(3) * (-t23 * t39 - t8 * t25), -g(1) * (-t12 * t25 - t6 * t28) - g(2) * (-t10 * t25 - t4 * t28) - g(3) * (t23 * t40 - t8 * t28), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t3;
