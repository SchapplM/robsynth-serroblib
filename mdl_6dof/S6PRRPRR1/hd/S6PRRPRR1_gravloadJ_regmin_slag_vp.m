% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRR1
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
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:28:23
% EndTime: 2021-01-16 03:28:24
% DurationCPUTime: 0.33s
% Computational Cost: add. (352->78), mult. (510->132), div. (0->0), fcn. (624->14), ass. (0->46)
t26 = cos(pkin(6));
t45 = cos(pkin(11));
t52 = cos(qJ(2));
t39 = t45 * t52;
t24 = sin(pkin(11));
t30 = sin(qJ(2));
t48 = t24 * t30;
t11 = t26 * t48 - t39;
t40 = t45 * t30;
t44 = t24 * t52;
t13 = t26 * t40 + t44;
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t25 = sin(pkin(6));
t41 = t25 * t45;
t46 = t25 * t32;
t47 = t25 * t30;
t54 = -g(1) * (t11 * t29 + t24 * t46) - g(2) * (-t13 * t29 - t32 * t41) - g(3) * (t26 * t32 - t29 * t47);
t53 = g(3) * t25;
t23 = qJ(3) + pkin(12);
t22 = qJ(5) + t23;
t18 = cos(t22);
t28 = sin(qJ(6));
t51 = t18 * t28;
t31 = cos(qJ(6));
t50 = t18 * t31;
t49 = t24 * t25;
t43 = t28 * t52;
t42 = t31 * t52;
t17 = sin(t22);
t38 = g(1) * (t11 * t17 + t18 * t49) + g(2) * (-t13 * t17 - t18 * t41) + g(3) * (-t17 * t47 + t26 * t18);
t35 = -g(1) * t11 + g(2) * t13 + g(3) * t47;
t12 = -t26 * t39 + t48;
t14 = t26 * t44 + t40;
t33 = g(1) * t14 + g(2) * t12 - t52 * t53;
t27 = qJ(4) + pkin(8);
t21 = cos(t23);
t20 = sin(t23);
t19 = t32 * pkin(3) + pkin(2);
t10 = t26 * t17 + t18 * t47;
t8 = -t11 * t18 + t17 * t49;
t6 = t13 * t18 - t17 * t41;
t4 = g(1) * t8 + g(2) * t6 + g(3) * t10;
t2 = t38 * t31;
t1 = t38 * t28;
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t33, t35, 0, 0, 0, 0, 0, t33 * t32, -t33 * t29, t33 * t21, -t33 * t20, -t35, -g(1) * (-t11 * t27 - t14 * t19) - g(2) * (-t12 * t19 + t13 * t27) - (t52 * t19 + t27 * t30) * t53, 0, 0, 0, 0, 0, t33 * t18, -t33 * t17, 0, 0, 0, 0, 0, -g(1) * (-t11 * t28 - t14 * t50) - g(2) * (-t12 * t50 + t13 * t28) - (t18 * t42 + t28 * t30) * t53, -g(1) * (-t11 * t31 + t14 * t51) - g(2) * (t12 * t51 + t13 * t31) - (-t18 * t43 + t30 * t31) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -g(1) * (t11 * t32 - t29 * t49) - g(2) * (-t13 * t32 + t29 * t41) - g(3) * (-t26 * t29 - t30 * t46), -g(1) * (t11 * t20 + t21 * t49) - g(2) * (-t13 * t20 - t21 * t41) - g(3) * (-t20 * t47 + t26 * t21), -g(1) * (t11 * t21 - t20 * t49) - g(2) * (-t13 * t21 + t20 * t41) - g(3) * (-t26 * t20 - t21 * t47), 0, t54 * pkin(3), 0, 0, 0, 0, 0, -t38, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t31 - t8 * t28) - g(2) * (t12 * t31 - t6 * t28) - g(3) * (-t10 * t28 - t25 * t42), -g(1) * (-t14 * t28 - t8 * t31) - g(2) * (-t12 * t28 - t6 * t31) - g(3) * (-t10 * t31 + t25 * t43);];
taug_reg = t3;
