% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:41:23
% EndTime: 2021-01-15 20:41:25
% DurationCPUTime: 0.21s
% Computational Cost: add. (187->56), mult. (278->76), div. (0->0), fcn. (291->10), ass. (0->42)
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t49 = -g(1) * t34 - g(2) * t31;
t25 = qJ(2) + pkin(8);
t21 = sin(t25);
t46 = g(3) * t21;
t28 = -qJ(3) - pkin(6);
t45 = t31 * t28;
t29 = sin(qJ(4));
t44 = t31 * t29;
t32 = cos(qJ(4));
t43 = t31 * t32;
t42 = t34 * t29;
t41 = t34 * t32;
t22 = cos(t25);
t11 = t22 * t42 - t43;
t9 = t22 * t44 + t41;
t40 = g(1) * t9 - g(2) * t11;
t15 = g(1) * t31 - g(2) * t34;
t39 = -pkin(4) * t32 - qJ(5) * t29;
t1 = g(1) * t11 + g(2) * t9 + t29 * t46;
t10 = t22 * t43 - t42;
t12 = t22 * t41 + t44;
t37 = g(1) * t12 + g(2) * t10 + t32 * t46;
t36 = -g(3) * t22 - t49 * t21;
t30 = sin(qJ(2));
t33 = cos(qJ(2));
t35 = -g(3) * t33 - t30 * t49;
t27 = cos(pkin(8));
t26 = sin(pkin(8));
t23 = t33 * pkin(2);
t20 = t23 + pkin(1);
t17 = t28 * t34;
t14 = -t26 * pkin(3) + t27 * pkin(7);
t13 = t27 * pkin(3) + t26 * pkin(7) + pkin(2);
t8 = t15 * t21;
t7 = -t22 * t49 + t46;
t5 = t13 * t33 + t14 * t30 + pkin(1);
t4 = t36 * t32;
t3 = t36 * t29;
t2 = g(1) * t10 - g(2) * t12;
t6 = [0, t15, -t49, 0, 0, 0, 0, 0, t15 * t33, -t15 * t30, t15 * t22, -t8, t49, -g(1) * (-t31 * t20 - t17) - g(2) * (t34 * t20 - t45), 0, 0, 0, 0, 0, t2, -t40, t2, t8, t40, -g(1) * (-t10 * pkin(4) - t9 * qJ(5) - t5 * t31 - t17) - g(2) * (t12 * pkin(4) + t11 * qJ(5) + t5 * t34 - t45); 0, 0, 0, 0, 0, 0, 0, 0, t35, g(3) * t30 - t33 * t49, t36, t7, 0, t35 * pkin(2), 0, 0, 0, 0, 0, t4, -t3, t4, -t7, t3, -g(3) * (t21 * pkin(7) + t23 + (pkin(3) - t39) * t22) + t49 * (-t13 * t30 + t14 * t33 + t21 * t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t37, t1, 0, -t37, -g(1) * (-t11 * pkin(4) + t12 * qJ(5)) - g(2) * (-t9 * pkin(4) + t10 * qJ(5)) - (-pkin(4) * t29 + qJ(5) * t32) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t6;
