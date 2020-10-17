% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR14_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:37
% EndTime: 2019-12-31 20:38:39
% DurationCPUTime: 0.30s
% Computational Cost: add. (235->72), mult. (463->131), div. (0->0), fcn. (574->12), ass. (0->41)
t24 = sin(qJ(2));
t25 = sin(qJ(1));
t27 = cos(qJ(2));
t28 = cos(qJ(1));
t36 = cos(pkin(5));
t34 = t28 * t36;
t10 = t25 * t24 - t27 * t34;
t23 = sin(qJ(5));
t26 = cos(qJ(5));
t11 = t24 * t34 + t25 * t27;
t19 = pkin(10) + qJ(4);
t17 = sin(t19);
t18 = cos(t19);
t21 = sin(pkin(5));
t39 = t21 * t28;
t4 = t11 * t18 - t17 * t39;
t48 = -t10 * t26 + t4 * t23;
t47 = t10 * t23 + t4 * t26;
t46 = g(3) * t21;
t43 = t18 * t23;
t42 = t18 * t26;
t41 = t21 * t24;
t40 = t21 * t25;
t38 = t23 * t27;
t37 = t26 * t27;
t35 = t25 * t36;
t12 = t28 * t24 + t27 * t35;
t33 = g(1) * t10 - g(2) * t12;
t32 = t11 * t17 + t18 * t39;
t13 = -t24 * t35 + t28 * t27;
t6 = -t13 * t17 + t18 * t40;
t31 = g(1) * t6 - g(2) * t32 + g(3) * (-t17 * t41 + t36 * t18);
t30 = -g(1) * t12 - g(2) * t10 + t27 * t46;
t29 = g(1) * t13 + g(2) * t11 + g(3) * t41;
t22 = cos(pkin(10));
t20 = sin(pkin(10));
t9 = t36 * t17 + t18 * t41;
t7 = t13 * t18 + t17 * t40;
t2 = t12 * t23 + t7 * t26;
t1 = t12 * t26 - t7 * t23;
t3 = [0, g(1) * t25 - g(2) * t28, g(1) * t28 + g(2) * t25, 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t13, -t33, -g(1) * (-t11 * t22 + t20 * t39) - g(2) * (t13 * t22 + t20 * t40), -g(1) * (t11 * t20 + t22 * t39) - g(2) * (-t13 * t20 + t22 * t40), t33, -g(1) * (-t25 * pkin(1) - t11 * pkin(2) + pkin(7) * t39 - t10 * qJ(3)) - g(2) * (t28 * pkin(1) + t13 * pkin(2) + pkin(7) * t40 + t12 * qJ(3)), 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t7, -g(1) * t32 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t47 - g(2) * t2, -g(1) * t48 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, -t30 * t22, t30 * t20, -t29, -g(1) * (-t12 * pkin(2) + t13 * qJ(3)) - g(2) * (-t10 * pkin(2) + t11 * qJ(3)) - (pkin(2) * t27 + qJ(3) * t24) * t46, 0, 0, 0, 0, 0, -t30 * t18, t30 * t17, 0, 0, 0, 0, 0, -g(1) * (-t12 * t42 + t13 * t23) - g(2) * (-t10 * t42 + t11 * t23) - (t18 * t37 + t23 * t24) * t46, -g(1) * (t12 * t43 + t13 * t26) - g(2) * (t10 * t43 + t11 * t26) - (-t18 * t38 + t24 * t26) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, g(1) * t7 + g(2) * t4 + g(3) * t9, 0, 0, 0, 0, 0, -t31 * t26, t31 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t48 - g(3) * (-t21 * t37 - t9 * t23), g(1) * t2 + g(2) * t47 - g(3) * (t21 * t38 - t9 * t26);];
taug_reg = t3;
