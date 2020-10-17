% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:05
% EndTime: 2019-12-31 19:33:06
% DurationCPUTime: 0.20s
% Computational Cost: add. (140->46), mult. (176->71), div. (0->0), fcn. (178->10), ass. (0->35)
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t6 = g(1) * t23 + g(2) * t21;
t16 = qJ(2) + pkin(8);
t10 = sin(t16);
t12 = cos(t16);
t26 = -g(3) * t12 + t6 * t10;
t38 = g(3) * t10;
t15 = pkin(9) + qJ(5);
t9 = sin(t15);
t36 = t21 * t9;
t35 = t23 * t9;
t11 = cos(t15);
t34 = t21 * t11;
t17 = sin(pkin(9));
t33 = t21 * t17;
t18 = cos(pkin(9));
t32 = t21 * t18;
t31 = t23 * t11;
t30 = t23 * t17;
t29 = t23 * t18;
t5 = g(1) * t21 - g(2) * t23;
t28 = t12 * pkin(3) + t10 * qJ(4);
t20 = sin(qJ(2));
t22 = cos(qJ(2));
t24 = -g(3) * t22 + t6 * t20;
t19 = -qJ(3) - pkin(6);
t13 = t22 * pkin(2);
t8 = t13 + pkin(1);
t7 = t23 * t8;
t4 = t12 * t31 + t36;
t3 = -t12 * t35 + t34;
t2 = -t12 * t34 + t35;
t1 = t12 * t36 + t31;
t14 = [0, t5, t6, 0, 0, 0, 0, 0, t5 * t22, -t5 * t20, -t6, -g(1) * (-t23 * t19 - t21 * t8) - g(2) * (-t21 * t19 + t7), -g(1) * (-t12 * t32 + t30) - g(2) * (t12 * t29 + t33), -g(1) * (t12 * t33 + t29) - g(2) * (-t12 * t30 + t32), t5 * t10, -g(2) * t7 + (g(1) * t19 - g(2) * t28) * t23 + (-g(1) * (-t28 - t8) + g(2) * t19) * t21, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t24, g(3) * t20 + t6 * t22, 0, t24 * pkin(2), t26 * t18, -t26 * t17, -t6 * t12 - t38, -g(3) * (t13 + t28) + t6 * (pkin(2) * t20 + pkin(3) * t10 - qJ(4) * t12), 0, 0, 0, 0, 0, t26 * t11, -t26 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t9 * t38, g(1) * t4 - g(2) * t2 + t11 * t38;];
taug_reg = t14;
