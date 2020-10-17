% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:58
% EndTime: 2019-12-31 20:21:59
% DurationCPUTime: 0.15s
% Computational Cost: add. (135->37), mult. (166->61), div. (0->0), fcn. (181->10), ass. (0->38)
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t12 = g(1) * t26 + g(2) * t23;
t18 = qJ(2) + pkin(9);
t14 = sin(t18);
t15 = cos(t18);
t29 = -g(3) * t15 + t12 * t14;
t39 = g(3) * t14;
t19 = qJ(4) + qJ(5);
t16 = sin(t19);
t37 = t23 * t16;
t17 = cos(t19);
t36 = t23 * t17;
t21 = sin(qJ(4));
t35 = t23 * t21;
t24 = cos(qJ(4));
t34 = t23 * t24;
t33 = t26 * t16;
t32 = t26 * t17;
t31 = t26 * t21;
t30 = t26 * t24;
t11 = g(1) * t23 - g(2) * t26;
t22 = sin(qJ(2));
t25 = cos(qJ(2));
t27 = -g(3) * t25 + t12 * t22;
t20 = -qJ(3) - pkin(6);
t13 = t25 * pkin(2) + pkin(1);
t10 = t15 * t30 + t35;
t9 = -t15 * t31 + t34;
t8 = -t15 * t34 + t31;
t7 = t15 * t35 + t30;
t6 = t15 * t32 + t37;
t5 = -t15 * t33 + t36;
t4 = -t15 * t36 + t33;
t3 = t15 * t37 + t32;
t2 = g(1) * t6 - g(2) * t4 + t17 * t39;
t1 = -g(1) * t5 + g(2) * t3 + t16 * t39;
t28 = [0, t11, t12, 0, 0, 0, 0, 0, t11 * t25, -t11 * t22, -t12, -g(1) * (-t23 * t13 - t26 * t20) - g(2) * (t26 * t13 - t23 * t20), 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t27, g(3) * t22 + t12 * t25, 0, t27 * pkin(2), 0, 0, 0, 0, 0, t29 * t24, -t29 * t21, 0, 0, 0, 0, 0, t29 * t17, -t29 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t21 * t39, g(1) * t10 - g(2) * t8 + t24 * t39, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t28;
