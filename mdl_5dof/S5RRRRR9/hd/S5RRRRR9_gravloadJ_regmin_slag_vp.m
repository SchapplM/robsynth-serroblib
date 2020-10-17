% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:52
% EndTime: 2019-12-31 22:29:53
% DurationCPUTime: 0.16s
% Computational Cost: add. (208->40), mult. (238->71), div. (0->0), fcn. (274->10), ass. (0->40)
t24 = sin(qJ(2));
t27 = cos(qJ(2));
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t32 = g(1) * t28 + g(2) * t25;
t30 = -g(3) * t27 + t32 * t24;
t41 = g(3) * t24;
t39 = t25 * t27;
t22 = qJ(3) + qJ(4);
t21 = qJ(5) + t22;
t17 = sin(t21);
t38 = t28 * t17;
t18 = cos(t21);
t37 = t28 * t18;
t19 = sin(t22);
t36 = t28 * t19;
t20 = cos(t22);
t35 = t28 * t20;
t23 = sin(qJ(3));
t34 = t28 * t23;
t26 = cos(qJ(3));
t33 = t28 * t26;
t31 = g(1) * t25 - g(2) * t28;
t16 = t25 * t23 + t27 * t33;
t15 = t25 * t26 - t27 * t34;
t14 = -t26 * t39 + t34;
t13 = t23 * t39 + t33;
t12 = t25 * t19 + t27 * t35;
t11 = t25 * t20 - t27 * t36;
t10 = -t20 * t39 + t36;
t9 = t19 * t39 + t35;
t8 = t25 * t17 + t27 * t37;
t7 = t25 * t18 - t27 * t38;
t6 = -t18 * t39 + t38;
t5 = t17 * t39 + t37;
t4 = g(1) * t12 - g(2) * t10 + t20 * t41;
t3 = -g(1) * t11 + g(2) * t9 + t19 * t41;
t2 = g(1) * t8 - g(2) * t6 + t18 * t41;
t1 = -g(1) * t7 + g(2) * t5 + t17 * t41;
t29 = [0, t31, t32, 0, 0, 0, 0, 0, t31 * t27, -t31 * t24, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, t30, t32 * t27 + t41, 0, 0, 0, 0, 0, t30 * t26, -t30 * t23, 0, 0, 0, 0, 0, t30 * t20, -t30 * t19, 0, 0, 0, 0, 0, t30 * t18, -t30 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t23 * t41, g(1) * t16 - g(2) * t14 + t26 * t41, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t29;
