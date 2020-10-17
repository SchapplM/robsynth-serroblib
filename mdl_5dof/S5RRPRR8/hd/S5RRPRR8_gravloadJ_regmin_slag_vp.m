% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR8
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:13
% EndTime: 2019-12-31 20:18:13
% DurationCPUTime: 0.11s
% Computational Cost: add. (127->28), mult. (138->44), div. (0->0), fcn. (141->8), ass. (0->28)
t18 = sin(qJ(1));
t21 = cos(qJ(1));
t10 = g(1) * t21 + g(2) * t18;
t14 = qJ(2) + pkin(9) + qJ(4);
t11 = sin(t14);
t12 = cos(t14);
t3 = -g(3) * t12 + t10 * t11;
t28 = g(3) * t11;
t16 = sin(qJ(5));
t26 = t18 * t16;
t19 = cos(qJ(5));
t25 = t18 * t19;
t24 = t21 * t16;
t23 = t21 * t19;
t9 = g(1) * t18 - g(2) * t21;
t17 = sin(qJ(2));
t20 = cos(qJ(2));
t22 = -g(3) * t20 + t10 * t17;
t15 = -qJ(3) - pkin(6);
t13 = t20 * pkin(2) + pkin(1);
t8 = t12 * t23 + t26;
t7 = -t12 * t24 + t25;
t6 = -t12 * t25 + t24;
t5 = t12 * t26 + t23;
t4 = t10 * t12 + t28;
t2 = t3 * t19;
t1 = t3 * t16;
t27 = [0, t9, t10, 0, 0, 0, 0, 0, t9 * t20, -t9 * t17, -t10, -g(1) * (-t18 * t13 - t21 * t15) - g(2) * (t21 * t13 - t18 * t15), 0, 0, 0, 0, 0, t9 * t12, -t9 * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, t22, g(3) * t17 + t10 * t20, 0, t22 * pkin(2), 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 + t16 * t28, g(1) * t8 - g(2) * t6 + t19 * t28;];
taug_reg = t27;
