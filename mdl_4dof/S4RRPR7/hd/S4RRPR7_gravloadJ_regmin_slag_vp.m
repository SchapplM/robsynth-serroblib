% Calculate minimal parameter regressor of gravitation load for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:44
% EndTime: 2019-12-31 17:06:44
% DurationCPUTime: 0.11s
% Computational Cost: add. (53->25), mult. (94->41), div. (0->0), fcn. (97->8), ass. (0->25)
t14 = sin(qJ(1));
t17 = cos(qJ(1));
t6 = g(1) * t17 + g(2) * t14;
t10 = qJ(2) + pkin(7);
t8 = sin(t10);
t9 = cos(t10);
t25 = -g(3) * t9 + t6 * t8;
t24 = g(3) * t8;
t12 = sin(qJ(4));
t22 = t14 * t12;
t15 = cos(qJ(4));
t21 = t14 * t15;
t20 = t17 * t12;
t19 = t17 * t15;
t5 = g(1) * t14 - g(2) * t17;
t13 = sin(qJ(2));
t16 = cos(qJ(2));
t18 = -g(3) * t16 + t6 * t13;
t11 = -qJ(3) - pkin(5);
t7 = t16 * pkin(2) + pkin(1);
t4 = t9 * t19 + t22;
t3 = -t9 * t20 + t21;
t2 = -t9 * t21 + t20;
t1 = t9 * t22 + t19;
t23 = [0, t5, t6, 0, 0, 0, 0, 0, t5 * t16, -t5 * t13, -t6, -g(1) * (-t17 * t11 - t14 * t7) - g(2) * (-t14 * t11 + t17 * t7), 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t18, g(3) * t13 + t6 * t16, 0, t18 * pkin(2), 0, 0, 0, 0, 0, t25 * t15, -t25 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t12 * t24, g(1) * t4 - g(2) * t2 + t15 * t24;];
taug_reg = t23;
