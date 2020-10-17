% Calculate minimal parameter regressor of gravitation load for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:04
% EndTime: 2019-12-31 17:21:05
% DurationCPUTime: 0.14s
% Computational Cost: add. (86->39), mult. (227->62), div. (0->0), fcn. (245->6), ass. (0->28)
t15 = sin(qJ(1));
t18 = cos(qJ(1));
t24 = g(1) * t18 + g(2) * t15;
t14 = sin(qJ(2));
t29 = g(3) * t14;
t17 = cos(qJ(2));
t28 = t15 * t17;
t13 = sin(qJ(3));
t27 = t18 * t13;
t16 = cos(qJ(3));
t26 = t18 * t16;
t6 = t13 * t28 + t26;
t8 = -t15 * t16 + t17 * t27;
t25 = g(1) * t6 - g(2) * t8;
t23 = g(1) * t15 - g(2) * t18;
t22 = t17 * pkin(2) + t14 * pkin(6) + pkin(1);
t21 = pkin(3) * t16 + qJ(4) * t13 + pkin(2);
t1 = g(1) * t8 + g(2) * t6 + t13 * t29;
t7 = t16 * t28 - t27;
t9 = t15 * t13 + t17 * t26;
t20 = g(1) * t9 + g(2) * t7 + t16 * t29;
t19 = -g(3) * t17 + t24 * t14;
t10 = t23 * t14;
t5 = t24 * t17 + t29;
t4 = t19 * t16;
t3 = t19 * t13;
t2 = g(1) * t7 - g(2) * t9;
t11 = [0, t23, t24, 0, 0, 0, 0, 0, t23 * t17, -t10, 0, 0, 0, 0, 0, t2, -t25, t2, t10, t25, -g(1) * (-t7 * pkin(3) - t6 * qJ(4)) - g(2) * (t9 * pkin(3) + t8 * qJ(4)) + (-g(1) * pkin(5) - g(2) * t22) * t18 + (-g(2) * pkin(5) + g(1) * t22) * t15; 0, 0, 0, 0, 0, 0, 0, 0, t19, t5, 0, 0, 0, 0, 0, t4, -t3, t4, -t5, t3, (-t24 * pkin(6) - g(3) * t21) * t17 + (-g(3) * pkin(6) + t24 * t21) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t20, t1, 0, -t20, -g(1) * (-t8 * pkin(3) + t9 * qJ(4)) - g(2) * (-t6 * pkin(3) + t7 * qJ(4)) - (-pkin(3) * t13 + qJ(4) * t16) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t11;
