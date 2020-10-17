% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:16
% EndTime: 2019-12-31 18:43:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (116->34), mult. (145->52), div. (0->0), fcn. (143->8), ass. (0->30)
t11 = qJ(1) + pkin(8);
t10 = cos(t11);
t9 = sin(t11);
t24 = g(1) * t10 + g(2) * t9;
t13 = sin(qJ(4));
t16 = cos(qJ(4));
t17 = cos(qJ(3));
t29 = t13 * t17;
t3 = t10 * t16 + t9 * t29;
t14 = sin(qJ(3));
t31 = g(3) * t14;
t5 = -t10 * t29 + t9 * t16;
t36 = -g(1) * t5 + g(2) * t3 + t13 * t31;
t1 = -g(3) * t17 + t24 * t14;
t28 = t16 * t17;
t26 = pkin(4) * t13 + pkin(6);
t25 = g(1) * t9 - g(2) * t10;
t15 = sin(qJ(1));
t18 = cos(qJ(1));
t23 = g(1) * t15 - g(2) * t18;
t12 = -qJ(5) - pkin(7);
t8 = t16 * pkin(4) + pkin(3);
t21 = -t14 * t12 + t17 * t8;
t20 = -pkin(2) - t21;
t19 = t23 * pkin(1);
t7 = t25 * t14;
t6 = t10 * t28 + t9 * t13;
t4 = t10 * t13 - t9 * t28;
t2 = t24 * t17 + t31;
t22 = [0, t23, g(1) * t18 + g(2) * t15, t19, 0, 0, 0, 0, 0, t25 * t17, -t7, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t7, t19 + (-g(1) * t20 - g(2) * t26) * t9 + (-g(1) * t26 + g(2) * t20) * t10; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, t1 * t16, -t1 * t13, -t2, -g(3) * t21 + t24 * (t12 * t17 + t14 * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, g(1) * t6 - g(2) * t4 + t16 * t31, 0, t36 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t22;
