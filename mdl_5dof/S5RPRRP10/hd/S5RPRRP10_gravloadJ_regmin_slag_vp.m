% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP10
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
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:04
% EndTime: 2019-12-31 18:52:04
% DurationCPUTime: 0.14s
% Computational Cost: add. (118->38), mult. (161->54), div. (0->0), fcn. (159->8), ass. (0->29)
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t9 = g(1) * t23 + g(2) * t21;
t20 = sin(qJ(4));
t15 = pkin(8) + qJ(3);
t12 = sin(t15);
t34 = g(3) * t12;
t13 = cos(t15);
t22 = cos(qJ(4));
t29 = t23 * t22;
t32 = t21 * t20;
t4 = t13 * t32 + t29;
t30 = t23 * t20;
t31 = t21 * t22;
t6 = -t13 * t30 + t31;
t39 = -g(1) * t6 + g(2) * t4 + t20 * t34;
t1 = -g(3) * t13 + t9 * t12;
t27 = pkin(4) * t20 + pkin(6) + qJ(2);
t8 = g(1) * t21 - g(2) * t23;
t11 = t22 * pkin(4) + pkin(3);
t18 = -qJ(5) - pkin(7);
t26 = t13 * t11 - t12 * t18;
t17 = cos(pkin(8));
t24 = t17 * pkin(2) + pkin(1) + t26;
t7 = t13 * t29 + t32;
t5 = -t13 * t31 + t30;
t3 = t8 * t12;
t2 = t9 * t13 + t34;
t10 = [0, t8, t9, t8 * t17, -t8 * sin(pkin(8)), -t9, -g(1) * (-t21 * pkin(1) + t23 * qJ(2)) - g(2) * (t23 * pkin(1) + t21 * qJ(2)), 0, 0, 0, 0, 0, t8 * t13, -t3, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, -g(1) * t4 - g(2) * t6, t3, (-g(1) * t27 - g(2) * t24) * t23 + (g(1) * t24 - g(2) * t27) * t21; 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, t1 * t22, -t1 * t20, -t2, -g(3) * t26 + t9 * (t11 * t12 + t13 * t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, g(1) * t7 - g(2) * t5 + t22 * t34, 0, t39 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t10;
