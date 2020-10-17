% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP11
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
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:32
% EndTime: 2019-12-31 18:54:33
% DurationCPUTime: 0.18s
% Computational Cost: add. (174->49), mult. (249->70), div. (0->0), fcn. (265->8), ass. (0->32)
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t12 = g(1) * t26 + g(2) * t24;
t19 = pkin(8) + qJ(3);
t16 = sin(t19);
t36 = g(3) * t16;
t23 = sin(qJ(4));
t35 = t24 * t23;
t25 = cos(qJ(4));
t34 = t24 * t25;
t33 = t26 * t23;
t32 = t26 * t25;
t17 = cos(t19);
t7 = t17 * t35 + t32;
t9 = t17 * t33 - t34;
t31 = g(1) * t7 - g(2) * t9;
t11 = g(1) * t24 - g(2) * t26;
t21 = cos(pkin(8));
t30 = t21 * pkin(2) + t17 * pkin(3) + t16 * pkin(7) + pkin(1);
t29 = pkin(4) * t25 + qJ(5) * t23 + pkin(3);
t1 = g(1) * t9 + g(2) * t7 + t23 * t36;
t10 = t17 * t32 + t35;
t8 = t17 * t34 - t33;
t28 = g(1) * t10 + g(2) * t8 + t25 * t36;
t27 = -g(3) * t17 + t12 * t16;
t22 = -pkin(6) - qJ(2);
t6 = t11 * t16;
t5 = t12 * t17 + t36;
t4 = t27 * t25;
t3 = t27 * t23;
t2 = g(1) * t8 - g(2) * t10;
t13 = [0, t11, t12, t11 * t21, -t11 * sin(pkin(8)), -t12, -g(1) * (-t24 * pkin(1) + t26 * qJ(2)) - g(2) * (t26 * pkin(1) + t24 * qJ(2)), 0, 0, 0, 0, 0, t11 * t17, -t6, 0, 0, 0, 0, 0, t2, -t31, t2, t6, t31, -g(1) * (-t8 * pkin(4) - t7 * qJ(5)) - g(2) * (t10 * pkin(4) + t9 * qJ(5)) + (g(1) * t22 - g(2) * t30) * t26 + (g(1) * t30 + g(2) * t22) * t24; 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t5, 0, 0, 0, 0, 0, t4, -t3, t4, -t5, t3, (-t12 * pkin(7) - g(3) * t29) * t17 + (-g(3) * pkin(7) + t12 * t29) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t28, t1, 0, -t28, -g(1) * (-t9 * pkin(4) + t10 * qJ(5)) - g(2) * (-t7 * pkin(4) + t8 * qJ(5)) - (-pkin(4) * t23 + qJ(5) * t25) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t13;
