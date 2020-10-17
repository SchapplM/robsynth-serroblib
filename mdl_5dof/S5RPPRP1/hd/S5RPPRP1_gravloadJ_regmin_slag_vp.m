% Calculate minimal parameter regressor of gravitation load for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:48
% EndTime: 2020-01-03 11:25:49
% DurationCPUTime: 0.15s
% Computational Cost: add. (100->35), mult. (112->54), div. (0->0), fcn. (111->8), ass. (0->26)
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t21 = cos(qJ(4));
t37 = -t16 * (-qJ(5) - pkin(6)) + (t21 * pkin(4) + pkin(3)) * t17;
t15 = qJ(1) + pkin(7);
t11 = sin(t15);
t12 = cos(t15);
t19 = sin(qJ(4));
t28 = t17 * t19;
t1 = -t11 * t28 - t12 * t21;
t3 = -t11 * t21 + t12 * t28;
t33 = g(1) * t16;
t36 = -g(2) * t1 - g(3) * t3 + t19 * t33;
t20 = sin(qJ(1));
t32 = t20 * pkin(1) + t11 * pkin(2);
t30 = t11 * t19;
t27 = t17 * t21;
t22 = cos(qJ(1));
t26 = t22 * pkin(1) + t12 * pkin(2) + t11 * qJ(3);
t6 = g(2) * t12 + g(3) * t11;
t24 = -g(2) * t11 + g(3) * t12;
t23 = -g(2) * t22 - g(3) * t20;
t5 = t6 * t16;
t4 = t12 * t27 + t30;
t2 = t11 * t27 - t12 * t19;
t7 = [0, t23, g(2) * t20 - g(3) * t22, t23 * pkin(1), -t6 * t17, t5, t24, -g(2) * t26 - g(3) * (-t12 * qJ(3) + t32), 0, 0, 0, 0, 0, -g(2) * t4 - g(3) * t2, g(2) * t3 - g(3) * t1, -t5, -g(2) * (pkin(4) * t30 + t26) - g(3) * (t37 * t11 + t32) + (-g(2) * t37 - g(3) * (-pkin(4) * t19 - qJ(3))) * t12; 0, 0, 0, -g(1), 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, g(2) * t2 - g(3) * t4 + t21 * t33, 0, t36 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t17 + t24 * t16;];
taug_reg = t7;
