% Calculate minimal parameter regressor of gravitation load for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:01
% EndTime: 2019-12-31 17:44:01
% DurationCPUTime: 0.12s
% Computational Cost: add. (92->30), mult. (107->45), div. (0->0), fcn. (111->8), ass. (0->25)
t15 = qJ(1) + pkin(7);
t12 = sin(t15);
t29 = g(1) * t12;
t13 = cos(t15);
t21 = cos(qJ(1));
t28 = t21 * pkin(1) + t13 * pkin(2) + t12 * qJ(3);
t19 = sin(qJ(1));
t27 = -t19 * pkin(1) + t13 * qJ(3);
t8 = -g(1) * t13 - g(2) * t12;
t26 = -g(2) * t13 + t29;
t25 = g(1) * t19 - g(2) * t21;
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t24 = pkin(3) * t17 + qJ(4) * t16;
t18 = sin(qJ(5));
t20 = cos(qJ(5));
t23 = t16 * t20 - t17 * t18;
t22 = t16 * t18 + t17 * t20;
t6 = t26 * t17;
t5 = t26 * t16;
t4 = t22 * t13;
t3 = t23 * t13;
t2 = t22 * t12;
t1 = t23 * t12;
t7 = [0, t25, g(1) * t21 + g(2) * t19, t25 * pkin(1), t6, -t5, t8, -g(1) * (-t12 * pkin(2) + t27) - g(2) * t28, t6, t8, t5, -g(1) * t27 - g(2) * (t24 * t13 + t28) - (-pkin(2) - t24) * t29, 0, 0, 0, 0, 0, g(1) * t2 - g(2) * t4, g(1) * t1 - g(2) * t3; 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t17 + t8 * t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t1 + g(3) * t22, g(1) * t4 + g(2) * t2 + g(3) * t23;];
taug_reg = t7;
