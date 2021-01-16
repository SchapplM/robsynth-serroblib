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
% taug_reg [4x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
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
% StartTime: 2021-01-15 10:56:38
% EndTime: 2021-01-15 10:56:39
% DurationCPUTime: 0.10s
% Computational Cost: add. (69->27), mult. (112->44), div. (0->0), fcn. (115->8), ass. (0->25)
t14 = sin(qJ(1));
t17 = cos(qJ(1));
t6 = g(1) * t17 + g(2) * t14;
t10 = qJ(2) + pkin(7);
t8 = sin(t10);
t9 = cos(t10);
t19 = -g(3) * t9 + t6 * t8;
t25 = g(3) * t8;
t12 = sin(qJ(4));
t23 = t14 * t12;
t15 = cos(qJ(4));
t22 = t14 * t15;
t21 = t17 * t12;
t20 = t17 * t15;
t5 = g(1) * t14 - g(2) * t17;
t13 = sin(qJ(2));
t16 = cos(qJ(2));
t18 = -g(3) * t16 + t6 * t13;
t11 = -qJ(3) - pkin(5);
t7 = t16 * pkin(2) + pkin(1);
t4 = t9 * t20 + t23;
t3 = -t9 * t21 + t22;
t2 = -t9 * t22 + t21;
t1 = t9 * t23 + t20;
t24 = [0, t5, t6, 0, 0, 0, 0, 0, t5 * t16, -t5 * t13, t5 * t9, -t5 * t8, -t6, -g(1) * (-t11 * t17 - t14 * t7) - g(2) * (-t14 * t11 + t17 * t7), 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t18, g(3) * t13 + t6 * t16, t19, t6 * t9 + t25, 0, t18 * pkin(2), 0, 0, 0, 0, 0, t19 * t15, -t19 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t12 * t25, g(1) * t4 - g(2) * t2 + t15 * t25;];
taug_reg = t24;
