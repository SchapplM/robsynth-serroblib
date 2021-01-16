% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:45
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:44:14
% EndTime: 2021-01-15 11:44:15
% DurationCPUTime: 0.09s
% Computational Cost: add. (120->26), mult. (94->37), div. (0->0), fcn. (87->10), ass. (0->22)
t11 = qJ(3) + pkin(9);
t12 = qJ(1) + pkin(8);
t7 = sin(t12);
t9 = cos(t12);
t21 = g(2) * t9 + g(3) * t7;
t20 = g(2) * t7 - g(3) * t9;
t15 = sin(qJ(1));
t17 = cos(qJ(1));
t19 = -g(2) * t17 - g(3) * t15;
t14 = sin(qJ(3));
t16 = cos(qJ(3));
t18 = -g(1) * t16 + t20 * t14;
t13 = -qJ(4) - pkin(6);
t10 = qJ(5) + t11;
t8 = cos(t11);
t6 = sin(t11);
t5 = t16 * pkin(3) + pkin(2);
t4 = cos(t10);
t3 = sin(t10);
t2 = -g(1) * t4 + t20 * t3;
t1 = g(1) * t3 + t20 * t4;
t22 = [0, t19, g(2) * t15 - g(3) * t17, t19 * pkin(1), 0, 0, 0, 0, 0, -t21 * t16, t21 * t14, -t21 * t8, t21 * t6, -t20, -g(2) * (t17 * pkin(1) - t7 * t13 + t9 * t5) - g(3) * (t15 * pkin(1) + t9 * t13 + t7 * t5), 0, 0, 0, 0, 0, -t21 * t4, t21 * t3; 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, g(1) * t14 + t20 * t16, -g(1) * t8 + t20 * t6, g(1) * t6 + t20 * t8, 0, t18 * pkin(3), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t22;
