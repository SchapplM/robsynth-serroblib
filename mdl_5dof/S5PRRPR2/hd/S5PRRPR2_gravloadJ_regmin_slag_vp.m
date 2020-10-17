% Calculate minimal parameter regressor of gravitation load for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:35
% EndTime: 2019-12-05 16:17:35
% DurationCPUTime: 0.09s
% Computational Cost: add. (157->29), mult. (100->41), div. (0->0), fcn. (108->8), ass. (0->26)
t20 = sin(pkin(9));
t28 = g(3) * t20;
t21 = cos(pkin(9));
t22 = sin(qJ(5));
t27 = t21 * t22;
t23 = cos(qJ(5));
t26 = t21 * t23;
t19 = pkin(8) + qJ(2);
t18 = qJ(3) + t19;
t14 = sin(t18);
t15 = cos(t18);
t25 = t15 * pkin(3) + t14 * qJ(4);
t24 = -t14 * pkin(3) + t15 * qJ(4);
t9 = g(1) * t14 - g(2) * t15;
t17 = cos(t19);
t16 = sin(t19);
t10 = g(1) * t15 + g(2) * t14;
t8 = t9 * t21;
t7 = t9 * t20;
t6 = t14 * t22 + t15 * t26;
t5 = t14 * t23 - t15 * t27;
t4 = -t14 * t26 + t15 * t22;
t3 = t14 * t27 + t15 * t23;
t2 = -g(1) * t4 - g(2) * t6;
t1 = -g(1) * t3 - g(2) * t5;
t11 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t16 - g(2) * t17, g(1) * t17 + g(2) * t16, 0, t9, t10, t8, -t7, -t10, -g(1) * (-pkin(2) * t16 + t24) - g(2) * (pkin(2) * t17 + t25), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, t9, t10, t8, -t7, -t10, -g(1) * t24 - g(2) * t25, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t22 * t28, g(1) * t6 - g(2) * t4 + t23 * t28;];
taug_reg = t11;
