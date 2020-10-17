% Calculate minimal parameter regressor of gravitation load for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:06
% EndTime: 2019-12-05 15:36:07
% DurationCPUTime: 0.13s
% Computational Cost: add. (112->33), mult. (166->48), div. (0->0), fcn. (171->8), ass. (0->28)
t14 = sin(pkin(7));
t15 = cos(pkin(7));
t25 = g(1) * t15 + g(2) * t14;
t13 = qJ(2) + pkin(8);
t11 = sin(t13);
t30 = g(3) * t11;
t16 = sin(qJ(4));
t29 = t14 * t16;
t18 = cos(qJ(4));
t28 = t14 * t18;
t27 = t15 * t16;
t26 = t15 * t18;
t24 = pkin(4) * t18 + qJ(5) * t16 + pkin(3);
t12 = cos(t13);
t4 = t12 * t29 + t26;
t6 = t12 * t27 - t28;
t1 = g(1) * t6 + g(2) * t4 + t16 * t30;
t5 = t12 * t28 - t27;
t7 = t12 * t26 + t29;
t23 = g(1) * t7 + g(2) * t5 + t18 * t30;
t22 = -g(3) * t12 + t25 * t11;
t17 = sin(qJ(2));
t19 = cos(qJ(2));
t21 = -g(3) * t19 + t25 * t17;
t8 = -g(1) * t14 + g(2) * t15;
t3 = t22 * t18;
t2 = t22 * t16;
t9 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t21, g(3) * t17 + t25 * t19, t21 * pkin(2), 0, 0, 0, 0, 0, t3, -t2, t3, -t25 * t12 - t30, t2, -g(3) * (t19 * pkin(2) + t11 * pkin(6) + t24 * t12) + t25 * (pkin(2) * t17 - pkin(6) * t12 + t24 * t11); 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t23, t1, 0, -t23, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (-t4 * pkin(4) + t5 * qJ(5)) - (-pkin(4) * t16 + qJ(5) * t18) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t9;
