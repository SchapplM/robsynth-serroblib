% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:08:52
% EndTime: 2021-01-15 20:08:53
% DurationCPUTime: 0.10s
% Computational Cost: add. (149->31), mult. (110->35), div. (0->0), fcn. (99->8), ass. (0->27)
t21 = qJ(1) + qJ(2);
t16 = pkin(8) + t21;
t11 = sin(t16);
t12 = cos(t16);
t17 = sin(t21);
t13 = pkin(2) * t17;
t25 = cos(qJ(4));
t15 = t25 * pkin(4) + pkin(3);
t22 = -qJ(5) - pkin(7);
t30 = t11 * t15 + t12 * t22 + t13;
t18 = cos(t21);
t14 = pkin(2) * t18;
t29 = -t11 * t22 + t12 * t15 + t14;
t28 = g(2) * t12 + g(3) * t11;
t27 = g(2) * t11 - g(3) * t12;
t7 = -g(2) * t18 - g(3) * t17;
t23 = sin(qJ(4));
t2 = -g(1) * t25 + t27 * t23;
t26 = cos(qJ(1));
t24 = sin(qJ(1));
t20 = t26 * pkin(1);
t19 = t24 * pkin(1);
t6 = g(2) * t17 - g(3) * t18;
t4 = t28 * t25;
t3 = t28 * t23;
t1 = g(1) * t23 + t27 * t25;
t5 = [0, -g(2) * t26 - g(3) * t24, g(2) * t24 - g(3) * t26, 0, t7, t6, -g(2) * (t14 + t20) - g(3) * (t13 + t19), 0, 0, 0, 0, 0, -t4, t3, -t4, t3, -t27, -g(2) * (t20 + t29) - g(3) * (t19 + t30); 0, 0, 0, 0, t7, t6, t7 * pkin(2), 0, 0, 0, 0, 0, -t4, t3, -t4, t3, -t27, -g(2) * t29 - g(3) * t30; 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28;];
taug_reg = t5;
