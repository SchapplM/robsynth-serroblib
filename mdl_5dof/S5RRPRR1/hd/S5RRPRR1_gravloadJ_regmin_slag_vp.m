% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_gravloadJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:13:18
% EndTime: 2021-01-15 21:13:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (91->27), mult. (156->44), div. (0->0), fcn. (159->8), ass. (0->30)
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t14 = g(1) * t23 + g(2) * t20;
t17 = qJ(2) + qJ(4);
t15 = sin(t17);
t16 = cos(t17);
t3 = -g(3) * t16 + t14 * t15;
t22 = cos(qJ(2));
t30 = pkin(1) * t22;
t29 = g(3) * t15;
t18 = sin(qJ(5));
t27 = t20 * t18;
t21 = cos(qJ(5));
t26 = t20 * t21;
t25 = t23 * t18;
t24 = t23 * t21;
t13 = g(1) * t20 - g(2) * t23;
t19 = sin(qJ(2));
t5 = -g(3) * t22 + t14 * t19;
t12 = t13 * t22;
t11 = t13 * t19;
t10 = t16 * t24 + t27;
t9 = -t16 * t25 + t26;
t8 = -t16 * t26 + t25;
t7 = t16 * t27 + t24;
t6 = g(3) * t19 + t14 * t22;
t4 = t14 * t16 + t29;
t2 = t3 * t21;
t1 = t3 * t18;
t28 = [0, t13, t14, 0, 0, 0, 0, 0, t12, -t11, t12, -t11, -t14, -g(1) * (t23 * qJ(3) - t20 * t30) - g(2) * (t20 * qJ(3) + t23 * t30), 0, 0, 0, 0, 0, t13 * t16, -t13 * t15, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t5, t6, 0, t5 * pkin(1), 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t18 * t29, g(1) * t10 - g(2) * t8 + t21 * t29;];
taug_reg = t28;
