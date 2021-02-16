% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:35:11
% EndTime: 2021-01-15 21:35:12
% DurationCPUTime: 0.12s
% Computational Cost: add. (143->32), mult. (156->50), div. (0->0), fcn. (159->10), ass. (0->31)
t21 = sin(qJ(1));
t24 = cos(qJ(1));
t10 = g(1) * t24 + g(2) * t21;
t17 = qJ(2) + pkin(9);
t16 = qJ(4) + t17;
t11 = sin(t16);
t12 = cos(t16);
t3 = -g(3) * t12 + t10 * t11;
t31 = g(3) * t11;
t19 = sin(qJ(5));
t29 = t21 * t19;
t22 = cos(qJ(5));
t28 = t21 * t22;
t27 = t24 * t19;
t26 = t24 * t22;
t9 = g(1) * t21 - g(2) * t24;
t20 = sin(qJ(2));
t23 = cos(qJ(2));
t25 = -g(3) * t23 + t10 * t20;
t18 = -qJ(3) - pkin(6);
t15 = cos(t17);
t14 = sin(t17);
t13 = t23 * pkin(2) + pkin(1);
t8 = t12 * t26 + t29;
t7 = -t12 * t27 + t28;
t6 = -t12 * t28 + t27;
t5 = t12 * t29 + t26;
t4 = t10 * t12 + t31;
t2 = t3 * t22;
t1 = t3 * t19;
t30 = [0, t9, t10, 0, 0, 0, 0, 0, t9 * t23, -t9 * t20, t9 * t15, -t9 * t14, -t10, -g(1) * (-t21 * t13 - t18 * t24) - g(2) * (t24 * t13 - t21 * t18), 0, 0, 0, 0, 0, t9 * t12, -t9 * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, t25, g(3) * t20 + t10 * t23, -g(3) * t15 + t10 * t14, g(3) * t14 + t10 * t15, 0, t25 * pkin(2), 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 + t19 * t31, g(1) * t8 - g(2) * t6 + t22 * t31;];
taug_reg = t30;
