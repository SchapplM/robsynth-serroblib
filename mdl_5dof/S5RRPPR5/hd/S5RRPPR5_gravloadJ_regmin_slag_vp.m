% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:36:30
% EndTime: 2021-01-15 19:36:31
% DurationCPUTime: 0.13s
% Computational Cost: add. (145->40), mult. (189->53), div. (0->0), fcn. (197->10), ass. (0->32)
t28 = sin(qJ(1));
t30 = cos(qJ(1));
t13 = g(1) * t30 + g(2) * t28;
t23 = qJ(2) + pkin(8);
t19 = sin(t23);
t20 = cos(t23);
t33 = sin(qJ(5));
t34 = cos(qJ(5));
t7 = -t19 * t33 - t20 * t34;
t8 = t19 * t34 - t20 * t33;
t39 = -g(3) * t8 + t13 * t7;
t38 = g(3) * t7 + t13 * t8;
t35 = g(3) * t20;
t26 = -qJ(3) - pkin(6);
t32 = t28 * t26;
t12 = g(1) * t28 - g(2) * t30;
t27 = sin(qJ(2));
t29 = cos(qJ(2));
t31 = -g(3) * t29 + t13 * t27;
t25 = cos(pkin(8));
t24 = sin(pkin(8));
t21 = t29 * pkin(2);
t18 = t21 + pkin(1);
t17 = t26 * t30;
t11 = -t24 * pkin(3) + qJ(4) * t25;
t10 = pkin(3) * t25 + qJ(4) * t24 + pkin(2);
t6 = t12 * t20;
t5 = t12 * t19;
t4 = g(3) * t19 + t13 * t20;
t3 = t13 * t19 - t35;
t1 = t10 * t29 + t11 * t27 + pkin(1);
t2 = [0, t12, t13, 0, 0, 0, 0, 0, t12 * t29, -t12 * t27, t6, -t5, -t13, -g(1) * (-t28 * t18 - t17) - g(2) * (t30 * t18 - t32), t6, -t13, t5, -g(1) * (-t1 * t28 - t17) - g(2) * (t1 * t30 - t32), 0, 0, 0, 0, 0, -t12 * t7, t12 * t8; 0, 0, 0, 0, 0, 0, 0, 0, t31, g(3) * t27 + t13 * t29, t3, t4, 0, t31 * pkin(2), t3, 0, -t4, -g(3) * (t20 * pkin(3) + t19 * qJ(4) + t21) - t13 * (-t10 * t27 + t11 * t29), 0, 0, 0, 0, 0, t38, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 - t13 * (t24 * t29 + t25 * t27), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t39;];
taug_reg = t2;
