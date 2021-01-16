% Calculate minimal parameter regressor of gravitation load for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:15
% EndTime: 2021-01-15 14:39:15
% DurationCPUTime: 0.11s
% Computational Cost: add. (77->29), mult. (183->46), div. (0->0), fcn. (191->6), ass. (0->31)
t20 = sin(qJ(3));
t23 = cos(qJ(3));
t25 = cos(qJ(1));
t30 = t25 * t23;
t22 = sin(qJ(1));
t24 = cos(qJ(2));
t32 = t22 * t24;
t11 = t20 * t32 + t30;
t31 = t25 * t20;
t13 = t22 * t23 - t24 * t31;
t21 = sin(qJ(2));
t34 = g(3) * t21;
t1 = -g(1) * t13 + g(2) * t11 + t20 * t34;
t27 = g(1) * t25 + g(2) * t22;
t7 = -g(3) * t24 + t27 * t21;
t18 = t23 * pkin(3) + pkin(2);
t19 = qJ(4) + pkin(6);
t28 = t18 * t24 + t19 * t21;
t26 = g(1) * t22 - g(2) * t25;
t17 = t20 * pkin(3) + pkin(5);
t15 = t26 * t21;
t14 = t22 * t20 + t24 * t30;
t12 = -t23 * t32 + t31;
t9 = pkin(1) + t28;
t8 = t27 * t24 + t34;
t6 = t7 * t23;
t5 = t7 * t20;
t4 = -g(1) * t12 - g(2) * t14;
t3 = -g(1) * t11 - g(2) * t13;
t2 = g(1) * t14 - g(2) * t12 + t23 * t34;
t10 = [0, t26, t27, 0, 0, 0, 0, 0, t26 * t24, -t15, 0, 0, 0, 0, 0, t4, t3, t4, t3, t15, -g(1) * (t17 * t25 - t9 * t22) - g(2) * (t17 * t22 + t9 * t25); 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, t6, -t5, -t8, -g(3) * t28 - t27 * (-t21 * t18 + t19 * t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t10;
