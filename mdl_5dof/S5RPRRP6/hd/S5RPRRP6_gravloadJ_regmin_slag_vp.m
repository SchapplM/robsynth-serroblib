% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 18:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 18:08:27
% EndTime: 2021-01-15 18:08:28
% DurationCPUTime: 0.14s
% Computational Cost: add. (154->35), mult. (197->52), div. (0->0), fcn. (203->8), ass. (0->35)
t17 = qJ(1) + pkin(8);
t15 = sin(t17);
t16 = cos(t17);
t31 = g(1) * t16 + g(2) * t15;
t22 = cos(qJ(4));
t19 = sin(qJ(4));
t23 = cos(qJ(3));
t35 = t19 * t23;
t11 = t15 * t22 - t16 * t35;
t20 = sin(qJ(3));
t37 = g(3) * t20;
t9 = t15 * t35 + t16 * t22;
t1 = -g(1) * t11 + g(2) * t9 + t19 * t37;
t7 = -g(3) * t23 + t31 * t20;
t34 = t22 * t23;
t32 = pkin(4) * t19 + pkin(6);
t30 = g(1) * t15 - g(2) * t16;
t21 = sin(qJ(1));
t24 = cos(qJ(1));
t29 = g(1) * t21 - g(2) * t24;
t14 = t22 * pkin(4) + pkin(3);
t18 = -qJ(5) - pkin(7);
t28 = t23 * t14 - t20 * t18;
t26 = pkin(2) + t28;
t25 = t29 * pkin(1);
t13 = t30 * t20;
t12 = t15 * t19 + t16 * t34;
t10 = -t15 * t34 + t16 * t19;
t8 = t31 * t23 + t37;
t6 = t7 * t22;
t5 = t7 * t19;
t4 = -g(1) * t10 - g(2) * t12;
t3 = -g(1) * t9 - g(2) * t11;
t2 = g(1) * t12 - g(2) * t10 + t22 * t37;
t27 = [0, t29, g(1) * t24 + g(2) * t21, t25, 0, 0, 0, 0, 0, t30 * t23, -t13, 0, 0, 0, 0, 0, t4, t3, t4, t3, t13, t25 + (-g(1) * t32 - g(2) * t26) * t16 + (g(1) * t26 - g(2) * t32) * t15; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, t6, -t5, -t8, -g(3) * t28 + t31 * (t14 * t20 + t18 * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t27;
