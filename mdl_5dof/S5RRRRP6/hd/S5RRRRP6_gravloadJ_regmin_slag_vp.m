% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:34
% EndTime: 2019-12-31 21:54:34
% DurationCPUTime: 0.18s
% Computational Cost: add. (161->38), mult. (211->55), div. (0->0), fcn. (208->8), ass. (0->34)
t22 = cos(qJ(2));
t21 = cos(qJ(4));
t11 = t21 * pkin(4) + pkin(3);
t16 = qJ(2) + qJ(3);
t13 = sin(t16);
t14 = cos(t16);
t17 = -qJ(5) - pkin(8);
t30 = t14 * t11 - t13 * t17;
t45 = t22 * pkin(2) + t30;
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t29 = g(1) * t23 + g(2) * t20;
t18 = sin(qJ(4));
t39 = g(3) * t13;
t33 = t23 * t21;
t36 = t20 * t18;
t6 = t14 * t36 + t33;
t34 = t23 * t18;
t35 = t20 * t21;
t8 = -t14 * t34 + t35;
t44 = -g(1) * t8 + g(2) * t6 + t18 * t39;
t3 = -g(3) * t14 + t29 * t13;
t31 = pkin(4) * t18 + pkin(6) + pkin(7);
t28 = g(1) * t20 - g(2) * t23;
t27 = t11 * t13 + t14 * t17;
t26 = pkin(1) + t45;
t19 = sin(qJ(2));
t9 = t14 * t33 + t36;
t7 = -t14 * t35 + t34;
t5 = t28 * t13;
t4 = t29 * t14 + t39;
t2 = t3 * t21;
t1 = t3 * t18;
t10 = [0, t28, t29, 0, 0, 0, 0, 0, t28 * t22, -t28 * t19, 0, 0, 0, 0, 0, t28 * t14, -t5, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t5, (-g(1) * t31 - g(2) * t26) * t23 + (g(1) * t26 - g(2) * t31) * t20; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t22 + t29 * t19, g(3) * t19 + t29 * t22, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(3) * t45 + t29 * (pkin(2) * t19 + t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(3) * t30 + t27 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, g(1) * t9 - g(2) * t7 + t21 * t39, 0, t44 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3;];
taug_reg = t10;
