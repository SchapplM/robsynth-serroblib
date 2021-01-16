% Calculate minimal parameter regressor of gravitation load for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:24:48
% EndTime: 2021-01-15 22:24:49
% DurationCPUTime: 0.17s
% Computational Cost: add. (230->46), mult. (205->58), div. (0->0), fcn. (189->8), ass. (0->35)
t25 = qJ(2) + qJ(3);
t19 = pkin(8) + t25;
t16 = sin(t19);
t17 = cos(t19);
t30 = t17 * pkin(4) + t16 * qJ(5);
t20 = sin(t25);
t37 = pkin(3) * t20;
t36 = pkin(4) * t16;
t21 = cos(t25);
t18 = pkin(3) * t21;
t28 = cos(qJ(2));
t22 = t28 * pkin(2);
t35 = t18 + t22;
t34 = qJ(5) * t17;
t33 = t18 + t30;
t26 = sin(qJ(2));
t9 = -t26 * pkin(2) - t37;
t32 = t9 - t36;
t31 = -t36 - t37;
t27 = sin(qJ(1));
t29 = cos(qJ(1));
t13 = g(1) * t29 + g(2) * t27;
t12 = g(1) * t27 - g(2) * t29;
t3 = -g(3) * t21 + t13 * t20;
t24 = -qJ(4) - pkin(7) - pkin(6);
t11 = t29 * t34;
t10 = t27 * t34;
t8 = pkin(1) + t35;
t7 = t29 * t8;
t6 = t12 * t17;
t5 = t12 * t16;
t4 = g(3) * t20 + t13 * t21;
t2 = g(3) * t16 + t13 * t17;
t1 = -g(3) * t17 + t13 * t16;
t14 = [0, t12, t13, 0, 0, 0, 0, 0, t12 * t28, -t12 * t26, 0, 0, 0, 0, 0, t12 * t21, -t12 * t20, t6, -t5, -t13, -g(1) * (-t29 * t24 - t27 * t8) - g(2) * (-t27 * t24 + t7), t6, -t13, t5, -g(2) * t7 + (g(1) * t24 - g(2) * t30) * t29 + (-g(1) * (-t30 - t8) + g(2) * t24) * t27; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t28 + t13 * t26, g(3) * t26 + t13 * t28, 0, 0, 0, 0, 0, t3, t4, t1, t2, 0, -g(3) * t35 - t13 * t9, t1, 0, -t2, -g(1) * (t32 * t29 + t11) - g(2) * (t32 * t27 + t10) - g(3) * (t22 + t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, t1, t2, 0, t3 * pkin(3), t1, 0, -t2, -g(1) * (t31 * t29 + t11) - g(2) * (t31 * t27 + t10) - g(3) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t14;
