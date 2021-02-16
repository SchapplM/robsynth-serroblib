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
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:11
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-16 00:10:23
% EndTime: 2021-01-16 00:10:24
% DurationCPUTime: 0.20s
% Computational Cost: add. (205->40), mult. (279->55), div. (0->0), fcn. (284->8), ass. (0->37)
t26 = cos(qJ(2));
t25 = cos(qJ(4));
t15 = t25 * pkin(4) + pkin(3);
t20 = qJ(2) + qJ(3);
t17 = sin(t20);
t18 = cos(t20);
t21 = -qJ(5) - pkin(8);
t34 = t18 * t15 - t17 * t21;
t48 = t26 * pkin(2) + t34;
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t33 = g(1) * t27 + g(2) * t24;
t37 = t27 * t25;
t22 = sin(qJ(4));
t40 = t24 * t22;
t10 = t18 * t40 + t37;
t38 = t27 * t22;
t39 = t24 * t25;
t12 = -t18 * t38 + t39;
t43 = g(3) * t17;
t1 = -g(1) * t12 + g(2) * t10 + t22 * t43;
t7 = -g(3) * t18 + t33 * t17;
t35 = pkin(4) * t22 + pkin(6) + pkin(7);
t32 = g(1) * t24 - g(2) * t27;
t31 = t15 * t17 + t18 * t21;
t30 = pkin(1) + t48;
t23 = sin(qJ(2));
t13 = t18 * t37 + t40;
t11 = -t18 * t39 + t38;
t9 = t32 * t17;
t8 = t33 * t18 + t43;
t6 = t7 * t25;
t5 = t7 * t22;
t4 = -g(1) * t11 - g(2) * t13;
t3 = -g(1) * t10 - g(2) * t12;
t2 = g(1) * t13 - g(2) * t11 + t25 * t43;
t14 = [0, t32, t33, 0, 0, 0, 0, 0, t32 * t26, -t32 * t23, 0, 0, 0, 0, 0, t32 * t18, -t9, 0, 0, 0, 0, 0, t4, t3, t4, t3, t9, (-g(1) * t35 - g(2) * t30) * t27 + (g(1) * t30 - g(2) * t35) * t24; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t26 + t33 * t23, g(3) * t23 + t33 * t26, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, t6, -t5, -t8, -g(3) * t48 + t33 * (pkin(2) * t23 + t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, t6, -t5, -t8, -g(3) * t34 + t33 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t14;
