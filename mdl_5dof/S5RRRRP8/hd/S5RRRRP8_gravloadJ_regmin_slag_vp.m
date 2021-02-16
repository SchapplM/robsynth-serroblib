% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP8
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
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:21:20
% EndTime: 2021-01-16 00:21:21
% DurationCPUTime: 0.22s
% Computational Cost: add. (217->47), mult. (310->76), div. (0->0), fcn. (331->8), ass. (0->42)
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t38 = g(1) * t32 + g(2) * t29;
t26 = qJ(3) + qJ(4);
t22 = sin(t26);
t28 = sin(qJ(2));
t47 = g(3) * t28;
t23 = cos(t26);
t42 = t32 * t23;
t31 = cos(qJ(2));
t44 = t29 * t31;
t7 = t22 * t44 + t42;
t43 = t32 * t22;
t9 = t29 * t23 - t31 * t43;
t1 = -g(1) * t9 + g(2) * t7 + t22 * t47;
t11 = -g(3) * t31 + t38 * t28;
t27 = sin(qJ(3));
t19 = t27 * pkin(3) + pkin(4) * t22;
t45 = pkin(6) + t19;
t41 = t32 * t27;
t30 = cos(qJ(3));
t40 = t32 * t30;
t20 = t30 * pkin(3) + pkin(4) * t23;
t37 = g(1) * t29 - g(2) * t32;
t18 = pkin(2) + t20;
t25 = -qJ(5) - pkin(8) - pkin(7);
t36 = t31 * t18 - t28 * t25;
t34 = pkin(1) + t36;
t17 = t37 * t28;
t16 = t29 * t27 + t31 * t40;
t15 = t29 * t30 - t31 * t41;
t14 = -t30 * t44 + t41;
t13 = t27 * t44 + t40;
t12 = t38 * t31 + t47;
t10 = t29 * t22 + t31 * t42;
t8 = -t23 * t44 + t43;
t6 = t11 * t23;
t5 = t11 * t22;
t4 = -g(1) * t8 - g(2) * t10;
t3 = -g(1) * t7 - g(2) * t9;
t2 = g(1) * t10 - g(2) * t8 + t23 * t47;
t21 = [0, t37, t38, 0, 0, 0, 0, 0, t37 * t31, -t17, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, t4, t3, t4, t3, t17, (-g(1) * t45 - g(2) * t34) * t32 + (g(1) * t34 - g(2) * t45) * t29; 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, t11 * t30, -t11 * t27, 0, 0, 0, 0, 0, t6, -t5, t6, -t5, -t12, -g(3) * t36 + t38 * (t18 * t28 + t25 * t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t27 * t47, g(1) * t16 - g(2) * t14 + t30 * t47, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, -g(1) * (-t32 * t31 * t19 + t29 * t20) - g(2) * (-t19 * t44 - t32 * t20) + t19 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11;];
taug_reg = t21;
