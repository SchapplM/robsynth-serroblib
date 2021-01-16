% Calculate minimal parameter regressor of gravitation load for
% S5RRRPP1
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
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:14:48
% EndTime: 2021-01-15 22:14:49
% DurationCPUTime: 0.17s
% Computational Cost: add. (223->44), mult. (185->47), div. (0->0), fcn. (171->8), ass. (0->32)
t21 = qJ(3) + pkin(8);
t15 = sin(t21);
t16 = cos(t21);
t43 = t16 * pkin(4) + t15 * qJ(5);
t22 = qJ(1) + qJ(2);
t17 = sin(t22);
t18 = cos(t22);
t7 = g(1) * t17 - g(2) * t18;
t8 = g(1) * t18 + g(2) * t17;
t25 = sin(qJ(1));
t37 = t25 * pkin(1);
t23 = -qJ(4) - pkin(7);
t36 = t18 * t23;
t26 = cos(qJ(3));
t19 = t26 * pkin(3);
t14 = t19 + pkin(2);
t12 = t18 * t14;
t34 = t43 * t18 + t12;
t33 = -t17 * t23 + t12;
t31 = -t17 * t14 - t36;
t24 = sin(qJ(3));
t29 = -g(3) * t26 + t8 * t24;
t28 = (-g(1) * (-t14 - t43) + g(2) * t23) * t17;
t27 = cos(qJ(1));
t20 = t27 * pkin(1);
t6 = t7 * t26;
t5 = t7 * t24;
t4 = t7 * t16;
t3 = t7 * t15;
t2 = g(3) * t15 + t8 * t16;
t1 = -g(3) * t16 + t8 * t15;
t9 = [0, g(1) * t25 - g(2) * t27, g(1) * t27 + g(2) * t25, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, t4, -t3, -t8, -g(1) * (t31 - t37) - g(2) * (t20 + t33), t4, -t8, t3, -g(1) * (-t36 - t37) - g(2) * (t20 + t34) + t28; 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, t4, -t3, -t8, -g(1) * t31 - g(2) * t33, t4, -t8, t3, g(1) * t36 - g(2) * t34 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, g(3) * t24 + t8 * t26, t1, t2, 0, t29 * pkin(3), t1, 0, -t2, -g(3) * (t19 + t43) + t8 * (pkin(3) * t24 + pkin(4) * t15 - qJ(5) * t16); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t9;
