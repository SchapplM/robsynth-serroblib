% Calculate minimal parameter regressor of gravitation load for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:14:03
% EndTime: 2021-01-15 11:14:04
% DurationCPUTime: 0.12s
% Computational Cost: add. (149->38), mult. (129->44), div. (0->0), fcn. (117->8), ass. (0->25)
t15 = qJ(3) + pkin(8);
t11 = cos(t15);
t9 = sin(t15);
t24 = t11 * pkin(4) + t9 * qJ(5);
t16 = qJ(1) + pkin(7);
t10 = sin(t16);
t12 = cos(t16);
t28 = g(1) * t12 + g(2) * t10;
t21 = cos(qJ(1));
t20 = cos(qJ(3));
t13 = t20 * pkin(3);
t8 = t13 + pkin(2);
t30 = t21 * pkin(1) + t12 * t8;
t27 = g(1) * t10 - g(2) * t12;
t19 = sin(qJ(1));
t26 = g(1) * t19 - g(2) * t21;
t17 = -qJ(4) - pkin(6);
t25 = -t19 * pkin(1) - t12 * t17;
t18 = sin(qJ(3));
t22 = -g(3) * t20 + t28 * t18;
t4 = t27 * t11;
t3 = t27 * t9;
t2 = g(3) * t9 + t28 * t11;
t1 = -g(3) * t11 + t28 * t9;
t5 = [0, t26, g(1) * t21 + g(2) * t19, t26 * pkin(1), 0, 0, 0, 0, 0, t27 * t20, -t27 * t18, t4, -t3, -t28, -g(1) * (-t10 * t8 + t25) - g(2) * (-t10 * t17 + t30), t4, -t28, t3, -g(1) * t25 - g(2) * (t24 * t12 + t30) + (-g(1) * (-t24 - t8) + g(2) * t17) * t10; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, g(3) * t18 + t28 * t20, t1, t2, 0, t22 * pkin(3), t1, 0, -t2, -g(3) * (t13 + t24) + t28 * (pkin(3) * t18 + pkin(4) * t9 - qJ(5) * t11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, 0, 0, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t5;
