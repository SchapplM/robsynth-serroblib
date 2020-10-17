% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:11
% EndTime: 2019-12-31 18:22:12
% DurationCPUTime: 0.17s
% Computational Cost: add. (140->39), mult. (153->64), div. (0->0), fcn. (158->10), ass. (0->32)
t13 = qJ(1) + pkin(8);
t11 = cos(t13);
t9 = sin(t13);
t26 = g(1) * t11 + g(2) * t9;
t16 = sin(qJ(3));
t18 = cos(qJ(3));
t5 = -g(3) * t18 + t16 * t26;
t12 = pkin(9) + qJ(5);
t8 = sin(t12);
t35 = t9 * t8;
t33 = g(3) * t16;
t10 = cos(t12);
t31 = t9 * t10;
t30 = t11 * t18;
t14 = sin(pkin(9));
t29 = t14 * t18;
t15 = cos(pkin(9));
t28 = t15 * t18;
t27 = g(1) * t9 - g(2) * t11;
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t25 = g(1) * t17 - g(2) * t19;
t24 = t18 * pkin(3) + t16 * qJ(4);
t22 = pkin(2) + t24;
t21 = t25 * pkin(1);
t7 = t27 * t16;
t6 = t18 * t26 + t33;
t4 = t10 * t30 + t35;
t3 = -t30 * t8 + t31;
t2 = t11 * t8 - t18 * t31;
t1 = t11 * t10 + t18 * t35;
t20 = [0, t25, g(1) * t19 + g(2) * t17, t21, 0, 0, 0, 0, 0, t27 * t18, -t7, -g(1) * (t11 * t14 - t28 * t9) - g(2) * (t11 * t28 + t9 * t14), -g(1) * (t11 * t15 + t29 * t9) - g(2) * (-t11 * t29 + t9 * t15), t7, t21 + (-g(2) * pkin(6) + g(1) * t22) * t9 + (-g(1) * pkin(6) - g(2) * t22) * t11, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t5 * t15, -t5 * t14, -t6, -g(3) * t24 + t26 * (pkin(3) * t16 - qJ(4) * t18), 0, 0, 0, 0, 0, t5 * t10, -t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t33 * t8, g(1) * t4 - g(2) * t2 + t10 * t33;];
taug_reg = t20;
