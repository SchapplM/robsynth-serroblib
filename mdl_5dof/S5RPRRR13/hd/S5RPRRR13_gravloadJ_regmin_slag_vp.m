% Calculate minimal parameter regressor of gravitation load for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% taug_reg [5x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:27
% EndTime: 2019-12-31 19:15:28
% DurationCPUTime: 0.14s
% Computational Cost: add. (92->36), mult. (158->56), div. (0->0), fcn. (176->8), ass. (0->32)
t19 = sin(qJ(1));
t22 = cos(qJ(1));
t11 = g(1) * t19 - g(2) * t22;
t18 = sin(qJ(3));
t21 = cos(qJ(3));
t24 = -g(3) * t18 + t11 * t21;
t33 = g(3) * t21;
t16 = qJ(4) + qJ(5);
t13 = sin(t16);
t32 = t19 * t13;
t14 = cos(t16);
t31 = t19 * t14;
t17 = sin(qJ(4));
t30 = t19 * t17;
t20 = cos(qJ(4));
t29 = t19 * t20;
t28 = t22 * t13;
t27 = t22 * t14;
t26 = t22 * t17;
t25 = t22 * t20;
t12 = g(1) * t22 + g(2) * t19;
t10 = t18 * t25 - t30;
t9 = t18 * t26 + t29;
t8 = t18 * t29 + t26;
t7 = -t18 * t30 + t25;
t6 = t18 * t27 - t32;
t5 = t18 * t28 + t31;
t4 = t18 * t31 + t28;
t3 = -t18 * t32 + t27;
t2 = g(1) * t4 - g(2) * t6 + t14 * t33;
t1 = -g(1) * t3 - g(2) * t5 + t13 * t33;
t15 = [0, t11, t12, -t11, -t12, -g(1) * (-t19 * pkin(1) + t22 * qJ(2)) - g(2) * (t22 * pkin(1) + t19 * qJ(2)), 0, 0, 0, 0, 0, -t12 * t18, -t12 * t21, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t11 * t18 + t33, 0, 0, 0, 0, 0, -t24 * t20, t24 * t17, 0, 0, 0, 0, 0, -t24 * t14, t24 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t17 * t33, g(1) * t8 - g(2) * t10 + t20 * t33, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t15;
