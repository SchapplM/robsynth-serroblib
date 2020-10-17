% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:17
% EndTime: 2019-12-31 18:04:18
% DurationCPUTime: 0.11s
% Computational Cost: add. (96->33), mult. (171->56), div. (0->0), fcn. (193->8), ass. (0->32)
t25 = sin(qJ(1));
t34 = g(1) * t25;
t27 = cos(qJ(1));
t33 = t27 * pkin(1) + t25 * qJ(2);
t14 = g(1) * t27 + g(2) * t25;
t13 = -g(2) * t27 + t34;
t22 = sin(pkin(8));
t23 = cos(pkin(8));
t32 = pkin(2) * t23 + qJ(3) * t22;
t21 = qJ(4) + qJ(5);
t15 = sin(t21);
t16 = cos(t21);
t31 = t23 * t15 - t22 * t16;
t30 = t22 * t15 + t23 * t16;
t24 = sin(qJ(4));
t26 = cos(qJ(4));
t29 = t22 * t26 - t23 * t24;
t28 = t22 * t24 + t23 * t26;
t18 = t27 * qJ(2);
t12 = t13 * t23;
t11 = t13 * t22;
t10 = t28 * t27;
t9 = t29 * t27;
t8 = t28 * t25;
t7 = t29 * t25;
t6 = t30 * t27;
t5 = t31 * t27;
t4 = t30 * t25;
t3 = t31 * t25;
t2 = g(1) * t6 + g(2) * t4 - g(3) * t31;
t1 = g(1) * t5 + g(2) * t3 + g(3) * t30;
t17 = [0, t13, t14, t12, -t11, -t14, -g(1) * (-t25 * pkin(1) + t18) - g(2) * t33, t12, -t14, t11, -g(1) * t18 - g(2) * (t32 * t27 + t33) - (-pkin(1) - t32) * t34, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t10, g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t6, -g(1) * t3 + g(2) * t5; 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t23 - t14 * t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7 + g(3) * t28, g(1) * t10 + g(2) * t8 + g(3) * t29, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t17;
