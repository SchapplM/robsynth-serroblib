% Calculate minimal parameter regressor of gravitation load for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:24
% EndTime: 2019-12-31 17:49:24
% DurationCPUTime: 0.10s
% Computational Cost: add. (128->35), mult. (104->40), div. (0->0), fcn. (95->8), ass. (0->22)
t14 = qJ(1) + pkin(7);
t11 = cos(t14);
t9 = sin(t14);
t24 = g(1) * t11 + g(2) * t9;
t18 = sin(qJ(1));
t26 = t18 * pkin(1);
t25 = g(1) * t9 - g(2) * t11;
t19 = cos(qJ(1));
t23 = g(1) * t18 - g(2) * t19;
t13 = pkin(8) + qJ(4);
t10 = cos(t13);
t8 = sin(t13);
t21 = t10 * pkin(4) + t8 * qJ(5);
t16 = cos(pkin(8));
t20 = t16 * pkin(3) + pkin(2) + t21;
t17 = -pkin(6) - qJ(3);
t12 = t19 * pkin(1);
t4 = t25 * t10;
t3 = t25 * t8;
t2 = g(3) * t8 + t24 * t10;
t1 = -g(3) * t10 + t24 * t8;
t5 = [0, t23, g(1) * t19 + g(2) * t18, t23 * pkin(1), t25 * t16, -t25 * sin(pkin(8)), -t24, -g(1) * (-t9 * pkin(2) + t11 * qJ(3) - t26) - g(2) * (t11 * pkin(2) + t9 * qJ(3) + t12), 0, 0, 0, 0, 0, t4, -t3, t4, -t24, t3, g(1) * t26 - g(2) * t12 + (g(1) * t20 + g(2) * t17) * t9 + (g(1) * t17 - g(2) * t20) * t11; 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, -t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t21 + t24 * (pkin(4) * t8 - qJ(5) * t10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t5;
