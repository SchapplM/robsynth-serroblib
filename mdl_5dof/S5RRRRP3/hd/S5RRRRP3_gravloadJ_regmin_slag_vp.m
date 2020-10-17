% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP3
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
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:25
% EndTime: 2019-12-31 21:49:26
% DurationCPUTime: 0.12s
% Computational Cost: add. (241->33), mult. (164->38), div. (0->0), fcn. (155->8), ass. (0->27)
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t28 = t24 * pkin(4) + t22 * qJ(5);
t38 = -pkin(3) - t28;
t21 = qJ(1) + qJ(2);
t20 = qJ(3) + t21;
t17 = cos(t20);
t16 = sin(t20);
t37 = g(1) * t16;
t5 = -g(2) * t17 + t37;
t6 = g(1) * t17 + g(2) * t16;
t13 = t17 * pkin(8);
t18 = sin(t21);
t31 = -pkin(2) * t18 + t13;
t30 = t16 * pkin(8) - t38 * t17;
t19 = cos(t21);
t29 = pkin(2) * t19 + t30;
t26 = t38 * t37;
t25 = cos(qJ(1));
t23 = sin(qJ(1));
t8 = g(1) * t19 + g(2) * t18;
t7 = g(1) * t18 - g(2) * t19;
t4 = t5 * t24;
t3 = t5 * t22;
t2 = g(3) * t22 + t6 * t24;
t1 = -g(3) * t24 + t6 * t22;
t9 = [0, g(1) * t23 - g(2) * t25, g(1) * t25 + g(2) * t23, 0, t7, t8, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * (-t23 * pkin(1) + t31) - g(2) * (t25 * pkin(1) + t29) - t26; 0, 0, 0, 0, t7, t8, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * t31 - g(2) * t29 - t26; 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * t13 - g(2) * t30 - t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t28 + t6 * (pkin(4) * t22 - qJ(5) * t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t9;
