% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:43
% EndTime: 2019-12-31 19:49:43
% DurationCPUTime: 0.11s
% Computational Cost: add. (173->34), mult. (130->39), div. (0->0), fcn. (119->8), ass. (0->29)
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t28 = t24 * pkin(4) + t22 * qJ(5);
t40 = -pkin(3) - t28;
t21 = qJ(1) + qJ(2);
t17 = pkin(8) + t21;
t15 = cos(t17);
t14 = sin(t17);
t37 = g(1) * t14;
t39 = -g(2) * t15 + t37;
t29 = g(1) * t15 + g(2) * t14;
t18 = sin(t21);
t38 = pkin(2) * t18;
t19 = cos(t21);
t16 = pkin(2) * t19;
t31 = t14 * pkin(7) - t40 * t15 + t16;
t23 = sin(qJ(1));
t30 = -t23 * pkin(1) - t38;
t6 = g(1) * t18 - g(2) * t19;
t26 = t40 * t37;
t25 = cos(qJ(1));
t20 = t25 * pkin(1);
t12 = t15 * pkin(7);
t7 = g(1) * t19 + g(2) * t18;
t4 = t39 * t24;
t3 = t39 * t22;
t2 = g(3) * t22 + t29 * t24;
t1 = -g(3) * t24 + t29 * t22;
t5 = [0, g(1) * t23 - g(2) * t25, g(1) * t25 + g(2) * t23, 0, t6, t7, -g(1) * t30 - g(2) * (t16 + t20), 0, 0, 0, 0, 0, t4, -t3, t4, -t29, t3, -g(1) * (t12 + t30) - g(2) * (t20 + t31) - t26; 0, 0, 0, 0, t6, t7, t6 * pkin(2), 0, 0, 0, 0, 0, t4, -t3, t4, -t29, t3, -g(1) * (t12 - t38) - g(2) * t31 - t26; 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t28 + t29 * (pkin(4) * t22 - qJ(5) * t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t5;
