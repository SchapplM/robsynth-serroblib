% Calculate minimal parameter regressor of gravitation load for
% S5RPPRP1
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
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:45
% EndTime: 2022-01-23 09:12:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (127->35), mult. (144->54), div. (0->0), fcn. (151->8), ass. (0->28)
t18 = sin(pkin(8));
t19 = cos(pkin(8));
t23 = cos(qJ(4));
t39 = -t18 * (-qJ(5) - pkin(6)) + (t23 * pkin(4) + pkin(3)) * t19;
t21 = sin(qJ(4));
t36 = g(3) * t18;
t17 = qJ(1) + pkin(7);
t14 = sin(t17);
t15 = cos(t17);
t32 = t19 * t21;
t5 = t14 * t32 + t15 * t23;
t7 = t14 * t23 - t15 * t32;
t1 = -g(1) * t7 + g(2) * t5 + t21 * t36;
t34 = t15 * t21;
t31 = t19 * t23;
t24 = cos(qJ(1));
t29 = t24 * pkin(1) + t15 * pkin(2) + t14 * qJ(3);
t22 = sin(qJ(1));
t28 = -t22 * pkin(1) + t15 * qJ(3);
t27 = -g(1) * t15 - g(2) * t14;
t26 = g(1) * t14 - g(2) * t15;
t25 = g(1) * t22 - g(2) * t24;
t8 = t14 * t21 + t15 * t31;
t6 = -t14 * t31 + t34;
t4 = -g(1) * t6 - g(2) * t8;
t3 = -g(1) * t5 - g(2) * t7;
t2 = g(1) * t8 - g(2) * t6 + t23 * t36;
t9 = [0, t25, g(1) * t24 + g(2) * t22, t25 * pkin(1), t26 * t19, t27, -g(1) * (-t14 * pkin(2) + t28) - g(2) * t29, 0, 0, 0, 0, 0, t4, t3, t4, t3, t26 * t18, -g(1) * (pkin(4) * t34 + t28) - g(2) * (t39 * t15 + t29) + (-g(1) * (-pkin(2) - t39) - g(2) * pkin(4) * t21) * t14; 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t19 + t18 * t27;];
taug_reg = t9;
