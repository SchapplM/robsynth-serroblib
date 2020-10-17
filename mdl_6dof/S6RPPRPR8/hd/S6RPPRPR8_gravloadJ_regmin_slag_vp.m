% Calculate minimal parameter regressor of gravitation load for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:38:34
% EndTime: 2019-05-05 14:38:34
% DurationCPUTime: 0.18s
% Computational Cost: add. (120->47), mult. (170->56), div. (0->0), fcn. (169->8), ass. (0->30)
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t39 = -g(1) * t22 + g(2) * t24;
t17 = pkin(9) + qJ(4);
t11 = sin(t17);
t12 = cos(t17);
t1 = g(3) * t12 - t11 * t39;
t36 = g(3) * t11;
t21 = sin(qJ(6));
t34 = t22 * t21;
t23 = cos(qJ(6));
t33 = t22 * t23;
t32 = t24 * t21;
t31 = t24 * t23;
t30 = t24 * pkin(1) + t22 * qJ(2);
t29 = g(2) * t30;
t10 = g(1) * t24 + g(2) * t22;
t27 = -t11 * pkin(4) + t12 * qJ(5);
t18 = sin(pkin(9));
t26 = pkin(3) * t18 - t27;
t20 = -pkin(7) - qJ(3);
t14 = t24 * qJ(2);
t8 = -t12 * t34 + t31;
t7 = -t12 * t33 - t32;
t6 = -t12 * t32 - t33;
t5 = -t12 * t31 + t34;
t4 = t10 * t12;
t3 = t10 * t11;
t2 = -t12 * t39 - t36;
t9 = [0, -t39, t10, t39, -t10, -g(1) * (-t22 * pkin(1) + t14) - t29, -t10 * t18, -t10 * cos(pkin(9)) -t39, -g(1) * (t14 + (-pkin(1) - qJ(3)) * t22) - g(2) * (t24 * qJ(3) + t30) 0, 0, 0, 0, 0, -t3, -t4, -t39, t3, t4, -g(1) * t14 - t29 + (-g(1) * t26 + g(2) * t20) * t24 + (-g(1) * (-pkin(1) + t20) - g(2) * t26) * t22, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, t39, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, t2, -t1, -g(3) * t27 + t39 * (pkin(4) * t12 + qJ(5) * t11) 0, 0, 0, 0, 0, -t1 * t21, -t1 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 - t23 * t36, g(1) * t8 - g(2) * t6 + t21 * t36;];
taug_reg  = t9;
