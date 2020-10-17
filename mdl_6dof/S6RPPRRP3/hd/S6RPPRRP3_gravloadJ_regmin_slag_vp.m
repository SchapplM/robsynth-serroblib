% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:51:54
% EndTime: 2019-05-05 14:51:55
% DurationCPUTime: 0.23s
% Computational Cost: add. (205->57), mult. (251->74), div. (0->0), fcn. (265->8), ass. (0->34)
t19 = qJ(1) + pkin(9);
t16 = sin(t19);
t17 = cos(t19);
t47 = -g(1) * t16 + g(2) * t17;
t21 = sin(qJ(4));
t24 = cos(qJ(4));
t37 = t24 * pkin(8);
t46 = t21 * pkin(4) - t37;
t44 = -g(3) * t21 - t47 * t24;
t39 = g(3) * t24;
t20 = sin(qJ(5));
t36 = t20 * t21;
t23 = cos(qJ(5));
t35 = t21 * t23;
t25 = cos(qJ(1));
t34 = t25 * pkin(1) + t17 * pkin(2) + t16 * qJ(3);
t22 = sin(qJ(1));
t33 = -t22 * pkin(1) + t17 * qJ(3);
t6 = t16 * t36 - t17 * t23;
t8 = t16 * t23 + t17 * t36;
t32 = g(1) * t8 + g(2) * t6;
t31 = g(1) * t17 + g(2) * t16;
t29 = g(1) * t22 - g(2) * t25;
t28 = pkin(5) * t23 + qJ(6) * t20 + pkin(4);
t1 = g(1) * t6 - g(2) * t8 + t20 * t39;
t7 = t16 * t35 + t17 * t20;
t9 = -t16 * t20 + t17 * t35;
t27 = g(1) * t7 - g(2) * t9 + t23 * t39;
t10 = t31 * t24;
t5 = -t21 * t47 + t39;
t4 = t44 * t23;
t3 = t44 * t20;
t2 = -g(1) * t9 - g(2) * t7;
t11 = [0, t29, g(1) * t25 + g(2) * t22, t29 * pkin(1), t47, -t31, -g(1) * (-t16 * pkin(2) + t33) - g(2) * t34, 0, 0, 0, 0, 0, -t31 * t21, -t10, 0, 0, 0, 0, 0, t2, t32, t2, t10, -t32, -g(1) * (t9 * pkin(5) + t8 * qJ(6) + t46 * t17 + t33) - g(2) * (t7 * pkin(5) + t17 * pkin(7) + t6 * qJ(6) + t34) + (-g(1) * (-pkin(2) - pkin(7)) - g(2) * t46) * t16; 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t5, 0, 0, 0, 0, 0, -t4, t3, -t4, -t5, -t3, -g(3) * (-t28 * t21 + t37) + t47 * (pkin(8) * t21 + t28 * t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t27, t1, 0, -t27, -g(1) * (-t6 * pkin(5) + t7 * qJ(6)) - g(2) * (t8 * pkin(5) - t9 * qJ(6)) - (-pkin(5) * t20 + qJ(6) * t23) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t11;
