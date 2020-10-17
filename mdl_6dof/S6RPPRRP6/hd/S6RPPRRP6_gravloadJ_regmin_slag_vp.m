% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:01:23
% EndTime: 2019-05-05 15:01:23
% DurationCPUTime: 0.20s
% Computational Cost: add. (107->56), mult. (263->69), div. (0->0), fcn. (277->6), ass. (0->34)
t21 = sin(qJ(4));
t24 = cos(qJ(4));
t38 = t24 * pkin(8);
t46 = t21 * pkin(4) - t38;
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t45 = -g(1) * t25 - g(2) * t22;
t44 = -g(3) * t21 - t45 * t24;
t40 = g(3) * t24;
t20 = sin(qJ(5));
t37 = t22 * t20;
t23 = cos(qJ(5));
t36 = t22 * t23;
t35 = t25 * t20;
t34 = t25 * t23;
t33 = -pkin(1) - qJ(3);
t32 = t25 * pkin(1) + t22 * qJ(2);
t31 = t25 * qJ(3) + t32;
t6 = t21 * t37 - t34;
t8 = t21 * t35 + t36;
t30 = g(1) * t6 - g(2) * t8;
t11 = g(1) * t22 - g(2) * t25;
t29 = pkin(5) * t23 + qJ(6) * t20 + pkin(4);
t1 = g(1) * t8 + g(2) * t6 + t20 * t40;
t7 = t21 * t36 + t35;
t9 = t21 * t34 - t37;
t28 = g(1) * t9 + g(2) * t7 + t23 * t40;
t17 = t25 * qJ(2);
t10 = t11 * t24;
t5 = -t21 * t45 + t40;
t4 = t44 * t23;
t3 = t44 * t20;
t2 = g(1) * t7 - g(2) * t9;
t12 = [0, t11, -t45, -t11, t45, -g(1) * (-t22 * pkin(1) + t17) - g(2) * t32, t45, t11, -g(1) * (t33 * t22 + t17) - g(2) * t31, 0, 0, 0, 0, 0, t11 * t21, t10, 0, 0, 0, 0, 0, t2, -t30, t2, -t10, t30, -g(1) * (-t7 * pkin(5) - t25 * pkin(7) - t6 * qJ(6) + t17) - g(2) * (t9 * pkin(5) + t8 * qJ(6) + t46 * t25 + t31) + (-g(1) * (t33 - t46) + g(2) * pkin(7)) * t22; 0, 0, 0, 0, 0, -t11, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t5, 0, 0, 0, 0, 0, -t4, t3, -t4, -t5, -t3, -g(3) * (-t29 * t21 + t38) + t45 * (pkin(8) * t21 + t29 * t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t28, t1, 0, -t28, -g(1) * (-t8 * pkin(5) + t9 * qJ(6)) - g(2) * (-t6 * pkin(5) + t7 * qJ(6)) - (-pkin(5) * t20 + qJ(6) * t23) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t12;
