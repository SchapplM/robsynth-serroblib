% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:34:22
% EndTime: 2019-05-05 17:34:22
% DurationCPUTime: 0.20s
% Computational Cost: add. (281->60), mult. (264->79), div. (0->0), fcn. (275->10), ass. (0->41)
t20 = qJ(3) + pkin(10);
t16 = cos(t20);
t14 = sin(t20);
t45 = t14 * pkin(8);
t49 = t16 * pkin(4) + t45;
t21 = qJ(1) + pkin(9);
t15 = sin(t21);
t17 = cos(t21);
t37 = g(1) * t17 + g(2) * t15;
t46 = g(3) * t14;
t23 = sin(qJ(5));
t43 = t15 * t23;
t26 = cos(qJ(5));
t42 = t15 * t26;
t41 = t17 * t23;
t40 = t17 * t26;
t27 = cos(qJ(3));
t18 = t27 * pkin(3);
t13 = t18 + pkin(2);
t28 = cos(qJ(1));
t39 = t28 * pkin(1) + t17 * t13;
t5 = t16 * t43 + t40;
t7 = t16 * t41 - t42;
t38 = g(1) * t5 - g(2) * t7;
t36 = g(1) * t15 - g(2) * t17;
t25 = sin(qJ(1));
t35 = g(1) * t25 - g(2) * t28;
t22 = -qJ(4) - pkin(7);
t34 = -t25 * pkin(1) - t17 * t22;
t33 = pkin(5) * t26 + qJ(6) * t23 + pkin(4);
t1 = g(1) * t7 + g(2) * t5 + t23 * t46;
t6 = t16 * t42 - t41;
t8 = t16 * t40 + t43;
t32 = g(1) * t8 + g(2) * t6 + t26 * t46;
t31 = -g(3) * t16 + t14 * t37;
t24 = sin(qJ(3));
t30 = -g(3) * t27 + t24 * t37;
t4 = t31 * t26;
t3 = t31 * t23;
t2 = g(1) * t6 - g(2) * t8;
t9 = [0, t35, g(1) * t28 + g(2) * t25, t35 * pkin(1), 0, 0, 0, 0, 0, t36 * t27, -t36 * t24, -t37, -g(1) * (-t15 * t13 + t34) - g(2) * (-t15 * t22 + t39) 0, 0, 0, 0, 0, t2, -t38, t2, t36 * t14, t38, -g(1) * (-t6 * pkin(5) - t5 * qJ(6) + t34) - g(2) * (t8 * pkin(5) + t7 * qJ(6) + t49 * t17 + t39) + (-g(1) * (-t13 - t49) + g(2) * t22) * t15; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, g(3) * t24 + t27 * t37, 0, t30 * pkin(3), 0, 0, 0, 0, 0, t4, -t3, t4, -t16 * t37 - t46, t3, -g(3) * (t16 * t33 + t18 + t45) + t37 * (pkin(3) * t24 - pkin(8) * t16 + t14 * t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t32, t1, 0, -t32, -g(1) * (-t7 * pkin(5) + t8 * qJ(6)) - g(2) * (-t5 * pkin(5) + t6 * qJ(6)) - (-pkin(5) * t23 + qJ(6) * t26) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t9;
