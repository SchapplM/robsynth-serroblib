% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRP1
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
% taug_reg [6x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:30:13
% EndTime: 2019-05-05 17:30:14
% DurationCPUTime: 0.19s
% Computational Cost: add. (189->50), mult. (176->68), div. (0->0), fcn. (169->10), ass. (0->37)
t16 = qJ(1) + pkin(9);
t10 = sin(t16);
t12 = cos(t16);
t31 = g(1) * t12 + g(2) * t10;
t15 = qJ(3) + pkin(10);
t11 = cos(t15);
t22 = cos(qJ(5));
t34 = t12 * t22;
t19 = sin(qJ(5));
t37 = t10 * t19;
t1 = t11 * t37 + t34;
t35 = t12 * t19;
t36 = t10 * t22;
t3 = -t11 * t35 + t36;
t9 = sin(t15);
t43 = g(3) * t9;
t47 = -g(1) * t3 + g(2) * t1 + t19 * t43;
t46 = -g(3) * t11 + t31 * t9;
t21 = sin(qJ(1));
t39 = t21 * pkin(1);
t24 = cos(qJ(1));
t23 = cos(qJ(3));
t13 = t23 * pkin(3);
t8 = t13 + pkin(2);
t38 = t24 * pkin(1) + t12 * t8;
t18 = -qJ(4) - pkin(7);
t32 = pkin(5) * t19 - t18;
t30 = g(1) * t10 - g(2) * t12;
t29 = g(1) * t21 - g(2) * t24;
t17 = -qJ(6) - pkin(8);
t7 = t22 * pkin(5) + pkin(4);
t28 = t11 * t7 - t9 * t17;
t20 = sin(qJ(3));
t25 = -g(3) * t23 + t31 * t20;
t4 = t11 * t34 + t37;
t2 = -t11 * t36 + t35;
t5 = [0, t29, g(1) * t24 + g(2) * t21, t29 * pkin(1), 0, 0, 0, 0, 0, t30 * t23, -t30 * t20, -t31, -g(1) * (-t10 * t8 - t12 * t18 - t39) - g(2) * (-t10 * t18 + t38) 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3, t30 * t9, g(1) * t39 - g(2) * t38 + (-g(1) * t32 - g(2) * t28) * t12 + (-g(1) * (-t28 - t8) - g(2) * t32) * t10; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, g(3) * t20 + t31 * t23, 0, t25 * pkin(3), 0, 0, 0, 0, 0, t46 * t22, -t46 * t19, -t31 * t11 - t43, -g(3) * (t13 + t28) + t31 * (pkin(3) * t20 + t11 * t17 + t7 * t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, 0, 0, 0, 0, 0, 0, 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, g(1) * t4 - g(2) * t2 + t22 * t43, 0, t47 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46;];
taug_reg  = t5;
