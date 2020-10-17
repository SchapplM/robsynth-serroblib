% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:39:57
% EndTime: 2019-05-05 22:39:58
% DurationCPUTime: 0.20s
% Computational Cost: add. (250->55), mult. (222->69), div. (0->0), fcn. (220->10), ass. (0->38)
t27 = pkin(10) + qJ(3);
t23 = cos(t27);
t24 = qJ(4) + t27;
t20 = sin(t24);
t21 = cos(t24);
t38 = t21 * pkin(4) + t20 * qJ(5);
t46 = pkin(3) * t23 + t38;
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t15 = g(1) * t33 + g(2) * t31;
t4 = g(3) * t20 + t15 * t21;
t45 = pkin(4) * t20;
t43 = g(3) * t21;
t30 = sin(qJ(6));
t42 = t31 * t30;
t32 = cos(qJ(6));
t41 = t31 * t32;
t40 = t33 * t30;
t39 = t33 * t32;
t37 = qJ(5) * t21;
t22 = sin(t27);
t36 = -pkin(3) * t22 - t45;
t14 = g(1) * t31 - g(2) * t33;
t29 = cos(pkin(10));
t35 = t29 * pkin(2) + pkin(1) + t46;
t26 = -pkin(8) - pkin(7) - qJ(2);
t13 = t33 * t37;
t12 = t31 * t37;
t10 = -t20 * t42 + t39;
t9 = t20 * t41 + t40;
t8 = t20 * t40 + t41;
t7 = t20 * t39 - t42;
t6 = t14 * t21;
t5 = t14 * t20;
t3 = t15 * t20 - t43;
t2 = t4 * t32;
t1 = t4 * t30;
t11 = [0, t14, t15, t14 * t29, -t14 * sin(pkin(10)) -t15, -g(1) * (-t31 * pkin(1) + t33 * qJ(2)) - g(2) * (t33 * pkin(1) + t31 * qJ(2)) 0, 0, 0, 0, 0, t14 * t23, -t14 * t22, 0, 0, 0, 0, 0, t6, -t5, -t15, -t6, t5 (g(1) * t26 - g(2) * t35) * t33 + (g(1) * t35 + g(2) * t26) * t31, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7; 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t23 + t15 * t22, g(3) * t22 + t15 * t23, 0, 0, 0, 0, 0, t3, t4, 0, -t3, -t4, -g(1) * (t36 * t33 + t13) - g(2) * (t36 * t31 + t12) - g(3) * t46, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, -t3, -t4, -g(1) * (-t33 * t45 + t13) - g(2) * (-t31 * t45 + t12) - g(3) * t38, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t32 * t43, g(1) * t8 - g(2) * t10 - t30 * t43;];
taug_reg  = t11;
