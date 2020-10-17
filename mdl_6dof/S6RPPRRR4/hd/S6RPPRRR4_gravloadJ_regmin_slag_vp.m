% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:41:33
% EndTime: 2019-05-05 15:41:34
% DurationCPUTime: 0.19s
% Computational Cost: add. (160->43), mult. (292->72), div. (0->0), fcn. (368->10), ass. (0->36)
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t34 = sin(pkin(10));
t35 = cos(pkin(10));
t41 = sin(qJ(1));
t42 = cos(qJ(1));
t7 = -t41 * t34 - t42 * t35;
t8 = t42 * t34 - t41 * t35;
t32 = g(1) * t7 + g(2) * t8;
t26 = -g(3) * t24 + t32 * t22;
t44 = g(3) * t22;
t20 = qJ(5) + qJ(6);
t14 = sin(t20);
t40 = t14 * t24;
t15 = cos(t20);
t39 = t15 * t24;
t21 = sin(qJ(5));
t38 = t21 * t24;
t23 = cos(qJ(5));
t37 = t23 * t24;
t36 = t42 * pkin(1) + t41 * qJ(2);
t33 = g(1) * t8 - g(2) * t7;
t31 = -t41 * pkin(1) + t42 * qJ(2);
t30 = t7 * t14 + t8 * t39;
t29 = -t7 * t15 + t8 * t40;
t28 = t7 * t21 + t8 * t37;
t27 = -t7 * t23 + t8 * t38;
t10 = g(1) * t42 + g(2) * t41;
t9 = g(1) * t41 - g(2) * t42;
t6 = t8 * t21 - t7 * t37;
t5 = t8 * t23 + t7 * t38;
t4 = t8 * t14 - t7 * t39;
t3 = t8 * t15 + t7 * t40;
t2 = g(1) * t4 - g(2) * t30 - t15 * t44;
t1 = -g(1) * t3 - g(2) * t29 - t14 * t44;
t11 = [0, t9, t10, t9, -t10, -g(1) * t31 - g(2) * t36, -t33, t32, -g(1) * (-t41 * pkin(2) + t31) - g(2) * (t42 * pkin(2) + t36) 0, 0, 0, 0, 0, -t33 * t24, t33 * t22, 0, 0, 0, 0, 0, -g(1) * t28 - g(2) * t6, g(1) * t27 - g(2) * t5, 0, 0, 0, 0, 0, -g(1) * t30 - g(2) * t4, g(1) * t29 - g(2) * t3; 0, 0, 0, 0, 0, -t9, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t32 * t24 - t44, 0, 0, 0, 0, 0, -t26 * t23, t26 * t21, 0, 0, 0, 0, 0, -t26 * t15, t26 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t27 - t21 * t44, g(1) * t6 - g(2) * t28 - t23 * t44, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t11;
