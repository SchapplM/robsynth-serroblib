% Calculate minimal parameter regressor of gravitation load for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPPRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:41:26
% EndTime: 2019-05-05 13:41:26
% DurationCPUTime: 0.14s
% Computational Cost: add. (124->41), mult. (214->58), div. (0->0), fcn. (260->10), ass. (0->28)
t18 = pkin(10) + qJ(5);
t11 = sin(t18);
t12 = cos(t18);
t31 = sin(pkin(9));
t32 = cos(pkin(9));
t36 = sin(qJ(1));
t37 = cos(qJ(1));
t4 = -t36 * t31 - t37 * t32;
t5 = t37 * t31 - t36 * t32;
t28 = g(1) * t4 + g(2) * t5;
t40 = -g(3) * t12 + t28 * t11;
t39 = g(3) * t11;
t21 = sin(qJ(6));
t35 = t12 * t21;
t22 = cos(qJ(6));
t34 = t12 * t22;
t33 = t37 * pkin(1) + t36 * qJ(2);
t30 = t37 * pkin(2) + t33;
t29 = g(1) * t5 - g(2) * t4;
t27 = -t36 * pkin(1) + t37 * qJ(2);
t26 = t4 * t21 + t5 * t34;
t25 = -t4 * t22 + t5 * t35;
t24 = -t36 * pkin(2) + t27;
t7 = g(1) * t37 + g(2) * t36;
t6 = g(1) * t36 - g(2) * t37;
t2 = t5 * t21 - t4 * t34;
t1 = t5 * t22 + t4 * t35;
t3 = [0, t6, t7, t6, -t7, -g(1) * t27 - g(2) * t33, -t29, t28, -g(1) * t24 - g(2) * t30, -t29 * cos(pkin(10)) t29 * sin(pkin(10)) -t28, -g(1) * (t5 * pkin(3) + t4 * qJ(4) + t24) - g(2) * (-t4 * pkin(3) + t5 * qJ(4) + t30) 0, 0, 0, 0, 0, -t29 * t12, t29 * t11, 0, 0, 0, 0, 0, -g(1) * t26 - g(2) * t2, g(1) * t25 - g(2) * t1; 0, 0, 0, 0, 0, -t6, 0, 0, -t6, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t28 * t12 - t39, 0, 0, 0, 0, 0, -t40 * t22, t40 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t25 - t21 * t39, g(1) * t2 - g(2) * t26 - t22 * t39;];
taug_reg  = t3;
