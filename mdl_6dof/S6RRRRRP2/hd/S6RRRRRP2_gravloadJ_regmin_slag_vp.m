% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:32:44
% EndTime: 2019-05-08 04:32:45
% DurationCPUTime: 0.28s
% Computational Cost: add. (477->67), mult. (442->91), div. (0->0), fcn. (454->10), ass. (0->51)
t34 = sin(qJ(5));
t37 = cos(qJ(5));
t63 = pkin(5) * t37 + qJ(6) * t34;
t33 = qJ(2) + qJ(3);
t30 = qJ(4) + t33;
t26 = sin(t30);
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t46 = g(1) * t39 + g(2) * t36;
t42 = t46 * t26;
t27 = cos(t30);
t62 = t27 * pkin(4) + t26 * pkin(10);
t28 = sin(t33);
t61 = pkin(3) * t28;
t59 = pkin(10) * t27;
t56 = g(3) * t26;
t55 = g(3) * t34;
t54 = t36 * t34;
t53 = t36 * t37;
t52 = t39 * t34;
t51 = t39 * t37;
t49 = t63 * t27 + t62;
t29 = cos(t33);
t23 = pkin(3) * t29;
t48 = t23 + t49;
t10 = t27 * t54 + t51;
t12 = t27 * t52 - t53;
t47 = g(1) * t10 - g(2) * t12;
t45 = g(1) * t36 - g(2) * t39;
t38 = cos(qJ(2));
t31 = t38 * pkin(2);
t44 = t23 + t31 + pkin(1) + t62;
t1 = g(1) * t12 + g(2) * t10 + t26 * t55;
t11 = t27 * t53 - t52;
t13 = t27 * t51 + t54;
t41 = g(1) * t13 + g(2) * t11 + t37 * t56;
t5 = -g(3) * t27 + t42;
t40 = (pkin(4) + t63) * t42;
t35 = sin(qJ(2));
t32 = -pkin(9) - pkin(8) - pkin(7);
t19 = t39 * t59;
t17 = t36 * t59;
t15 = -t35 * pkin(2) - t61;
t9 = t45 * t26;
t8 = g(3) * t28 + t46 * t29;
t7 = -g(3) * t29 + t46 * t28;
t6 = t46 * t27 + t56;
t4 = t5 * t37;
t3 = -t27 * t55 + t34 * t42;
t2 = g(1) * t11 - g(2) * t13;
t14 = [0, t45, t46, 0, 0, 0, 0, 0, t45 * t38, -t45 * t35, 0, 0, 0, 0, 0, t45 * t29, -t45 * t28, 0, 0, 0, 0, 0, t45 * t27, -t9, 0, 0, 0, 0, 0, t2, -t47, t2, t9, t47, -g(1) * (-t11 * pkin(5) - t10 * qJ(6)) - g(2) * (t13 * pkin(5) + t12 * qJ(6)) + (g(1) * t32 - g(2) * t44) * t39 + (g(1) * t44 + g(2) * t32) * t36; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t38 + t46 * t35, g(3) * t35 + t46 * t38, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * (t39 * t15 + t19) - g(2) * (t36 * t15 + t17) - g(3) * (t31 + t48) + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * (-t39 * t61 + t19) - g(2) * (-t36 * t61 + t17) - g(3) * t48 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * t19 - g(2) * t17 - g(3) * t49 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t41, t1, 0, -t41, -g(1) * (-t12 * pkin(5) + t13 * qJ(6)) - g(2) * (-t10 * pkin(5) + t11 * qJ(6)) - (-pkin(5) * t34 + qJ(6) * t37) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t14;
