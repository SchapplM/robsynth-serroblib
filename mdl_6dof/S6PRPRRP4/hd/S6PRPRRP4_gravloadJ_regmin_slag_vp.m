% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:50:28
% EndTime: 2019-05-04 23:50:29
% DurationCPUTime: 0.29s
% Computational Cost: add. (439->82), mult. (785->125), div. (0->0), fcn. (972->12), ass. (0->56)
t43 = sin(pkin(10));
t48 = sin(qJ(2));
t50 = cos(qJ(2));
t62 = cos(pkin(10));
t63 = cos(pkin(6));
t58 = t63 * t62;
t29 = t43 * t48 - t50 * t58;
t60 = t43 * t63;
t31 = t62 * t48 + t50 * t60;
t76 = -g(1) * t31 - g(2) * t29;
t30 = t43 * t50 + t48 * t58;
t32 = -t48 * t60 + t62 * t50;
t41 = pkin(11) + qJ(4);
t39 = sin(t41);
t40 = cos(t41);
t44 = sin(pkin(6));
t59 = t44 * t62;
t66 = t44 * t48;
t67 = t43 * t44;
t53 = g(3) * (-t39 * t66 + t63 * t40) + g(2) * (-t30 * t39 - t40 * t59) + g(1) * (-t32 * t39 + t40 * t67);
t70 = g(3) * t44;
t47 = sin(qJ(5));
t69 = t40 * t47;
t49 = cos(qJ(5));
t68 = t40 * t49;
t65 = t44 * t50;
t64 = t49 * t50;
t61 = t47 * t65;
t24 = t63 * t39 + t40 * t66;
t19 = t24 * t47 + t44 * t64;
t16 = t30 * t40 - t39 * t59;
t6 = t16 * t47 - t29 * t49;
t18 = t32 * t40 + t39 * t67;
t8 = t18 * t47 - t31 * t49;
t1 = g(1) * t8 + g(2) * t6 + g(3) * t19;
t20 = t24 * t49 - t61;
t7 = t16 * t49 + t29 * t47;
t9 = t18 * t49 + t31 * t47;
t55 = g(1) * t9 + g(2) * t7 + g(3) * t20;
t10 = -t29 * t69 - t30 * t49;
t12 = -t31 * t69 - t32 * t49;
t21 = t40 * t61 - t49 * t66;
t54 = g(1) * t12 + g(2) * t10 + g(3) * t21;
t52 = g(1) * t18 + g(2) * t16 + g(3) * t24;
t14 = g(3) * t65 + t76;
t51 = g(1) * t32 + g(2) * t30 + g(3) * t66;
t46 = -pkin(8) - qJ(3);
t45 = cos(pkin(11));
t22 = (t40 * t64 + t47 * t48) * t44;
t13 = -t31 * t68 + t32 * t47;
t11 = -t29 * t68 + t30 * t47;
t5 = t14 * t39;
t4 = t53 * t49;
t3 = t53 * t47;
t2 = -g(1) * t13 - g(2) * t11 - g(3) * t22;
t15 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t14, t51, -t14 * t45, t14 * sin(pkin(11)) -t51, -g(1) * (-t31 * pkin(2) + t32 * qJ(3)) - g(2) * (-t29 * pkin(2) + t30 * qJ(3)) - (pkin(2) * t50 + qJ(3) * t48) * t70, 0, 0, 0, 0, 0, -t14 * t40, t5, 0, 0, 0, 0, 0, t2, t54, t2, -t5, -t54, -g(1) * (t13 * pkin(5) + t12 * qJ(6) - t32 * t46) - g(2) * (t11 * pkin(5) + t10 * qJ(6) - t30 * t46) - g(3) * (t22 * pkin(5) + t21 * qJ(6)) + t46 * t48 * t70 + (-t50 * t70 - t76) * (t45 * pkin(3) + pkin(4) * t40 + pkin(9) * t39 + pkin(2)); 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t52, 0, 0, 0, 0, 0, -t4, t3, -t4, -t52, -t3, -t52 * pkin(9) - t53 * (pkin(5) * t49 + qJ(6) * t47 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t55, t1, 0, -t55, -g(1) * (-t8 * pkin(5) + t9 * qJ(6)) - g(2) * (-t6 * pkin(5) + t7 * qJ(6)) - g(3) * (-t19 * pkin(5) + t20 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t15;
