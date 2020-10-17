% Calculate minimal parameter regressor of gravitation load for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [6x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:29:09
% EndTime: 2019-05-04 20:29:09
% DurationCPUTime: 0.30s
% Computational Cost: add. (453->70), mult. (1256->122), div. (0->0), fcn. (1627->14), ass. (0->56)
t58 = sin(pkin(12));
t59 = sin(pkin(11));
t48 = t59 * t58;
t62 = cos(pkin(12));
t63 = cos(pkin(11));
t55 = t63 * t62;
t65 = cos(pkin(6));
t39 = -t65 * t55 + t48;
t60 = sin(pkin(7));
t61 = sin(pkin(6));
t52 = t61 * t60;
t64 = cos(pkin(7));
t76 = t39 * t64 + t63 * t52;
t50 = t59 * t62;
t53 = t63 * t58;
t40 = t65 * t50 + t53;
t49 = t59 * t61;
t75 = t40 * t64 - t60 * t49;
t74 = t62 * t64 * t61 + t65 * t60;
t30 = sin(qJ(3));
t51 = t61 * t58;
t68 = cos(qJ(3));
t15 = t74 * t30 + t68 * t51;
t20 = t65 * t53 + t50;
t7 = t20 * t68 - t76 * t30;
t21 = -t65 * t48 + t55;
t9 = t21 * t68 - t75 * t30;
t73 = g(1) * t9 + g(2) * t7 + g(3) * t15;
t14 = t30 * t51 - t74 * t68;
t6 = t20 * t30 + t76 * t68;
t8 = t21 * t30 + t75 * t68;
t44 = g(1) * t8 + g(2) * t6 + g(3) * t14;
t28 = sin(qJ(5));
t32 = cos(qJ(4));
t67 = t28 * t32;
t31 = cos(qJ(5));
t66 = t31 * t32;
t54 = t63 * t61;
t29 = sin(qJ(4));
t38 = -t62 * t52 + t65 * t64;
t10 = t15 * t29 - t38 * t32;
t33 = t39 * t60 - t64 * t54;
t2 = t7 * t29 - t33 * t32;
t34 = t40 * t60 + t64 * t49;
t4 = t9 * t29 - t34 * t32;
t46 = g(1) * t4 + g(2) * t2 + g(3) * t10;
t11 = t15 * t32 + t38 * t29;
t3 = t33 * t29 + t7 * t32;
t5 = t34 * t29 + t9 * t32;
t45 = g(1) * t5 + g(2) * t3 + g(3) * t11;
t35 = -g(1) * (-t5 * t28 + t8 * t31) - g(2) * (-t3 * t28 + t6 * t31) - g(3) * (-t11 * t28 + t14 * t31);
t27 = -qJ(6) - pkin(10);
t26 = t31 * pkin(5) + pkin(4);
t18 = -g(1) * t49 + g(2) * t54 - g(3) * t65;
t1 = t44 * t29;
t12 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, t44, t73, 0, 0, 0, 0, 0, t44 * t32, -t1, 0, 0, 0, 0, 0, -g(1) * (t9 * t28 - t8 * t66) - g(2) * (t7 * t28 - t6 * t66) - g(3) * (-t14 * t66 + t15 * t28) -g(1) * (t9 * t31 + t8 * t67) - g(2) * (t7 * t31 + t6 * t67) - g(3) * (t14 * t67 + t15 * t31) t1, -t73 * (pkin(5) * t28 + pkin(9)) + t44 * (t26 * t32 - t27 * t29 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t45, 0, 0, 0, 0, 0, t46 * t31, -t46 * t28, -t45, -g(1) * (-t4 * t26 - t5 * t27) - g(2) * (-t2 * t26 - t3 * t27) - g(3) * (-t10 * t26 - t11 * t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -g(1) * (-t8 * t28 - t5 * t31) - g(2) * (-t6 * t28 - t3 * t31) - g(3) * (-t11 * t31 - t14 * t28) 0, t35 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46;];
taug_reg  = t12;
