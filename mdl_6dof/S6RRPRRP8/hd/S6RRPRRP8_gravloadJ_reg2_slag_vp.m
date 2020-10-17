% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:22:47
% EndTime: 2019-05-06 18:22:49
% DurationCPUTime: 0.51s
% Computational Cost: add. (527->109), mult. (578->152), div. (0->0), fcn. (597->10), ass. (0->67)
t46 = cos(pkin(10));
t34 = t46 * pkin(3) + pkin(2);
t44 = pkin(10) + qJ(4);
t36 = cos(t44);
t24 = pkin(4) * t36 + t34;
t50 = cos(qJ(2));
t20 = t50 * t24;
t47 = -pkin(8) - qJ(3);
t43 = -pkin(9) + t47;
t48 = sin(qJ(2));
t88 = -t48 * t43 + t20;
t49 = sin(qJ(1));
t51 = cos(qJ(1));
t30 = g(1) * t51 + g(2) * t49;
t17 = -g(3) * t50 + t30 * t48;
t35 = sin(t44);
t87 = pkin(4) * t35;
t86 = g(1) * t49;
t83 = g(3) * t48;
t45 = sin(pkin(10));
t81 = t45 * pkin(3);
t37 = qJ(5) + t44;
t32 = sin(t37);
t80 = t32 * t48;
t33 = cos(t37);
t79 = t33 * t48;
t77 = t48 * t51;
t76 = t49 * t36;
t75 = t49 * t50;
t74 = t50 * t51;
t73 = t51 * t32;
t72 = t51 * t33;
t71 = t51 * t35;
t70 = t51 * t36;
t69 = t51 * t45;
t68 = t51 * t46;
t67 = t51 * pkin(1) + t49 * pkin(7);
t66 = t50 * t71;
t10 = t33 * t75 - t73;
t9 = t32 * t75 + t72;
t65 = -t9 * pkin(5) + t10 * qJ(6);
t11 = -t49 * t33 + t50 * t73;
t12 = t49 * t32 + t50 * t72;
t64 = -t11 * pkin(5) + t12 * qJ(6);
t63 = g(1) * t9 - g(2) * t11;
t62 = -g(2) * t51 + t86;
t61 = t50 * pkin(2) + t48 * qJ(3);
t59 = pkin(5) * t33 + qJ(6) * t32;
t57 = t50 * t34 - t48 * t47;
t13 = t35 * t75 + t70;
t25 = t81 + t87;
t54 = t24 * t74 + t49 * t25 - t43 * t77 + t67;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t80;
t3 = g(1) * t12 + g(2) * t10 + g(3) * t79;
t40 = t51 * pkin(7);
t53 = t51 * t25 + t40 + (-pkin(1) - t88) * t49;
t18 = t30 * t50 + t83;
t31 = pkin(4) * t76;
t26 = qJ(6) * t79;
t23 = t62 * t48;
t16 = t49 * t35 + t50 * t70;
t15 = -t66 + t76;
t14 = -t36 * t75 + t71;
t6 = t17 * t33;
t5 = t17 * t32;
t4 = g(1) * t10 - g(2) * t12;
t2 = [0, 0, 0, 0, 0, 0, t62, t30, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t50, -t23, -t30, -g(1) * (-t49 * pkin(1) + t40) - g(2) * t67, 0, 0, 0, 0, 0, 0, -g(1) * (-t46 * t75 + t69) - g(2) * (t49 * t45 + t50 * t68) -g(1) * (t45 * t75 + t68) - g(2) * (t49 * t46 - t50 * t69) t23, -g(1) * t40 - g(2) * (t61 * t51 + t67) - (-pkin(1) - t61) * t86, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, t23, -g(1) * (pkin(3) * t69 + t40) - g(2) * (t34 * t74 - t47 * t77 + t67) + (-g(1) * (-pkin(1) - t57) - g(2) * t81) * t49, 0, 0, 0, 0, 0, 0, t4, -t63, t23, -g(1) * t53 - g(2) * t54, 0, 0, 0, 0, 0, 0, t4, t23, t63, -g(1) * (-t10 * pkin(5) - t9 * qJ(6) + t53) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t54); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t46, -t17 * t45, -t18, -g(3) * t61 + t30 * (pkin(2) * t48 - qJ(3) * t50) 0, 0, 0, 0, 0, 0, t17 * t36, -t17 * t35, -t18, -g(3) * t57 + t30 * (t34 * t48 + t47 * t50) 0, 0, 0, 0, 0, 0, t6, -t5, -t18, -g(3) * t88 + t30 * (t24 * t48 + t43 * t50) 0, 0, 0, 0, 0, 0, t6, -t18, t5, -g(3) * t20 + (-g(3) * t59 + t30 * t43) * t50 + (g(3) * t43 + t30 * (t24 + t59)) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t35 * t83, g(1) * t16 - g(2) * t14 + t36 * t83, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t31 + (g(2) * t70 + t18 * t35) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (-pkin(4) * t66 + t31 + t64) - g(2) * (-t13 * pkin(4) + t65) - g(3) * (t26 + (-pkin(5) * t32 - t87) * t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t64 - g(2) * t65 - g(3) * (-pkin(5) * t80 + t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
