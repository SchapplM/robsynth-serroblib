% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:22:13
% EndTime: 2019-05-06 02:22:16
% DurationCPUTime: 0.88s
% Computational Cost: add. (1014->110), mult. (2849->183), div. (0->0), fcn. (3733->14), ass. (0->73)
t100 = sin(pkin(7));
t67 = cos(qJ(1));
t108 = sin(qJ(1));
t102 = cos(pkin(12));
t104 = cos(pkin(6));
t92 = t104 * t102;
t99 = sin(pkin(12));
t76 = t108 * t99 - t67 * t92;
t101 = sin(pkin(6));
t103 = cos(pkin(7));
t90 = t101 * t103;
t116 = t76 * t100 - t67 * t90;
t109 = cos(qJ(3));
t89 = t101 * t100;
t122 = t76 * t103 + t67 * t89;
t91 = t104 * t99;
t52 = t108 * t102 + t67 * t91;
t64 = sin(qJ(3));
t37 = -t52 * t109 + t122 * t64;
t63 = sin(qJ(4));
t66 = cos(qJ(4));
t22 = -t116 * t63 + t37 * t66;
t34 = t122 * t109 + t52 * t64;
t62 = sin(qJ(5));
t65 = cos(qJ(5));
t9 = t22 * t62 + t34 * t65;
t10 = t22 * t65 - t34 * t62;
t21 = t116 * t66 + t37 * t63;
t71 = t108 * t92 + t67 * t99;
t119 = t71 * t103 - t108 * t89;
t118 = t100 * t104 + t102 * t90;
t53 = t67 * t102 - t108 * t91;
t38 = t119 * t109 + t53 * t64;
t88 = t101 * t99;
t44 = -t118 * t109 + t64 * t88;
t78 = g(1) * t38 + g(2) * t34 + g(3) * t44;
t39 = t53 * t109 - t119 * t64;
t68 = -t71 * t100 - t108 * t90;
t23 = t39 * t63 + t68 * t66;
t45 = t109 * t88 + t118 * t64;
t51 = -t102 * t89 + t104 * t103;
t80 = -g(3) * (-t45 * t63 + t51 * t66) - g(2) * t21 + g(1) * t23;
t107 = t62 * t66;
t106 = t65 * t66;
t93 = t108 * t101;
t105 = t67 * pkin(1) + qJ(2) * t93;
t98 = t67 * t101;
t24 = t39 * t66 - t68 * t63;
t11 = t24 * t62 - t38 * t65;
t97 = g(1) * t9 + g(2) * t11;
t96 = -t108 * pkin(1) + qJ(2) * t98;
t95 = g(1) * t21 + g(2) * t23;
t33 = t45 * t66 + t51 * t63;
t17 = t33 * t62 - t44 * t65;
t1 = g(1) * t11 - g(2) * t9 + g(3) * t17;
t12 = t24 * t65 + t38 * t62;
t18 = t33 * t65 + t44 * t62;
t82 = g(1) * t12 - g(2) * t10 + g(3) * t18;
t13 = -t34 * t107 + t37 * t65;
t15 = -t38 * t107 - t39 * t65;
t25 = -t44 * t107 - t45 * t65;
t81 = g(1) * t15 + g(2) * t13 + g(3) * t25;
t79 = g(1) * t24 - g(2) * t22 + g(3) * t33;
t50 = -g(1) * t93 + g(2) * t98 - g(3) * t104;
t26 = -t44 * t106 + t45 * t62;
t16 = -t38 * t106 + t39 * t62;
t14 = -t34 * t106 - t37 * t62;
t6 = t78 * t63;
t5 = t80 * t65;
t4 = t80 * t62;
t3 = -g(1) * t10 - g(2) * t12;
t2 = -g(1) * t16 - g(2) * t14 - g(3) * t26;
t7 = [0, g(1) * t108 - g(2) * t67, g(1) * t67 + g(2) * t108, g(1) * t52 - g(2) * t53, -g(1) * t76 + g(2) * t71, -g(1) * t98 - g(2) * t93, -g(1) * t96 - g(2) * t105, 0, 0, 0, 0, 0, -g(1) * t37 - g(2) * t39, -g(1) * t34 + g(2) * t38, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t24, t95, 0, 0, 0, 0, 0, t3, t97, t3, -t95, -t97, -g(1) * (-t52 * pkin(2) + t37 * pkin(3) + t22 * pkin(4) + t10 * pkin(5) - pkin(10) * t34 + t21 * pkin(11) + t9 * qJ(6) + t96) - g(2) * (t53 * pkin(2) + t39 * pkin(3) + t24 * pkin(4) + t12 * pkin(5) + t38 * pkin(10) + t23 * pkin(11) + t11 * qJ(6) + t105) + (g(1) * t116 + g(2) * t68) * pkin(9); 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, g(1) * t39 - g(2) * t37 + g(3) * t45, 0, 0, 0, 0, 0, t78 * t66, -t6, 0, 0, 0, 0, 0, t2, t81, t2, t6, -t81, -g(1) * (t16 * pkin(5) + t39 * pkin(10) + t15 * qJ(6)) - g(2) * (t14 * pkin(5) - pkin(10) * t37 + t13 * qJ(6)) - g(3) * (t26 * pkin(5) + t45 * pkin(10) + t25 * qJ(6)) + t78 * (pkin(4) * t66 + pkin(11) * t63 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t79, 0, 0, 0, 0, 0, t5, -t4, t5, -t79, t4, -t79 * pkin(11) + t80 * (pkin(5) * t65 + qJ(6) * t62 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t82, t1, 0, -t82, -g(1) * (-t11 * pkin(5) + t12 * qJ(6)) - g(2) * (pkin(5) * t9 - qJ(6) * t10) - g(3) * (-t17 * pkin(5) + t18 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t7;
