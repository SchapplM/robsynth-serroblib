% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:44
% EndTime: 2019-12-05 15:19:51
% DurationCPUTime: 1.87s
% Computational Cost: add. (2473->247), mult. (6914->397), div. (0->0), fcn. (6025->12), ass. (0->137)
t69 = sin(pkin(6));
t76 = sin(qJ(3));
t139 = t69 * t76;
t73 = cos(pkin(5));
t61 = qJD(1) * t73 + qJD(2);
t116 = t61 * t139;
t70 = sin(pkin(5));
t71 = cos(pkin(11));
t72 = cos(pkin(6));
t136 = t71 * t72;
t68 = sin(pkin(11));
t79 = cos(qJ(3));
t91 = t76 * t136 + t68 * t79;
t85 = t91 * t70;
t35 = qJD(1) * t85 + t116;
t75 = sin(qJ(4));
t78 = cos(qJ(4));
t99 = pkin(4) * t75 - pkin(9) * t78;
t57 = t99 * qJD(4);
t24 = (t57 + t35) * qJD(3);
t138 = t69 * t79;
t140 = t68 * t76;
t160 = (t79 * t136 - t140) * t70;
t162 = qJD(1) * t160 + t61 * t138;
t30 = t162 * qJD(3);
t33 = qJD(3) * pkin(8) + t35;
t129 = qJD(1) * t70;
t111 = t71 * t129;
t46 = -t69 * t111 + t61 * t72;
t95 = t33 * t75 - t46 * t78;
t7 = -t95 * qJD(4) + t78 * t30;
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t22 = t33 * t78 + t46 * t75;
t20 = qJD(4) * pkin(9) + t22;
t34 = -t129 * t140 + t79 * (t72 * t111 + t61 * t69);
t59 = -pkin(4) * t78 - pkin(9) * t75 - pkin(3);
t27 = t59 * qJD(3) - t34;
t96 = t20 * t74 - t27 * t77;
t1 = -t96 * qJD(5) + t74 * t24 + t77 * t7;
t126 = qJD(3) * t78;
t62 = -qJD(5) + t126;
t165 = -t96 * t62 + t1;
t122 = qJD(5) * t74;
t107 = t75 * t122;
t119 = t77 * qJD(4);
t164 = -t78 * t119 + t107;
t163 = t73 * t138 + t160;
t161 = pkin(8) * qJD(5) * t78 + t35 - t57;
t6 = t20 * t77 + t27 * t74;
t2 = -qJD(5) * t6 + t77 * t24 - t74 * t7;
t159 = t6 * t62 - t2;
t117 = qJD(4) * qJD(5);
t121 = qJD(5) * t77;
t123 = qJD(4) * t78;
t84 = t75 * t121 + t74 * t123;
t45 = qJD(3) * t84 + t74 * t117;
t39 = t73 * t139 + t85;
t50 = -t69 * t70 * t71 + t72 * t73;
t25 = t39 * t75 - t50 * t78;
t8 = qJD(4) * t22 + t75 * t30;
t158 = t25 * t8;
t51 = t75 * t139 - t78 * t72;
t156 = t51 * t8;
t154 = t74 * t8;
t153 = t75 * t8;
t152 = t77 * t8;
t19 = -qJD(4) * pkin(4) + t95;
t151 = t19 * t74;
t150 = t19 * t77;
t31 = qJD(3) * t35;
t149 = t31 * t163;
t148 = t31 * t79;
t44 = t164 * qJD(3) - t77 * t117;
t147 = t44 * t74;
t146 = t45 * t77;
t128 = qJD(3) * t75;
t53 = t74 * t128 - t119;
t145 = t53 * t62;
t125 = qJD(4) * t74;
t55 = t77 * t128 + t125;
t144 = t55 * t53;
t143 = t55 * t62;
t142 = t62 * t74;
t141 = t62 * t77;
t81 = qJD(3) ^ 2;
t137 = t69 * t81;
t135 = t74 * t78;
t134 = t77 * t78;
t124 = qJD(4) * t75;
t133 = -t74 * t124 * pkin(8) + t59 * t122 - t34 * t135 + t161 * t77;
t132 = t75 * t119 * pkin(8) - t59 * t121 + t34 * t134 + t161 * t74;
t66 = t75 ^ 2;
t67 = t78 ^ 2;
t131 = t66 - t67;
t130 = qJD(3) * pkin(3);
t127 = qJD(3) * t76;
t118 = qJD(3) * qJD(4);
t114 = t76 * t137;
t112 = t75 * t81 * t78;
t110 = t69 * t127;
t109 = qJD(3) * t138;
t106 = t62 * t121;
t104 = t75 * t118;
t32 = -t34 - t130;
t103 = -qJD(3) * t32 - t30;
t102 = t78 * t109;
t101 = t75 * t109;
t100 = t78 * t104;
t98 = -t6 * t74 + t77 * t96;
t26 = t39 * t78 + t50 * t75;
t12 = -t163 * t74 + t26 * t77;
t11 = -t163 * t77 - t26 * t74;
t94 = qJD(3) * t66 - t62 * t78;
t52 = t78 * t139 + t72 * t75;
t42 = -t77 * t138 - t52 * t74;
t93 = t74 * t138 - t52 * t77;
t80 = qJD(4) ^ 2;
t90 = pkin(8) * t80;
t89 = qJD(4) * (t32 + t34 - t130);
t82 = t7 * t78 + t153 + (-t22 * t75 + t78 * t95) * qJD(4);
t56 = t99 * qJD(3);
t48 = pkin(8) * t134 + t59 * t74;
t47 = -pkin(8) * t135 + t59 * t77;
t41 = t52 * qJD(4) + t101;
t40 = -t51 * qJD(4) + t102;
t37 = t39 * qJD(3);
t36 = t163 * qJD(3);
t18 = t93 * qJD(5) + t77 * t110 - t74 * t40;
t17 = t42 * qJD(5) + t74 * t110 + t77 * t40;
t14 = t56 * t74 - t77 * t95;
t13 = t56 * t77 + t74 * t95;
t10 = -t25 * qJD(4) + t36 * t78;
t9 = t26 * qJD(4) + t36 * t75;
t4 = t11 * qJD(5) + t10 * t77 + t37 * t74;
t3 = -t12 * qJD(5) - t10 * t74 + t37 * t77;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * qJD(3), -t36 * qJD(3), 0, t30 * t39 - t34 * t37 + t35 * t36 - t149, 0, 0, 0, 0, 0, 0, -t9 * qJD(4) + (-t124 * t163 - t37 * t78) * qJD(3), -t10 * qJD(4) + (-t123 * t163 + t37 * t75) * qJD(3), (t10 * t78 + t75 * t9 + (t25 * t78 - t26 * t75) * qJD(4)) * qJD(3), t10 * t22 + t26 * t7 + t32 * t37 + t9 * t95 - t149 + t158, 0, 0, 0, 0, 0, 0, t11 * t104 + t25 * t45 - t3 * t62 + t53 * t9, -t104 * t12 - t25 * t44 + t4 * t62 + t55 * t9, t11 * t44 - t12 * t45 - t3 * t55 - t4 * t53, t1 * t12 + t11 * t2 + t19 * t9 - t3 * t96 + t4 * t6 + t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, -t79 * t137, 0, (t30 * t76 - t148 + (-t34 * t76 + t35 * t79) * qJD(3)) * t69, 0, 0, 0, 0, 0, 0, -t78 * t114 + (-t41 - t101) * qJD(4), t75 * t114 + (-t40 - t102) * qJD(4), (t40 * t78 + t41 * t75 + (t51 * t78 - t52 * t75) * qJD(4)) * qJD(3), t95 * t41 + t22 * t40 + t156 + t52 * t7 + (t32 * t127 - t148) * t69, 0, 0, 0, 0, 0, 0, t42 * t104 - t18 * t62 + t41 * t53 + t45 * t51, t104 * t93 + t17 * t62 + t41 * t55 - t44 * t51, -t17 * t53 - t18 * t55 + t42 * t44 + t45 * t93, -t1 * t93 + t17 * t6 - t18 * t96 + t19 * t41 + t2 * t42 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t91 * t129 - t116 + t35) * qJD(3), (-t162 + t34) * qJD(3), 0, 0, 0.2e1 * t100, -0.2e1 * t131 * t118, t80 * t78, -0.2e1 * t100, -t80 * t75, 0, t75 * t89 - t78 * t90, t75 * t90 + t78 * t89, (-t66 - t67) * t34 * qJD(3) + t82, -pkin(3) * t31 - t32 * t35 + (-t22 * t78 - t75 * t95) * t34 + t82 * pkin(8), -t44 * t77 * t75 - t164 * t55, (-t53 * t77 - t55 * t74) * t123 + (t147 - t146 + (t53 * t74 - t55 * t77) * qJD(5)) * t75, t62 * t107 + t44 * t78 + (t55 * t75 + t94 * t77) * qJD(4), t45 * t74 * t75 + t53 * t84, t75 * t106 + t45 * t78 + (-t53 * t75 - t94 * t74) * qJD(4), (-t62 - t126) * t124, t133 * t62 + (-t2 + (pkin(8) * t53 + t151) * qJD(4)) * t78 + (t19 * t121 + pkin(8) * t45 - t34 * t53 + t154 + (qJD(3) * t47 - t96) * qJD(4)) * t75, -t132 * t62 + (t1 + (pkin(8) * t55 + t150) * qJD(4)) * t78 + (-t19 * t122 - pkin(8) * t44 - t34 * t55 + t152 + (-qJD(3) * t48 - t6) * qJD(4)) * t75, t44 * t47 - t45 * t48 + t133 * t55 + t132 * t53 + t98 * t123 + (-t1 * t74 - t2 * t77 + (-t6 * t77 - t74 * t96) * qJD(5)) * t75, -t19 * t34 * t75 + t1 * t48 + t2 * t47 - t132 * t6 + t133 * t96 + (t123 * t19 + t153) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, t131 * t81, 0, t112, 0, 0, t103 * t75, t103 * t78, 0, 0, -t55 * t141 - t147, (-t44 + t145) * t77 + (-t45 + t143) * t74, -t106 + (t62 * t134 + (-t55 + t125) * t75) * qJD(3), -t53 * t142 - t146, t62 * t122 + (-t62 * t135 + (t53 + t119) * t75) * qJD(3), t62 * t128, -pkin(4) * t45 + t13 * t62 - t22 * t53 - t152 + (pkin(9) * t141 + t151) * qJD(5) + (t96 * t75 + (-pkin(9) * t124 - t19 * t78) * t74) * qJD(3), pkin(4) * t44 - t14 * t62 - t22 * t55 + t154 + (-pkin(9) * t142 + t150) * qJD(5) + (-t19 * t134 + (-pkin(9) * t119 + t6) * t75) * qJD(3), t13 * t55 + t14 * t53 + ((qJD(5) * t55 - t45) * pkin(9) + t165) * t77 + ((qJD(5) * t53 - t44) * pkin(9) + t159) * t74, -pkin(4) * t8 + t13 * t96 - t14 * t6 - t19 * t22 + (qJD(5) * t98 + t1 * t77 - t2 * t74) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, -t53 ^ 2 + t55 ^ 2, -t44 - t145, -t144, -t143 - t45, t104, -t19 * t55 - t159, t19 * t53 - t165, 0, 0;];
tauc_reg = t5;
