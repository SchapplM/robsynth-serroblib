% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:06:52
% EndTime: 2019-12-05 18:06:58
% DurationCPUTime: 0.97s
% Computational Cost: add. (1563->175), mult. (4157->264), div. (0->0), fcn. (2876->6), ass. (0->130)
t131 = qJD(3) + qJD(4);
t88 = sin(pkin(8));
t84 = t88 ^ 2;
t168 = 0.2e1 * t84;
t90 = sin(qJ(4));
t91 = sin(qJ(3));
t92 = cos(qJ(4));
t93 = cos(qJ(3));
t63 = t90 * t93 + t92 * t91;
t54 = t63 * t88;
t167 = t91 * t93;
t89 = cos(pkin(8));
t85 = t89 ^ 2;
t166 = t168 + t85;
t66 = -t89 * pkin(2) - t88 * pkin(6) - pkin(1);
t57 = t66 * qJD(1) + qJD(2);
t53 = t93 * t57;
t144 = qJ(2) * t91;
t161 = pkin(7) * t88;
t98 = -t89 * t144 - t93 * t161;
t38 = t98 * qJD(1) + t53;
t135 = t89 * qJD(1);
t77 = -qJD(3) + t135;
t24 = -t77 * pkin(3) + t38;
t141 = qJD(1) * t88;
t127 = t91 * t141;
t143 = qJ(2) * t93;
t128 = t89 * t143;
t154 = t91 * t57;
t39 = -pkin(7) * t127 + qJD(1) * t128 + t154;
t31 = t90 * t39;
t116 = t92 * t24 - t31;
t153 = t92 * t93;
t129 = t88 * t153;
t110 = qJD(1) * t129;
t111 = t90 * t127;
t50 = t110 - t111;
t44 = t50 * qJ(5);
t165 = t44 - t116;
t121 = -t66 + t161;
t164 = t121 * t91 - t128;
t163 = t50 ^ 2;
t72 = -qJD(4) + t77;
t5 = -t72 * pkin(4) - t165;
t162 = t5 + t165;
t99 = qJD(1) * t63;
t47 = t88 * t99;
t58 = pkin(3) * t127 + qJ(2) * t141;
t36 = t47 * pkin(4) + qJD(5) + t58;
t160 = t36 * t50;
t159 = t50 * t47;
t158 = t58 * t50;
t94 = qJD(1) ^ 2;
t157 = t84 * t94;
t156 = t88 * t91;
t155 = t90 * t91;
t33 = t92 * t39;
t152 = t92 * t38 - t31;
t62 = t153 - t155;
t151 = (t131 - t135) * t62;
t42 = t131 * t63;
t150 = t89 * t99 - t42;
t132 = qJD(1) * qJD(2);
t119 = t89 * t132;
t138 = qJD(3) * t93;
t149 = t93 * t119 + t57 * t138;
t139 = qJD(2) * t93;
t148 = t66 * t138 + t89 * t139;
t147 = t131 * t111;
t75 = t88 * pkin(3) * t138;
t56 = qJD(1) * t75 + t88 * t132;
t60 = t88 * qJD(2) + t75;
t61 = pkin(3) * t156 + t88 * qJ(2);
t146 = t84 + t85;
t145 = t91 ^ 2 - t93 ^ 2;
t142 = t47 * qJ(5);
t140 = qJD(2) * t91;
t137 = qJD(4) * t90;
t136 = qJD(4) * t92;
t134 = qJD(3) + t77;
t133 = qJ(2) * qJD(3);
t126 = t36 * t141;
t125 = t47 * t141;
t124 = t50 * t141;
t123 = t89 * t140;
t122 = t91 * t133;
t120 = t146 * t94;
t118 = qJD(1) * qJD(3) * t84;
t97 = t98 * qJD(3);
t19 = qJD(1) * t97 + t149;
t20 = -qJD(3) * t154 + (-t123 + (pkin(7) * t156 - t128) * qJD(3)) * qJD(1);
t117 = -t90 * t19 + t92 * t20;
t34 = t97 + t148;
t35 = t164 * qJD(3) - t123;
t115 = -t90 * t34 + t92 * t35;
t114 = -t90 * t38 - t33;
t113 = -t24 * t136 + t39 * t137 - t92 * t19 - t90 * t20;
t112 = qJD(1) * t134;
t104 = t131 * t129;
t23 = qJD(1) * t104 - t147;
t109 = t23 * pkin(4) + t56;
t108 = t88 * t112;
t107 = -t90 * t24 - t33;
t40 = -t121 * t93 + (-pkin(3) - t144) * t89;
t106 = t164 * t92 - t90 * t40;
t105 = 0.2e1 * t146 * t132;
t103 = t58 * t47 + t113;
t102 = qJD(3) * t88 * (t77 + t135);
t101 = t40 * t136 + t137 * t164 + t92 * t34 + t90 * t35;
t96 = -t77 ^ 2 - t157;
t95 = t107 * qJD(4) + t117;
t29 = t131 * t54;
t81 = t92 * pkin(3) + pkin(4);
t55 = t62 * t88;
t45 = t47 ^ 2;
t30 = -t131 * t88 * t155 + t104;
t22 = qJD(1) * t29;
t14 = -t45 + t163;
t13 = -t131 * t110 - t50 * t72 + t147;
t12 = -t42 * t141 - t47 * t72;
t11 = -t54 * qJ(5) - t106;
t10 = -t89 * pkin(4) - t55 * qJ(5) + t164 * t90 + t92 * t40;
t9 = -t44 + t152;
t8 = t114 + t142;
t7 = -t107 - t142;
t4 = t29 * qJ(5) + t106 * qJD(4) - t55 * qJD(5) + t115;
t3 = -t30 * qJ(5) - t54 * qJD(5) + t101;
t2 = t22 * qJ(5) - t50 * qJD(5) + t95;
t1 = -t23 * qJ(5) - t47 * qJD(5) - t113;
t6 = [0, 0, 0, 0, 0, t105, qJ(2) * t105, -0.2e1 * t118 * t167, 0.2e1 * t145 * t118, t91 * t102, t93 * t102, 0, t77 * t123 + (-(-t91 * t66 - t128) * t77 + t89 * t154) * qJD(3) + t166 * qJD(1) * (t93 * t133 + t140), (-t89 * t122 + t148) * t77 + t149 * t89 + (-t166 * t122 + t139 * t168) * qJD(1), -t22 * t55 - t50 * t29, t22 * t54 - t55 * t23 + t29 * t47 - t50 * t30, t22 * t89 + t29 * t72, t23 * t89 + t30 * t72, 0, -t115 * t72 - t117 * t89 + t60 * t47 + t61 * t23 + t56 * t54 + t58 * t30 + (-t106 * t72 - t107 * t89) * qJD(4), t101 * t72 - t113 * t89 - t61 * t22 - t58 * t29 + t60 * t50 + t56 * t55, -t1 * t54 + t10 * t22 - t11 * t23 - t2 * t55 + t5 * t29 - t3 * t47 - t7 * t30 - t4 * t50, t1 * t11 + t7 * t3 + t2 * t10 + t5 * t4 + t109 * (t54 * pkin(4) + t61) + t36 * (t30 * pkin(4) + t60); 0, 0, 0, 0, 0, -t120, -qJ(2) * t120, 0, 0, 0, 0, 0, t96 * t91, t96 * t93, 0, 0, 0, 0, 0, -t150 * t72 - t125, t151 * t72 - t124, -t150 * t50 - t151 * t47 + t62 * t22 - t63 * t23, t1 * t63 + t150 * t5 + t151 * t7 + t2 * t62 - t126; 0, 0, 0, 0, 0, 0, 0, t157 * t167, -t145 * t157, -t91 * t108, -t93 * t108, 0, (-t134 * t57 - t119) * t91 + (-t89 * t112 - t157) * t143, -t53 * t77 + (t134 * t135 + t157) * t144 - t149, t159, t14, t12, t13, 0, t114 * t72 - t93 * pkin(3) * t125 - t158 + (-t33 + (pkin(3) * t72 - t24) * t90) * qJD(4) + t117, -t152 * t72 + (-t93 * t124 + t72 * t136) * pkin(3) + t103, t81 * t22 + (t7 + t8) * t50 + (-t5 + t9) * t47 + (-t23 * t90 + (-t47 * t92 + t50 * t90) * qJD(4)) * pkin(3), -pkin(4) * t160 + t2 * t81 - t5 * t8 - t7 * t9 + (-t93 * t126 + t1 * t90 + (-t5 * t90 + t7 * t92) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t14, t12, t13, 0, t107 * t72 - t158 + t95, -t116 * t72 + t103, pkin(4) * t22 - t162 * t47, t162 * t7 + (t2 - t160) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45 - t163, t7 * t47 + t5 * t50 + t109;];
tauc_reg = t6;
