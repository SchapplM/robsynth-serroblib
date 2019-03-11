% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:13
% EndTime: 2019-03-09 03:06:17
% DurationCPUTime: 1.71s
% Computational Cost: add. (3883->304), mult. (9141->397), div. (0->0), fcn. (6285->8), ass. (0->152)
t106 = sin(pkin(9)) * pkin(1) + pkin(7);
t173 = qJ(4) + t106;
t120 = cos(qJ(3));
t178 = cos(pkin(10));
t154 = t178 * t120;
t102 = qJD(1) * t154;
t114 = sin(pkin(10));
t118 = sin(qJ(3));
t167 = t118 * qJD(1);
t89 = -t114 * t167 + t102;
t86 = qJD(5) - t89;
t130 = -t114 * t118 + t154;
t117 = sin(qJ(5));
t164 = qJD(1) * qJD(3);
t156 = t118 * t164;
t128 = qJD(3) * t102 - t114 * t156;
t119 = cos(qJ(5));
t98 = t114 * t120 + t178 * t118;
t91 = t98 * qJD(1);
t74 = t117 * qJD(3) + t119 * t91;
t45 = t74 * qJD(5) + t117 * t128;
t165 = t119 * qJD(3);
t72 = t117 * t91 - t165;
t90 = t98 * qJD(3);
t206 = t130 * t45 - t90 * t72;
t169 = qJD(5) * t117;
t44 = -qJD(5) * t165 - t119 * t128 + t91 * t169;
t186 = t130 * t44 + t74 * t90;
t151 = t173 * qJD(1);
t76 = t120 * qJD(2) - t151 * t118;
t163 = qJD(1) * qJD(4);
t77 = t118 * qJD(2) + t151 * t120;
t205 = -t77 * qJD(3) - t118 * t163;
t67 = t114 * t77;
t184 = qJD(3) * pkin(3);
t70 = t76 + t184;
t30 = t178 * t70 - t67;
t27 = -qJD(3) * pkin(4) - t30;
t13 = t72 * pkin(5) - t74 * qJ(6) + t27;
t105 = t114 * pkin(3) + pkin(8);
t84 = qJD(1) * t90;
t183 = t105 * t84;
t204 = t86 * t13 - t183;
t203 = t74 ^ 2;
t202 = t86 ^ 2;
t158 = t178 * t77;
t31 = t114 * t70 + t158;
t28 = qJD(3) * pkin(8) + t31;
t108 = -cos(pkin(9)) * pkin(1) - pkin(2);
t138 = -t120 * pkin(3) + t108;
t129 = t138 * qJD(1);
t88 = qJD(4) + t129;
t41 = -t89 * pkin(4) - t91 * pkin(8) + t88;
t12 = t117 * t41 + t119 * t28;
t8 = t86 * qJ(6) + t12;
t201 = t8 * t86;
t200 = t84 * pkin(5);
t199 = t12 * t86;
t64 = t76 * qJD(3) + t120 * t163;
t20 = t114 * t64 - t178 * t205;
t198 = t20 * t98;
t93 = t130 * qJD(3);
t197 = t27 * t93;
t196 = t72 * t89;
t195 = t74 * t72;
t155 = t74 * t86;
t194 = t74 * t91;
t193 = t74 * t93;
t192 = t91 * t72;
t191 = t98 * t84;
t141 = pkin(5) * t117 - qJ(6) * t119;
t32 = t114 * t76 + t158;
t190 = t117 * qJD(6) - t86 * t141 + t32;
t180 = t119 * t98;
t181 = t119 * t93;
t189 = -t45 * t180 - t72 * t181;
t33 = t178 * t76 - t67;
t52 = pkin(3) * t167 + t91 * pkin(4) - t89 * pkin(8);
t188 = t117 * t52 + t119 * t33;
t168 = qJD(5) * t119;
t187 = -t117 * t45 - t72 * t168;
t57 = -pkin(4) * t130 - t98 * pkin(8) + t138;
t95 = t173 * t118;
t96 = t173 * t120;
t60 = -t114 * t95 + t178 * t96;
t185 = t117 * t57 + t119 * t60;
t81 = t117 * t84;
t182 = t117 * t86;
t82 = t119 * t84;
t179 = t84 * qJ(6);
t177 = qJD(5) * t72;
t176 = qJD(5) * t98;
t121 = qJD(3) ^ 2;
t175 = t121 * t118;
t174 = t121 * t120;
t11 = -t117 * t28 + t119 * t41;
t172 = qJD(6) - t11;
t171 = t118 ^ 2 - t120 ^ 2;
t100 = qJD(1) * t108;
t170 = qJD(5) * t105;
t162 = t118 * t184;
t161 = t86 * t170;
t160 = t98 * t169;
t159 = t98 * t168;
t150 = qJD(3) * t173;
t78 = t120 * qJD(4) - t118 * t150;
t79 = -t118 * qJD(4) - t120 * t150;
t39 = t114 * t78 - t178 * t79;
t59 = t114 * t96 + t178 * t95;
t21 = t114 * t205 + t178 * t64;
t103 = pkin(3) * t156;
t46 = t84 * pkin(4) - t128 * pkin(8) + t103;
t153 = t117 * t21 - t119 * t46 + t28 * t168 + t41 * t169;
t149 = t86 * t160;
t148 = t74 * t159;
t4 = t45 * pkin(5) + t44 * qJ(6) - t74 * qJD(6) + t20;
t147 = -t4 - t161;
t107 = -t178 * pkin(3) - pkin(4);
t146 = t13 * t93 + t4 * t98;
t7 = -t86 * pkin(5) + t172;
t145 = -t117 * t8 + t119 * t7;
t144 = t117 * t7 + t119 * t8;
t143 = -t60 * t84 + t198;
t142 = -t86 * t93 - t191;
t140 = t81 + (-t119 * t89 + t168) * t86;
t137 = -t86 * t169 + t89 * t182 + t82;
t136 = 0.2e1 * qJD(3) * t100;
t135 = t13 * t74 + t153;
t134 = t84 * t180 + t86 * t181 - t149;
t133 = t117 * t46 + t119 * t21 + t41 * t168 - t28 * t169;
t40 = t114 * t79 + t178 * t78;
t53 = t90 * pkin(4) - t93 * pkin(8) + t162;
t132 = t117 * t53 + t119 * t40 + t57 * t168 - t60 * t169;
t131 = t86 * t27 - t183;
t125 = t142 * t117 - t86 * t159;
t1 = t86 * qJD(6) + t133 + t179;
t2 = t153 - t200;
t124 = t145 * qJD(5) + t1 * t119 + t2 * t117;
t123 = t125 - t206;
t122 = qJD(1) ^ 2;
t94 = -t119 * pkin(5) - t117 * qJ(6) + t107;
t37 = t74 * pkin(5) + t72 * qJ(6);
t22 = t141 * t98 + t59;
t17 = t72 * t86 - t44;
t15 = pkin(5) * t130 + t117 * t60 - t119 * t57;
t14 = -qJ(6) * t130 + t185;
t10 = -t91 * pkin(5) + t117 * t33 - t119 * t52;
t9 = t91 * qJ(6) + t188;
t6 = (pkin(5) * t93 + qJ(6) * t176) * t117 + (-qJ(6) * t93 + (pkin(5) * qJD(5) - qJD(6)) * t98) * t119 + t39;
t5 = -t90 * pkin(5) + t185 * qJD(5) + t117 * t40 - t119 * t53;
t3 = t90 * qJ(6) - qJD(6) * t130 + t132;
t16 = [0, 0, 0, 0, 0.2e1 * t120 * t156, -0.2e1 * t171 * t164, t174, -t175, 0, -t106 * t174 + t118 * t136, t106 * t175 + t120 * t136, t128 * t59 + t130 * t21 - t30 * t93 - t31 * t90 + t39 * t91 + t40 * t89 + t143, t20 * t59 + t21 * t60 - t30 * t39 + t31 * t40 + (t88 + t129) * t162, -t74 * t160 + (-t44 * t98 + t193) * t119, -t148 + (-t193 + (t44 + t177) * t98) * t117 + t189, t134 + t186, t125 + t206, -t130 * t84 + t86 * t90, t153 * t130 + t11 * t90 + t39 * t72 + t59 * t45 + ((-qJD(5) * t60 + t53) * t86 + t57 * t84 + t27 * t176) * t119 + ((-qJD(5) * t57 - t40) * t86 + t197 + t143) * t117, -t132 * t86 - t185 * t84 + t133 * t130 - t12 * t90 + t39 * t74 - t59 * t44 - t27 * t160 + (t197 + t198) * t119, t117 * t146 + t13 * t159 + t130 * t2 - t15 * t84 + t22 * t45 - t5 * t86 + t6 * t72 - t7 * t90, -t14 * t45 - t15 * t44 - t3 * t72 + t5 * t74 + t145 * t93 + (-qJD(5) * t144 - t1 * t117 + t2 * t119) * t98, -t1 * t130 - t119 * t146 + t13 * t160 + t14 * t84 + t22 * t44 + t3 * t86 - t6 * t74 + t8 * t90, t1 * t14 + t13 * t6 + t2 * t15 + t4 * t22 + t8 * t3 + t7 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, -t174, -t128 * t130 + t93 * t89 + t90 * t91 - t191, -t130 * t20 + t21 * t98 - t30 * t90 + t31 * t93, 0, 0, 0, 0, 0, t123, t119 * t142 + t149 + t186, t123, t148 + (t193 + (-t44 + t177) * t98) * t117 + t189, t134 - t186, t124 * t98 + t13 * t90 - t130 * t4 + t144 * t93; 0, 0, 0, 0, -t118 * t122 * t120, t171 * t122, 0, 0, 0, -t100 * t167, -t100 * t120 * qJD(1) (t31 - t32) * t91 + (-t33 + t30) * t89 + (-t114 * t84 - t178 * t128) * pkin(3), t30 * t32 - t31 * t33 + (t114 * t21 - t88 * t167 - t178 * t20) * pkin(3), -t44 * t117 + t119 * t155 (-t44 + t196) * t119 - t74 * t182 + t187, t140 - t194, t137 + t192, -t86 * t91, t107 * t45 - t11 * t91 - t32 * t72 + (-t20 + (-t52 - t170) * t86) * t119 + (t33 * t86 + t131) * t117, -t107 * t44 + t188 * t86 + t12 * t91 - t32 * t74 + (t20 + t161) * t117 + t131 * t119, t10 * t86 + t117 * t204 + t147 * t119 - t190 * t72 + t94 * t45 + t7 * t91, -t10 * t74 + t9 * t72 + (-t105 * t45 - t7 * t89 + t1 + (t105 * t74 + t7) * qJD(5)) * t119 + (-t105 * t44 + t8 * t89 + t2 + (t105 * t72 - t8) * qJD(5)) * t117, t147 * t117 - t119 * t204 + t190 * t74 + t94 * t44 - t8 * t91 - t9 * t86, -t7 * t10 + t124 * t105 - t190 * t13 + t4 * t94 - t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89 ^ 2 - t91 ^ 2, t30 * t91 - t31 * t89 + t103, 0, 0, 0, 0, 0, t137 - t192, -t202 * t119 - t194 - t81, -t182 * t86 - t192 + t82 (t44 + t196) * t119 + t117 * t155 + t187, t140 + t194, -t13 * t91 + (-t2 + t201) * t119 + (t7 * t86 + t1) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, -t72 ^ 2 + t203, t17, t155 - t45, t84, -t27 * t74 - t153 + t199, t11 * t86 + t27 * t72 - t133, -t37 * t72 - t135 + t199 + 0.2e1 * t200, pkin(5) * t44 - t45 * qJ(6) + (-t12 + t8) * t74 + (t7 - t172) * t72, 0.2e1 * t179 - t13 * t72 + t37 * t74 + (0.2e1 * qJD(6) - t11) * t86 + t133, -t2 * pkin(5) + t1 * qJ(6) - t7 * t12 - t13 * t37 + t172 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t91 + t195, t17, -t202 - t203, t135 - t200 - t201;];
tauc_reg  = t16;
