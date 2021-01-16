% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:57
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:54:48
% EndTime: 2021-01-16 02:55:01
% DurationCPUTime: 2.94s
% Computational Cost: add. (4254->363), mult. (10950->501), div. (0->0), fcn. (8308->10), ass. (0->179)
t125 = sin(pkin(11));
t129 = sin(qJ(3));
t132 = cos(qJ(3));
t218 = qJ(4) + pkin(8);
t170 = qJD(3) * t218;
t143 = -t129 * qJD(4) - t132 * t170;
t204 = cos(pkin(11));
t94 = t132 * qJD(4) - t129 * t170;
t55 = t125 * t143 + t204 * t94;
t168 = t204 * t132;
t145 = -t125 * t129 + t168;
t133 = cos(qJ(2));
t126 = sin(pkin(6));
t194 = qJD(1) * t126;
t175 = t133 * t194;
t78 = t145 * t175;
t213 = t55 - t78;
t114 = qJD(2) * t168;
t192 = qJD(2) * t129;
t97 = t125 * t192 - t114;
t93 = qJD(5) + t97;
t106 = t125 * t132 + t204 * t129;
t214 = -t106 * t175 + t125 * t94 - t204 * t143;
t131 = cos(qJ(5));
t188 = qJD(5) * t131;
t128 = sin(qJ(5));
t99 = t106 * qJD(3);
t91 = qJD(2) * t99;
t86 = t128 * t91;
t235 = -t93 * t188 - t86;
t130 = sin(qJ(2));
t176 = t130 * t194;
t189 = qJD(5) * t128;
t102 = t145 * qJD(3);
t187 = t129 * qJD(3);
t184 = pkin(3) * t187;
t57 = t99 * pkin(4) - t102 * pkin(9) + t184;
t121 = -t132 * pkin(3) - pkin(2);
t66 = -pkin(4) * t145 - t106 * pkin(9) + t121;
t110 = t218 * t129;
t111 = t218 * t132;
t80 = -t125 * t110 + t204 * t111;
t234 = -t66 * t188 + t80 * t189 - t213 * t131 + (t176 - t57) * t128;
t87 = t131 * t91;
t233 = -t93 * t189 + t87;
t108 = qJD(2) * pkin(8) + t176;
t127 = cos(pkin(6));
t193 = qJD(1) * t127;
t174 = t129 * t193;
t149 = -t132 * t108 - t174;
t156 = qJD(4) + t175;
t177 = qJD(3) * t132 * qJ(4);
t232 = t149 * qJD(3) + (-t156 * t129 - t177) * qJD(2);
t100 = t106 * qJD(2);
t231 = qJD(3) * t100;
t165 = qJ(4) * qJD(2) + t108;
t76 = t165 * t132 + t174;
t67 = t125 * t76;
t115 = t132 * t193;
t75 = -t165 * t129 + t115;
t71 = qJD(3) * pkin(3) + t75;
t30 = t204 * t71 - t67;
t27 = -qJD(3) * pkin(4) - t30;
t186 = t131 * qJD(3);
t81 = t128 * t100 - t186;
t83 = t128 * qJD(3) + t131 * t100;
t15 = t81 * pkin(5) - t83 * qJ(6) + t27;
t118 = t125 * pkin(3) + pkin(9);
t209 = t118 * t91;
t230 = t93 * t15 - t209;
t185 = qJD(2) * qJD(3);
t171 = t129 * t185;
t112 = t125 * t171;
t144 = qJD(3) * t114 - t112;
t50 = t83 * qJD(5) + t128 * t144;
t229 = t83 ^ 2;
t228 = t93 ^ 2;
t172 = t204 * t76;
t31 = t125 * t71 + t172;
t28 = qJD(3) * pkin(9) + t31;
t92 = t121 * qJD(2) + qJD(4) - t175;
t44 = t97 * pkin(4) - t100 * pkin(9) + t92;
t13 = t128 * t44 + t131 * t28;
t7 = t93 * qJ(6) + t13;
t227 = t7 * t93;
t226 = t91 * pkin(5);
t225 = t99 * qJ(6) - qJD(6) * t145 - t234;
t212 = t128 * t66 + t131 * t80;
t58 = t128 * t78 - t131 * t176;
t224 = -t99 * pkin(5) + t212 * qJD(5) + t128 * t55 - t131 * t57 - t58;
t157 = pkin(5) * t128 - qJ(6) * t131;
t158 = t131 * pkin(5) + t128 * qJ(6);
t223 = -t157 * t102 - (t158 * qJD(5) - qJD(6) * t131) * t106 - t214;
t222 = pkin(3) * t129;
t221 = t13 * t93;
t220 = t81 * t97;
t219 = t83 * t81;
t169 = t83 * t93;
t33 = t125 * t75 + t172;
t217 = t128 * qJD(6) - t93 * t157 + t33;
t35 = t204 * t75 - t67;
t183 = pkin(3) * t192;
t56 = t100 * pkin(4) + t97 * pkin(9) + t183;
t216 = t128 * t56 + t131 * t35;
t215 = -t128 * t50 - t81 * t188;
t211 = qJD(2) * pkin(2);
t210 = t100 * t81;
t208 = t128 * t93;
t49 = -qJD(5) * t186 + t100 * t189 - t131 * t144;
t207 = t49 * t128;
t206 = t83 * t100;
t205 = t91 * qJ(6);
t203 = t102 * t128;
t202 = t102 * t131;
t201 = t126 * t130;
t200 = t126 * t133;
t135 = qJD(2) ^ 2;
t199 = t126 * t135;
t134 = qJD(3) ^ 2;
t198 = t134 * t129;
t197 = t134 * t132;
t12 = -t128 * t28 + t131 * t44;
t196 = qJD(6) - t12;
t96 = pkin(3) * t171 + qJD(2) * t176;
t195 = t129 ^ 2 - t132 ^ 2;
t191 = qJD(2) * t130;
t190 = qJD(5) * t118;
t182 = t93 * t190;
t180 = t130 * t199;
t179 = t126 * t191;
t178 = qJD(2) * t200;
t48 = (-t129 * t108 + t115) * qJD(3) + (-qJ(4) * t187 + t156 * t132) * qJD(2);
t173 = t204 * t48;
t18 = t125 * t48 - t204 * t232;
t19 = t232 * t125 + t173;
t41 = t91 * pkin(4) - t144 * pkin(9) + t96;
t167 = t128 * t19 - t131 * t41 + t28 * t188 + t44 * t189;
t79 = t204 * t110 + t125 * t111;
t164 = t129 * t178;
t163 = t132 * t178;
t3 = t50 * pkin(5) + t49 * qJ(6) - t83 * qJD(6) + t18;
t162 = -t3 - t182;
t119 = -t204 * pkin(3) - pkin(4);
t6 = -t93 * pkin(5) + t196;
t161 = -t128 * t7 + t131 * t6;
t159 = t18 * t106 - t80 * t91;
t155 = t131 * t97 * t93 - t235;
t154 = -t97 * t208 + t233;
t153 = t15 * t83 + t167;
t103 = t127 * t129 + t132 * t201;
t152 = t127 * t132 - t129 * t201;
t62 = t204 * t103 + t125 * t152;
t46 = t128 * t62 + t131 * t200;
t47 = -t128 * t200 + t131 * t62;
t151 = t128 * t41 + t131 * t19 + t44 * t188 - t28 * t189;
t142 = t103 * qJD(3);
t136 = -t142 - t164;
t74 = t152 * qJD(3) + t163;
t34 = t125 * t136 + t204 * t74;
t10 = t47 * qJD(5) + t128 * t34 - t131 * t179;
t32 = t125 * t74 - t204 * t136;
t61 = t125 * t103 - t204 * t152;
t148 = -t10 * t93 + t32 * t81 - t46 * t91 + t61 * t50;
t147 = t211 * qJD(2);
t146 = t93 * t27 - t209;
t140 = -0.2e1 * qJD(3) * t211;
t9 = -t46 * qJD(5) + t128 * t179 + t131 * t34;
t138 = t32 * t83 - t47 * t91 - t61 * t49 - t9 * t93;
t104 = -t158 + t119;
t43 = t83 * pkin(5) + t81 * qJ(6);
t36 = t157 * t106 + t79;
t25 = pkin(5) * t145 + t128 * t80 - t131 * t66;
t24 = -qJ(6) * t145 + t212;
t21 = t81 * t93 - t49;
t14 = -t100 * pkin(5) + t128 * t35 - t131 * t56;
t11 = t100 * qJ(6) + t216;
t2 = t167 - t226;
t1 = t93 * qJD(6) + t151 + t205;
t4 = [0, 0, -t180, -t133 * t199, 0, 0, 0, 0, 0, -t132 * t180 + (-t142 - 0.2e1 * t164) * qJD(3), t129 * t180 + (-t74 - t163) * qJD(3), -t32 * qJD(3) + (-t133 * t91 + t97 * t191) * t126, -t34 * qJD(3) + (t100 * t191 - t133 * t144) * t126, t32 * t100 + t144 * t61 - t34 * t97 - t62 * t91, t18 * t61 + t19 * t62 - t30 * t32 + t31 * t34 + (-t133 * t96 + t92 * t191) * t126, 0, 0, 0, 0, 0, t148, t138, t148, t10 * t83 - t46 * t49 - t47 * t50 - t9 * t81, -t138, t1 * t47 + t6 * t10 + t15 * t32 + t2 * t46 + t3 * t61 + t7 * t9; 0, 0, 0, 0, 0.2e1 * t132 * t171, -0.2e1 * t195 * t185, t197, -t198, 0, -pkin(8) * t197 + t129 * t140, pkin(8) * t198 + t132 * t140, -t97 * t176 - t96 * t145 + t121 * t91 + t92 * t99 + (t97 * t222 - t214) * qJD(3), -t100 * t176 + t92 * t102 + t96 * t106 - t121 * t112 + (t100 * t222 + t114 * t121 - t213) * qJD(3), t214 * t100 - t30 * t102 + t79 * t144 + t145 * t19 - t213 * t97 - t31 * t99 + t159, t96 * t121 + t18 * t79 + t19 * t80 + (-t176 + t184) * t92 + t213 * t31 - t214 * t30, t83 * t202 + (-t49 * t131 - t83 * t189) * t106, (-t128 * t83 - t131 * t81) * t102 + (t207 - t131 * t50 + (t128 * t81 - t131 * t83) * qJD(5)) * t106, t233 * t106 + t145 * t49 + t93 * t202 + t83 * t99, t235 * t106 + t145 * t50 - t93 * t203 - t81 * t99, -t145 * t91 + t93 * t99, t167 * t145 + t12 * t99 + t79 * t50 + t58 * t93 + t214 * t81 + ((-qJD(5) * t80 + t57) * t93 + t66 * t91 + t27 * qJD(5) * t106) * t131 + ((-qJD(5) * t66 - t55) * t93 + t27 * t102 + t159) * t128, -t212 * t91 + t151 * t145 - t13 * t99 - t79 * t49 + t27 * t202 + t234 * t93 + t214 * t83 + (t18 * t131 - t27 * t189) * t106, t15 * t203 + t2 * t145 - t25 * t91 + t36 * t50 - t6 * t99 - t224 * t93 - t223 * t81 + (t3 * t128 + t15 * t188) * t106, -t24 * t50 - t25 * t49 + t224 * t83 - t225 * t81 + t161 * t102 + (-t1 * t128 + t2 * t131 + (-t128 * t6 - t131 * t7) * qJD(5)) * t106, -t15 * t202 - t1 * t145 + t24 * t91 + t36 * t49 + t7 * t99 + t225 * t93 + t223 * t83 + (-t3 * t131 + t15 * t189) * t106, t1 * t24 - t223 * t15 + t2 * t25 + t224 * t6 + t225 * t7 + t3 * t36; 0, 0, 0, 0, -t129 * t135 * t132, t195 * t135, 0, 0, 0, t129 * t147, t132 * t147, t33 * qJD(3) - t92 * t100 - t97 * t183 - t18, -t173 + t92 * t97 + (-t125 * t149 + t35) * qJD(3) + (t125 * t177 + (-pkin(3) * t100 + t125 * t156) * t129) * qJD(2), (-t30 + t35) * t97 + (-t33 + t31) * t100 + (-t125 * t91 - t204 * t144) * pkin(3), t30 * t33 - t31 * t35 + (t125 * t19 - t204 * t18 - t92 * t192) * pkin(3), t131 * t169 - t207, (-t49 - t220) * t131 - t83 * t208 + t215, t155 - t206, t154 + t210, -t93 * t100, -t12 * t100 + t119 * t50 - t33 * t81 + (-t18 + (-t56 - t190) * t93) * t131 + (t35 * t93 + t146) * t128, -t119 * t49 + t216 * t93 + t13 * t100 - t33 * t83 + (t18 + t182) * t128 + t146 * t131, t6 * t100 + t104 * t50 + t230 * t128 + t162 * t131 + t14 * t93 - t217 * t81, t11 * t81 - t14 * t83 + (-t118 * t50 + t6 * t97 + t1 + (t118 * t83 + t6) * qJD(5)) * t131 + (-t118 * t49 - t7 * t97 + t2 + (t118 * t81 - t7) * qJD(5)) * t128, -t7 * t100 + t104 * t49 - t11 * t93 + t162 * t128 - t230 * t131 + t217 * t83, t3 * t104 - t7 * t11 - t6 * t14 - t217 * t15 + (qJD(5) * t161 + t1 * t131 + t2 * t128) * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t231, -t112 + (t114 - t97) * qJD(3), -t100 ^ 2 - t97 ^ 2, t30 * t100 + t31 * t97 + t96, 0, 0, 0, 0, 0, t154 - t210, -t228 * t131 - t206 - t86, -t208 * t93 - t210 + t87, (t49 - t220) * t131 + t128 * t169 + t215, t155 + t206, -t15 * t100 + (-t2 + t227) * t131 + (t93 * t6 + t1) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219, -t81 ^ 2 + t229, t21, t169 - t50, t91, -t27 * t83 - t167 + t221, t12 * t93 + t27 * t81 - t151, -t43 * t81 - t153 + t221 + 0.2e1 * t226, pkin(5) * t49 - t50 * qJ(6) + (-t13 + t7) * t83 + (t6 - t196) * t81, 0.2e1 * t205 - t15 * t81 + t43 * t83 + (0.2e1 * qJD(6) - t12) * t93 + t151, -t2 * pkin(5) + t1 * qJ(6) - t6 * t13 - t15 * t43 + t196 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219 - t231, t21, -t228 - t229, t153 - t226 - t227;];
tauc_reg = t4;
