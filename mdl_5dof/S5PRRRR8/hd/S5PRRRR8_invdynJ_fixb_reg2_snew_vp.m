% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:44
% EndTime: 2019-12-05 17:16:54
% DurationCPUTime: 2.60s
% Computational Cost: add. (8947->319), mult. (17846->460), div. (0->0), fcn. (13069->12), ass. (0->195)
t172 = sin(qJ(5));
t174 = sin(qJ(3));
t177 = cos(qJ(4));
t173 = sin(qJ(4));
t178 = cos(qJ(3));
t208 = t178 * t173;
t137 = (t174 * t177 + t208) * qJD(2);
t203 = qJD(2) * qJD(3);
t195 = t178 * t203;
t202 = t174 * qJDD(2);
t141 = t195 + t202;
t161 = t178 * qJDD(2);
t196 = t174 * t203;
t142 = t161 - t196;
t190 = t173 * t141 - t177 * t142;
t101 = -t137 * qJD(4) - t190;
t100 = qJDD(5) - t101;
t165 = qJD(3) + qJD(4);
t176 = cos(qJ(5));
t119 = t172 * t137 - t176 * t165;
t121 = t176 * t137 + t172 * t165;
t99 = t121 * t119;
t232 = t100 - t99;
t237 = t172 * t232;
t206 = qJD(2) * t174;
t135 = -t177 * t178 * qJD(2) + t173 * t206;
t116 = t137 * t135;
t164 = qJDD(3) + qJDD(4);
t231 = -t116 + t164;
t236 = t173 * t231;
t235 = t176 * t232;
t234 = t177 * t231;
t169 = sin(pkin(5));
t170 = cos(pkin(5));
t219 = sin(pkin(10));
t220 = cos(pkin(10));
t187 = t219 * g(1) - t220 * g(2);
t207 = -g(3) + qJDD(1);
t233 = t169 * t207 + t170 * t187;
t114 = t135 * pkin(4) - t137 * pkin(9);
t230 = t165 ^ 2;
t125 = -t169 * t187 + t170 * t207;
t146 = -t220 * g(1) - t219 * g(2);
t175 = sin(qJ(2));
t179 = cos(qJ(2));
t109 = t179 * t146 + t233 * t175;
t181 = qJD(2) ^ 2;
t184 = -t181 * pkin(2) + qJDD(2) * pkin(7) + t109;
t183 = t174 * t184;
t211 = t174 * t181;
t218 = qJDD(3) * pkin(3);
t150 = qJD(3) * pkin(3) - pkin(8) * t206;
t167 = t178 ^ 2;
t163 = t167 * t181;
t87 = t174 * t125 + t178 * t184;
t72 = -pkin(3) * t163 + t142 * pkin(8) - qJD(3) * t150 + t87;
t222 = t177 * t72;
t227 = t141 * pkin(8);
t43 = t222 + t173 * (-t183 + t218 - t227) + (pkin(3) * t211 + pkin(8) * t203 + t125) * t208;
t33 = -t230 * pkin(4) + t164 * pkin(9) - t135 * t114 + t43;
t189 = t175 * t146 - t233 * t179;
t104 = -qJDD(2) * pkin(2) - t181 * pkin(7) + t189;
t84 = -t142 * pkin(3) - pkin(8) * t163 + t150 * t206 + t104;
t102 = -t135 * qJD(4) + t177 * t141 + t173 * t142;
t129 = t165 * t135;
t93 = t102 - t129;
t41 = -t93 * pkin(9) + (t137 * t165 - t101) * pkin(4) + t84;
t16 = t172 * t33 - t176 * t41;
t17 = t172 * t41 + t176 * t33;
t7 = t172 * t16 + t176 * t17;
t131 = qJD(5) + t135;
t191 = t172 * t102 - t176 * t164;
t65 = (qJD(5) - t131) * t121 + t191;
t117 = t119 ^ 2;
t118 = t121 ^ 2;
t130 = t131 ^ 2;
t133 = t135 ^ 2;
t134 = t137 ^ 2;
t153 = t178 * t211;
t86 = -t178 * t125 + t183;
t182 = pkin(8) * t195 - t227 - t86;
t42 = t173 * t72 - t177 * (pkin(3) * t153 + t182 + t218);
t32 = -t164 * pkin(4) - t230 * pkin(9) + t137 * t114 + t42;
t229 = -pkin(4) * t32 + pkin(9) * t7;
t228 = pkin(4) * t173;
t29 = t172 * t32;
t74 = t100 + t99;
t226 = t172 * t74;
t225 = t173 * t84;
t20 = t173 * t43 - t177 * t42;
t224 = t174 * t20;
t30 = t176 * t32;
t223 = t176 * t74;
t221 = t177 * t84;
t217 = t131 * t172;
t216 = t131 * t176;
t215 = t165 * t173;
t214 = t165 * t177;
t112 = t116 + t164;
t213 = t173 * t112;
t147 = qJDD(3) + t153;
t212 = t174 * t147;
t210 = t177 * t112;
t148 = qJDD(3) - t153;
t209 = t178 * t148;
t204 = qJD(5) + t131;
t97 = -t118 - t130;
t49 = -t172 * t97 - t223;
t188 = -t176 * t102 - t172 * t164;
t70 = t119 * t204 + t188;
t201 = pkin(4) * t70 + pkin(9) * t49 + t29;
t88 = -t130 - t117;
t46 = t176 * t88 - t237;
t67 = -t121 * t204 - t191;
t200 = pkin(4) * t67 + pkin(9) * t46 - t30;
t199 = t173 * t99;
t198 = t177 * t99;
t197 = -pkin(4) * t177 - pkin(3);
t21 = t173 * t42 + t177 * t43;
t55 = t174 * t86 + t178 * t87;
t107 = t131 * t119;
t79 = -t119 * qJD(5) - t188;
t69 = t107 + t79;
t37 = t172 * t69 - t176 * t65;
t81 = t117 + t118;
t194 = pkin(4) * t81 + pkin(9) * t37 + t7;
t2 = t173 * t7 - t177 * t32;
t3 = t173 * t32 + t177 * t7;
t1 = -t174 * t2 + t178 * t3;
t6 = -t176 * t16 + t172 * t17;
t186 = (-qJD(4) + t165) * t137 - t190;
t180 = qJD(3) ^ 2;
t166 = t174 ^ 2;
t162 = t166 * t181;
t152 = -t163 - t180;
t151 = -t162 - t180;
t145 = t162 + t163;
t144 = (t166 + t167) * qJDD(2);
t143 = t161 - 0.2e1 * t196;
t140 = 0.2e1 * t195 + t202;
t127 = -t134 + t230;
t126 = t133 - t230;
t124 = -t134 - t230;
t123 = -t174 * t151 - t209;
t122 = t178 * t152 - t212;
t115 = t134 - t133;
t110 = -t230 - t133;
t106 = -t118 + t130;
t105 = t117 - t130;
t103 = -t133 - t134;
t98 = t118 - t117;
t96 = -t173 * t124 - t210;
t95 = t177 * t124 - t213;
t94 = t102 + t129;
t89 = (qJD(4) + t165) * t137 + t190;
t83 = t177 * t110 - t236;
t82 = t173 * t110 + t234;
t78 = -t121 * qJD(5) - t191;
t77 = (-t119 * t176 + t121 * t172) * t131;
t76 = (-t119 * t172 - t121 * t176) * t131;
t68 = -t107 + t79;
t62 = -t121 * t217 + t176 * t79;
t61 = t121 * t216 + t172 * t79;
t60 = t119 * t216 - t172 * t78;
t59 = t119 * t217 + t176 * t78;
t58 = -t174 * t95 + t178 * t96;
t57 = t173 * t94 + t177 * t186;
t56 = t173 * t186 - t177 * t94;
t54 = t176 * t105 - t226;
t53 = -t172 * t106 + t235;
t52 = t172 * t105 + t223;
t51 = t176 * t106 + t237;
t50 = -t174 * t82 + t178 * t83;
t48 = t176 * t97 - t226;
t45 = t172 * t88 + t235;
t38 = -t172 * t68 + t176 * t67;
t36 = t172 * t67 + t176 * t68;
t35 = -t172 * t65 - t176 * t69;
t28 = -t174 * t56 + t178 * t57;
t27 = -t173 * t70 + t177 * t49;
t26 = t173 * t49 + t177 * t70;
t25 = -t173 * t67 + t177 * t46;
t24 = t173 * t46 + t177 * t67;
t23 = -t173 * t81 + t177 * t37;
t22 = t173 * t37 + t177 * t81;
t19 = -pkin(9) * t48 + t30;
t18 = -pkin(9) * t45 + t29;
t13 = -t174 * t26 + t178 * t27;
t12 = -t174 * t24 + t178 * t25;
t11 = -pkin(4) * t48 + t17;
t10 = -pkin(4) * t45 + t16;
t9 = -t174 * t22 + t178 * t23;
t8 = t178 * t21 - t224;
t4 = -pkin(9) * t35 - t6;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t207, 0, 0, 0, 0, 0, 0, (qJDD(2) * t179 - t175 * t181) * t169, (-qJDD(2) * t175 - t179 * t181) * t169, 0, t170 * t125 + (t109 * t175 - t179 * t189) * t169, 0, 0, 0, 0, 0, 0, t170 * (t178 * t147 + t174 * t152) + (t175 * t122 + t179 * t143) * t169, t170 * (-t174 * t148 + t178 * t151) + (t175 * t123 - t179 * t140) * t169, (t144 * t175 + t145 * t179) * t169, t170 * (t174 * t87 - t178 * t86) + (-t179 * t104 + t175 * t55) * t169, 0, 0, 0, 0, 0, 0, t170 * (t174 * t83 + t178 * t82) + (t175 * t50 - t179 * t89) * t169, t170 * (t174 * t96 + t178 * t95) + (t175 * t58 - t179 * t93) * t169, t170 * (t174 * t57 + t178 * t56) + (-t179 * t103 + t175 * t28) * t169, t170 * (t174 * t21 + t178 * t20) + (t175 * t8 - t179 * t84) * t169, 0, 0, 0, 0, 0, 0, t170 * (t174 * t25 + t178 * t24) + (t175 * t12 - t179 * t45) * t169, t170 * (t174 * t27 + t178 * t26) + (t175 * t13 - t179 * t48) * t169, t170 * (t174 * t23 + t178 * t22) + (t175 * t9 - t179 * t35) * t169, t170 * (t174 * t3 + t178 * t2) + (t175 * t1 - t179 * t6) * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t189, -t109, 0, 0, (t141 + t195) * t174, t178 * t140 + t174 * t143, t212 + t178 * (-t162 + t180), (t142 - t196) * t178, t174 * (t163 - t180) + t209, 0, pkin(2) * t143 + pkin(7) * t122 - t178 * t104, -pkin(2) * t140 + pkin(7) * t123 + t174 * t104, pkin(2) * t145 + pkin(7) * t144 + t55, -pkin(2) * t104 + pkin(7) * t55, t174 * (t177 * t102 - t137 * t215) + t178 * (t173 * t102 + t137 * t214), t174 * (-t173 * t93 - t177 * t89) + t178 * (-t173 * t89 + t177 * t93), t174 * (-t173 * t127 + t234) + t178 * (t177 * t127 + t236), t174 * (-t173 * t101 + t135 * t214) + t178 * (t177 * t101 + t135 * t215), t174 * (t177 * t126 - t213) + t178 * (t173 * t126 + t210), (t174 * (-t135 * t177 + t137 * t173) + t178 * (-t135 * t173 - t137 * t177)) * t165, t174 * (-pkin(8) * t82 + t225) + t178 * (-pkin(3) * t89 + pkin(8) * t83 - t221) - pkin(2) * t89 + pkin(7) * t50, t174 * (-pkin(8) * t95 + t221) + t178 * (-pkin(3) * t93 + pkin(8) * t96 + t225) - pkin(2) * t93 + pkin(7) * t58, t174 * (-pkin(8) * t56 - t20) + t178 * (-pkin(3) * t103 + pkin(8) * t57 + t21) - pkin(2) * t103 + pkin(7) * t28, -pkin(8) * t224 + t178 * (-pkin(3) * t84 + pkin(8) * t21) - pkin(2) * t84 + pkin(7) * t8, t174 * (t177 * t62 + t199) + t178 * (t173 * t62 - t198), t174 * (t173 * t98 + t177 * t38) + t178 * (t173 * t38 - t177 * t98), t174 * (t173 * t69 + t177 * t53) + t178 * (t173 * t53 - t177 * t69), t174 * (t177 * t60 - t199) + t178 * (t173 * t60 + t198), t174 * (-t173 * t65 + t177 * t54) + t178 * (t173 * t54 + t177 * t65), t174 * (t173 * t100 + t177 * t77) + t178 * (-t177 * t100 + t173 * t77), t174 * (-pkin(8) * t24 - t173 * t10 + t177 * t18) + t178 * (-pkin(3) * t45 + pkin(8) * t25 + t177 * t10 + t173 * t18) - pkin(2) * t45 + pkin(7) * t12, t174 * (-pkin(8) * t26 - t173 * t11 + t177 * t19) + t178 * (-pkin(3) * t48 + pkin(8) * t27 + t177 * t11 + t173 * t19) - pkin(2) * t48 + pkin(7) * t13, t174 * (-pkin(8) * t22 + t177 * t4) + t178 * (pkin(8) * t23 + t173 * t4) + pkin(7) * t9 + (t174 * t228 + t178 * t197 - pkin(2)) * t35, (t174 * (-pkin(9) * t177 + t228) + t178 * (-pkin(9) * t173 + t197) - pkin(2)) * t6 + (pkin(7) + pkin(8)) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t162 - t163, t202, t153, t161, qJDD(3), -t86, -t87, 0, 0, t116, t115, t94, -t116, t186, t164, pkin(3) * t82 - t42, -t222 - t173 * t182 + (-t147 * t173 + t95) * pkin(3), pkin(3) * t56, pkin(3) * t20, t61, t36, t51, t59, t52, t76, pkin(3) * t24 + t200, pkin(3) * t26 + t201, pkin(3) * t22 + t194, pkin(3) * t2 + t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t115, t94, -t116, t186, t164, -t42, -t43, 0, 0, t61, t36, t51, t59, t52, t76, t200, t201, t194, t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t98, t69, -t99, -t65, t100, -t16, -t17, 0, 0;];
tauJ_reg = t5;
