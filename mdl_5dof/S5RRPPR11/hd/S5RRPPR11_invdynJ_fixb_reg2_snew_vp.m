% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPR11
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPR11_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:47:59
% EndTime: 2019-12-31 19:48:07
% DurationCPUTime: 2.54s
% Computational Cost: add. (9028->315), mult. (20094->410), div. (0->0), fcn. (12062->8), ass. (0->193)
t249 = -2 * qJD(4);
t166 = sin(qJ(2));
t159 = t166 ^ 2;
t172 = qJD(1) ^ 2;
t154 = t159 * t172;
t171 = qJD(2) ^ 2;
t142 = -t154 - t171;
t169 = cos(qJ(2));
t210 = t166 * t172;
t198 = t169 * t210;
t139 = qJDD(2) - t198;
t209 = t169 * t139;
t248 = pkin(6) * (t166 * t142 + t209);
t162 = sin(pkin(8));
t202 = qJD(1) * qJD(2);
t150 = t169 * t202;
t152 = t166 * qJDD(1);
t131 = t152 + t150;
t163 = cos(pkin(8));
t207 = qJD(1) * t169;
t121 = t162 * qJD(2) + t163 * t207;
t123 = t163 * qJD(2) - t162 * t207;
t95 = t123 * t121;
t242 = -t95 + t131;
t247 = t162 * t242;
t246 = t163 * t242;
t165 = sin(qJ(5));
t124 = qJDD(5) + t131;
t168 = cos(qJ(5));
t90 = t168 * t121 + t165 * t123;
t92 = -t165 * t121 + t168 * t123;
t70 = t92 * t90;
t243 = -t70 + t124;
t245 = t165 * t243;
t244 = t168 * t243;
t187 = t131 + t150;
t241 = t187 * qJ(3);
t153 = t169 * qJDD(1);
t195 = t166 * t202;
t132 = t153 - t195;
t104 = t163 * qJDD(2) - t162 * t132;
t204 = t166 * qJD(1);
t109 = t121 * t204;
t84 = t104 + t109;
t138 = pkin(3) * t204 - qJD(2) * qJ(4);
t146 = pkin(2) * t195;
t196 = qJD(3) * t204;
t149 = -0.2e1 * t196;
t160 = t169 ^ 2;
t167 = sin(qJ(1));
t170 = cos(qJ(1));
t194 = t167 * g(1) - t170 * g(2);
t181 = -qJDD(1) * pkin(1) - t194;
t231 = -pkin(2) - qJ(4);
t52 = -t138 * t204 + t146 + t149 + (-pkin(3) * t160 - pkin(6)) * t172 + t231 * t132 - t241 + t181;
t216 = t166 * qJ(3);
t233 = t169 * pkin(2);
t188 = -t216 - t233;
t128 = t188 * qJD(1);
t189 = t170 * g(1) + t167 * g(2);
t220 = qJDD(1) * pkin(6);
t116 = -t172 * pkin(1) - t189 + t220;
t215 = t166 * t116;
t180 = -qJDD(2) * pkin(2) - t171 * qJ(3) + t128 * t204 + qJDD(3) + t215;
t65 = t131 * pkin(3) - qJDD(2) * qJ(4) + (-pkin(3) * t202 - qJ(4) * t210 + g(3)) * t169 + t180;
t190 = t123 * t249 - t162 * t52 + t163 * t65;
t30 = t121 * t249 + t162 * t65 + t163 * t52;
t217 = t160 * t172;
t240 = t209 + t166 * (-t171 + t217);
t133 = t153 - 0.2e1 * t195;
t144 = -t171 - t217;
t137 = qJDD(2) + t198;
t213 = t166 * t137;
t238 = pkin(6) * (-t169 * t144 + t213) - pkin(1) * t133;
t88 = t90 ^ 2;
t89 = t92 ^ 2;
t119 = t121 ^ 2;
t120 = t123 ^ 2;
t147 = qJD(5) + t204;
t145 = t147 ^ 2;
t237 = 0.2e1 * qJD(3);
t173 = t242 * pkin(4) - t84 * pkin(7) + t190;
t103 = -t162 * qJDD(2) - t163 * t132;
t105 = pkin(4) * t204 - t123 * pkin(7);
t24 = -t119 * pkin(4) + t103 * pkin(7) - t105 * t204 + t30;
t10 = t165 * t173 + t168 * t24;
t9 = t165 * t24 - t168 * t173;
t6 = t165 * t10 - t168 * t9;
t235 = t162 * t6;
t234 = t163 * t6;
t156 = t166 * g(3);
t232 = t169 * g(3);
t177 = (qJD(1) * t128 + t116) * t169 - t156 - t171 * pkin(2);
t201 = qJDD(2) * qJ(3);
t61 = t201 + qJDD(4) + t132 * pkin(3) - qJ(4) * t217 + (t237 + t138) * qJD(2) + t177;
t229 = t162 * t61;
t86 = t95 + t131;
t228 = t162 * t86;
t227 = t163 * t61;
t226 = t163 * t86;
t33 = -t103 * pkin(4) - t119 * pkin(7) + t123 * t105 + t61;
t225 = t165 * t33;
t63 = t70 + t124;
t224 = t165 * t63;
t197 = t123 * t204;
t82 = t103 + t197;
t55 = t162 * t82 - t163 * t84;
t223 = t166 * t55;
t222 = t168 * t33;
t221 = t168 * t63;
t219 = t147 * t165;
t218 = t147 * t168;
t214 = t166 * t133;
t135 = t154 + t217;
t208 = pkin(1) * t135 + (t159 + t160) * t220;
t203 = qJD(5) + t147;
t200 = t166 * t70;
t199 = t166 * t95;
t7 = t168 * t10 + t165 * t9;
t100 = t215 + t232;
t101 = t169 * t116 - t156;
t192 = t166 * t100 + t169 * t101;
t191 = -t168 * t103 + t165 * t104;
t15 = t162 * t30 + t163 * t190;
t185 = -t162 * t190 + t163 * t30;
t184 = t165 * t103 + t168 * t104;
t183 = t169 * (-t154 + t171) + t213;
t179 = t169 * t231 - pkin(1) - t216;
t178 = (-qJD(5) + t147) * t92 - t191;
t58 = -t90 * qJD(5) + t184;
t115 = t172 * pkin(6) - t181;
t76 = t180 + t232;
t176 = qJD(2) * t237 + t177;
t175 = t132 * pkin(2) + t115 - t146;
t74 = t176 + t201;
t174 = t175 + 0.2e1 * t196;
t136 = t154 - t217;
t130 = t152 + 0.2e1 * t150;
t117 = t166 * t131;
t108 = -t120 - t154;
t107 = -t120 + t154;
t106 = t119 - t154;
t99 = t166 * t150 + t117;
t98 = (t132 - t195) * t169;
t94 = t169 * t130 + t214;
t93 = -t154 - t119;
t83 = t104 - t109;
t81 = -t103 + t197;
t80 = t147 * t90;
t79 = -t119 - t120;
t78 = -t89 + t145;
t77 = t88 - t145;
t75 = -t89 - t145;
t71 = t163 * t108 - t228;
t69 = t89 - t88;
t68 = -t145 - t88;
t66 = t162 * t93 + t246;
t57 = -t92 * qJD(5) - t191;
t54 = (t165 * t92 - t168 * t90) * t147;
t53 = (-t165 * t90 - t168 * t92) * t147;
t50 = -t88 - t89;
t49 = t58 + t80;
t48 = t58 - t80;
t47 = -t203 * t90 + t184;
t44 = t203 * t92 + t191;
t43 = t168 * t77 - t224;
t42 = -t165 * t78 + t244;
t41 = t165 * t77 + t221;
t40 = t168 * t78 + t245;
t39 = t168 * t58 - t92 * t219;
t38 = t165 * t58 + t92 * t218;
t37 = -t165 * t57 + t90 * t218;
t36 = t168 * t57 + t90 * t219;
t35 = -t165 * t75 - t221;
t34 = t168 * t75 - t224;
t32 = t168 * t68 - t245;
t31 = t165 * t68 + t244;
t28 = t165 * t49 + t168 * t178;
t27 = -t165 * t48 - t168 * t44;
t26 = t165 * t178 - t168 * t49;
t25 = -t165 * t44 + t168 * t48;
t21 = t162 * t35 + t163 * t34;
t20 = -pkin(7) * t34 + t222;
t19 = -pkin(7) * t31 + t225;
t17 = t162 * t32 + t163 * t31;
t14 = -pkin(4) * t47 + pkin(7) * t35 + t225;
t13 = -pkin(4) * t44 + pkin(7) * t32 - t222;
t11 = t162 * t28 + t163 * t26;
t5 = -pkin(4) * t33 + pkin(7) * t7;
t4 = -pkin(7) * t26 - t6;
t3 = -pkin(4) * t50 + pkin(7) * t28 + t7;
t1 = t162 * t7 + t234;
t2 = [0, 0, 0, 0, 0, qJDD(1), t194, t189, 0, 0, t99, t94, t183, t98, t240, 0, t169 * t115 - t238, -pkin(1) * t130 - t166 * t115 - t248, t192 + t208, pkin(1) * t115 + pkin(6) * t192, 0, -t183, -t240, t99, t94, t98, t166 * (qJ(3) * t135 + t180) + (pkin(2) * t135 + t156 + t74) * t169 + t208, t169 * (-pkin(2) * t133 + t149 - t175) + (-t169 * t187 - t214) * qJ(3) + t238, t166 * t174 + t248 + (pkin(1) + t233) * t130 + (t130 + t187) * t216, pkin(6) * (t166 * t76 + t169 * t74) + (pkin(1) - t188) * (t174 + t241), t199 + t169 * (-t162 * t104 - t163 * t197), t166 * (t120 - t119) + t169 * (t162 * t81 - t163 * t83), t166 * t84 + t169 * (-t163 * t107 - t247), -t199 + t169 * (-t163 * t103 - t162 * t109), t166 * t82 + t169 * (-t162 * t106 - t226), t117 + t169 * (t121 * t162 + t123 * t163) * t204, t166 * (pkin(3) * t66 + t190) + t169 * (pkin(3) * t81 + t227) + pkin(6) * (t166 * t66 + t169 * t81) + t179 * (t163 * t93 - t247), t166 * (pkin(3) * t71 - t30) + t169 * (pkin(3) * t83 - t229) + pkin(6) * (t166 * t71 + t169 * t83) + t179 * (-t162 * t108 - t226), pkin(3) * t223 + t169 * (pkin(3) * t79 - t185) + pkin(6) * (t169 * t79 + t223) + t179 * (t162 * t84 + t163 * t82), t179 * t185 + (pkin(3) + pkin(6)) * (t166 * t15 + t169 * t61), t200 + t169 * (-t162 * t39 - t163 * t38), t166 * t69 + t169 * (-t162 * t27 - t163 * t25), t166 * t49 + t169 * (-t162 * t42 - t163 * t40), -t200 + t169 * (-t162 * t37 - t163 * t36), t166 * t178 + t169 * (-t162 * t43 - t163 * t41), t166 * t124 + t169 * (-t162 * t54 - t163 * t53), t166 * (pkin(3) * t17 + pkin(4) * t31 - t9) + t169 * (pkin(3) * t44 - t163 * t13 - t162 * t19) + pkin(6) * (t166 * t17 + t169 * t44) + t179 * (-t162 * t31 + t163 * t32), t166 * (pkin(3) * t21 + pkin(4) * t34 - t10) + t169 * (pkin(3) * t47 - t163 * t14 - t162 * t20) + pkin(6) * (t166 * t21 + t169 * t47) + t179 * (-t162 * t34 + t163 * t35), t166 * (pkin(3) * t11 + pkin(4) * t26) + t169 * (pkin(3) * t50 - t162 * t4 - t163 * t3) + pkin(6) * (t166 * t11 + t169 * t50) + t179 * (-t162 * t26 + t163 * t28), t166 * (pkin(3) * t1 + pkin(4) * t6) + t169 * (pkin(3) * t33 + pkin(7) * t235 - t163 * t5) + pkin(6) * (t166 * t1 + t169 * t33) + t179 * (t163 * t7 - t235); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, t136, t152, t198, t153, qJDD(2), -t100, -t101, 0, 0, qJDD(2), -t152, -t153, -t198, t136, t198, (-pkin(2) * t166 + qJ(3) * t169) * qJDD(1), -pkin(2) * t137 - qJ(3) * t144 + t76, -pkin(2) * t142 + (qJDD(2) + t139) * qJ(3) + t176, -pkin(2) * t76 + qJ(3) * t74, t163 * t104 - t162 * t197, -t162 * t83 - t163 * t81, -t162 * t107 + t246, -t162 * t103 + t163 * t109, t163 * t106 - t228, (-t121 * t163 + t123 * t162) * t204, qJ(3) * t81 + t231 * t66 + t229, qJ(3) * t83 + t231 * t71 + t227, qJ(3) * t79 + t231 * t55 - t15, qJ(3) * t61 + t231 * t15, -t162 * t38 + t163 * t39, -t162 * t25 + t163 * t27, -t162 * t40 + t163 * t42, -t162 * t36 + t163 * t37, -t162 * t41 + t163 * t43, -t162 * t53 + t163 * t54, qJ(3) * t44 - t162 * t13 + t163 * t19 + t231 * t17, qJ(3) * t47 - t162 * t14 + t163 * t20 + t231 * t21, qJ(3) * t50 + t231 * t11 - t162 * t3 + t163 * t4, -pkin(7) * t234 + qJ(3) * t33 + t231 * t1 - t162 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, t137, t142, t76, 0, 0, 0, 0, 0, 0, t66, t71, t55, t15, 0, 0, 0, 0, 0, 0, t17, t21, t11, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t83, t79, t61, 0, 0, 0, 0, 0, 0, t44, t47, t50, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t69, t49, -t70, t178, t124, -t9, -t10, 0, 0;];
tauJ_reg = t2;
