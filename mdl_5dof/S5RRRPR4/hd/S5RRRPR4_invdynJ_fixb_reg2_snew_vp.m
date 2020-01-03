% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:39
% EndTime: 2019-12-31 21:11:44
% DurationCPUTime: 1.59s
% Computational Cost: add. (5816->226), mult. (7607->276), div. (0->0), fcn. (4282->8), ass. (0->149)
t157 = cos(qJ(3));
t147 = qJD(1) + qJD(2);
t184 = qJD(3) * t147;
t179 = t157 * t184;
t144 = qJDD(1) + qJDD(2);
t154 = sin(qJ(3));
t191 = t154 * t144;
t115 = 0.2e1 * t179 + t191;
t155 = sin(qJ(2));
t158 = cos(qJ(2));
t159 = qJD(3) ^ 2;
t142 = t147 ^ 2;
t150 = t154 ^ 2;
t197 = t150 * t142;
t128 = t159 + t197;
t202 = t142 * t157;
t133 = t154 * t202;
t127 = qJDD(3) - t133;
t189 = t157 * t127;
t93 = -t154 * t128 + t189;
t225 = pkin(1) * (t158 * t115 + t155 * t93);
t153 = sin(qJ(5));
t143 = qJDD(3) - qJDD(5);
t156 = cos(qJ(5));
t103 = (-t154 * t153 - t157 * t156) * t147;
t198 = t147 * t157;
t199 = t147 * t154;
t105 = -t153 * t198 + t156 * t199;
t203 = t105 * t103;
t218 = -t143 + t203;
t224 = t153 * t218;
t223 = t156 * t218;
t222 = pkin(2) * t115 + pkin(7) * t93;
t145 = qJD(3) - qJD(5);
t204 = t103 * t145;
t116 = t179 + t191;
t137 = t154 * t184;
t188 = t157 * t144;
t173 = -t137 + t188;
t60 = t103 * qJD(5) + t156 * t116 - t153 * t173;
t220 = t60 - t204;
t219 = t116 + t179;
t212 = t157 * g(3);
t214 = sin(qJ(1));
t215 = cos(qJ(1));
t164 = t215 * g(1) + t214 * g(2);
t125 = -qJD(1) ^ 2 * pkin(1) - t164;
t163 = t214 * g(1) - t215 * g(2);
t161 = qJDD(1) * pkin(1) + t163;
t83 = t158 * t125 + t155 * t161;
t76 = -t142 * pkin(2) + t144 * pkin(7) + t83;
t65 = t154 * t76 + t212;
t66 = -t154 * g(3) + t157 * t76;
t39 = t154 * t65 + t157 * t66;
t151 = t157 ^ 2;
t196 = t151 * t142;
t89 = t189 + t154 * (-t159 + t196);
t174 = t153 * t116 + t156 * t173;
t49 = (qJD(5) + t145) * t105 + t174;
t169 = -qJDD(3) * pkin(3) - t159 * qJ(4) + qJDD(4) + t212;
t170 = -t157 * pkin(3) - t154 * qJ(4);
t114 = t170 * t147;
t175 = t114 * t147 + t76;
t55 = t175 * t154 + t169;
t101 = t103 ^ 2;
t102 = t105 ^ 2;
t141 = t145 ^ 2;
t217 = 2 * qJD(4);
t216 = pkin(3) + pkin(4);
t166 = qJDD(3) * qJ(4) + (qJD(3) * t217) + t114 * t198 + t66;
t54 = -t159 * pkin(3) + t166;
t25 = t154 * t55 + t157 * t54;
t82 = -t155 * t125 + t158 * t161;
t75 = -t144 * pkin(2) - t142 * pkin(7) - t82;
t162 = -t173 * pkin(3) - t219 * qJ(4) + t75;
t41 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t199 + t162;
t211 = -pkin(2) * t41 + pkin(7) * t25;
t210 = -pkin(2) * t75 + pkin(7) * t39;
t168 = -qJD(3) * pkin(4) - pkin(8) * t199;
t31 = -t173 * pkin(4) + pkin(8) * t196 - t168 * t199 + t41;
t209 = t153 * t31;
t71 = t203 + t143;
t208 = t153 * t71;
t207 = t156 * t31;
t206 = t156 * t71;
t117 = -0.2e1 * t137 + t188;
t130 = -t159 - t196;
t126 = qJDD(3) + t133;
t194 = t154 * t126;
t92 = t157 * t130 - t194;
t205 = pkin(2) * t117 + pkin(7) * t92;
t201 = t145 * t153;
t200 = t145 * t156;
t195 = t154 * t117;
t190 = t157 * t115;
t185 = t150 + t151;
t120 = t185 * t144;
t123 = t185 * t142;
t186 = pkin(2) * t123 + pkin(7) * t120;
t182 = t154 * t75 - t222;
t181 = -t157 * t75 + t205;
t35 = -pkin(4) * t196 - t173 * pkin(8) + qJD(3) * t168 + t54;
t36 = -qJDD(3) * pkin(4) + (-t116 + t179) * pkin(8) + (-pkin(4) * t202 + t175) * t154 + t169;
t16 = t153 * t35 - t156 * t36;
t17 = t153 * t36 + t156 * t35;
t7 = t153 * t17 - t156 * t16;
t8 = t153 * t16 + t156 * t17;
t2 = t154 * t7 + t157 * t8;
t180 = pkin(7) * t2 - pkin(2) * t31 + t157 * (-pkin(8) * t8 - t216 * t31) + t154 * (-pkin(8) * t7 - qJ(4) * t31);
t53 = t204 + t60;
t22 = -t153 * t49 - t156 * t53;
t23 = t153 * t53 - t156 * t49;
t13 = t154 * t22 + t157 * t23;
t61 = -t101 - t102;
t178 = t154 * (-pkin(8) * t22 + qJ(4) * t61 - t7) + t157 * (-pkin(8) * t23 + t216 * t61 - t8) + pkin(2) * t61 + pkin(7) * t13;
t70 = -t141 - t101;
t44 = t153 * t70 + t223;
t45 = t156 * t70 - t224;
t19 = t154 * t44 + t157 * t45;
t48 = (qJD(5) - t145) * t105 + t174;
t177 = t154 * (-pkin(8) * t44 + qJ(4) * t48 - t209) + t157 * (-pkin(8) * t45 + t216 * t48 - t207) + pkin(2) * t48 + pkin(7) * t19;
t95 = -t102 - t141;
t56 = t156 * t95 + t208;
t57 = -t153 * t95 + t206;
t27 = t154 * t56 + t157 * t57;
t176 = t154 * (-pkin(8) * t56 + qJ(4) * t220 - t207) + t157 * (-pkin(8) * t57 + t216 * t220 + t209) + pkin(2) * t220 + pkin(7) * t27;
t172 = t154 * (qJ(4) * t123 + t55) + t157 * ((t123 - t159) * pkin(3) + t166) + t186;
t171 = t186 + t39;
t79 = t190 + t195;
t160 = t199 * t217 - t162;
t167 = pkin(3) * t190 + t154 * (-pkin(3) * t137 + qJ(4) * t115 + t160) + t222;
t165 = qJ(4) * t195 + t205 + t157 * ((t117 - t137) * pkin(3) + t160);
t124 = (t150 - t151) * t142;
t97 = -t102 + t141;
t96 = t101 - t141;
t91 = t194 + t157 * (t159 - t197);
t85 = t219 * t154;
t84 = (t173 - t137) * t157;
t81 = pkin(1) * (t155 * t120 + t158 * t123);
t77 = t102 - t101;
t64 = pkin(1) * (t158 * t117 + t155 * t92);
t59 = -t105 * qJD(5) - t174;
t37 = (t154 * (-t103 * t156 - t105 * t153) + t157 * (t103 * t153 - t105 * t156)) * t145;
t29 = t154 * (t156 * t96 + t208) + t157 * (-t153 * t96 + t206);
t28 = t154 * (-t153 * t97 + t223) + t157 * (-t156 * t97 - t224);
t21 = t154 * (t105 * t201 + t156 * t60) + t157 * (t105 * t200 - t153 * t60);
t20 = t154 * (t103 * t200 - t153 * t59) + t157 * (-t103 * t201 - t156 * t59);
t12 = t154 * (-t153 * t220 - t156 * t48) + t157 * (t153 * t48 - t156 * t220);
t1 = [0, 0, 0, 0, 0, qJDD(1), t163, t164, 0, 0, 0, 0, 0, 0, 0, t144, pkin(1) * (-t155 * t142 + t158 * t144) + t82, pkin(1) * (-t158 * t142 - t155 * t144) - t83, 0, pkin(1) * (t155 * t83 + t158 * t82), t85, t79, t91, t84, t89, 0, t64 + t181, t182 - t225, t81 + t171, pkin(1) * (t155 * t39 - t158 * t75) + t210, t85, t91, -t79, 0, -t89, t84, t165 + t64, t81 + t172, t167 + t225, pkin(1) * t155 * t25 + (-pkin(1) * t158 + t170) * t41 + t211, t21, t12, t28, t20, t29, t37, pkin(1) * (t155 * t19 + t158 * t48) + t177, pkin(1) * (t155 * t27 + t158 * t220) + t176, pkin(1) * (t155 * t13 + t158 * t61) + t178, pkin(1) * (t155 * t2 - t158 * t31) + t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, t82, -t83, 0, 0, t85, t79, t91, t84, t89, 0, t181, t182, t171, t210, t85, t91, -t79, 0, -t89, t84, t165, t172, t167, t170 * t41 + t211, t21, t12, t28, t20, t29, t37, t177, t176, t178, t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, t124, t191, t133, t188, qJDD(3), -t65, -t66, 0, 0, -t133, t191, -t124, qJDD(3), -t188, t133, pkin(3) * t126 + qJ(4) * t130 - t55, (-pkin(3) * t154 + qJ(4) * t157) * t144, qJ(4) * t127 + (t128 - t159) * pkin(3) + t166, -pkin(3) * t55 + qJ(4) * t54, t203, -t77, -t53, -t203, t49, t143, qJ(4) * t45 - t216 * t44 + t16, qJ(4) * t57 - t216 * t56 + t17, qJ(4) * t23 - t216 * t22, qJ(4) * t8 - t216 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t191, -t128, t55, 0, 0, 0, 0, 0, 0, t44, t56, t22, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t203, t77, t53, t203, -t49, -t143, -t16, -t17, 0, 0;];
tauJ_reg = t1;
