% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:26
% EndTime: 2022-01-20 11:17:32
% DurationCPUTime: 2.41s
% Computational Cost: add. (12687->254), mult. (17686->362), div. (0->0), fcn. (11446->10), ass. (0->177)
t175 = sin(qJ(5));
t167 = qJDD(1) + qJDD(2);
t174 = cos(pkin(9));
t220 = t174 * t167;
t154 = -qJDD(4) + t220;
t146 = -qJDD(5) + t154;
t176 = sin(qJ(4));
t179 = cos(qJ(5));
t180 = cos(qJ(4));
t170 = qJD(1) + qJD(2);
t173 = sin(pkin(9));
t224 = t170 * t173;
t118 = (-t180 * t175 - t176 * t179) * t224;
t221 = t173 * t180;
t206 = t170 * t221;
t223 = t170 * t176;
t207 = t173 * t223;
t120 = -t175 * t207 + t179 * t206;
t230 = t120 * t118;
t244 = -t146 + t230;
t247 = t175 * t244;
t246 = t179 * t244;
t178 = sin(qJ(1));
t182 = cos(qJ(1));
t199 = t178 * g(1) - t182 * g(2);
t144 = qJDD(1) * pkin(1) + t199;
t192 = t182 * g(1) + t178 * g(2);
t145 = -qJD(1) ^ 2 * pkin(1) - t192;
t177 = sin(qJ(2));
t181 = cos(qJ(2));
t114 = t177 * t144 + t181 * t145;
t166 = t170 ^ 2;
t106 = -t166 * pkin(2) + t167 * qJ(3) + t114;
t242 = 2 * qJD(3);
t245 = t170 * t242 + t106;
t193 = -t174 * pkin(3) - t173 * pkin(7);
t194 = -t173 * g(3) + t245 * t174;
t238 = t174 * g(3);
t56 = t174 * t194 + t173 * (t245 * t173 + t238);
t126 = (qJD(4) * t170 * t180 + t167 * t176) * t173;
t219 = t174 * t170;
t156 = -qJD(4) + t219;
t132 = t156 * t206;
t108 = t132 - t126;
t243 = t173 * t108;
t168 = t173 ^ 2;
t169 = t174 ^ 2;
t140 = (t168 + t169) * t166;
t150 = -qJD(5) + t156;
t149 = qJD(4) * t207;
t127 = t167 * t221 - t149;
t196 = t179 * t126 + t175 * t127;
t63 = (qJD(5) + t150) * t120 + t196;
t116 = t118 ^ 2;
t117 = t120 ^ 2;
t147 = t150 ^ 2;
t153 = t156 ^ 2;
t226 = t168 * t166;
t209 = t180 * t226;
t136 = t193 * t170;
t78 = t136 * t219 + t194;
t113 = t181 * t144 - t177 * t145;
t105 = -t167 * pkin(2) - t166 * qJ(3) + qJDD(3) - t113;
t92 = t193 * t167 + t105;
t85 = t180 * t92;
t43 = -t154 * pkin(4) - t127 * pkin(8) + t85 + (pkin(8) * t156 * t224 - pkin(4) * t209 - t78) * t176;
t125 = -t156 * pkin(4) - pkin(8) * t206;
t155 = t176 ^ 2 * t226;
t50 = t176 * t92 + t180 * t78;
t44 = -pkin(4) * t155 - t126 * pkin(8) + t156 * t125 + t50;
t21 = t175 * t44 - t179 * t43;
t22 = t175 * t43 + t179 * t44;
t10 = t175 * t22 - t179 * t21;
t241 = pkin(4) * t10;
t112 = t118 * t150;
t191 = -t175 * t126 + t179 * t127;
t80 = t118 * qJD(5) + t191;
t66 = t112 + t80;
t40 = -t175 * t63 - t179 * t66;
t240 = pkin(4) * t40;
t49 = t176 * t78 - t85;
t215 = t242 + t136;
t77 = t238 + (t215 * t170 + t106) * t173;
t25 = t174 * (t176 * t49 + t180 * t50) + t173 * t77;
t35 = t176 * t50 - t180 * t49;
t237 = -pkin(2) * t35 + qJ(3) * t25;
t51 = t238 + t126 * pkin(4) - pkin(8) * t155 + (t106 + (t125 * t180 + t215) * t170) * t173;
t236 = t175 * t51;
t89 = t146 + t230;
t235 = t175 * t89;
t234 = t179 * t51;
t233 = t179 * t89;
t232 = t180 * t10;
t231 = -pkin(2) * t105 + qJ(3) * t56;
t229 = t150 * t175;
t228 = t150 * t179;
t227 = t167 * t180;
t225 = t168 * t174;
t222 = t173 * t167;
t143 = t176 * t209;
t123 = -t143 + t154;
t218 = t176 * t123;
t124 = -t143 - t154;
t217 = t180 * t124;
t216 = t181 * t167;
t214 = qJD(5) - t150;
t11 = t175 * t21 + t179 * t22;
t5 = t174 * (-t176 * t10 + t180 * t11) + t173 * t51;
t7 = t176 * t11 + t232;
t211 = t173 * (-pkin(8) * t232 - t176 * (-pkin(4) * t51 + pkin(8) * t11) - pkin(7) * t7) + t174 * (-pkin(3) * t7 - t241) - pkin(2) * t7 + qJ(3) * t5;
t210 = t156 * t223;
t208 = t166 * t225;
t205 = t180 ^ 2 * t226;
t204 = t174 * t230;
t134 = (-t174 * t169 - t225) * t166;
t203 = pkin(2) * t220 + qJ(3) * t134 - t174 * t105;
t41 = t175 * t66 - t179 * t63;
t81 = -t116 - t117;
t16 = t174 * (-t176 * t40 + t180 * t41) + t173 * t81;
t19 = t176 * t41 + t180 * t40;
t202 = t173 * (t180 * (-pkin(8) * t40 - t10) - t176 * (-pkin(4) * t81 + pkin(8) * t41 + t11) - pkin(7) * t19) + t174 * (-pkin(3) * t19 - t240) - pkin(2) * t19 + qJ(3) * t16;
t86 = -t147 - t116;
t53 = t175 * t86 + t246;
t189 = pkin(4) * t53 - t21;
t54 = t179 * t86 - t247;
t62 = t214 * t120 + t196;
t27 = t174 * (-t176 * t53 + t180 * t54) + t173 * t62;
t37 = t176 * t54 + t180 * t53;
t201 = t173 * (t180 * (-pkin(8) * t53 + t236) - t176 * (-pkin(4) * t62 + pkin(8) * t54 - t234) - pkin(7) * t37) + t174 * (-pkin(3) * t37 - t189) - pkin(2) * t37 + qJ(3) * t27;
t100 = -t117 - t147;
t58 = t179 * t100 + t235;
t186 = pkin(4) * t58 - t22;
t59 = -t175 * t100 + t233;
t67 = t214 * t118 + t191;
t29 = t174 * (-t176 * t58 + t180 * t59) + t173 * t67;
t39 = t176 * t59 + t180 * t58;
t200 = t173 * (t180 * (-pkin(8) * t58 + t234) - t176 * (-pkin(4) * t67 + pkin(8) * t59 + t236) - pkin(7) * t39) + t174 * (-pkin(3) * t39 - t186) - pkin(2) * t39 + qJ(3) * t29;
t115 = -t205 - t153;
t185 = -t149 + (t210 + t227) * t173;
t69 = t174 * (-t176 * t115 + t180 * t123) + t185 * t173;
t95 = t180 * t115 + t218;
t198 = t173 * (-pkin(7) * t95 + t180 * t77) + t174 * (-pkin(3) * t95 + t50) - pkin(2) * t95 + qJ(3) * t69;
t128 = -t155 - t153;
t73 = t174 * (-t176 * t124 + t180 * t128) - t243;
t99 = t176 * t128 + t217;
t197 = t173 * (-pkin(7) * t99 + t176 * t77) + t174 * (-pkin(3) * t99 + t49) - pkin(2) * t99 + qJ(3) * t73;
t161 = t168 * t167;
t162 = t169 * t167;
t138 = t162 + t161;
t195 = pkin(2) * t140 + qJ(3) * t138 + t56;
t133 = t173 * t140;
t190 = -pkin(2) * t222 + qJ(3) * t133 + t173 * t105;
t107 = t132 + t126;
t109 = -t149 + (-t210 + t227) * t173;
t61 = t174 * (-t180 * t107 + t176 * t109) - t173 * (t155 + t205);
t75 = -t176 * t107 - t180 * t109;
t188 = qJ(3) * t61 - t173 * t35 + (-pkin(2) + t193) * t75;
t148 = 0.2e1 * t173 * t220;
t141 = t174 * t154;
t129 = -t155 + t205;
t111 = -t117 + t147;
t110 = t116 - t147;
t96 = t117 - t116;
t83 = (t173 * (t156 * t207 + t127) - t176 * t208) * t180;
t82 = (t180 * t208 - t243) * t176;
t79 = -t120 * qJD(5) - t196;
t72 = t173 * (t180 * (t155 - t153) + t218) + t174 * t107;
t71 = t173 * (t217 - t176 * (t153 - t205)) - t174 * t109;
t65 = -t112 + t80;
t60 = t173 * (t180 * t108 - t176 * t185) - t174 * t129;
t47 = t174 * t146 + t173 * (t180 * (-t118 * t179 - t120 * t175) - t176 * (-t118 * t175 + t120 * t179)) * t150;
t34 = t173 * (t180 * (t120 * t229 + t179 * t80) - t176 * (-t120 * t228 + t175 * t80)) + t204;
t33 = t173 * (t180 * (t118 * t228 - t175 * t79) - t176 * (t118 * t229 + t179 * t79)) - t204;
t31 = t173 * (t180 * (t179 * t110 + t235) - t176 * (t175 * t110 - t233)) + t174 * t63;
t30 = t173 * (t180 * (-t175 * t111 + t246) - t176 * (t179 * t111 + t247)) - t174 * t66;
t17 = t173 * (t180 * (-t175 * t65 - t179 * t62) - t176 * (-t175 * t62 + t179 * t65)) - t174 * t96;
t1 = [0, 0, 0, 0, 0, qJDD(1), t199, t192, 0, 0, 0, 0, 0, 0, 0, t167, pkin(1) * (-t177 * t166 + t216) + t113, pkin(1) * (-t181 * t166 - t177 * t167) - t114, 0, pkin(1) * (t181 * t113 + t177 * t114), t161, t148, 0, t162, 0, 0, pkin(1) * (t177 * t134 + t174 * t216) + t203, pkin(1) * (t177 * t133 - t173 * t216) + t190, pkin(1) * (t177 * t138 + t181 * t140) + t195, pkin(1) * (-t181 * t105 + t177 * t56) + t231, t83, t60, t71, t82, t72, t141, pkin(1) * (t177 * t73 - t181 * t99) + t197, pkin(1) * (t177 * t69 - t181 * t95) + t198, pkin(1) * (t177 * t61 - t181 * t75) + t188, pkin(1) * t177 * t25 + (-pkin(1) * t181 + t193) * t35 + t237, t34, t17, t30, t33, t31, t47, pkin(1) * (t177 * t27 - t181 * t37) + t201, pkin(1) * (t177 * t29 - t181 * t39) + t200, pkin(1) * (t177 * t16 - t181 * t19) + t202, pkin(1) * (t177 * t5 - t181 * t7) + t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, t113, -t114, 0, 0, t161, t148, 0, t162, 0, 0, t203, t190, t195, t231, t83, t60, t71, t82, t72, t141, t197, t198, t188, t193 * t35 + t237, t34, t17, t30, t33, t31, t47, t201, t200, t202, t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t222, -t140, t105, 0, 0, 0, 0, 0, 0, t99, t95, t75, t35, 0, 0, 0, 0, 0, 0, t37, t39, t19, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t129, t109, -t143, -t107, -t154, -t49, -t50, 0, 0, -t230, t96, t66, t230, -t63, -t146, t189, t186, t240, t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230, t96, t66, t230, -t63, -t146, -t21, -t22, 0, 0;];
tauJ_reg = t1;
