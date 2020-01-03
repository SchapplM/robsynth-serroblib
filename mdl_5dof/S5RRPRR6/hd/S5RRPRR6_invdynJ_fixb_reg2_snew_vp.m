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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:05:55
% EndTime: 2020-01-03 12:06:02
% DurationCPUTime: 2.29s
% Computational Cost: add. (12687->254), mult. (17686->362), div. (0->0), fcn. (11446->10), ass. (0->177)
t174 = sin(qJ(5));
t166 = qJDD(1) + qJDD(2);
t173 = cos(pkin(9));
t219 = t173 * t166;
t154 = -qJDD(4) + t219;
t146 = -qJDD(5) + t154;
t175 = sin(qJ(4));
t178 = cos(qJ(5));
t179 = cos(qJ(4));
t169 = qJD(1) + qJD(2);
t172 = sin(pkin(9));
t223 = t169 * t172;
t118 = (-t179 * t174 - t175 * t178) * t223;
t220 = t172 * t179;
t205 = t169 * t220;
t222 = t169 * t175;
t206 = t172 * t222;
t120 = -t174 * t206 + t178 * t205;
t229 = t120 * t118;
t243 = -t146 + t229;
t246 = t174 * t243;
t245 = t178 * t243;
t177 = sin(qJ(1));
t181 = cos(qJ(1));
t192 = -t181 * g(2) - t177 * g(3);
t144 = qJDD(1) * pkin(1) + t192;
t191 = t177 * g(2) - t181 * g(3);
t145 = -qJD(1) ^ 2 * pkin(1) - t191;
t176 = sin(qJ(2));
t180 = cos(qJ(2));
t114 = t176 * t144 + t180 * t145;
t165 = t169 ^ 2;
t106 = -t165 * pkin(2) + t166 * qJ(3) + t114;
t241 = 2 * qJD(3);
t244 = t169 * t241 + t106;
t193 = -t173 * pkin(3) - t172 * pkin(7);
t194 = -t172 * g(1) + t244 * t173;
t237 = t173 * g(1);
t56 = t173 * t194 + t172 * (t244 * t172 + t237);
t126 = (qJD(4) * t169 * t179 + t166 * t175) * t172;
t218 = t173 * t169;
t156 = -qJD(4) + t218;
t132 = t156 * t205;
t108 = t132 - t126;
t242 = t172 * t108;
t167 = t172 ^ 2;
t168 = t173 ^ 2;
t140 = (t167 + t168) * t165;
t150 = -qJD(5) + t156;
t149 = qJD(4) * t206;
t127 = t166 * t220 - t149;
t196 = t178 * t126 + t174 * t127;
t63 = (qJD(5) + t150) * t120 + t196;
t116 = t118 ^ 2;
t117 = t120 ^ 2;
t147 = t150 ^ 2;
t153 = t156 ^ 2;
t225 = t167 * t165;
t208 = t179 * t225;
t136 = t193 * t169;
t78 = t136 * t218 + t194;
t113 = t180 * t144 - t176 * t145;
t105 = -t166 * pkin(2) - t165 * qJ(3) + qJDD(3) - t113;
t92 = t193 * t166 + t105;
t85 = t179 * t92;
t43 = -t154 * pkin(4) - t127 * pkin(8) + t85 + (pkin(8) * t156 * t223 - pkin(4) * t208 - t78) * t175;
t125 = -t156 * pkin(4) - pkin(8) * t205;
t155 = t175 ^ 2 * t225;
t50 = t175 * t92 + t179 * t78;
t44 = -pkin(4) * t155 - t126 * pkin(8) + t156 * t125 + t50;
t21 = t174 * t44 - t178 * t43;
t22 = t174 * t43 + t178 * t44;
t10 = t174 * t22 - t178 * t21;
t240 = pkin(4) * t10;
t112 = t118 * t150;
t190 = -t174 * t126 + t178 * t127;
t80 = t118 * qJD(5) + t190;
t66 = t112 + t80;
t40 = -t174 * t63 - t178 * t66;
t239 = pkin(4) * t40;
t49 = t175 * t78 - t85;
t214 = t241 + t136;
t77 = t237 + (t214 * t169 + t106) * t172;
t25 = t173 * (t175 * t49 + t179 * t50) + t172 * t77;
t35 = t175 * t50 - t179 * t49;
t236 = -pkin(2) * t35 + qJ(3) * t25;
t51 = t237 + t126 * pkin(4) - pkin(8) * t155 + (t106 + (t125 * t179 + t214) * t169) * t172;
t235 = t174 * t51;
t89 = t146 + t229;
t234 = t174 * t89;
t233 = t178 * t51;
t232 = t178 * t89;
t231 = t179 * t10;
t230 = -pkin(2) * t105 + qJ(3) * t56;
t228 = t150 * t174;
t227 = t150 * t178;
t226 = t166 * t179;
t224 = t167 * t173;
t221 = t172 * t166;
t143 = t175 * t208;
t123 = -t143 + t154;
t217 = t175 * t123;
t124 = -t143 - t154;
t216 = t179 * t124;
t215 = t180 * t166;
t213 = qJD(5) - t150;
t11 = t174 * t21 + t178 * t22;
t5 = t173 * (-t175 * t10 + t179 * t11) + t172 * t51;
t7 = t175 * t11 + t231;
t210 = t172 * (-pkin(8) * t231 - t175 * (-pkin(4) * t51 + pkin(8) * t11) - pkin(7) * t7) + t173 * (-pkin(3) * t7 - t240) - pkin(2) * t7 + qJ(3) * t5;
t209 = t156 * t222;
t207 = t165 * t224;
t204 = t179 ^ 2 * t225;
t203 = t173 * t229;
t134 = (-t173 * t168 - t224) * t165;
t202 = pkin(2) * t219 + qJ(3) * t134 - t173 * t105;
t41 = t174 * t66 - t178 * t63;
t81 = -t116 - t117;
t16 = t173 * (-t175 * t40 + t179 * t41) + t172 * t81;
t19 = t175 * t41 + t179 * t40;
t201 = t172 * (t179 * (-pkin(8) * t40 - t10) - t175 * (-pkin(4) * t81 + pkin(8) * t41 + t11) - pkin(7) * t19) + t173 * (-pkin(3) * t19 - t239) - pkin(2) * t19 + qJ(3) * t16;
t86 = -t147 - t116;
t53 = t174 * t86 + t245;
t188 = pkin(4) * t53 - t21;
t54 = t178 * t86 - t246;
t62 = t213 * t120 + t196;
t27 = t173 * (-t175 * t53 + t179 * t54) + t172 * t62;
t37 = t175 * t54 + t179 * t53;
t200 = t172 * (t179 * (-pkin(8) * t53 + t235) - t175 * (-pkin(4) * t62 + pkin(8) * t54 - t233) - pkin(7) * t37) + t173 * (-pkin(3) * t37 - t188) - pkin(2) * t37 + qJ(3) * t27;
t100 = -t117 - t147;
t58 = t178 * t100 + t234;
t185 = pkin(4) * t58 - t22;
t59 = -t174 * t100 + t232;
t67 = t213 * t118 + t190;
t29 = t173 * (-t175 * t58 + t179 * t59) + t172 * t67;
t39 = t175 * t59 + t179 * t58;
t199 = t172 * (t179 * (-pkin(8) * t58 + t233) - t175 * (-pkin(4) * t67 + pkin(8) * t59 + t235) - pkin(7) * t39) + t173 * (-pkin(3) * t39 - t185) - pkin(2) * t39 + qJ(3) * t29;
t115 = -t204 - t153;
t184 = -t149 + (t209 + t226) * t172;
t69 = t173 * (-t175 * t115 + t179 * t123) + t184 * t172;
t95 = t179 * t115 + t217;
t198 = t172 * (-pkin(7) * t95 + t179 * t77) + t173 * (-pkin(3) * t95 + t50) - pkin(2) * t95 + qJ(3) * t69;
t128 = -t155 - t153;
t73 = t173 * (-t175 * t124 + t179 * t128) - t242;
t99 = t175 * t128 + t216;
t197 = t172 * (-pkin(7) * t99 + t175 * t77) + t173 * (-pkin(3) * t99 + t49) - pkin(2) * t99 + qJ(3) * t73;
t161 = t167 * t166;
t162 = t168 * t166;
t138 = t162 + t161;
t195 = pkin(2) * t140 + qJ(3) * t138 + t56;
t133 = t172 * t140;
t189 = -pkin(2) * t221 + qJ(3) * t133 + t172 * t105;
t107 = t132 + t126;
t109 = -t149 + (-t209 + t226) * t172;
t61 = t173 * (-t179 * t107 + t175 * t109) - t172 * (t155 + t204);
t75 = -t175 * t107 - t179 * t109;
t187 = qJ(3) * t61 - t172 * t35 + (-pkin(2) + t193) * t75;
t148 = 0.2e1 * t172 * t219;
t141 = t173 * t154;
t129 = -t155 + t204;
t111 = -t117 + t147;
t110 = t116 - t147;
t96 = t117 - t116;
t83 = (t172 * (t156 * t206 + t127) - t175 * t207) * t179;
t82 = (t179 * t207 - t242) * t175;
t79 = -t120 * qJD(5) - t196;
t72 = t172 * (t179 * (t155 - t153) + t217) + t173 * t107;
t71 = t172 * (t216 - t175 * (t153 - t204)) - t173 * t109;
t65 = -t112 + t80;
t60 = t172 * (t179 * t108 - t175 * t184) - t173 * t129;
t47 = t173 * t146 + t172 * (t179 * (-t118 * t178 - t120 * t174) - t175 * (-t118 * t174 + t120 * t178)) * t150;
t34 = t172 * (t179 * (t120 * t228 + t178 * t80) - t175 * (-t120 * t227 + t174 * t80)) + t203;
t33 = t172 * (t179 * (t118 * t227 - t174 * t79) - t175 * (t118 * t228 + t178 * t79)) - t203;
t31 = t172 * (t179 * (t178 * t110 + t234) - t175 * (t174 * t110 - t232)) + t173 * t63;
t30 = t172 * (t179 * (-t174 * t111 + t245) - t175 * (t178 * t111 + t246)) - t173 * t66;
t17 = t172 * (t179 * (-t174 * t65 - t178 * t62) - t175 * (-t174 * t62 + t178 * t65)) - t173 * t96;
t1 = [0, 0, 0, 0, 0, qJDD(1), t192, t191, 0, 0, 0, 0, 0, 0, 0, t166, pkin(1) * (-t176 * t165 + t215) + t113, pkin(1) * (-t180 * t165 - t176 * t166) - t114, 0, pkin(1) * (t180 * t113 + t176 * t114), t161, t148, 0, t162, 0, 0, pkin(1) * (t176 * t134 + t173 * t215) + t202, pkin(1) * (t176 * t133 - t172 * t215) + t189, pkin(1) * (t176 * t138 + t180 * t140) + t195, pkin(1) * (-t180 * t105 + t176 * t56) + t230, t83, t60, t71, t82, t72, t141, pkin(1) * (t176 * t73 - t180 * t99) + t197, pkin(1) * (t176 * t69 - t180 * t95) + t198, pkin(1) * (t176 * t61 - t180 * t75) + t187, pkin(1) * t176 * t25 + (-pkin(1) * t180 + t193) * t35 + t236, t34, t17, t30, t33, t31, t47, pkin(1) * (t176 * t27 - t180 * t37) + t200, pkin(1) * (t176 * t29 - t180 * t39) + t199, pkin(1) * (t176 * t16 - t180 * t19) + t201, pkin(1) * (t176 * t5 - t180 * t7) + t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t113, -t114, 0, 0, t161, t148, 0, t162, 0, 0, t202, t189, t195, t230, t83, t60, t71, t82, t72, t141, t197, t198, t187, t193 * t35 + t236, t34, t17, t30, t33, t31, t47, t200, t199, t201, t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, t221, -t140, t105, 0, 0, 0, 0, 0, 0, t99, t95, t75, t35, 0, 0, 0, 0, 0, 0, t37, t39, t19, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t129, t109, -t143, -t107, -t154, -t49, -t50, 0, 0, -t229, t96, t66, t229, -t63, -t146, t188, t185, t239, t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t229, t96, t66, t229, -t63, -t146, -t21, -t22, 0, 0;];
tauJ_reg = t1;
