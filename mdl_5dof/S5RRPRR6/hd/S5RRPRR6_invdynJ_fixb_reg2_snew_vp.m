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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:36:11
% EndTime: 2019-12-05 18:36:19
% DurationCPUTime: 2.30s
% Computational Cost: add. (12687->254), mult. (17686->362), div. (0->0), fcn. (11446->10), ass. (0->177)
t176 = sin(qJ(5));
t168 = qJDD(1) + qJDD(2);
t175 = cos(pkin(9));
t221 = t175 * t168;
t154 = -qJDD(4) + t221;
t146 = -qJDD(5) + t154;
t177 = sin(qJ(4));
t180 = cos(qJ(5));
t181 = cos(qJ(4));
t171 = qJD(1) + qJD(2);
t174 = sin(pkin(9));
t225 = t171 * t174;
t118 = (-t176 * t181 - t177 * t180) * t225;
t222 = t174 * t181;
t206 = t171 * t222;
t224 = t171 * t177;
t207 = t174 * t224;
t120 = -t176 * t207 + t180 * t206;
t231 = t120 * t118;
t245 = -t146 + t231;
t248 = t176 * t245;
t247 = t180 * t245;
t179 = sin(qJ(1));
t183 = cos(qJ(1));
t216 = g(2) * t183 + g(3) * t179;
t144 = qJDD(1) * pkin(1) + t216;
t193 = -g(2) * t179 + g(3) * t183;
t145 = -qJD(1) ^ 2 * pkin(1) - t193;
t178 = sin(qJ(2));
t182 = cos(qJ(2));
t114 = t178 * t144 + t182 * t145;
t167 = t171 ^ 2;
t106 = -pkin(2) * t167 + qJ(3) * t168 + t114;
t243 = 2 * qJD(3);
t246 = t171 * t243 + t106;
t194 = -t175 * pkin(3) - pkin(7) * t174;
t195 = -t174 * g(1) + t175 * t246;
t239 = t175 * g(1);
t56 = t175 * t195 + t174 * (t174 * t246 + t239);
t126 = (qJD(4) * t171 * t181 + t168 * t177) * t174;
t220 = t175 * t171;
t156 = -qJD(4) + t220;
t132 = t156 * t206;
t108 = t132 - t126;
t244 = t174 * t108;
t169 = t174 ^ 2;
t170 = t175 ^ 2;
t140 = (t169 + t170) * t167;
t150 = -qJD(5) + t156;
t149 = qJD(4) * t207;
t127 = t168 * t222 - t149;
t197 = t126 * t180 + t176 * t127;
t63 = (qJD(5) + t150) * t120 + t197;
t116 = t118 ^ 2;
t117 = t120 ^ 2;
t147 = t150 ^ 2;
t153 = t156 ^ 2;
t227 = t169 * t167;
t209 = t181 * t227;
t136 = t194 * t171;
t78 = t136 * t220 + t195;
t113 = t144 * t182 - t178 * t145;
t105 = -pkin(2) * t168 - qJ(3) * t167 + qJDD(3) - t113;
t92 = t168 * t194 + t105;
t85 = t181 * t92;
t43 = -t154 * pkin(4) - t127 * pkin(8) + t85 + (pkin(8) * t156 * t225 - pkin(4) * t209 - t78) * t177;
t125 = -pkin(4) * t156 - pkin(8) * t206;
t155 = t177 ^ 2 * t227;
t50 = t177 * t92 + t181 * t78;
t44 = -pkin(4) * t155 - pkin(8) * t126 + t125 * t156 + t50;
t21 = t176 * t44 - t180 * t43;
t22 = t176 * t43 + t180 * t44;
t10 = t176 * t22 - t180 * t21;
t242 = pkin(4) * t10;
t112 = t118 * t150;
t192 = -t176 * t126 + t180 * t127;
t80 = qJD(5) * t118 + t192;
t66 = t112 + t80;
t40 = -t176 * t63 - t180 * t66;
t241 = pkin(4) * t40;
t49 = t177 * t78 - t85;
t215 = t243 + t136;
t77 = t239 + (t171 * t215 + t106) * t174;
t25 = t175 * (t177 * t49 + t181 * t50) + t174 * t77;
t35 = t177 * t50 - t181 * t49;
t238 = -pkin(2) * t35 + qJ(3) * t25;
t51 = t239 + t126 * pkin(4) - pkin(8) * t155 + (t106 + (t125 * t181 + t215) * t171) * t174;
t237 = t176 * t51;
t89 = t146 + t231;
t236 = t176 * t89;
t235 = t180 * t51;
t234 = t180 * t89;
t233 = t181 * t10;
t232 = -pkin(2) * t105 + qJ(3) * t56;
t230 = t150 * t176;
t229 = t150 * t180;
t228 = t168 * t181;
t226 = t169 * t175;
t223 = t174 * t168;
t143 = t177 * t209;
t123 = -t143 + t154;
t219 = t177 * t123;
t124 = -t143 - t154;
t218 = t181 * t124;
t217 = t182 * t168;
t214 = qJD(5) - t150;
t11 = t176 * t21 + t180 * t22;
t5 = t175 * (-t10 * t177 + t11 * t181) + t174 * t51;
t7 = t11 * t177 + t233;
t211 = t174 * (-pkin(8) * t233 - t177 * (-pkin(4) * t51 + pkin(8) * t11) - pkin(7) * t7) + t175 * (-pkin(3) * t7 - t242) - pkin(2) * t7 + qJ(3) * t5;
t210 = t156 * t224;
t208 = t167 * t226;
t205 = t181 ^ 2 * t227;
t204 = t175 * t231;
t134 = (-t170 * t175 - t226) * t167;
t203 = pkin(2) * t221 + qJ(3) * t134 - t105 * t175;
t41 = t176 * t66 - t180 * t63;
t81 = -t116 - t117;
t16 = t175 * (-t177 * t40 + t181 * t41) + t174 * t81;
t19 = t177 * t41 + t181 * t40;
t202 = t174 * (t181 * (-pkin(8) * t40 - t10) - t177 * (-pkin(4) * t81 + pkin(8) * t41 + t11) - pkin(7) * t19) + t175 * (-pkin(3) * t19 - t241) - pkin(2) * t19 + qJ(3) * t16;
t86 = -t147 - t116;
t53 = t176 * t86 + t247;
t190 = pkin(4) * t53 - t21;
t54 = t180 * t86 - t248;
t62 = t120 * t214 + t197;
t27 = t175 * (-t177 * t53 + t181 * t54) + t174 * t62;
t37 = t177 * t54 + t181 * t53;
t201 = t174 * (t181 * (-pkin(8) * t53 + t237) - t177 * (-pkin(4) * t62 + pkin(8) * t54 - t235) - pkin(7) * t37) + t175 * (-pkin(3) * t37 - t190) - pkin(2) * t37 + qJ(3) * t27;
t100 = -t117 - t147;
t58 = t100 * t180 + t236;
t187 = pkin(4) * t58 - t22;
t59 = -t100 * t176 + t234;
t67 = t118 * t214 + t192;
t29 = t175 * (-t177 * t58 + t181 * t59) + t174 * t67;
t39 = t177 * t59 + t181 * t58;
t200 = t174 * (t181 * (-pkin(8) * t58 + t235) - t177 * (-pkin(4) * t67 + pkin(8) * t59 + t237) - pkin(7) * t39) + t175 * (-pkin(3) * t39 - t187) - pkin(2) * t39 + qJ(3) * t29;
t115 = -t205 - t153;
t186 = -t149 + (t210 + t228) * t174;
t69 = t175 * (-t115 * t177 + t123 * t181) + t186 * t174;
t95 = t115 * t181 + t219;
t199 = t174 * (-pkin(7) * t95 + t181 * t77) + t175 * (-pkin(3) * t95 + t50) - pkin(2) * t95 + qJ(3) * t69;
t128 = -t155 - t153;
t73 = t175 * (-t124 * t177 + t128 * t181) - t244;
t99 = t128 * t177 + t218;
t198 = t174 * (-pkin(7) * t99 + t177 * t77) + t175 * (-pkin(3) * t99 + t49) - pkin(2) * t99 + qJ(3) * t73;
t161 = t169 * t168;
t162 = t170 * t168;
t138 = t162 + t161;
t196 = pkin(2) * t140 + qJ(3) * t138 + t56;
t133 = t174 * t140;
t191 = -pkin(2) * t223 + qJ(3) * t133 + t105 * t174;
t107 = t132 + t126;
t109 = -t149 + (-t210 + t228) * t174;
t61 = t175 * (-t107 * t181 + t109 * t177) - t174 * (t155 + t205);
t75 = -t107 * t177 - t109 * t181;
t189 = qJ(3) * t61 - t174 * t35 + (-pkin(2) + t194) * t75;
t148 = 0.2e1 * t174 * t221;
t141 = t175 * t154;
t129 = -t155 + t205;
t111 = -t117 + t147;
t110 = t116 - t147;
t96 = t117 - t116;
t83 = (t174 * (t156 * t207 + t127) - t177 * t208) * t181;
t82 = (t181 * t208 - t244) * t177;
t79 = -qJD(5) * t120 - t197;
t72 = t174 * (t181 * (t155 - t153) + t219) + t175 * t107;
t71 = t174 * (t218 - t177 * (t153 - t205)) - t175 * t109;
t65 = -t112 + t80;
t60 = t174 * (t181 * t108 - t177 * t186) - t175 * t129;
t47 = t175 * t146 + t174 * (t181 * (-t118 * t180 - t120 * t176) - t177 * (-t118 * t176 + t120 * t180)) * t150;
t34 = t174 * (t181 * (t120 * t230 + t180 * t80) - t177 * (-t120 * t229 + t176 * t80)) + t204;
t33 = t174 * (t181 * (t118 * t229 - t176 * t79) - t177 * (t118 * t230 + t180 * t79)) - t204;
t31 = t174 * (t181 * (t110 * t180 + t236) - t177 * (t110 * t176 - t234)) + t175 * t63;
t30 = t174 * (t181 * (-t111 * t176 + t247) - t177 * (t111 * t180 + t248)) - t175 * t66;
t17 = t174 * (t181 * (-t176 * t65 - t180 * t62) - t177 * (-t176 * t62 + t180 * t65)) - t175 * t96;
t1 = [0, 0, 0, 0, 0, qJDD(1), t216, t193, 0, 0, 0, 0, 0, 0, 0, t168, pkin(1) * (-t167 * t178 + t217) + t113, pkin(1) * (-t167 * t182 - t168 * t178) - t114, 0, pkin(1) * (t113 * t182 + t114 * t178), t161, t148, 0, t162, 0, 0, pkin(1) * (t134 * t178 + t175 * t217) + t203, pkin(1) * (t133 * t178 - t174 * t217) + t191, pkin(1) * (t138 * t178 + t140 * t182) + t196, pkin(1) * (-t105 * t182 + t178 * t56) + t232, t83, t60, t71, t82, t72, t141, pkin(1) * (t178 * t73 - t182 * t99) + t198, pkin(1) * (t178 * t69 - t182 * t95) + t199, pkin(1) * (t178 * t61 - t182 * t75) + t189, pkin(1) * t178 * t25 + (-pkin(1) * t182 + t194) * t35 + t238, t34, t17, t30, t33, t31, t47, pkin(1) * (t178 * t27 - t182 * t37) + t201, pkin(1) * (t178 * t29 - t182 * t39) + t200, pkin(1) * (t16 * t178 - t182 * t19) + t202, pkin(1) * (t178 * t5 - t182 * t7) + t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t113, -t114, 0, 0, t161, t148, 0, t162, 0, 0, t203, t191, t196, t232, t83, t60, t71, t82, t72, t141, t198, t199, t189, t194 * t35 + t238, t34, t17, t30, t33, t31, t47, t201, t200, t202, t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t221, t223, -t140, t105, 0, 0, 0, 0, 0, 0, t99, t95, t75, t35, 0, 0, 0, 0, 0, 0, t37, t39, t19, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t129, t109, -t143, -t107, -t154, -t49, -t50, 0, 0, -t231, t96, t66, t231, -t63, -t146, t190, t187, t241, t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231, t96, t66, t231, -t63, -t146, -t21, -t22, 0, 0;];
tauJ_reg = t1;
