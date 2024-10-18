% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRR13
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
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRR13_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:16
% EndTime: 2024-09-27 17:32:19
% DurationCPUTime: 1.25s
% Computational Cost: add. (24076->248), mult. (28602->374), div. (0->0), fcn. (19064->12), ass. (0->183)
t182 = sin(qJ(5));
t187 = cos(qJ(5));
t177 = qJD(1) + qJD(2);
t172 = qJD(3) + t177;
t180 = sin(pkin(5));
t188 = cos(qJ(4));
t231 = t180 * t188;
t223 = t172 * t231;
t183 = sin(qJ(4));
t232 = t180 * t183;
t224 = t172 * t232;
t135 = t182 * t224 - t187 * t223;
t234 = t172 * t180;
t137 = (t188 * t182 + t183 * t187) * t234;
t116 = t137 * t135;
t175 = qJDD(1) + qJDD(2);
t171 = qJDD(3) + t175;
t181 = cos(pkin(5));
t165 = t181 * t171 + qJDD(4);
t157 = qJDD(5) + t165;
t250 = -t116 + t157;
t252 = t182 * t250;
t251 = t187 * t250;
t218 = qJD(4) * t234;
t143 = t171 * t232 + t188 * t218;
t144 = t171 * t231 - t183 * t218;
t102 = -t135 * qJD(5) + t187 * t143 + t182 * t144;
t167 = t181 * t172 + qJD(4);
t162 = qJD(5) + t167;
t124 = t162 * t135;
t249 = -t124 + t102;
t170 = t172 ^ 2;
t176 = t180 ^ 2;
t233 = t176 * t183;
t225 = t170 * t233;
t155 = t188 * t225;
t141 = t165 + t155;
t186 = sin(qJ(1));
t191 = cos(qJ(1));
t217 = t186 * g(1) - t191 * g(2);
t160 = qJDD(1) * pkin(1) + t217;
t205 = t191 * g(1) + t186 * g(2);
t161 = -qJD(1) ^ 2 * pkin(1) - t205;
t185 = sin(qJ(2));
t190 = cos(qJ(2));
t130 = t190 * t160 - t185 * t161;
t126 = t175 * pkin(2) + t130;
t131 = t185 * t160 + t190 * t161;
t174 = t177 ^ 2;
t127 = -t174 * pkin(2) + t131;
t184 = sin(qJ(3));
t189 = cos(qJ(3));
t105 = t184 * t126 + t189 * t127;
t246 = pkin(9) * t180;
t100 = -t170 * pkin(3) + t171 * t246 + t105;
t104 = t189 * t126 - t184 * t127;
t235 = t170 * t180;
t99 = t171 * pkin(3) + pkin(9) * t235 + t104;
t243 = t181 * t99;
t215 = t183 * t100 - t188 * t243;
t237 = t167 * t172;
t245 = t143 * pkin(10);
t248 = (pkin(4) * t225 + (pkin(10) * t237 - g(3)) * t180) * t188 + t165 * pkin(4) - t215 - t245;
t133 = t135 ^ 2;
t134 = t137 ^ 2;
t158 = t162 ^ 2;
t166 = t167 ^ 2;
t247 = pkin(4) * t181;
t64 = g(3) * t231 + t215;
t65 = -g(3) * t232 + t188 * t100 + t183 * t243;
t95 = t181 * g(3) + t180 * t99;
t32 = t180 * t95 + (t183 * t65 - t188 * t64) * t181;
t40 = t183 * t64 + t188 * t65;
t244 = pkin(3) * t32 + t40 * t246;
t197 = -t167 * pkin(4) + pkin(10) * t224;
t179 = t188 ^ 2;
t236 = t170 * t176;
t221 = t179 * t236;
t69 = t144 * pkin(4) + pkin(10) * t221 + t197 * t224 + t95;
t242 = t182 * t69;
t52 = -pkin(4) * t221 + t144 * pkin(10) + t167 * t197 + t65;
t241 = t187 * t52;
t240 = t187 * t69;
t239 = t162 * t182;
t238 = t162 * t187;
t111 = t116 + t157;
t230 = t182 * t111;
t229 = t183 * t141;
t228 = t187 * t111;
t142 = -t155 + t165;
t227 = t188 * t142;
t16 = t184 * t40 + t189 * t32;
t226 = pkin(2) * t16 + t244;
t178 = t183 ^ 2;
t222 = t178 * t236;
t220 = t181 * t116;
t148 = t167 * t223;
t219 = t148 + t143;
t29 = t182 * t52 - t187 * t248;
t30 = t248 * t182 + t241;
t14 = t182 * t29 + t187 * t30;
t119 = -t148 + t143;
t147 = t167 * t224;
t120 = t144 + t147;
t77 = (-t188 * t119 + t183 * t120) * t181 - (-t178 - t179) * t176 * t235;
t98 = t183 * t119 + t188 * t120;
t216 = pkin(3) * t77 + t65 * t231 + t64 * t232 + t98 * t246;
t129 = -t222 - t166;
t113 = -t183 * t129 - t227;
t88 = -t180 * t219 + (t188 * t129 - t183 * t142) * t181;
t214 = pkin(3) * t88 + t113 * t246 - t181 * t65 - t95 * t232;
t145 = -t166 - t221;
t117 = t188 * t145 - t229;
t121 = t144 - t147;
t92 = t180 * t121 + (t141 * t188 + t145 * t183) * t181;
t213 = pkin(3) * t92 + t117 * t246 - t181 * t64 + t95 * t231;
t212 = t182 * t143 - t187 * t144;
t103 = -t133 - t134;
t13 = t182 * t30 - t187 * t29;
t195 = (-qJD(5) + t162) * t137 - t212;
t85 = t124 + t102;
t56 = t182 * t195 - t187 * t85;
t57 = t182 * t85 + t187 * t195;
t26 = -t180 * t103 + (t183 * t57 + t188 * t56) * t181;
t34 = -t183 * t56 + t188 * t57;
t211 = (-pkin(10) * t56 - t13) * t232 + (-pkin(4) * t103 + pkin(10) * t57 + t14) * t231 + t56 * t247 + pkin(3) * t26 + t34 * t246;
t109 = -t158 - t133;
t74 = t182 * t109 + t251;
t20 = pkin(4) * t74 - t29;
t75 = t187 * t109 - t252;
t81 = (qJD(5) + t162) * t137 + t212;
t38 = -t180 * t81 + (t183 * t75 + t188 * t74) * t181;
t50 = -t183 * t74 + t188 * t75;
t210 = t181 * t20 + (-pkin(4) * t81 + pkin(10) * t75 + t240) * t231 + pkin(3) * t38 + (-pkin(10) * t74 - t242) * t232 + t50 * t246;
t118 = -t134 - t158;
t79 = t187 * t118 - t230;
t22 = -t241 - t182 * (pkin(10) * t148 - t245 - t64) + (-t182 * t141 + t79) * pkin(4);
t80 = -t182 * t118 - t228;
t42 = -t180 * t249 + (t183 * t80 + t188 * t79) * t181;
t54 = -t183 * t79 + t188 * t80;
t209 = t181 * t22 + (-pkin(4) * t249 + pkin(10) * t80 - t242) * t231 + pkin(3) * t42 + (-pkin(10) * t79 - t240) * t232 + t54 * t246;
t59 = t184 * t98 + t189 * t77;
t208 = pkin(2) * t59 + t216;
t68 = t184 * t113 + t189 * t88;
t207 = pkin(2) * t68 + t214;
t71 = t184 * t117 + t189 * t92;
t206 = pkin(2) * t71 + t213;
t199 = t184 * t170 - t189 * t171;
t204 = -pkin(2) * t199 + t104;
t12 = t184 * t34 + t189 * t26;
t203 = pkin(2) * t12 + t211;
t18 = t184 * t50 + t189 * t38;
t202 = pkin(2) * t18 + t210;
t24 = t184 * t54 + t189 * t42;
t201 = pkin(2) * t24 + t209;
t200 = -t170 * t181 + t237;
t151 = -t189 * t170 - t184 * t171;
t198 = pkin(2) * t151 - t105;
t4 = t180 * t69 + (t13 * t188 + t14 * t183) * t181;
t6 = -t183 * t13 + t188 * t14;
t196 = pkin(3) * t4 + t6 * t246 + (pkin(4) * t69 + pkin(10) * t14) * t231 + (-pkin(10) * t232 + t247) * t13;
t2 = t184 * t6 + t189 * t4;
t194 = pkin(2) * t2 + t196;
t154 = t181 * t165;
t146 = (t178 - t179) * t236;
t123 = -t134 + t158;
t122 = t133 - t158;
t115 = t134 - t133;
t107 = (t200 * t188 * t176 + t143 * t180) * t183;
t106 = (t144 * t180 - t200 * t233) * t188;
t101 = -t137 * qJD(5) - t212;
t94 = t181 * t120 + (t183 * (-t166 + t221) + t227) * t180;
t93 = t181 * t119 + (t229 + t188 * (t166 - t222)) * t180;
t78 = t181 * t146 + (t183 * t121 + t188 * t219) * t180;
t73 = t189 * t104 + t184 * t105;
t72 = pkin(2) * t73;
t66 = t181 * t157 + (t183 * (-t135 * t187 + t137 * t182) + t188 * (-t135 * t182 - t137 * t187)) * t180 * t162;
t46 = t220 + (t183 * (t187 * t102 - t137 * t239) + t188 * (t182 * t102 + t137 * t238)) * t180;
t45 = -t220 + (t183 * (-t182 * t101 + t135 * t238) + t188 * (t187 * t101 + t135 * t239)) * t180;
t44 = t181 * t85 + (t183 * (-t182 * t123 + t251) + t188 * (t187 * t123 + t252)) * t180;
t43 = t181 * t195 + (t183 * (t187 * t122 - t230) + t188 * (t182 * t122 + t228)) * t180;
t27 = t181 * t115 + (t183 * (-t182 * t249 - t187 * t81) + t188 * (-t182 * t81 + t187 * t249)) * t180;
t1 = [0, 0, 0, 0, 0, qJDD(1), t217, t205, 0, 0, 0, 0, 0, 0, 0, t175, pkin(1) * (-t185 * t174 + t190 * t175) + t130, pkin(1) * (-t190 * t174 - t185 * t175) - t131, 0, pkin(1) * (t190 * t130 + t185 * t131), 0, 0, 0, 0, 0, t171, pkin(1) * (t185 * t151 - t190 * t199) + t204, pkin(1) * (t190 * t151 + t185 * t199) + t198, 0, pkin(1) * (t185 * (-t184 * t104 + t189 * t105) + t190 * t73) + t72, t107, t78, t93, t106, t94, t154, pkin(1) * (t185 * (t189 * t117 - t184 * t92) + t190 * t71) + t206, pkin(1) * (t185 * (t189 * t113 - t184 * t88) + t190 * t68) + t207, pkin(1) * (t185 * (-t184 * t77 + t189 * t98) + t190 * t59) + t208, pkin(1) * (t185 * (-t184 * t32 + t189 * t40) + t190 * t16) + t226, t46, t27, t44, t45, t43, t66, pkin(1) * (t185 * (-t184 * t38 + t189 * t50) + t190 * t18) + t202, pkin(1) * (t185 * (-t184 * t42 + t189 * t54) + t190 * t24) + t201, pkin(1) * (t185 * (-t184 * t26 + t189 * t34) + t190 * t12) + t203, pkin(1) * (t185 * (-t184 * t4 + t189 * t6) + t190 * t2) + t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, t130, -t131, 0, 0, 0, 0, 0, 0, 0, t171, t204, t198, 0, t72, t107, t78, t93, t106, t94, t154, t206, t207, t208, t226, t46, t27, t44, t45, t43, t66, t202, t201, t203, t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, t104, -t105, 0, 0, t107, t78, t93, t106, t94, t154, t213, t214, t216, t244, t46, t27, t44, t45, t43, t66, t210, t209, t211, t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, t146, t119, t155, t120, t165, -t64, -t65, 0, 0, t116, t115, t85, -t116, t195, t157, t20, t22, pkin(4) * t56, pkin(4) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t115, t85, -t116, t195, t157, -t29, -t30, 0, 0;];
tauJ_reg = t1;
