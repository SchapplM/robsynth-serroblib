% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:40:46
% EndTime: 2019-03-09 04:40:55
% DurationCPUTime: 3.19s
% Computational Cost: add. (7353->366), mult. (19038->476), div. (0->0), fcn. (14372->8), ass. (0->177)
t179 = sin(qJ(4));
t215 = qJD(4) * t179;
t178 = cos(pkin(9));
t182 = cos(qJ(3));
t220 = t182 * t178;
t165 = qJD(1) * t220;
t176 = sin(pkin(9));
t180 = sin(qJ(3));
t227 = t180 * t176;
t205 = qJD(1) * t227;
t140 = t165 - t205;
t229 = t179 * t140;
t265 = t215 - t229;
t151 = t182 * t176 + t180 * t178;
t175 = sin(pkin(10));
t177 = cos(pkin(10));
t181 = cos(qJ(4));
t255 = -t175 * t179 + t177 * t181;
t92 = t255 * t151;
t133 = qJD(4) - t140;
t213 = t181 * qJD(3);
t254 = t151 * qJD(1);
t116 = t179 * t254 - t213;
t118 = t179 * qJD(3) + t181 * t254;
t66 = t177 * t116 + t175 * t118;
t264 = t133 * t66;
t150 = t175 * t181 + t177 * t179;
t139 = t150 * qJD(4);
t241 = -t150 * t140 + t139;
t214 = qJD(4) * t181;
t225 = t181 * t140;
t240 = (t214 - t225) * t177 - t265 * t175;
t263 = t151 * qJD(2);
t206 = t151 * t214;
t149 = -t220 + t227;
t143 = t149 * qJD(3);
t228 = t179 * t143;
t262 = t206 - t228;
t261 = t254 * qJD(3);
t197 = -t175 * t116 + t177 * t118;
t260 = t197 ^ 2;
t250 = -qJ(5) - pkin(8);
t203 = qJD(4) * t250;
t135 = t181 * qJD(5) + t179 * t203;
t186 = -t179 * qJD(5) + t181 * t203;
t251 = pkin(7) + qJ(2);
t156 = t251 * t176;
t153 = qJD(1) * t156;
t157 = t251 * t178;
t154 = qJD(1) * t157;
t258 = -t182 * t153 - t180 * t154;
t103 = pkin(3) * t254 - t140 * pkin(8);
t95 = t181 * t103;
t48 = pkin(4) * t254 - qJ(5) * t225 - t179 * t258 + t95;
t243 = t179 * t103 + t181 * t258;
t56 = -qJ(5) * t229 + t243;
t246 = (t186 - t48) * t177 + (-t135 + t56) * t175;
t101 = -qJD(3) * pkin(3) - t258;
t64 = t116 * pkin(4) + qJD(5) + t101;
t35 = t66 * pkin(5) - qJ(6) * t197 + t64;
t259 = t35 * t197;
t108 = -t180 * t153 + t182 * t154;
t257 = pkin(4) * t265 - t108;
t256 = -t182 * t156 - t180 * t157;
t161 = qJD(3) * t165;
t128 = -qJD(3) * t205 + t161;
t77 = qJD(4) * t118 + t179 * t128;
t253 = t133 ^ 2;
t102 = qJD(3) * pkin(8) + t108;
t170 = -t178 * pkin(2) - pkin(1);
t155 = t170 * qJD(1) + qJD(2);
t80 = -t140 * pkin(3) - pkin(8) * t254 + t155;
t59 = t181 * t102 + t179 * t80;
t47 = -t116 * qJ(5) + t59;
t43 = t177 * t47;
t58 = -t179 * t102 + t181 * t80;
t46 = -t118 * qJ(5) + t58;
t21 = t175 * t46 + t43;
t252 = t21 * t197;
t144 = t151 * qJD(3);
t129 = qJD(1) * t144;
t187 = t149 * qJD(2);
t72 = -qJD(1) * t187 + qJD(3) * t258;
t89 = t129 * pkin(3) - t128 * pkin(8);
t86 = t181 * t89;
t184 = -qJD(4) * t59 - t179 * t72 + t86;
t76 = qJD(4) * t213 + t181 * t128 - t215 * t254;
t13 = t129 * pkin(4) - t76 * qJ(5) - t118 * qJD(5) + t184;
t190 = -t102 * t215 + t179 * t89 + t181 * t72 + t80 * t214;
t18 = -t77 * qJ(5) - t116 * qJD(5) + t190;
t3 = t177 * t13 - t175 * t18;
t4 = t175 * t13 + t177 * t18;
t29 = t175 * t48 + t177 * t56;
t24 = qJ(6) * t254 + t29;
t88 = t177 * t135 + t175 * t186;
t249 = -t24 + t88;
t248 = -pkin(5) * t254 + t246;
t106 = t149 * pkin(3) - t151 * pkin(8) + t170;
t113 = -t180 * t156 + t182 * t157;
t109 = t181 * t113;
t195 = qJ(5) * t143 - qJD(5) * t151;
t81 = t256 * qJD(3) - t187;
t104 = t144 * pkin(3) + t143 * pkin(8);
t96 = t181 * t104;
t26 = t144 * pkin(4) - t179 * t81 + t96 + t195 * t181 + (-t109 + (qJ(5) * t151 - t106) * t179) * qJD(4);
t211 = t179 * t104 + t106 * t214 + t181 * t81;
t30 = -qJ(5) * t206 + (-qJD(4) * t113 + t195) * t179 + t211;
t8 = t175 * t26 + t177 * t30;
t247 = -t241 * pkin(5) + t240 * qJ(6) + t150 * qJD(6) - t257;
t39 = t133 * pkin(4) + t46;
t20 = t175 * t39 + t43;
t234 = t151 * t181;
t98 = t181 * t106;
t53 = t149 * pkin(4) - qJ(5) * t234 - t179 * t113 + t98;
t235 = t151 * t179;
t242 = t179 * t106 + t109;
t60 = -qJ(5) * t235 + t242;
t34 = t175 * t53 + t177 * t60;
t245 = t175 * t47;
t244 = t76 * t179;
t239 = t116 * t133;
t238 = t118 * t133;
t237 = t118 * t254;
t236 = t254 * t116;
t230 = t179 * t129;
t122 = t181 * t129;
t224 = t181 * t143;
t22 = t177 * t46 - t245;
t219 = qJD(6) - t22;
t218 = t176 ^ 2 + t178 ^ 2;
t217 = qJD(3) * t180;
t216 = qJD(3) * t182;
t212 = qJD(1) * qJD(2);
t210 = t129 * qJ(6) + t4;
t207 = -t181 * pkin(4) - pkin(3);
t204 = t250 * t179;
t49 = t175 * t76 + t177 * t77;
t202 = t218 * qJD(1) ^ 2;
t201 = t133 * t181;
t73 = t263 * qJD(1) - t153 * t217 + t154 * t216;
t82 = -t156 * t217 + t157 * t216 + t263;
t2 = -t129 * pkin(5) - t3;
t158 = t250 * t181;
t114 = -t175 * t158 - t177 * t204;
t115 = -t177 * t158 + t175 * t204;
t50 = -t175 * t77 + t177 * t76;
t200 = t114 * t50 - t115 * t49 - t88 * t66;
t199 = -t66 ^ 2 - t260;
t198 = pkin(4) * t235 - t256;
t7 = -t175 * t30 + t177 * t26;
t19 = t177 * t39 - t245;
t33 = -t175 * t60 + t177 * t53;
t194 = 0.2e1 * t218 * t212;
t193 = -t133 * t265 + t122;
t192 = t262 * pkin(4) + t82;
t52 = t77 * pkin(4) + t73;
t191 = -t151 * t215 - t224;
t189 = -pkin(8) * t129 + t133 * t101;
t185 = -t150 * t49 + t197 * t241 - t240 * t66 - t255 * t50;
t9 = t49 * pkin(5) - t50 * qJ(6) - qJD(6) * t197 + t52;
t169 = -t177 * pkin(4) - pkin(5);
t166 = t175 * pkin(4) + qJ(6);
t105 = -pkin(5) * t255 - t150 * qJ(6) + t207;
t91 = t150 * t151;
t62 = t139 * t151 - t175 * t228 + t177 * t224;
t61 = -qJD(4) * t92 + t143 * t150;
t41 = t91 * pkin(5) - t92 * qJ(6) + t198;
t37 = t118 * pkin(4) + pkin(5) * t197 + qJ(6) * t66;
t32 = -t149 * pkin(5) - t33;
t31 = t149 * qJ(6) + t34;
t15 = t133 * qJ(6) + t20;
t14 = -t133 * pkin(5) + qJD(6) - t19;
t10 = -t61 * pkin(5) + t62 * qJ(6) - t92 * qJD(6) + t192;
t6 = -t144 * pkin(5) - t7;
t5 = t144 * qJ(6) + t149 * qJD(6) + t8;
t1 = t133 * qJD(6) + t210;
t11 = [0, 0, 0, 0, 0, t194, qJ(2) * t194, t128 * t151 - t143 * t254, -t128 * t149 - t151 * t129 - t143 * t140 - t144 * t254, -t143 * qJD(3), -t144 * qJD(3), 0, -t82 * qJD(3) + t170 * t129 + t155 * t144, -t81 * qJD(3) + t170 * t128 - t155 * t143, t118 * t191 + t76 * t234 -(-t116 * t181 - t118 * t179) * t143 + (-t244 - t181 * t77 + (t116 * t179 - t118 * t181) * qJD(4)) * t151, t118 * t144 + t151 * t122 + t133 * t191 + t76 * t149, -t116 * t144 - t262 * t133 - t77 * t149 - t151 * t230, t129 * t149 + t133 * t144 (-t113 * t214 + t96) * t133 + t98 * t129 + (-t102 * t214 + t86) * t149 + t58 * t144 + t82 * t116 - t256 * t77 + t101 * t206 + ((-qJD(4) * t106 - t81) * t133 - t113 * t129 + (-qJD(4) * t80 - t72) * t149 + t73 * t151 - t101 * t143) * t179 -(-t113 * t215 + t211) * t133 - t242 * t129 - t190 * t149 - t59 * t144 + t82 * t118 - t256 * t76 + t73 * t234 + t191 * t101, t19 * t62 - t197 * t7 + t20 * t61 - t3 * t92 - t33 * t50 - t34 * t49 - t4 * t91 - t8 * t66, t19 * t7 + t192 * t64 + t198 * t52 + t20 * t8 + t3 * t33 + t4 * t34, t10 * t66 - t32 * t129 - t6 * t133 - t14 * t144 - t2 * t149 - t35 * t61 + t41 * t49 + t9 * t91, -t1 * t91 - t14 * t62 + t15 * t61 + t197 * t6 + t2 * t92 - t31 * t49 + t32 * t50 - t5 * t66, t1 * t149 - t10 * t197 + t31 * t129 + t5 * t133 + t15 * t144 + t35 * t62 - t41 * t50 - t9 * t92, t1 * t31 + t35 * t10 + t14 * t6 + t15 * t5 + t2 * t32 + t9 * t41; 0, 0, 0, 0, 0, -t202, -qJ(2) * t202, 0, 0, 0, 0, 0, 0.2e1 * t261, t161 + (t140 - t205) * qJD(3), 0, 0, 0, 0, 0, t193 - t236, -t253 * t181 - t230 - t237, t185, t4 * t150 - t241 * t19 + t240 * t20 - t254 * t64 + t255 * t3, t129 * t255 - t241 * t133 - t254 * t66, t185, t150 * t129 + t133 * t240 + t197 * t254, t1 * t150 + t14 * t241 + t15 * t240 - t2 * t255 - t254 * t35; 0, 0, 0, 0, 0, 0, 0, -t254 * t140, -t140 ^ 2 + t254 ^ 2, t161 + (-t140 - t205) * qJD(3), 0, 0, t108 * qJD(3) - t155 * t254 - t73, -t155 * t140 + t149 * t212, t118 * t201 + t244 (t76 - t239) * t181 + (-t77 - t238) * t179, t133 * t201 + t230 - t237, t193 + t236, -t133 * t254, -pkin(3) * t77 - t108 * t116 - t58 * t254 - t73 * t181 + (-pkin(8) * t214 - t95) * t133 + (t133 * t258 + t189) * t179, -pkin(3) * t76 - t108 * t118 + t59 * t254 + t73 * t179 + (pkin(8) * t215 + t243) * t133 + t189 * t181, -t3 * t150 - t240 * t19 - t197 * t246 - t241 * t20 + t255 * t4 + t29 * t66 + t200, t4 * t115 - t3 * t114 + t52 * t207 + t257 * t64 + (t88 - t29) * t20 + t246 * t19, t105 * t49 - t114 * t129 + t248 * t133 + t14 * t254 + t241 * t35 - t247 * t66 - t255 * t9, t1 * t255 + t240 * t14 - t241 * t15 + t2 * t150 - t197 * t248 + t24 * t66 + t200, -t105 * t50 + t115 * t129 + t249 * t133 - t15 * t254 - t9 * t150 + t197 * t247 - t240 * t35, t1 * t115 + t9 * t105 + t2 * t114 - t248 * t14 + t249 * t15 - t247 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118 * t116, -t116 ^ 2 + t118 ^ 2, t76 + t239, t238 - t77, t129, -t101 * t118 + t59 * t133 + t184, t101 * t116 + t58 * t133 - t190, t20 * t197 - t252 + (-t175 * t49 - t177 * t50) * pkin(4) + (-t19 + t22) * t66, t19 * t21 - t20 * t22 + (-t118 * t64 + t175 * t4 + t177 * t3) * pkin(4), t21 * t133 - t259 - t37 * t66 + (pkin(5) - t169) * t129 + t3, t15 * t197 - t166 * t49 + t169 * t50 - t252 + (t14 - t219) * t66, t166 * t129 - t35 * t66 + t37 * t197 + (0.2e1 * qJD(6) - t22) * t133 + t210, t1 * t166 - t14 * t21 + t15 * t219 + t2 * t169 - t35 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, t19 * t197 + t20 * t66 + t52, t133 * t197 + t49, t199, -t50 + t264, -t14 * t197 + t15 * t66 + t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197 * t66 - t261, t50 + t264, -t253 - t260, -t15 * t133 + t2 + t259;];
tauc_reg  = t11;
