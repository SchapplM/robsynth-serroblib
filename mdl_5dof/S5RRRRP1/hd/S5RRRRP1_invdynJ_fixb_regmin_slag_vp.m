% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:03
% EndTime: 2019-12-05 18:46:09
% DurationCPUTime: 2.30s
% Computational Cost: add. (3545->294), mult. (8497->390), div. (0->0), fcn. (6172->12), ass. (0->177)
t156 = sin(qJ(4));
t160 = cos(qJ(4));
t161 = cos(qJ(3));
t157 = sin(qJ(3));
t158 = sin(qJ(2));
t229 = qJD(1) * t158;
t215 = t157 * t229;
t162 = cos(qJ(2));
t228 = qJD(1) * t162;
t93 = -t161 * t228 + t215;
t95 = -t157 * t228 - t161 * t229;
t194 = t156 * t93 + t160 * t95;
t254 = t194 ^ 2;
t62 = t156 * t95 - t160 * t93;
t58 = t62 ^ 2;
t12 = -t58 + t254;
t155 = qJ(2) + qJ(3);
t149 = qJ(4) + t155;
t137 = sin(t149);
t138 = cos(t149);
t159 = sin(qJ(1));
t163 = cos(qJ(1));
t197 = g(1) * t163 + g(2) * t159;
t263 = -g(3) * t138 + t137 * t197;
t150 = t162 * pkin(2);
t259 = -pkin(1) - t150;
t237 = qJ(5) * t62;
t245 = t194 * t62;
t55 = t194 * qJ(5);
t223 = qJD(1) * qJD(2);
t212 = t162 * t223;
t221 = t158 * qJDD(1);
t187 = -t212 - t221;
t222 = qJD(1) * qJD(3);
t261 = -t162 * t222 + t187;
t117 = t259 * qJD(1);
t260 = qJDD(1) * t259;
t213 = t158 * t223;
t220 = t162 * qJDD(1);
t186 = -t213 + t220;
t253 = pkin(6) + pkin(7);
t118 = t253 * t158;
t109 = qJD(1) * t118;
t236 = qJD(2) * pkin(2);
t102 = -t109 + t236;
t76 = t253 * t186;
t203 = -qJD(3) * t102 - t76;
t75 = qJDD(2) * pkin(2) + t253 * t187;
t119 = t253 * t162;
t111 = qJD(1) * t119;
t227 = qJD(3) * t157;
t91 = t111 * t227;
t13 = -t91 + (t261 * pkin(8) + t75) * t157 + ((-t158 * t222 + t186) * pkin(8) - t203) * t161;
t225 = qJD(4) * t156;
t152 = qJD(2) + qJD(3);
t96 = t157 * t111;
t204 = t161 * t102 - t96;
t88 = t95 * pkin(8);
t45 = t204 + t88;
t34 = pkin(3) * t152 + t45;
t100 = t161 * t111;
t193 = -t157 * t102 - t100;
t251 = pkin(8) * t93;
t46 = -t193 - t251;
t151 = qJDD(2) + qJDD(3);
t176 = qJD(3) * t193 - t157 * t76 + t161 * t75;
t47 = -t152 * t215 + t157 * t220 - t261 * t161;
t9 = pkin(3) * t151 - pkin(8) * t47 + t176;
t255 = (qJD(4) * t34 + t13) * t160 + t156 * t9 - t46 * t225;
t77 = pkin(3) * t93 + t117;
t174 = g(3) * t137 + t197 * t138 - t77 * t62 - t255;
t145 = qJD(4) + t152;
t105 = t157 * t162 + t158 * t161;
t185 = t105 * qJD(2);
t171 = -qJD(3) * t105 - t185;
t170 = t171 * qJD(1);
t192 = t157 * t158 - t161 * t162;
t167 = -t192 * qJDD(1) + t170;
t224 = qJD(4) * t160;
t189 = t156 * t167 + t160 * t47 - t93 * t224 + t95 * t225;
t5 = -t145 * t62 + t189;
t38 = t160 * t46;
t195 = -t156 * t34 - t38;
t179 = qJD(4) * t195 - t156 * t13 + t160 * t9;
t169 = t194 * t77 + t179 + t263;
t178 = qJD(4) * t194 - t156 * t47 + t160 * t167;
t6 = -t145 * t194 + t178;
t140 = pkin(2) * t161 + pkin(3);
t233 = t157 * t160;
t202 = t109 * t157 - t100;
t49 = t202 + t251;
t238 = -t161 * t109 - t96;
t50 = t88 + t238;
t258 = t156 * t50 - t160 * t49 - t140 * t225 + (-t157 * t224 + (-t156 * t161 - t233) * qJD(3)) * pkin(2);
t235 = t156 * t157;
t257 = -t140 * t224 - (-t157 * t225 + (t160 * t161 - t235) * qJD(3)) * pkin(2) + t156 * t49 + t160 * t50;
t232 = -t157 * t118 + t161 * t119;
t80 = pkin(3) * t192 + t259;
t256 = -pkin(8) * t105 - t157 * t119;
t252 = pkin(4) * t194;
t36 = t156 * t46;
t210 = t160 * t34 - t36;
t10 = t210 + t55;
t7 = pkin(4) * t145 + t10;
t250 = t10 - t7;
t244 = t95 * t93;
t243 = t160 * t45 - t36;
t107 = t161 * t118;
t56 = -t107 + t256;
t57 = -pkin(8) * t192 + t232;
t241 = t156 * t56 + t160 * t57;
t240 = -t55 - t257;
t239 = t237 + t258;
t147 = cos(t155);
t231 = pkin(3) * t147 + pkin(4) * t138;
t153 = t158 ^ 2;
t230 = -t162 ^ 2 + t153;
t226 = qJD(3) * t161;
t216 = qJD(2) * t253;
t110 = t158 * t216;
t112 = t162 * t216;
t219 = -t161 * t110 - t157 * t112 - t118 * t226;
t143 = t158 * t236;
t218 = t150 + t231;
t78 = pkin(2) * t229 - pkin(3) * t95;
t209 = -t156 * t45 - t38;
t206 = -t156 * t57 + t160 * t56;
t134 = pkin(2) * t213;
t89 = t134 + t260;
t201 = t89 + t260;
t200 = -pkin(2) * t235 + t160 * t140;
t11 = -t195 + t237;
t199 = -t11 * t194 + t62 * t7;
t146 = sin(t155);
t198 = -pkin(3) * t146 - pkin(4) * t137;
t196 = g(1) * t159 - g(2) * t163;
t191 = -0.2e1 * pkin(1) * t223 - pkin(6) * qJDD(2);
t27 = -pkin(8) * t185 + t256 * qJD(3) + t219;
t101 = t161 * t112;
t74 = t152 * t192;
t28 = pkin(8) * t74 - t232 * qJD(3) + t110 * t157 - t101;
t190 = t156 * t28 + t160 * t27 + t56 * t224 - t57 * t225;
t188 = t160 * t192;
t164 = qJD(2) ^ 2;
t182 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t164 + t196;
t165 = qJD(1) ^ 2;
t181 = pkin(1) * t165 - pkin(6) * qJDD(1) + t197;
t73 = t160 * t105 - t156 * t192;
t180 = 0.2e1 * t117 * t152;
t177 = -t241 * qJD(4) - t156 * t27 + t160 * t28;
t173 = g(3) * t146 + t117 * t93 + t197 * t147 - t157 * t75 + t161 * t203 + t91;
t168 = -g(3) * t147 + t117 * t95 + t197 * t146 + t176;
t65 = -pkin(3) * t171 + t143;
t31 = -pkin(3) * t170 + t80 * qJDD(1) + t134;
t166 = -pkin(4) * t178 + qJDD(5) + t31;
t148 = -qJ(5) - pkin(8) - t253;
t144 = qJDD(4) + t151;
t139 = pkin(3) * t160 + pkin(4);
t90 = pkin(2) * t233 + t140 * t156;
t86 = pkin(4) + t200;
t81 = pkin(1) + t218;
t72 = t105 * t156 + t188;
t48 = -t93 ^ 2 + t95 ^ 2;
t32 = -pkin(4) * t62 + qJD(5) + t77;
t30 = -t152 * t95 + t167;
t29 = t152 * t93 + t47;
t23 = qJD(4) * t73 - t156 * t74 - t160 * t171;
t22 = qJD(4) * t188 + t105 * t225 - t156 * t171 + t160 * t74;
t21 = -qJ(5) * t72 + t241;
t20 = -qJ(5) * t73 + t206;
t15 = t55 + t243;
t14 = t209 - t237;
t4 = qJ(5) * t22 - qJD(5) * t73 + t177;
t3 = -qJ(5) * t23 - qJD(5) * t72 + t190;
t2 = qJ(5) * t178 + qJD(5) * t62 + t255;
t1 = pkin(4) * t144 - qJ(5) * t189 + qJD(5) * t194 + t179;
t8 = [qJDD(1), t196, t197, qJDD(1) * t153 + 0.2e1 * t158 * t212, 0.2e1 * t158 * t220 - 0.2e1 * t230 * t223, qJDD(2) * t158 + t162 * t164, qJDD(2) * t162 - t158 * t164, 0, t158 * t191 + t162 * t182, -t158 * t182 + t162 * t191, t105 * t47 + t74 * t95, t167 * t105 - t171 * t95 - t192 * t47 + t74 * t93, t105 * t151 - t152 * t74, -t151 * t192 + t152 * t171, 0, t93 * t143 - t101 * t152 - t107 * t151 + t196 * t147 + (-qJD(3) * t119 * t152 + t158 * t180 - t162 * t201) * t161 + ((qJD(3) * t118 + t110) * t152 - t119 * t151 + t201 * t158 + t180 * t162) * t157, -t95 * t143 + t259 * t47 + t89 * t105 - t117 * t74 - (-t119 * t227 + t219) * t152 - t232 * t151 - t196 * t146, t189 * t73 + t194 * t22, t178 * t73 - t189 * t72 + t194 * t23 - t22 * t62, t144 * t73 - t145 * t22, -t144 * t72 - t145 * t23, 0, t196 * t138 + t206 * t144 + t177 * t145 - t178 * t80 + t77 * t23 + t31 * t72 - t62 * t65, -t196 * t137 - t241 * t144 - t190 * t145 + t189 * t80 - t194 * t65 - t77 * t22 + t31 * t73, -t1 * t73 - t11 * t23 + t178 * t21 - t189 * t20 + t194 * t4 - t2 * t72 + t22 * t7 + t3 * t62 - t197, t2 * t21 + t11 * t3 + t1 * t20 + t7 * t4 + t166 * (t72 * pkin(4) + t80) + t32 * (t23 * pkin(4) + t65) - g(1) * (-t148 * t163 - t159 * t81) - g(2) * (-t148 * t159 + t163 * t81); 0, 0, 0, -t158 * t165 * t162, t230 * t165, t221, t220, qJDD(2), -g(3) * t162 + t158 * t181, g(3) * t158 + t162 * t181, -t244, t48, t29, t30, t151, -t202 * t152 + (t161 * t151 - t152 * t227 - t93 * t229) * pkin(2) + t168, t238 * t152 + (-t157 * t151 - t152 * t226 + t95 * t229) * pkin(2) + t173, t245, t12, t5, t6, t144, t200 * t144 + t145 * t258 + t62 * t78 + t169, -t90 * t144 + t145 * t257 + t194 * t78 + t174, t178 * t90 - t189 * t86 + t194 * t239 + t240 * t62 + t199, t2 * t90 + t1 * t86 - t32 * (t78 - t252) - g(3) * t218 - t197 * (-pkin(2) * t158 + t198) + t239 * t7 + t240 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, t48, t29, t30, t151, -t152 * t193 + t168, t152 * t204 + t173, t245, t12, t5, t6, t144, -t209 * t145 + (t160 * t144 - t145 * t225 - t62 * t95) * pkin(3) + t169, t243 * t145 + (-t156 * t144 - t145 * t224 - t194 * t95) * pkin(3) + t174, -t139 * t189 - t14 * t194 - t15 * t62 + (t156 * t178 + (-t156 * t194 + t160 * t62) * qJD(4)) * pkin(3) + t199, t1 * t139 - t11 * t15 - t7 * t14 + t32 * t252 - g(3) * t231 - t197 * t198 + (t2 * t156 + t32 * t95 + (t11 * t160 - t156 * t7) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, t12, t5, t6, t144, -t145 * t195 + t169, t145 * t210 + t174, -pkin(4) * t189 - t250 * t62, -t250 * t11 + (t194 * t32 + t1 + t263) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58 - t254, -t11 * t62 - t194 * t7 + t166 - t196;];
tau_reg = t8;
