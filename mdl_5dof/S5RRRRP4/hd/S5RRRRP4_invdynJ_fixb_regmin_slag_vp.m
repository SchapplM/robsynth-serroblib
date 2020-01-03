% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP4
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
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:16
% EndTime: 2019-12-31 21:51:21
% DurationCPUTime: 1.69s
% Computational Cost: add. (3198->297), mult. (4731->353), div. (0->0), fcn. (3006->12), ass. (0->188)
t148 = qJ(1) + qJ(2);
t138 = cos(t148);
t151 = sin(qJ(2));
t242 = t151 * pkin(1);
t204 = qJD(1) * t242;
t154 = cos(qJ(2));
t240 = t154 * pkin(1);
t216 = -qJD(2) * t204 + qJDD(1) * t240;
t257 = -g(2) * t138 + t216;
t142 = qJDD(1) + qJDD(2);
t244 = t142 * pkin(2);
t256 = -t244 - t257;
t153 = cos(qJ(3));
t246 = cos(qJ(4));
t192 = t246 * qJD(4);
t198 = t246 * t153;
t255 = -qJD(3) * t198 - t153 * t192;
t144 = qJD(1) + qJD(2);
t149 = sin(qJ(4));
t150 = sin(qJ(3));
t222 = t149 * t150;
t201 = t144 * t222;
t68 = -t144 * t198 + t201;
t199 = t246 * t150;
t81 = t149 * t153 + t199;
t70 = t81 * t144;
t241 = t153 * pkin(3);
t127 = pkin(2) + t241;
t214 = qJD(1) * t154;
t203 = pkin(1) * t214;
t72 = -t127 * t144 - t203;
t26 = t68 * pkin(4) - t70 * qJ(5) + t72;
t143 = qJD(3) + qJD(4);
t178 = t143 * t222;
t53 = t178 + t255;
t220 = t153 * t142;
t188 = -t142 * t199 + t255 * t144 - t149 * t220;
t20 = t144 * t178 + t188;
t221 = t150 * t142;
t179 = -t142 * t198 + t149 * t221;
t54 = t143 * t81;
t21 = t54 * t144 + t179;
t211 = qJD(3) * t150;
t196 = t144 * t211;
t47 = pkin(3) * t196 - t127 * t142 - t216;
t7 = t21 * pkin(4) + t20 * qJ(5) - t70 * qJD(5) + t47;
t254 = t26 * t53 - t7 * t81;
t173 = t198 - t222;
t253 = -t173 * t7 + t26 * t54;
t252 = -t173 * t47 + t72 * t54;
t251 = t47 * t81 - t72 * t53;
t147 = qJ(3) + qJ(4);
t135 = sin(t147);
t137 = cos(t147);
t218 = t137 * pkin(4) + t135 * qJ(5);
t136 = sin(t148);
t217 = g(1) * t138 + g(2) * t136;
t141 = qJDD(3) + qJDD(4);
t130 = t141 * qJ(5);
t133 = t143 * qJD(5);
t250 = t130 + t133;
t134 = t141 * pkin(4);
t249 = qJDD(5) - t134;
t227 = t135 * t138;
t228 = t135 * t136;
t233 = g(1) * t228 - g(2) * t227;
t156 = -pkin(8) - pkin(7);
t108 = t156 * t150;
t139 = t153 * pkin(8);
t109 = t153 * pkin(7) + t139;
t174 = t246 * t108 - t149 * t109;
t200 = qJD(3) * t156;
t82 = t150 * t200;
t83 = t153 * t200;
t235 = t174 * qJD(4) + t149 * t83 - t173 * t203 + t246 * t82;
t60 = t149 * t108 + t246 * t109;
t248 = -t60 * t141 - t235 * t143 - t233;
t247 = t70 ^ 2;
t122 = g(1) * t136;
t243 = t144 * pkin(2);
t239 = t26 * t68;
t238 = t70 * t68;
t237 = t72 * t68;
t125 = pkin(7) + t242;
t236 = -pkin(8) - t125;
t234 = t60 * qJD(4) + t149 * t82 - t81 * t203 - t246 * t83;
t85 = t144 * pkin(7) + t204;
t194 = pkin(8) * t144 + t85;
t64 = t194 * t153;
t232 = t149 * t64;
t202 = t246 * t64;
t63 = t194 * t150;
t58 = qJD(3) * pkin(3) - t63;
t34 = t149 * t58 + t202;
t231 = t34 * t143;
t86 = -t203 - t243;
t230 = t153 * t122 + t86 * t211;
t37 = -t246 * t63 - t232;
t229 = pkin(3) * t192 + qJD(5) - t37;
t226 = t136 * t137;
t225 = t137 * t138;
t224 = t138 * t156;
t223 = t144 * t150;
t33 = t246 * t58 - t232;
t219 = qJD(5) - t33;
t145 = t150 ^ 2;
t215 = -t153 ^ 2 + t145;
t213 = qJD(2) * t151;
t212 = qJD(2) * t154;
t210 = qJD(3) * t153;
t209 = qJD(4) * t149;
t208 = qJDD(1) * t151;
t207 = pkin(4) * t225 + qJ(5) * t227 + t138 * t127;
t206 = t256 * t150 + t86 * t210;
t205 = pkin(1) * t212;
t131 = pkin(3) * t211;
t197 = t144 * t213;
t195 = t144 * t210;
t191 = qJD(3) * t236;
t74 = t142 * pkin(7) + (qJD(1) * t212 + t208) * pkin(1);
t30 = -t85 * t210 + qJDD(3) * pkin(3) - t150 * t74 + (-t195 - t221) * pkin(8);
t35 = -t85 * t211 + t153 * t74 + (-t196 + t220) * pkin(8);
t190 = t149 * t30 + t58 * t192 - t64 * t209 + t246 * t35;
t189 = t149 * t35 + t64 * t192 + t58 * t209 - t246 * t30;
t187 = t144 * t204;
t186 = g(1) * t226 - g(2) * t225;
t19 = t54 * pkin(4) + t53 * qJ(5) - t81 * qJD(5) + t131;
t182 = -t19 + t204;
t36 = -t149 * t63 + t202;
t181 = pkin(3) * t209 - t36;
t180 = -pkin(3) * t150 - pkin(4) * t135;
t45 = t70 * pkin(4) + t68 * qJ(5);
t22 = -t143 * pkin(4) + t219;
t25 = t143 * qJ(5) + t34;
t6 = t190 + t250;
t8 = t189 + t249;
t177 = t173 * t6 - t22 * t53 - t25 * t54 + t8 * t81 - t217;
t176 = t122 + t257;
t78 = t236 * t150;
t79 = t153 * t125 + t139;
t175 = -t149 * t79 + t246 * t78;
t49 = t149 * t78 + t246 * t79;
t172 = -g(1) * t225 - g(2) * t226 - g(3) * t135 + t190;
t61 = t150 * t191 + t153 * t205;
t62 = -t150 * t205 + t153 * t191;
t10 = t175 * qJD(4) + t149 * t62 + t246 * t61;
t171 = -t10 * t143 - t49 * t141 - t233;
t170 = -t204 + t131;
t169 = g(1) * t227 + g(2) * t228 - g(3) * t137 - t189;
t50 = -pkin(4) * t173 - t81 * qJ(5) - t127;
t168 = -t144 * t86 + t217 - t74;
t167 = t33 * t143 - t172;
t157 = qJD(3) ^ 2;
t166 = pkin(7) * t157 - t187 - t244;
t11 = t49 * qJD(4) + t149 * t61 - t246 * t62;
t165 = -t11 * t143 + t141 * t175 + t186;
t128 = -pkin(2) - t240;
t164 = pkin(1) * t197 + t125 * t157 + t128 * t142;
t163 = -t72 * t70 + t169;
t162 = t141 * t174 - t234 * t143 + t186;
t161 = -pkin(7) * qJDD(3) + (t203 - t243) * qJD(3);
t160 = -qJDD(3) * t125 + (t128 * t144 - t205) * qJD(3);
t159 = (-g(1) * (-t127 - t218) + g(2) * t156) * t136;
t158 = t26 * t70 - t169 + t249;
t155 = cos(qJ(1));
t152 = sin(qJ(1));
t140 = t144 ^ 2;
t132 = pkin(1) * t213;
t126 = -t246 * pkin(3) - pkin(4);
t118 = t149 * pkin(3) + qJ(5);
t103 = -t127 - t240;
t102 = qJDD(3) * t153 - t157 * t150;
t101 = qJDD(3) * t150 + t157 * t153;
t89 = qJ(5) * t225;
t87 = qJ(5) * t226;
t84 = t132 + t131;
t75 = t145 * t142 + 0.2e1 * t150 * t195;
t55 = -0.2e1 * t215 * t144 * qJD(3) + 0.2e1 * t150 * t220;
t46 = t50 - t240;
t40 = pkin(3) * t223 + t45;
t39 = t141 * t173 - t54 * t143;
t38 = t81 * t141 - t53 * t143;
t27 = -t68 ^ 2 + t247;
t17 = -t188 + (-t201 + t68) * t143;
t16 = t132 + t19;
t9 = -t20 * t81 - t70 * t53;
t1 = -t173 * t20 - t81 * t21 + t53 * t68 - t70 * t54;
t2 = [qJDD(1), g(1) * t152 - g(2) * t155, g(1) * t155 + g(2) * t152, t142, (t142 * t154 - t197) * pkin(1) + t176, ((-qJDD(1) - t142) * t151 + (-qJD(1) - t144) * t212) * pkin(1) + t217, t75, t55, t101, t102, 0, t160 * t150 + (-t164 - t256) * t153 + t230, t160 * t153 + (t164 - t122) * t150 + t206, t9, t1, t38, t39, 0, t103 * t21 + t84 * t68 + t165 + t252, -t103 * t20 + t84 * t70 + t171 + t251, t16 * t68 + t46 * t21 + t165 + t253, -t10 * t68 + t11 * t70 + t175 * t20 - t49 * t21 + t177, -t16 * t70 + t46 * t20 - t171 + t254, t6 * t49 + t25 * t10 + t7 * t46 + t26 * t16 - t8 * t175 + t22 * t11 - g(1) * (-t152 * pkin(1) - t224) - g(2) * (t155 * pkin(1) + t207) + t159; 0, 0, 0, t142, t176 + t187, (-t208 + (-qJD(2) + t144) * t214) * pkin(1) + t217, t75, t55, t101, t102, 0, t161 * t150 + (-t166 - t256) * t153 + t230, t161 * t153 + (t166 - t122) * t150 + t206, t9, t1, t38, t39, 0, -t127 * t21 + t170 * t68 + t162 + t252, t127 * t20 + t170 * t70 + t248 + t251, -t182 * t68 + t50 * t21 + t162 + t253, t174 * t20 - t60 * t21 + t234 * t70 - t235 * t68 + t177, t182 * t70 + t50 * t20 - t248 + t254, g(1) * t224 - g(2) * t207 - t174 * t8 - t182 * t26 + t234 * t22 + t235 * t25 + t7 * t50 + t6 * t60 + t159; 0, 0, 0, 0, 0, 0, -t150 * t140 * t153, t215 * t140, t221, t220, qJDD(3), -g(3) * t153 + t150 * t168, g(3) * t150 + t153 * t168, t238, t27, t17, -t179, t141, t36 * t143 + (t246 * t141 - t143 * t209 - t68 * t223) * pkin(3) + t163, t37 * t143 + t237 + (-t141 * t149 - t143 * t192 - t70 * t223) * pkin(3) - t172, -t126 * t141 - t143 * t181 - t40 * t68 - t158, -t118 * t21 - t126 * t20 + (t181 + t25) * t70 + (t22 - t229) * t68, t118 * t141 + t229 * t143 + t40 * t70 + t172 - t239 + t250, t6 * t118 + t8 * t126 - t26 * t40 - g(1) * (t138 * t180 + t89) - g(2) * (t136 * t180 + t87) - g(3) * (t218 + t241) + t229 * t25 + t181 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t238, t27, t17, -t179, t141, t163 + t231, t167 + t237, -t45 * t68 + t134 - t158 + t231, pkin(4) * t20 - t21 * qJ(5) + (t25 - t34) * t70 + (t22 - t219) * t68, t45 * t70 + 0.2e1 * t130 + 0.2e1 * t133 - t167 - t239, t6 * qJ(5) - t8 * pkin(4) - t26 * t45 - t22 * t34 - g(1) * (-pkin(4) * t227 + t89) - g(2) * (-pkin(4) * t228 + t87) - g(3) * t218 + t219 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141 + t238, t17, -t143 ^ 2 - t247, -t25 * t143 + t158;];
tau_reg = t2;
