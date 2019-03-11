% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:08:18
% EndTime: 2019-03-09 08:08:26
% DurationCPUTime: 2.75s
% Computational Cost: add. (3824->355), mult. (10083->486), div. (0->0), fcn. (7398->8), ass. (0->183)
t173 = cos(qJ(2));
t239 = cos(pkin(9));
t212 = t239 * t173;
t156 = qJD(1) * t212;
t168 = sin(pkin(9));
t171 = sin(qJ(2));
t224 = qJD(1) * t171;
t128 = t168 * t224 - t156;
t220 = qJD(6) - t128;
t264 = qJD(6) - t220;
t170 = sin(qJ(6));
t172 = cos(qJ(6));
t147 = t168 * t173 + t239 * t171;
t131 = t147 * qJD(1);
t167 = sin(pkin(10));
t169 = cos(pkin(10));
t191 = t167 * qJD(2) + t169 * t131;
t97 = t169 * qJD(2) - t167 * t131;
t51 = t170 * t191 + t172 * t97;
t263 = t220 * t51;
t53 = -t170 * t97 + t172 * t191;
t262 = t220 * t53;
t261 = t97 * t128;
t148 = t172 * t167 - t170 * t169;
t221 = qJD(6) * t172;
t222 = qJD(6) * t170;
t240 = -t148 * t128 + t167 * t221 - t169 * t222;
t73 = t148 * t147;
t146 = t170 * t167 + t172 * t169;
t260 = t191 ^ 2;
t219 = qJD(1) * qJD(2);
t259 = -0.2e1 * t219;
t130 = t147 * qJD(2);
t116 = qJD(1) * t130;
t215 = t171 * t219;
t117 = qJD(2) * t156 - t168 * t215;
t158 = pkin(2) * t215;
t48 = t116 * pkin(3) - t117 * qJ(4) - t131 * qJD(4) + t158;
t246 = -qJ(3) - pkin(7);
t214 = qJD(2) * t246;
t125 = t173 * qJD(3) + t171 * t214;
t109 = t125 * qJD(1);
t126 = -t171 * qJD(3) + t173 * t214;
t110 = t126 * qJD(1);
t59 = t239 * t109 + t168 * t110;
t57 = qJD(2) * qJD(4) + t59;
t54 = t167 * t57;
t15 = t169 * t48 - t54;
t10 = -t116 * pkin(4) - t15;
t217 = -t173 * pkin(2) - pkin(1);
t201 = t217 * qJD(1);
t152 = qJD(3) + t201;
t61 = t128 * pkin(3) - t131 * qJ(4) + t152;
t153 = t246 * t171;
t150 = qJD(1) * t153;
t245 = qJD(2) * pkin(2);
t143 = t150 + t245;
t154 = t246 * t173;
t151 = qJD(1) * t154;
t213 = t239 * t151;
t81 = t168 * t143 - t213;
t76 = qJD(2) * qJ(4) + t81;
t36 = t167 * t61 + t169 * t76;
t26 = t128 * qJ(5) + t36;
t258 = -t26 * t128 + t10;
t25 = qJD(6) * t53 - t148 * t117;
t241 = t220 * t146;
t255 = t148 * t116 + t220 * t241;
t127 = t128 ^ 2;
t254 = t169 * t116 - t167 * t127 + t131 * t97;
t253 = t167 * t116 + t169 * t127 + t131 * t191;
t161 = -t239 * pkin(2) - pkin(3);
t184 = -t167 * qJ(5) + t161;
t139 = -t169 * pkin(4) + t184;
t159 = t168 * pkin(2) + qJ(4);
t236 = t116 * t159;
t136 = t168 * t151;
t80 = t239 * t143 + t136;
t190 = qJD(2) * pkin(3) - qJD(4) + t80;
t179 = qJ(5) * t191 + t190;
t34 = -pkin(4) * t97 - t179;
t252 = -t117 * t139 + (qJD(4) - t34) * t128 + t236;
t238 = qJ(5) * t169;
t250 = pkin(4) + pkin(5);
t251 = t250 * t167 - t238;
t249 = pkin(8) * t167;
t58 = t168 * t109 - t239 * t110;
t91 = -t239 * t153 - t168 * t154;
t248 = t58 * t91;
t247 = -pkin(8) + t159;
t16 = t167 * t48 + t169 * t57;
t185 = -t168 * t171 + t212;
t133 = t185 * qJD(2);
t218 = t171 * t245;
t56 = t130 * pkin(3) - t133 * qJ(4) - t147 * qJD(4) + t218;
t71 = t239 * t125 + t168 * t126;
t28 = t167 * t56 + t169 * t71;
t69 = pkin(2) * t224 + t131 * pkin(3) + t128 * qJ(4);
t85 = t239 * t150 + t136;
t42 = t167 * t69 + t169 * t85;
t79 = -pkin(3) * t185 - t147 * qJ(4) + t217;
t92 = t168 * t153 - t239 * t154;
t46 = t167 * t79 + t169 * t92;
t244 = t131 * t51;
t242 = t53 * t131;
t237 = qJD(4) * t191;
t105 = t167 * t117;
t106 = t169 * t117;
t175 = qJD(1) ^ 2;
t230 = t173 * t175;
t174 = qJD(2) ^ 2;
t229 = t174 * t171;
t228 = t174 * t173;
t227 = t191 * qJD(5);
t225 = t171 ^ 2 - t173 ^ 2;
t223 = qJD(4) * t169;
t31 = t131 * qJ(5) + t42;
t39 = -qJ(5) * t185 + t46;
t4 = t54 + (-pkin(8) * t117 - t48) * t169 - t250 * t116;
t6 = t116 * qJ(5) + t128 * qJD(5) + t16;
t5 = pkin(8) * t105 + t6;
t216 = -t170 * t5 + t172 * t4;
t63 = t167 * t71;
t27 = t169 * t56 - t63;
t35 = -t167 * t76 + t169 * t61;
t77 = t167 * t85;
t41 = t169 * t69 - t77;
t87 = t167 * t92;
t45 = t169 * t79 - t87;
t211 = t128 * t191 + t105;
t210 = t106 + t261;
t209 = pkin(1) * t259;
t84 = t168 * t150 - t213;
t208 = t167 * qJD(5) - t128 * t251 + t84;
t70 = t168 * t125 - t239 * t126;
t207 = t220 ^ 2;
t12 = t130 * qJ(5) - qJD(5) * t185 + t28;
t206 = t146 * t116 - t240 * t220;
t204 = qJD(5) - t35;
t203 = t170 * t4 + t172 * t5;
t202 = -t97 ^ 2 - t260;
t14 = -pkin(8) * t97 + t26;
t9 = -pkin(8) * t191 - t250 * t128 + t204;
t2 = t172 * t14 + t170 * t9;
t200 = t170 * t14 - t172 * t9;
t199 = pkin(4) * t167 - t238;
t198 = t91 * t117 + t58 * t147;
t197 = -t167 * t35 + t169 * t36;
t22 = t87 + (-pkin(8) * t147 - t79) * t169 + t250 * t185;
t29 = t147 * t249 + t39;
t196 = -t170 * t29 + t172 * t22;
t195 = t170 * t22 + t172 * t29;
t141 = t247 * t169;
t188 = qJD(4) * t167 - qJD(6) * t141 - t77 - (pkin(8) * t128 - t69) * t169 + t250 * t131;
t140 = t247 * t167;
t187 = qJD(6) * t140 + t128 * t249 + t223 - t31;
t186 = -t146 * t117 + t191 * t222 + t97 * t221;
t74 = t146 * t147;
t183 = -pkin(4) * t105 + t227 - t58;
t182 = -t169 * t147 * qJD(5) + t70;
t17 = -qJ(5) * t106 - t183;
t49 = t199 * t147 + t91;
t181 = t117 * t49 + t133 * t34 + t147 * t17;
t180 = -t133 * t190 + t198;
t178 = -t236 + t117 * t161 + (-qJD(4) - t190) * t128;
t177 = (t167 * t191 + t169 * t97) * t128 + (-t167 ^ 2 - t169 ^ 2) * t117;
t113 = t250 * t169 - t184;
t86 = t97 * t223;
t44 = -t199 * t128 + t84;
t43 = -t147 * t251 - t91;
t40 = pkin(4) * t185 - t45;
t38 = qJD(6) * t74 - t148 * t133;
t37 = qJD(6) * t73 + t133 * t146;
t32 = -t131 * pkin(4) - t41;
t30 = t199 * t133 + t182;
t23 = -t128 * pkin(4) + t204;
t20 = t250 * t97 + t179;
t19 = -t133 * t251 - t182;
t18 = -t130 * pkin(4) - t27;
t11 = (-pkin(5) * t167 + t238) * t117 + t183;
t8 = t133 * t249 + t12;
t7 = t63 + (-pkin(8) * t133 - t56) * t169 - t250 * t130;
t1 = [0, 0, 0, 0.2e1 * t173 * t215, t225 * t259, t228, -t229, 0, -pkin(7) * t228 + t171 * t209, pkin(7) * t229 + t173 * t209, -t92 * t116 - t71 * t128 - t81 * t130 + t70 * t131 - t80 * t133 + t185 * t59 + t198, t248 + t59 * t92 - t80 * t70 + t81 * t71 + (t152 + t201) * t218, t45 * t116 + t27 * t128 + t35 * t130 - t15 * t185 + t167 * t180 - t70 * t97, -t46 * t116 - t28 * t128 - t36 * t130 + t16 * t185 + t169 * t180 + t191 * t70, -t27 * t191 + t28 * t97 + (-t117 * t45 - t133 * t35 - t147 * t15) * t169 + (-t117 * t46 - t133 * t36 - t147 * t16) * t167, t15 * t45 + t16 * t46 - t190 * t70 + t35 * t27 + t36 * t28 + t248, t10 * t185 - t40 * t116 - t18 * t128 - t23 * t130 + t167 * t181 - t30 * t97, t12 * t97 + t18 * t191 + (t10 * t147 + t117 * t40 + t133 * t23) * t169 + (-t117 * t39 - t133 * t26 - t147 * t6) * t167, t39 * t116 + t12 * t128 + t26 * t130 - t169 * t181 - t185 * t6 - t191 * t30, t10 * t40 + t26 * t12 + t17 * t49 + t23 * t18 + t34 * t30 + t6 * t39, -t186 * t74 + t53 * t37, -t186 * t73 - t74 * t25 - t37 * t51 - t53 * t38, -t74 * t116 - t53 * t130 - t185 * t186 + t220 * t37, -t73 * t116 + t51 * t130 - t185 * t25 - t220 * t38, -t116 * t185 - t130 * t220 (-t170 * t8 + t172 * t7) * t220 - t196 * t116 + t216 * t185 + t200 * t130 + t19 * t51 + t43 * t25 - t11 * t73 + t20 * t38 + (-t185 * t2 - t195 * t220) * qJD(6) -(t170 * t7 + t172 * t8) * t220 + t195 * t116 - t203 * t185 + t2 * t130 + t19 * t53 - t43 * t186 + t11 * t74 + t20 * t37 + (t185 * t200 - t196 * t220) * qJD(6); 0, 0, 0, -t171 * t230, t225 * t175, 0, 0, 0, t175 * pkin(1) * t171, pkin(1) * t230 (t81 - t84) * t131 + (-t80 + t85) * t128 + (-t116 * t168 - t239 * t117) * pkin(2), t80 * t84 - t81 * t85 + (-t152 * t224 + t168 * t59 - t239 * t58) * pkin(2), -t41 * t128 - t35 * t131 + t167 * t178 - t58 * t169 + t84 * t97, t42 * t128 + t36 * t131 + t58 * t167 + t169 * t178 - t191 * t84, t41 * t191 - t42 * t97 + t86 + (-t128 * t35 + t16) * t169 + (-t128 * t36 - t15 + t237) * t167, t58 * t161 - t35 * t41 - t36 * t42 + t190 * t84 + (-t15 * t167 + t16 * t169) * t159 + t197 * qJD(4), t32 * t128 + t23 * t131 - t17 * t169 + t44 * t97 + (qJD(5) * t97 - t252) * t167, -t31 * t97 - t32 * t191 + t86 + (t128 * t23 + t6) * t169 + (t237 + t258) * t167, -t31 * t128 - t26 * t131 + t44 * t191 + (-t17 + t227) * t167 + t252 * t169, t17 * t139 - t23 * t32 - t26 * t31 - t34 * t44 + (qJD(4) * t26 + t159 * t6) * t169 + (qJD(4) * t23 - qJD(5) * t34 + t10 * t159) * t167, -t148 * t186 - t241 * t53, t146 * t186 - t148 * t25 - t240 * t53 + t241 * t51, t242 - t255, t206 - t244, t220 * t131 -(t172 * t140 - t170 * t141) * t116 + t113 * t25 + t11 * t146 - t200 * t131 + t208 * t51 + t240 * t20 - (t170 * t187 - t172 * t188) * t220 (t170 * t140 + t172 * t141) * t116 - t113 * t186 + t11 * t148 - t2 * t131 + t208 * t53 - t241 * t20 - (t170 * t188 + t172 * t187) * t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131 ^ 2 - t127, t81 * t128 + t80 * t131 + t158, t254, -t253, t177, t128 * t197 + t131 * t190 + t15 * t169 + t16 * t167, t254, t177, t253, -t10 * t169 - t34 * t131 + t6 * t167 + (t167 * t23 + t169 * t26) * t128, 0, 0, 0, 0, 0, t206 + t244, t242 + t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t210, t202, t191 * t35 - t36 * t97 + t58, t211, t202, -t210, -t191 * t23 - t26 * t97 + t17, 0, 0, 0, 0, 0, -t25 - t262, t186 + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t131 - t191 * t97, t106 - t261, -t127 - t260, t191 * t34 + t258, 0, 0, 0, 0, 0, -t172 * t116 - t170 * t207 - t191 * t51, t170 * t116 - t172 * t207 - t191 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t51, -t51 ^ 2 + t53 ^ 2, -t186 + t263, -t25 + t262, -t116, -t2 * t264 - t20 * t53 + t216, t20 * t51 + t200 * t264 - t203;];
tauc_reg  = t1;
