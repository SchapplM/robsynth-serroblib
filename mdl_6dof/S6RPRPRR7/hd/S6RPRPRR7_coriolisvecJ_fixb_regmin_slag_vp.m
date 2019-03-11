% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:03
% EndTime: 2019-03-09 03:56:11
% DurationCPUTime: 2.49s
% Computational Cost: add. (4066->274), mult. (9166->384), div. (0->0), fcn. (6918->8), ass. (0->168)
t155 = cos(qJ(6));
t197 = qJD(6) * t155;
t150 = sin(pkin(10));
t151 = cos(pkin(10));
t154 = sin(qJ(3));
t157 = cos(qJ(3));
t172 = t150 * t157 + t151 * t154;
t117 = t172 * qJD(1);
t156 = cos(qJ(5));
t104 = t156 * t117;
t202 = qJD(1) * t157;
t203 = qJD(1) * t154;
t120 = -t150 * t203 + t151 * t202;
t153 = sin(qJ(5));
t68 = t153 * t120 + t104;
t254 = t155 * t68;
t264 = t197 + t254;
t145 = qJD(3) + qJD(5);
t200 = qJD(3) * t157;
t201 = qJD(3) * t154;
t118 = t150 * t201 - t151 * t200;
t119 = t172 * qJD(3);
t171 = t150 * t154 - t151 * t157;
t173 = -t153 * t171 + t156 * t172;
t42 = qJD(5) * t173 - t153 * t118 + t156 * t119;
t263 = t42 * t145;
t219 = t68 * t145;
t192 = qJD(1) * qJD(3);
t108 = t171 * t192;
t109 = t172 * t192;
t199 = qJD(5) * t153;
t32 = -qJD(5) * t104 + t153 * t108 - t156 * t109 - t120 * t199;
t262 = t32 + t219;
t251 = -qJD(6) - t68;
t261 = qJD(6) + t251;
t174 = -t153 * t117 + t156 * t120;
t152 = sin(qJ(6));
t198 = qJD(6) * t152;
t18 = t145 * t197 + t155 * t32 - t174 * t198;
t56 = t152 * t145 + t155 * t174;
t74 = -t153 * t172 - t156 * t171;
t260 = t18 * t74 - t42 * t56;
t19 = qJD(6) * t56 + t152 * t32;
t54 = -t155 * t145 + t152 * t174;
t259 = -t152 * t19 + t18 * t155 - t264 * t54;
t14 = t18 * t152;
t258 = t264 * t56 + t14;
t33 = qJD(5) * t174 - t156 * t108 - t153 * t109;
t29 = t152 * t33;
t59 = t251 * t197;
t226 = t29 - t59;
t230 = t56 * t174;
t257 = -t251 * t254 + t226 - t230;
t233 = t120 * pkin(8);
t158 = -pkin(1) - pkin(7);
t132 = t158 * qJD(1) + qJD(2);
t114 = -qJ(4) * t203 + t154 * t132;
t95 = t150 * t114;
t115 = -qJ(4) * t202 + t157 * t132;
t99 = qJD(3) * pkin(3) + t115;
t57 = t151 * t99 - t95;
t44 = qJD(3) * pkin(4) - t233 + t57;
t234 = t117 * pkin(8);
t215 = t151 * t114;
t58 = t150 * t99 + t215;
t45 = t58 - t234;
t20 = -t153 * t45 + t156 * t44;
t16 = -t145 * pkin(5) - t20;
t256 = t16 * t68;
t253 = t174 * t68;
t252 = t152 * t251;
t220 = t174 * t145;
t250 = -t33 + t220;
t248 = t174 ^ 2 - t68 ^ 2;
t36 = pkin(5) * t174 + t68 * pkin(9);
t194 = t157 * qJD(4);
t83 = -t132 * t201 + (qJ(4) * t201 - t194) * qJD(1);
t195 = t154 * qJD(4);
t84 = t132 * t200 + (-qJ(4) * t200 - t195) * qJD(1);
t48 = -t150 * t84 + t151 * t83;
t37 = t109 * pkin(8) + t48;
t49 = t150 * t83 + t151 * t84;
t38 = t108 * pkin(8) + t49;
t2 = (qJD(5) * t44 + t38) * t156 + t153 * t37 - t45 * t199;
t128 = pkin(3) * t203 + qJD(1) * qJ(2) + qJD(4);
t82 = t117 * pkin(4) + t128;
t247 = t82 * t68 - t2;
t229 = t174 * t54;
t244 = t251 * t174;
t31 = t155 * t33;
t243 = -t198 * t251 - t31;
t161 = qJD(5) * t74 - t156 * t118 - t153 * t119;
t242 = t161 * t145;
t21 = t153 * t44 + t156 * t45;
t17 = t145 * pkin(9) + t21;
t24 = t68 * pkin(5) - pkin(9) * t174 + t82;
t175 = t152 * t17 - t155 * t24;
t241 = t16 * t198 + t174 * t175;
t3 = qJD(5) * t21 + t153 * t38 - t156 * t37;
t5 = t152 * t24 + t155 * t17;
t240 = t3 * t152 + t16 * t197 + t5 * t174;
t239 = -t174 * t82 - t3;
t209 = qJ(4) - t158;
t129 = t209 * t154;
t130 = t209 * t157;
t80 = t150 * t129 - t151 * t130;
t62 = pkin(8) * t171 + t80;
t81 = -t151 * t129 - t150 * t130;
t63 = -pkin(8) * t172 + t81;
t27 = t153 * t62 + t156 * t63;
t210 = t154 * pkin(3) + qJ(2);
t100 = pkin(4) * t172 + t210;
t28 = pkin(5) * t173 - t74 * pkin(9) + t100;
t26 = t153 * t63 - t156 * t62;
t110 = t209 * t201 - t194;
t111 = -qJD(3) * t130 - t195;
t60 = t151 * t110 - t150 * t111;
t46 = t119 * pkin(8) + t60;
t61 = t150 * t110 + t151 * t111;
t47 = t118 * pkin(8) + t61;
t6 = -qJD(5) * t26 + t153 * t46 + t156 * t47;
t238 = -t16 * t42 - t173 * (qJD(6) * t24 + t2) + (qJD(6) * t28 + t6) * t251 - t27 * t33 + t3 * t74;
t146 = qJD(1) * qJD(2);
t237 = 0.2e1 * t146;
t235 = pkin(3) * t150;
t232 = t16 * t74;
t231 = t28 * t33;
t228 = t74 * t33;
t65 = t151 * t115 - t95;
t224 = t152 * t56;
t139 = t151 * pkin(3) + pkin(4);
t168 = t156 * t139 - t153 * t235;
t64 = -t150 * t115 - t215;
t50 = t64 + t234;
t51 = t65 - t233;
t218 = -t168 * qJD(5) + t153 * t50 + t156 * t51;
t169 = t153 * t139 + t156 * t235;
t217 = t169 * qJD(5) - t153 * t51 + t156 * t50;
t159 = qJD(3) ^ 2;
t214 = t159 * t154;
t213 = t159 * t157;
t160 = qJD(1) ^ 2;
t212 = t160 * qJ(2);
t211 = t160 * t157;
t187 = t157 * t192;
t206 = pkin(3) * t187 + t146;
t205 = t154 ^ 2 - t157 ^ 2;
t204 = -t159 - t160;
t196 = t128 * qJD(1);
t193 = pkin(3) * t200 + qJD(2);
t190 = 0.2e1 * qJD(1);
t188 = t74 * t198;
t88 = pkin(3) * t202 + t120 * pkin(4);
t113 = pkin(9) + t169;
t179 = qJD(6) * t113 + t36 + t88;
t177 = qJD(6) * t173 + qJD(1);
t76 = -t108 * pkin(4) + t206;
t85 = -t118 * pkin(4) + t193;
t176 = t251 * t42 + t228;
t170 = t252 * t68 - t243;
t166 = -t113 * t33 - t218 * t251 + t256;
t164 = -t58 * t118 - t57 * t119 - t171 * t48 + t172 * t49;
t112 = -pkin(5) - t168;
t10 = pkin(5) * t161 + pkin(9) * t42 + t85;
t9 = t33 * pkin(5) - t32 * pkin(9) + t76;
t8 = t155 * t9;
t7 = qJD(5) * t27 + t153 * t47 - t156 * t46;
t1 = [0, 0, 0, 0, t237, qJ(2) * t237, -0.2e1 * t154 * t187, 0.2e1 * t205 * t192, -t214, -t213, 0, -t158 * t214 + (qJ(2) * t200 + qJD(2) * t154) * t190, -t158 * t213 + (-qJ(2) * t201 + qJD(2) * t157) * t190, t81 * t108 + t80 * t109 - t61 * t117 - t60 * t120 - t164, t128 * t193 + t206 * t210 + t48 * t80 + t49 * t81 + t57 * t60 + t58 * t61, -t174 * t42 + t32 * t74, -t161 * t174 - t173 * t32 + t42 * t68 - t228, -t263, -t242, 0, t100 * t33 - t7 * t145 + t161 * t82 + t173 * t76 + t85 * t68, t100 * t32 - t6 * t145 + t174 * t85 - t42 * t82 + t76 * t74, t155 * t260 - t56 * t188 -(-t155 * t54 - t224) * t42 + (-t14 - t155 * t19 + (t152 * t54 - t155 * t56) * qJD(6)) * t74, t155 * t176 + t161 * t56 + t173 * t18 + t188 * t251, -t152 * t176 - t161 * t54 - t173 * t19 + t59 * t74, -t161 * t251 + t173 * t33, t26 * t19 - t175 * t161 + t7 * t54 + t8 * t173 + (-t10 * t251 + t231 + (-t17 * t173 + t251 * t27 + t232) * qJD(6)) * t155 + t238 * t152, t26 * t18 - t5 * t161 + t7 * t56 + ((-qJD(6) * t27 + t10) * t251 - t231 - (-qJD(6) * t17 + t9) * t173 - qJD(6) * t232) * t152 + t238 * t155; 0, 0, 0, 0, -t160, -t212, 0, 0, 0, 0, 0, t204 * t154, t204 * t157, t108 * t172 - t109 * t171 + t118 * t117 + t119 * t120, t164 - t196, 0, 0, 0, 0, 0, -qJD(1) * t68 - t263, -qJD(1) * t174 - t242, 0, 0, 0, 0, 0, -t173 * t29 - t74 * t19 + t42 * t54 - (-t152 * t161 - t155 * t177) * t251, -t173 * t31 - (t152 * t177 - t155 * t161) * t251 - t260; 0, 0, 0, 0, 0, 0, t154 * t211, -t205 * t160, 0, 0, 0, -qJ(2) * t211, t154 * t212 (t58 + t64) * t120 - (t57 - t65) * t117 + (t108 * t150 + t109 * t151) * pkin(3), -t57 * t64 - t58 * t65 + (t150 * t49 + t151 * t48 - t157 * t196) * pkin(3), t253, t248, t262, t250, 0, -t145 * t217 - t88 * t68 + t239, t145 * t218 - t174 * t88 + t247, t258, t224 * t251 + t259, t257, t170 + t229, t244, t112 * t19 + t217 * t54 + (t179 * t251 - t3) * t155 + t166 * t152 + t241, t112 * t18 + t155 * t166 - t179 * t252 + t217 * t56 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117 ^ 2 - t120 ^ 2, t58 * t117 + t57 * t120 + t206, 0, 0, 0, 0, 0, t33 + t220, t32 - t219, 0, 0, 0, 0, 0, t170 - t229, -t155 * t251 ^ 2 - t230 - t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, t248, t262, t250, 0, t21 * t145 + t239, t20 * t145 + t247, t258, t252 * t56 + t259, t257, -t251 * t252 + t229 + t31, t244, -pkin(5) * t19 - t3 * t155 + (-t152 * t20 + t155 * t36) * t251 - t21 * t54 + t152 * t256 - t226 * pkin(9) + t241, -pkin(5) * t18 - (t152 * t36 + t155 * t20) * t251 - t21 * t56 + t16 * t254 + t243 * pkin(9) + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t54 ^ 2 + t56 ^ 2, -t251 * t54 + t18, -t251 * t56 - t19, t33, -t152 * t2 - t16 * t56 - t261 * t5 + t8, -t152 * t9 - t155 * t2 + t16 * t54 + t175 * t261;];
tauc_reg  = t1;
