% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:38:56
% EndTime: 2019-03-08 20:39:04
% DurationCPUTime: 3.03s
% Computational Cost: add. (3606->332), mult. (9449->477), div. (0->0), fcn. (7819->12), ass. (0->187)
t152 = cos(pkin(12));
t247 = cos(qJ(4));
t199 = t247 * t152;
t141 = qJD(2) * t199;
t150 = sin(pkin(12));
t156 = sin(qJ(4));
t218 = t156 * t150;
t193 = qJD(2) * t218;
t117 = -t141 + t193;
t255 = qJD(5) + qJD(6);
t270 = t117 + t255;
t170 = t199 - t218;
t246 = pkin(8) + qJ(3);
t133 = t246 * t150;
t134 = t246 * t152;
t256 = -t247 * t133 - t156 * t134;
t58 = t170 * qJD(3) + t256 * qJD(4);
t151 = sin(pkin(6));
t160 = cos(qJ(2));
t222 = t151 * t160;
t165 = t170 * t222;
t98 = qJD(1) * t165;
t269 = -t58 + t98;
t120 = t170 * qJD(4);
t126 = t247 * t150 + t156 * t152;
t121 = t126 * qJD(4);
t157 = sin(qJ(2));
t214 = qJD(1) * t151;
t195 = t157 * t214;
t268 = -t121 * pkin(4) + t120 * pkin(9) + t195;
t132 = qJD(2) * qJ(3) + t195;
t153 = cos(pkin(6));
t213 = qJD(1) * t153;
t140 = t152 * t213;
t239 = pkin(8) * qJD(2);
t90 = t140 + (-t132 - t239) * t150;
t103 = t152 * t132 + t150 * t213;
t91 = t152 * t239 + t103;
t37 = t156 * t90 + t247 * t91;
t267 = t37 * qJD(4);
t258 = -t156 * t91 + t247 * t90;
t32 = -qJD(4) * pkin(4) - t258;
t119 = qJD(2) * t126;
t155 = sin(qJ(5));
t159 = cos(qJ(5));
t207 = t159 * qJD(4);
t99 = t155 * t119 - t207;
t27 = t99 * pkin(5) + t32;
t101 = t155 * qJD(4) + t159 * t119;
t154 = sin(qJ(6));
t158 = cos(qJ(6));
t41 = t154 * t101 + t158 * t99;
t266 = t27 * t41;
t129 = t154 * t159 + t158 * t155;
t240 = t270 * t129;
t178 = t158 * t101 - t154 * t99;
t265 = t178 * t41;
t211 = qJD(5) * t155;
t228 = t117 * t155;
t264 = t211 + t228;
t263 = t178 ^ 2 - t41 ^ 2;
t114 = qJD(5) + t117;
t110 = qJD(6) + t114;
t208 = qJD(6) * t158;
t209 = qJD(6) * t154;
t138 = qJD(4) * t141;
t111 = -qJD(4) * t193 + t138;
t54 = qJD(5) * t207 + t159 * t111 - t119 * t211;
t55 = qJD(5) * t101 + t155 * t111;
t9 = -t101 * t209 - t154 * t55 + t158 * t54 - t99 * t208;
t262 = t41 * t110 + t9;
t166 = t126 * t222;
t96 = -t156 * t133 + t247 * t134;
t243 = -qJD(1) * t166 + qJD(3) * t126 + qJD(4) * t96;
t260 = t155 * t98 - t268 * t159;
t210 = qJD(5) * t159;
t143 = -t152 * pkin(3) - pkin(2);
t79 = -pkin(4) * t170 - t126 * pkin(9) + t143;
t259 = t268 * t155 + t269 * t159 - t79 * t210 + t96 * t211;
t72 = t129 * t126;
t194 = t160 * t214;
t124 = (qJD(3) + t194) * qJD(2);
t257 = t170 * t124;
t216 = t159 * t120;
t168 = -t126 * t211 + t216;
t182 = qJD(3) - t194;
t112 = qJD(2) * t121;
t33 = qJD(4) * pkin(9) + t37;
t113 = t143 * qJD(2) + t182;
t45 = t117 * pkin(4) - t119 * pkin(9) + t113;
t19 = t155 * t45 + t159 * t33;
t24 = qJD(4) * t258 + t257;
t137 = qJD(2) * t195;
t53 = t112 * pkin(4) - t111 * pkin(9) + t137;
t47 = t159 * t53;
t163 = -qJD(5) * t19 - t155 * t24 + t47;
t2 = t112 * pkin(5) - t54 * pkin(10) + t163;
t173 = t155 * t53 + t159 * t24 + t45 * t210 - t33 * t211;
t3 = -t55 * pkin(10) + t173;
t203 = -t154 * t3 + t158 * t2;
t14 = -t99 * pkin(10) + t19;
t235 = t158 * t14;
t18 = -t155 * t33 + t159 * t45;
t13 = -t101 * pkin(10) + t18;
t8 = t114 * pkin(5) + t13;
t5 = t154 * t8 + t235;
t254 = -qJD(6) * t5 - t27 * t178 + t203;
t10 = qJD(6) * t178 + t154 * t54 + t158 * t55;
t253 = t110 * t178 - t10;
t128 = t154 * t155 - t158 * t159;
t241 = t270 * t128;
t252 = t241 * t110 - t129 * t112;
t11 = t14 * t209;
t192 = qJD(6) * t8 + t3;
t251 = t154 * t2 + t158 * t192 - t11;
t250 = pkin(9) + pkin(10);
t87 = t159 * t96;
t249 = -pkin(10) * t216 + t121 * pkin(5) - t155 * t58 + (-t87 + (pkin(10) * t126 - t79) * t155) * qJD(5) + t260;
t197 = t126 * t210;
t227 = t120 * t155;
t169 = t197 + t227;
t248 = pkin(10) * t169 + t259;
t245 = pkin(5) * t169 + t243;
t77 = t119 * pkin(4) + t117 * pkin(9);
t244 = t155 * t77 + t159 * t258;
t242 = t155 * t79 + t87;
t238 = qJD(2) * pkin(2);
t237 = t119 * t41;
t236 = t119 * t99;
t25 = t126 * t124 + t267;
t234 = t25 * t159;
t233 = t178 * t119;
t232 = t54 * t155;
t231 = t99 * t114;
t230 = t101 * t114;
t229 = t101 * t119;
t226 = t126 * t155;
t225 = t126 * t159;
t223 = t151 * t157;
t161 = qJD(2) ^ 2;
t221 = t151 * t161;
t219 = t155 * t112;
t105 = t159 * t112;
t215 = t150 ^ 2 + t152 ^ 2;
t212 = qJD(2) * t157;
t205 = t157 * t221;
t204 = t160 * t221;
t202 = qJD(5) * t250;
t198 = t151 * t212;
t189 = t215 * t124;
t188 = t114 * t159;
t186 = -t240 * t110 - t128 * t112;
t185 = t264 * pkin(5) - t37;
t135 = t250 * t155;
t184 = pkin(10) * t228 + qJD(6) * t135 + t155 * t202 + t244;
t136 = t250 * t159;
t70 = t159 * t77;
t183 = t119 * pkin(5) + qJD(6) * t136 - t155 * t258 + t70 + (pkin(10) * t117 + t202) * t159;
t75 = t159 * t79;
t26 = -pkin(5) * t170 - pkin(10) * t225 - t155 * t96 + t75;
t28 = -pkin(10) * t226 + t242;
t181 = t154 * t26 + t158 * t28;
t115 = -t150 * t223 + t153 * t152;
t116 = t153 * t150 + t152 * t223;
t66 = t156 * t115 + t247 * t116;
t174 = t155 * t222 - t159 * t66;
t48 = -t155 * t66 - t159 * t222;
t180 = t154 * t174 + t158 * t48;
t179 = t154 * t48 - t158 * t174;
t177 = (-t150 * t132 + t140) * t150 - t103 * t152;
t176 = -t264 * t114 + t105;
t171 = t247 * t115 - t156 * t116;
t167 = -pkin(9) * t112 + t114 * t32;
t146 = -t159 * pkin(5) - pkin(4);
t127 = t182 - t238;
t85 = t112 * t170;
t73 = t128 * t126;
t60 = pkin(5) * t226 - t256;
t35 = qJD(2) * t166 + qJD(4) * t66;
t34 = qJD(2) * t165 + qJD(4) * t171;
t22 = -t209 * t226 + (t255 * t225 + t227) * t158 + t168 * t154;
t21 = -t128 * t120 - t255 * t72;
t17 = qJD(5) * t174 - t155 * t34 + t159 * t198;
t16 = qJD(5) * t48 + t155 * t198 + t159 * t34;
t12 = t55 * pkin(5) + t25;
t4 = -t154 * t14 + t158 * t8;
t1 = [0, 0, -t205, -t204, -t152 * t205, t150 * t205, t215 * t204 (-t115 * t150 + t116 * t152) * t124 + (t127 * t157 + (-t177 - t195) * t160) * t151 * qJD(2), 0, 0, 0, 0, 0, -t35 * qJD(4) + (-t112 * t160 + t117 * t212) * t151, -t34 * qJD(4) + (-t111 * t160 + t119 * t212) * t151, 0, 0, 0, 0, 0, t48 * t112 + t17 * t114 - t171 * t55 + t35 * t99, t35 * t101 + t112 * t174 - t16 * t114 - t171 * t54, 0, 0, 0, 0, 0 (-qJD(6) * t179 - t154 * t16 + t158 * t17) * t110 + t180 * t112 + t35 * t41 - t171 * t10 -(qJD(6) * t180 + t154 * t17 + t158 * t16) * t110 - t179 * t112 + t35 * t178 - t171 * t9; 0, 0, 0, 0, 0, 0, t182 * qJD(2) * t215 + t189, -t177 * qJD(3) + qJ(3) * t189 + (t177 * t160 + (-t127 - t238) * t157) * t214, t111 * t126 + t119 * t120, t111 * t170 - t126 * t112 - t120 * t117 - t119 * t121, t120 * qJD(4), -t121 * qJD(4), 0, t143 * t112 + t113 * t121 - t243 * qJD(4) + (-qJD(2) * t170 - t117) * t195, t269 * qJD(4) + t143 * t111 + t113 * t120, t101 * t168 + t54 * t225 (-t101 * t155 - t159 * t99) * t120 + (-t232 - t159 * t55 + (-t101 * t159 + t155 * t99) * qJD(5)) * t126, t101 * t121 + t126 * t105 + t114 * t168 - t170 * t54, -t114 * t169 - t99 * t121 - t126 * t219 + t170 * t55, t114 * t121 - t85, t75 * t112 - (-t210 * t33 + t47) * t170 + t18 * t121 - t256 * t55 + t32 * t197 + t243 * t99 + (-t210 * t96 + t260) * t114 + ((-qJD(5) * t79 - t58) * t114 - t96 * t112 - (-qJD(5) * t45 - t24) * t170 + t25 * t126 + t32 * t120) * t155, -t242 * t112 + t173 * t170 - t19 * t121 - t256 * t54 + t32 * t216 + (-t211 * t32 + t234) * t126 + t259 * t114 + t243 * t101, t178 * t21 - t9 * t73, t73 * t10 - t178 * t22 - t21 * t41 - t9 * t72, t21 * t110 - t73 * t112 + t121 * t178 - t170 * t9, t10 * t170 - t22 * t110 - t72 * t112 - t41 * t121, t110 * t121 - t85 (-t154 * t28 + t158 * t26) * t112 - t203 * t170 + t4 * t121 + t60 * t10 + t12 * t72 + t27 * t22 + t245 * t41 + (t248 * t154 + t249 * t158) * t110 + (-t110 * t181 + t170 * t5) * qJD(6), -t181 * t112 + t251 * t170 - t5 * t121 + t60 * t9 - t12 * t73 + t27 * t21 + t245 * t178 + ((-qJD(6) * t26 + t248) * t158 + (qJD(6) * t28 - t249) * t154) * t110; 0, 0, 0, 0, 0, 0, -t215 * t161, qJD(2) * t177 + t137, 0, 0, 0, 0, 0, 0.2e1 * t119 * qJD(4), t138 + (-t117 - t193) * qJD(4), 0, 0, 0, 0, 0, t176 - t236, -t114 ^ 2 * t159 - t219 - t229, 0, 0, 0, 0, 0, t186 - t237, -t233 + t252; 0, 0, 0, 0, 0, 0, 0, 0, t119 * t117, -t117 ^ 2 + t119 ^ 2, t138 + (t117 - t193) * qJD(4), 0, 0, -t113 * t119 - t25 + t267, t113 * t117 - t257, t101 * t188 + t232 (t54 - t231) * t159 + (-t55 - t230) * t155, t114 * t188 + t219 - t229, t176 + t236, -t114 * t119, -pkin(4) * t55 - t18 * t119 - t234 - t37 * t99 + (-pkin(9) * t210 - t70) * t114 + (t114 * t258 + t167) * t155, -pkin(4) * t54 - t37 * t101 + t19 * t119 + t25 * t155 + (pkin(9) * t211 + t244) * t114 + t167 * t159, t9 * t129 - t178 * t241, -t129 * t10 - t9 * t128 - t178 * t240 + t241 * t41, -t233 - t252, t186 + t237, -t110 * t119 (-t158 * t135 - t154 * t136) * t112 + t146 * t10 + t12 * t128 - t4 * t119 + t185 * t41 + t240 * t27 + (t154 * t184 - t158 * t183) * t110 -(-t154 * t135 + t158 * t136) * t112 + t146 * t9 + t12 * t129 + t5 * t119 + t185 * t178 - t241 * t27 + (t154 * t183 + t158 * t184) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101 * t99, t101 ^ 2 - t99 ^ 2, t54 + t231, t230 - t55, t112, -t32 * t101 + t19 * t114 + t163, t18 * t114 + t32 * t99 - t173, t265, t263, t262, t253, t112 -(-t154 * t13 - t235) * t110 + (-t101 * t41 - t110 * t209 + t158 * t112) * pkin(5) + t254, t266 + t11 + (-t14 * t110 - t2) * t154 + (t13 * t110 - t192) * t158 + (-t101 * t178 - t110 * t208 - t154 * t112) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, t263, t262, t253, t112, t5 * t110 + t254, t4 * t110 - t251 + t266;];
tauc_reg  = t1;
