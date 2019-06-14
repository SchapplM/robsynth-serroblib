% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 00:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:07:53
% EndTime: 2019-05-05 00:08:01
% DurationCPUTime: 3.26s
% Computational Cost: add. (15401->367), mult. (29446->546), div. (0->0), fcn. (21903->14), ass. (0->227)
t205 = sin(qJ(6));
t206 = sin(qJ(5));
t207 = sin(qJ(4));
t210 = cos(qJ(5));
t211 = cos(qJ(4));
t166 = (t206 * t211 + t207 * t210) * qJD(2);
t195 = qJD(4) + qJD(5);
t209 = cos(qJ(6));
t147 = t166 * t205 - t209 * t195;
t149 = t166 * t209 + t195 * t205;
t121 = t149 * t147;
t238 = qJD(2) * qJD(4);
t230 = t211 * t238;
t237 = t207 * qJDD(2);
t170 = t230 + t237;
t191 = t211 * qJDD(2);
t231 = t207 * t238;
t224 = t191 - t231;
t225 = t206 * t170 - t210 * t224;
t123 = -qJD(5) * t166 - t225;
t122 = qJDD(6) - t123;
t266 = -t121 + t122;
t270 = t205 * t266;
t241 = qJD(2) * t207;
t164 = -t210 * t211 * qJD(2) + t206 * t241;
t143 = t166 * t164;
t194 = qJDD(4) + qJDD(5);
t265 = -t143 + t194;
t269 = t206 * t265;
t268 = t209 * t266;
t267 = t210 * t265;
t141 = pkin(5) * t164 - pkin(10) * t166;
t263 = t195 ^ 2;
t179 = qJD(4) * pkin(4) - pkin(9) * t241;
t197 = t211 ^ 2;
t214 = qJD(2) ^ 2;
t193 = t197 * t214;
t252 = sin(pkin(11));
t253 = cos(pkin(11));
t176 = t252 * g(1) - t253 * g(2);
t198 = -g(3) + qJDD(1);
t203 = cos(pkin(6));
t187 = t203 * t198;
t201 = sin(pkin(6));
t155 = -t176 * t201 + qJDD(3) + t187;
t208 = sin(qJ(2));
t212 = cos(qJ(2));
t220 = -t253 * g(1) - t252 * g(2);
t247 = t176 * t203;
t222 = t198 * t201 + t247;
t136 = t222 * t208 + t212 * t220;
t134 = -t214 * pkin(2) + t136;
t200 = sin(pkin(12));
t202 = cos(pkin(12));
t135 = -t208 * t220 + t222 * t212;
t215 = qJDD(2) * pkin(2) + t135;
t103 = t202 * t134 + t200 * t215;
t221 = -pkin(3) * t214 + qJDD(2) * pkin(8) + t103;
t95 = t207 * t155 + t211 * t221;
t71 = -pkin(4) * t193 + t224 * pkin(9) - qJD(4) * t179 + t95;
t255 = t210 * t71;
t218 = t207 * t221;
t243 = t207 * t214;
t260 = t170 * pkin(9);
t264 = qJDD(4) * pkin(4) + (pkin(4) * t243 + pkin(9) * t238 + t155) * t211 - t218 - t260;
t43 = t264 * t206 + t255;
t36 = -t263 * pkin(5) + t194 * pkin(10) - t164 * t141 + t43;
t124 = -t164 * qJD(5) + t210 * t170 + t206 * t224;
t159 = t195 * t164;
t115 = -t159 + t124;
t226 = t200 * t134 - t202 * t215;
t99 = -qJDD(2) * pkin(3) - t214 * pkin(8) + t226;
t87 = -t224 * pkin(4) - pkin(9) * t193 + t179 * t241 + t99;
t51 = -t115 * pkin(10) + (t166 * t195 - t123) * pkin(5) + t87;
t22 = t205 * t36 - t209 * t51;
t23 = t205 * t51 + t209 * t36;
t10 = t205 * t22 + t209 * t23;
t161 = qJD(6) + t164;
t227 = t205 * t124 - t209 * t194;
t88 = (qJD(6) - t161) * t149 + t227;
t145 = t147 ^ 2;
t146 = t149 ^ 2;
t160 = t161 ^ 2;
t162 = t164 ^ 2;
t163 = t166 ^ 2;
t42 = t206 * t71 - t210 * t264;
t35 = -t194 * pkin(5) - t263 * pkin(10) + t141 * t166 + t42;
t262 = -pkin(5) * t35 + pkin(10) * t10;
t261 = pkin(5) * t206;
t32 = t205 * t35;
t97 = t121 + t122;
t259 = t205 * t97;
t258 = t206 * t87;
t24 = t206 * t43 - t210 * t42;
t257 = t207 * t24;
t33 = t209 * t35;
t256 = t209 * t97;
t254 = t210 * t87;
t139 = t143 + t194;
t251 = t139 * t206;
t250 = t139 * t210;
t249 = t161 * t205;
t248 = t161 * t209;
t246 = t195 * t206;
t245 = t195 * t210;
t182 = t211 * t243;
t177 = qJDD(4) + t182;
t244 = t207 * t177;
t178 = qJDD(4) - t182;
t242 = t211 * t178;
t239 = qJD(6) + t161;
t119 = -t146 - t160;
t67 = -t119 * t205 - t256;
t223 = -t209 * t124 - t205 * t194;
t93 = t239 * t147 + t223;
t236 = pkin(5) * t93 + pkin(10) * t67 + t32;
t110 = -t160 - t145;
t64 = t110 * t209 - t270;
t90 = -t239 * t149 - t227;
t235 = pkin(5) * t90 + pkin(10) * t64 - t33;
t234 = t206 * t121;
t233 = t210 * t121;
t232 = -pkin(5) * t210 - pkin(4);
t25 = t206 * t42 + t210 * t43;
t94 = -t211 * t155 + t218;
t58 = t207 * t94 + t211 * t95;
t107 = t145 + t146;
t105 = -qJD(6) * t147 - t223;
t133 = t161 * t147;
t92 = t105 + t133;
t55 = t205 * t92 - t209 * t88;
t229 = pkin(5) * t107 + pkin(10) * t55 + t10;
t4 = t10 * t206 - t210 * t35;
t5 = t10 * t210 + t206 * t35;
t3 = -t207 * t4 + t211 * t5;
t171 = t191 - 0.2e1 * t231;
t9 = t205 * t23 - t209 * t22;
t219 = (-qJD(5) + t195) * t166 - t225;
t213 = qJD(4) ^ 2;
t196 = t207 ^ 2;
t192 = t196 * t214;
t181 = -t193 - t213;
t180 = -t192 - t213;
t175 = t192 + t193;
t174 = (t196 + t197) * qJDD(2);
t173 = -qJDD(2) * t200 - t202 * t214;
t172 = qJDD(2) * t202 - t200 * t214;
t169 = 0.2e1 * t230 + t237;
t157 = -t163 + t263;
t156 = t162 - t263;
t154 = -t163 - t263;
t153 = -t180 * t207 - t242;
t152 = t181 * t211 - t244;
t151 = -t178 * t207 + t180 * t211;
t150 = t177 * t211 + t181 * t207;
t144 = t174 * t200 + t175 * t202;
t142 = t163 - t162;
t137 = -t263 - t162;
t132 = t153 * t200 - t169 * t202;
t131 = t152 * t200 + t171 * t202;
t130 = -t146 + t160;
t129 = t145 - t160;
t125 = -t162 - t163;
t120 = t146 - t145;
t118 = -t154 * t206 - t250;
t117 = t154 * t210 - t251;
t116 = t159 + t124;
t111 = (qJD(5) + t195) * t166 + t225;
t109 = t137 * t210 - t269;
t108 = t137 * t206 + t267;
t104 = -qJD(6) * t149 - t227;
t101 = (-t147 * t209 + t149 * t205) * t161;
t100 = (-t147 * t205 - t149 * t209) * t161;
t91 = t105 - t133;
t83 = t105 * t209 - t149 * t249;
t82 = t105 * t205 + t149 * t248;
t81 = -t104 * t205 + t147 * t248;
t80 = t104 * t209 + t147 * t249;
t79 = -t117 * t207 + t118 * t211;
t78 = t117 * t211 + t118 * t207;
t77 = t116 * t206 + t210 * t219;
t76 = -t116 * t210 + t206 * t219;
t75 = t129 * t209 - t259;
t74 = -t130 * t205 + t268;
t73 = t129 * t205 + t256;
t72 = t130 * t209 + t270;
t69 = -t108 * t207 + t109 * t211;
t68 = t108 * t211 + t109 * t207;
t66 = t119 * t209 - t259;
t63 = t110 * t205 + t268;
t61 = t103 * t200 - t202 * t226;
t60 = -t115 * t202 + t200 * t79;
t59 = -t111 * t202 + t200 * t69;
t57 = t207 * t95 - t211 * t94;
t56 = -t205 * t91 + t209 * t90;
t54 = t205 * t90 + t209 * t91;
t53 = -t205 * t88 - t209 * t92;
t49 = -t207 * t76 + t211 * t77;
t48 = t207 * t77 + t211 * t76;
t47 = -t206 * t93 + t210 * t67;
t46 = t206 * t67 + t210 * t93;
t45 = -t206 * t90 + t210 * t64;
t44 = t206 * t64 + t210 * t90;
t40 = -t125 * t202 + t200 * t49;
t39 = -t107 * t206 + t210 * t55;
t38 = t107 * t210 + t206 * t55;
t37 = t200 * t58 - t202 * t99;
t31 = -pkin(10) * t66 + t33;
t30 = -pkin(10) * t63 + t32;
t29 = -t207 * t46 + t211 * t47;
t28 = t207 * t47 + t211 * t46;
t27 = -t207 * t44 + t211 * t45;
t26 = t207 * t45 + t211 * t44;
t19 = -t207 * t38 + t211 * t39;
t18 = t207 * t39 + t211 * t38;
t17 = t200 * t29 - t202 * t66;
t16 = t200 * t27 - t202 * t63;
t15 = -pkin(5) * t66 + t23;
t14 = -pkin(5) * t63 + t22;
t13 = t19 * t200 - t202 * t53;
t12 = t211 * t25 - t257;
t11 = t207 * t25 + t211 * t24;
t7 = t12 * t200 - t202 * t87;
t6 = -pkin(10) * t53 - t9;
t2 = t207 * t5 + t211 * t4;
t1 = t200 * t3 - t202 * t9;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t198, 0, 0, 0, 0, 0, 0, (qJDD(2) * t212 - t208 * t214) * t201, (-qJDD(2) * t208 - t212 * t214) * t201, 0, t203 * t187 + (t135 * t212 + t136 * t208 - t247) * t201, 0, 0, 0, 0, 0, 0, (t172 * t212 + t173 * t208) * t201, (-t172 * t208 + t173 * t212) * t201, 0, t203 * t155 + (t208 * (t103 * t202 + t200 * t226) + t212 * t61) * t201, 0, 0, 0, 0, 0, 0, t203 * t150 + (t208 * (t152 * t202 - t171 * t200) + t212 * t131) * t201, t203 * t151 + (t208 * (t153 * t202 + t169 * t200) + t212 * t132) * t201, (t208 * (t174 * t202 - t175 * t200) + t212 * t144) * t201, t203 * t57 + (t208 * (t200 * t99 + t202 * t58) + t212 * t37) * t201, 0, 0, 0, 0, 0, 0, t203 * t68 + (t208 * (t111 * t200 + t202 * t69) + t212 * t59) * t201, t203 * t78 + (t208 * (t115 * t200 + t202 * t79) + t212 * t60) * t201, t203 * t48 + (t208 * (t125 * t200 + t202 * t49) + t212 * t40) * t201, t203 * t11 + (t208 * (t12 * t202 + t200 * t87) + t212 * t7) * t201, 0, 0, 0, 0, 0, 0, t203 * t26 + (t208 * (t200 * t63 + t202 * t27) + t212 * t16) * t201, t203 * t28 + (t208 * (t200 * t66 + t202 * t29) + t212 * t17) * t201, t203 * t18 + (t208 * (t19 * t202 + t200 * t53) + t212 * t13) * t201, t203 * t2 + (t208 * (t200 * t9 + t202 * t3) + t212 * t1) * t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t135, -t136, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t172 - t226, pkin(2) * t173 - t103, 0, pkin(2) * t61, (t170 + t230) * t207, t169 * t211 + t171 * t207, t244 + t211 * (-t192 + t213), t171 * t211, t207 * (t193 - t213) + t242, 0, pkin(2) * t131 + pkin(3) * t171 + pkin(8) * t152 - t211 * t99, pkin(2) * t132 - pkin(3) * t169 + pkin(8) * t153 + t207 * t99, pkin(2) * t144 + pkin(3) * t175 + pkin(8) * t174 + t58, pkin(2) * t37 - pkin(3) * t99 + pkin(8) * t58, t207 * (t124 * t210 - t166 * t246) + t211 * (t124 * t206 + t166 * t245), t207 * (-t111 * t210 - t115 * t206) + t211 * (-t111 * t206 + t115 * t210), t207 * (-t157 * t206 + t267) + t211 * (t157 * t210 + t269), t207 * (-t123 * t206 + t164 * t245) + t211 * (t123 * t210 + t164 * t246), t207 * (t156 * t210 - t251) + t211 * (t156 * t206 + t250), (t207 * (-t164 * t210 + t166 * t206) + t211 * (-t164 * t206 - t166 * t210)) * t195, t207 * (-pkin(9) * t108 + t258) + t211 * (-pkin(4) * t111 + pkin(9) * t109 - t254) - pkin(3) * t111 + pkin(8) * t69 + pkin(2) * t59, t207 * (-pkin(9) * t117 + t254) + t211 * (-pkin(4) * t115 + pkin(9) * t118 + t258) - pkin(3) * t115 + pkin(8) * t79 + pkin(2) * t60, t207 * (-pkin(9) * t76 - t24) + t211 * (-pkin(4) * t125 + pkin(9) * t77 + t25) - pkin(3) * t125 + pkin(8) * t49 + pkin(2) * t40, -pkin(9) * t257 + t211 * (-pkin(4) * t87 + pkin(9) * t25) - pkin(3) * t87 + pkin(8) * t12 + pkin(2) * t7, t207 * (t210 * t83 + t234) + t211 * (t206 * t83 - t233), t207 * (t120 * t206 + t210 * t56) + t211 * (-t120 * t210 + t206 * t56), t207 * (t206 * t92 + t210 * t74) + t211 * (t206 * t74 - t210 * t92), t207 * (t210 * t81 - t234) + t211 * (t206 * t81 + t233), t207 * (-t206 * t88 + t210 * t75) + t211 * (t206 * t75 + t210 * t88), t207 * (t101 * t210 + t122 * t206) + t211 * (t101 * t206 - t122 * t210), t207 * (-pkin(9) * t44 - t14 * t206 + t210 * t30) + t211 * (-pkin(4) * t63 + pkin(9) * t45 + t14 * t210 + t206 * t30) - pkin(3) * t63 + pkin(8) * t27 + pkin(2) * t16, t207 * (-pkin(9) * t46 - t15 * t206 + t210 * t31) + t211 * (-pkin(4) * t66 + pkin(9) * t47 + t15 * t210 + t206 * t31) - pkin(3) * t66 + pkin(8) * t29 + pkin(2) * t17, t207 * (-pkin(9) * t38 + t210 * t6) + t211 * (pkin(9) * t39 + t206 * t6) + pkin(8) * t19 + pkin(2) * t13 + (t207 * t261 + t211 * t232 - pkin(3)) * t53, pkin(2) * t1 + (t207 * (-pkin(10) * t210 + t261) + t211 * (-pkin(10) * t206 + t232) - pkin(3)) * t9 + (pkin(8) + pkin(9)) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, 0, 0, 0, 0, 0, 0, t150, t151, 0, t57, 0, 0, 0, 0, 0, 0, t68, t78, t48, t11, 0, 0, 0, 0, 0, 0, t26, t28, t18, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t192 - t193, t237, t182, t191, qJDD(4), -t94, -t95, 0, 0, t143, t142, t116, -t143, t219, t194, pkin(4) * t108 - t42, -t255 - t206 * (pkin(9) * t230 - t260 - t94) + (-t177 * t206 + t117) * pkin(4), pkin(4) * t76, pkin(4) * t24, t82, t54, t72, t80, t73, t100, pkin(4) * t44 + t235, pkin(4) * t46 + t236, pkin(4) * t38 + t229, pkin(4) * t4 + t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t142, t116, -t143, t219, t194, -t42, -t43, 0, 0, t82, t54, t72, t80, t73, t100, t235, t236, t229, t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t120, t92, -t121, -t88, t122, -t22, -t23, 0, 0;];
tauJ_reg  = t8;
