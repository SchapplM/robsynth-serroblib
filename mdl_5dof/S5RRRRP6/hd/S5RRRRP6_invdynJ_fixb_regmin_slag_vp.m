% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP6
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:37
% EndTime: 2019-12-31 21:54:44
% DurationCPUTime: 2.49s
% Computational Cost: add. (3810->355), mult. (8549->462), div. (0->0), fcn. (5990->10), ass. (0->205)
t266 = cos(qJ(3));
t204 = qJD(3) * t266;
t193 = pkin(2) * t204;
t154 = cos(qJ(2));
t156 = -pkin(7) - pkin(6);
t120 = t156 * t154;
t111 = qJD(1) * t120;
t150 = sin(qJ(3));
t102 = t150 * t111;
t151 = sin(qJ(2));
t119 = t156 * t151;
t109 = qJD(1) * t119;
t69 = t266 * t109 + t102;
t275 = -t193 + t69;
t147 = qJ(2) + qJ(3);
t142 = sin(t147);
t155 = cos(qJ(1));
t232 = t142 * t155;
t152 = sin(qJ(1));
t233 = t142 * t152;
t274 = g(1) * t232 + g(2) * t233;
t187 = g(1) * t155 + g(2) * t152;
t153 = cos(qJ(4));
t260 = t153 * pkin(4);
t138 = pkin(3) + t260;
t143 = cos(t147);
t148 = -qJ(5) - pkin(8);
t196 = t143 * t138 - t142 * t148;
t229 = t150 * t154;
t108 = t266 * t151 + t229;
t259 = t154 * pkin(2);
t140 = pkin(1) + t259;
t273 = -pkin(8) * t108 - t140;
t149 = sin(qJ(4));
t205 = qJD(1) * t266;
t224 = qJD(1) * t151;
t99 = t150 * t224 - t154 * t205;
t240 = t149 * t99;
t272 = -qJ(5) * t240 + t153 * qJD(5);
t208 = t266 * t154;
t173 = -t150 * t151 + t208;
t217 = qJD(2) + qJD(3);
t72 = t217 * t173;
t161 = t72 * qJD(1);
t160 = t108 * qJDD(1) + t161;
t227 = t153 * t155;
t231 = t149 * t152;
t88 = t143 * t231 + t227;
t228 = t152 * t153;
t230 = t149 * t155;
t90 = -t143 * t230 + t228;
t271 = -g(1) * t90 + g(2) * t88;
t262 = g(3) * t143;
t270 = -t262 + t274;
t101 = -qJD(1) * t229 - t151 * t205;
t170 = t153 * t101 - t149 * t217;
t216 = qJDD(2) + qJDD(3);
t29 = -qJD(4) * t170 + t160 * t149 - t153 * t216;
t269 = t170 ^ 2;
t135 = g(3) * t142;
t261 = g(3) * t149;
t246 = qJD(2) * pkin(2);
t104 = t109 + t246;
t65 = t266 * t104 + t102;
t54 = -t217 * pkin(3) - t65;
t258 = t54 * t72;
t257 = t54 * t99;
t194 = t153 * t217;
t79 = -t101 * t149 - t194;
t95 = qJD(4) + t99;
t256 = t79 * t95;
t255 = t170 * t95;
t118 = t140 * qJD(1);
t52 = pkin(3) * t99 + pkin(8) * t101 - t118;
t103 = t266 * t111;
t66 = t150 * t104 - t103;
t55 = t217 * pkin(8) + t66;
t30 = -t149 * t55 + t153 * t52;
t17 = qJ(5) * t170 + t30;
t13 = pkin(4) * t95 + t17;
t254 = -t17 + t13;
t62 = -pkin(3) * t101 + pkin(8) * t99;
t253 = t149 * t62 + t153 * t65;
t53 = pkin(2) * t224 + t62;
t252 = t149 * t53 + t153 * t69;
t137 = t150 * pkin(2) + pkin(8);
t226 = -qJ(5) - t137;
t195 = qJD(4) * t226;
t251 = t149 * t195 + t153 * t193 - t252 + t272;
t144 = t153 * qJ(5);
t190 = -t101 * pkin(4) + t99 * t144;
t51 = t153 * t53;
t250 = t153 * t195 - t190 - t51 + (-qJD(5) + t275) * t149;
t64 = -pkin(3) * t173 + t273;
t84 = t150 * t119 - t266 * t120;
t76 = t153 * t84;
t249 = t149 * t64 + t76;
t200 = qJD(4) * t148;
t248 = t149 * t200 - t253 + t272;
t58 = t153 * t62;
t247 = t153 * t200 - t190 - t58 + (-qJD(5) + t65) * t149;
t245 = t101 * t99;
t244 = t13 * t153;
t222 = qJD(4) * t149;
t28 = -qJD(4) * t194 - t101 * t222 - t149 * t216 - t153 * t160;
t243 = t149 * t28;
t219 = t151 * qJDD(1);
t184 = -qJDD(1) * t208 + t150 * t219;
t73 = t217 * t108;
t44 = qJD(1) * t73 + t184;
t43 = qJDD(4) + t44;
t242 = t149 * t43;
t241 = t149 * t72;
t239 = t153 * t43;
t238 = t153 * t170;
t197 = t153 * t95;
t237 = t95 * t101;
t236 = qJD(4) * t95;
t235 = t108 * t149;
t145 = t151 ^ 2;
t225 = -t154 ^ 2 + t145;
t223 = qJD(3) * t150;
t221 = qJD(4) * t153;
t220 = qJD(1) * qJD(2);
t218 = t154 * qJDD(1);
t213 = t151 * t246;
t37 = pkin(3) * t73 - pkin(8) * t72 + t213;
t209 = qJD(2) * t156;
t110 = t151 * t209;
t112 = t154 * t209;
t174 = t266 * t119 + t150 * t120;
t40 = t174 * qJD(3) + t266 * t110 + t150 * t112;
t215 = t149 * t37 + t153 * t40 + t64 * t221;
t214 = pkin(8) * t236;
t211 = t137 * t236;
t48 = t54 * t222;
t210 = t187 * t143 + t135;
t207 = t108 * t221;
t202 = t154 * t220;
t75 = qJDD(2) * pkin(2) - t156 * (-t202 - t219);
t203 = t151 * t220;
t78 = t156 * (-t203 + t218);
t199 = t104 * t223 - t111 * t204 - t150 * t78 - t266 * t75;
t24 = -t216 * pkin(3) + t199;
t206 = -t24 - t262;
t201 = pkin(4) * t149 - t156;
t164 = t104 * t204 + t111 * t223 + t150 * t75 - t266 * t78;
t23 = t216 * pkin(8) + t164;
t198 = -qJD(4) * t52 - t23;
t139 = -t266 * pkin(2) - pkin(3);
t68 = t150 * t109 - t103;
t192 = pkin(2) * t223 - t68;
t191 = (t222 + t240) * pkin(4);
t189 = -pkin(8) * t43 + t257;
t133 = pkin(2) * t203;
t16 = t44 * pkin(3) - pkin(8) * t161 + t273 * qJDD(1) + t133;
t15 = t153 * t16;
t188 = -t55 * t221 + t15;
t186 = g(1) * t152 - g(2) * t155;
t185 = -t137 * t43 + t257;
t31 = t149 * t52 + t153 * t55;
t18 = -qJ(5) * t79 + t31;
t183 = -t149 * t18 - t244;
t182 = -qJ(5) * t72 - qJD(5) * t108;
t181 = t138 * t142 + t143 * t148;
t180 = -t31 * t101 + t143 * t261 + t24 * t149 + t54 * t221;
t179 = t30 * t101 + t274 * t153 + t48;
t178 = t187 * t142;
t177 = -0.2e1 * pkin(1) * t220 - pkin(6) * qJDD(2);
t176 = t140 + t196;
t175 = t149 * t16 + t153 * t23 + t52 * t221 - t55 * t222;
t157 = qJD(2) ^ 2;
t167 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t157 + t186;
t158 = qJD(1) ^ 2;
t166 = pkin(1) * t158 - pkin(6) * qJDD(1) + t187;
t165 = -t118 * t101 - t199 + t270;
t7 = t29 * pkin(4) + qJDD(5) + t24;
t41 = t84 * qJD(3) + t150 * t110 - t266 * t112;
t1 = pkin(4) * t43 + qJ(5) * t28 - t31 * qJD(4) + qJD(5) * t170 - t149 * t23 + t15;
t3 = -qJ(5) * t29 - qJD(5) * t79 + t175;
t163 = t183 * qJD(4) - t1 * t149 + t3 * t153 - t18 * t240 - t99 * t244 - t210;
t162 = -t118 * t99 - t164 + t210;
t117 = pkin(8) * t153 + t144;
t116 = t148 * t149;
t106 = t137 * t153 + t144;
t105 = t226 * t149;
t96 = -t140 * qJDD(1) + t133;
t91 = t143 * t227 + t231;
t89 = -t143 * t228 + t230;
t77 = t79 ^ 2;
t60 = t153 * t64;
t45 = t101 ^ 2 - t99 ^ 2;
t38 = t79 * pkin(4) + qJD(5) + t54;
t36 = -t184 + (-qJD(1) * t108 - t101) * t217;
t35 = t99 * t217 + t160;
t34 = t153 * t37;
t32 = -qJ(5) * t235 + t249;
t26 = -pkin(4) * t173 - t108 * t144 - t149 * t84 + t60;
t10 = -t101 * t170 + t95 * t197 + t242;
t9 = -t95 ^ 2 * t149 - t101 * t79 + t239;
t8 = -t170 * t197 - t243;
t6 = -qJ(5) * t207 + (-qJD(4) * t84 + t182) * t149 + t215;
t5 = pkin(4) * t73 - t149 * t40 + t34 + t182 * t153 + (-t76 + (qJ(5) * t108 - t64) * t149) * qJD(4);
t4 = (-t28 - t256) * t153 + (-t29 + t255) * t149;
t2 = [qJDD(1), t186, t187, qJDD(1) * t145 + 0.2e1 * t151 * t202, 0.2e1 * t151 * t218 - 0.2e1 * t225 * t220, qJDD(2) * t151 + t154 * t157, qJDD(2) * t154 - t151 * t157, 0, t177 * t151 + t167 * t154, -t167 * t151 + t177 * t154, -t101 * t72 + t160 * t108, t101 * t73 - t108 * t44 + t160 * t173 - t72 * t99, t108 * t216 + t72 * t217, t173 * t216 - t73 * t217, 0, -t118 * t73 - t140 * t44 + t186 * t143 - t173 * t96 + t174 * t216 + t99 * t213 - t41 * t217, -g(1) * t233 + g(2) * t232 - t101 * t213 + t96 * t108 - t118 * t72 - t140 * t160 - t84 * t216 - t40 * t217, -t72 * t238 + (-t153 * t28 + t170 * t222) * t108, (t149 * t170 - t153 * t79) * t72 + (t243 - t153 * t29 + (t149 * t79 + t238) * qJD(4)) * t108, t72 * t197 + t173 * t28 - t73 * t170 + (-t95 * t222 + t239) * t108, -t95 * t241 + t173 * t29 - t73 * t79 + (-t95 * t221 - t242) * t108, -t173 * t43 + t73 * t95, (-t84 * t221 + t34) * t95 + t60 * t43 - t188 * t173 + t30 * t73 + t41 * t79 - t174 * t29 + t54 * t207 - g(1) * t89 - g(2) * t91 + ((-qJD(4) * t64 - t40) * t95 - t84 * t43 - t198 * t173 + t24 * t108 + t258) * t149, -(-t84 * t222 + t215) * t95 - t249 * t43 + t175 * t173 - t31 * t73 - t41 * t170 + t174 * t28 + t153 * t258 - g(1) * t88 - g(2) * t90 + (t24 * t153 - t48) * t108, t26 * t28 - t29 * t32 + t5 * t170 - t6 * t79 + t183 * t72 + t186 * t142 + (-t1 * t153 - t149 * t3 + (t13 * t149 - t153 * t18) * qJD(4)) * t108, t3 * t32 + t18 * t6 + t1 * t26 + t13 * t5 + t7 * (pkin(4) * t235 - t174) + t38 * ((t207 + t241) * pkin(4) + t41) - g(1) * (-t152 * t176 + t155 * t201) - g(2) * (t152 * t201 + t155 * t176); 0, 0, 0, -t151 * t158 * t154, t225 * t158, t219, t218, qJDD(2), -g(3) * t154 + t151 * t166, g(3) * t151 + t154 * t166, -t245, t45, t35, t36, t216, t68 * t217 + (t266 * t216 - t217 * t223 - t99 * t224) * pkin(2) + t165, t69 * t217 + (t101 * t224 - t150 * t216 - t217 * t204) * pkin(2) + t162, t8, t4, t10, t9, t237, t139 * t29 - t51 * t95 + t192 * t79 + (t206 - t211) * t153 + (t275 * t95 + t185) * t149 + t179, -t139 * t28 + t252 * t95 - t192 * t170 + (-t95 * t193 + t185) * t153 + (-t178 + t211) * t149 + t180, t105 * t28 - t106 * t29 + t170 * t250 - t251 * t79 + t163, t3 * t106 + t1 * t105 + t7 * (t139 - t260) - g(3) * (t196 + t259) + (t103 + (pkin(2) * qJD(3) - t109) * t150 + t191) * t38 + t251 * t18 + t250 * t13 + t187 * (pkin(2) * t151 + t181); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, t45, t35, t36, t216, t66 * t217 + t165, t65 * t217 + t162, t8, t4, t10, t9, t237, -pkin(3) * t29 - t58 * t95 - t66 * t79 + (t65 * t95 + t189) * t149 + (t206 - t214) * t153 + t179, pkin(3) * t28 + t253 * t95 + t66 * t170 + t189 * t153 + (-t178 + t214) * t149 + t180, t116 * t28 - t117 * t29 + t170 * t247 - t248 * t79 + t163, t3 * t117 + t1 * t116 - t7 * t138 - g(3) * t196 + (t191 - t66) * t38 + t248 * t18 + t247 * t13 + t187 * t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170 * t79, -t77 + t269, -t28 + t256, -t29 - t255, t43, t31 * t95 + t54 * t170 + (t198 + t135) * t149 + t188 + t271, g(1) * t91 - g(2) * t89 + t153 * t135 + t30 * t95 + t54 * t79 - t175, pkin(4) * t28 - t254 * t79, t254 * t18 + (t142 * t261 + t170 * t38 + t1 + t271) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77 - t269, -t13 * t170 + t18 * t79 - t270 + t7;];
tau_reg = t2;
