% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:49
% EndTime: 2019-12-31 20:24:18
% DurationCPUTime: 12.77s
% Computational Cost: add. (8178->600), mult. (25033->887), div. (0->0), fcn. (19488->10), ass. (0->273)
t189 = sin(pkin(5));
t188 = sin(pkin(10));
t193 = sin(qJ(2));
t196 = cos(qJ(2));
t266 = cos(pkin(10));
t235 = t266 * t196;
t205 = -t188 * t193 + t235;
t201 = t189 * t205;
t153 = qJD(1) * t201;
t328 = -t153 + qJD(4);
t292 = pkin(2) * t188;
t186 = pkin(8) + t292;
t192 = sin(qJ(4));
t253 = qJD(4) * t192;
t236 = t266 * t193;
t256 = qJD(1) * t189;
t247 = t196 * t256;
t154 = -t188 * t247 - t236 * t256;
t248 = t193 * t256;
t234 = pkin(2) * t248;
t108 = -pkin(3) * t154 - pkin(8) * t153 + t234;
t195 = cos(qJ(4));
t190 = cos(pkin(5));
t293 = pkin(1) * t190;
t183 = t193 * t293;
t262 = t189 * t196;
t286 = pkin(7) + qJ(3);
t326 = t262 * t286 + t183;
t141 = t326 * qJD(1);
t132 = t188 * t141;
t184 = t196 * t293;
t179 = qJD(1) * t184;
t242 = t286 * t193;
t230 = t189 * t242;
t140 = -qJD(1) * t230 + t179;
t97 = t140 * t266 - t132;
t52 = t192 * t108 + t195 * t97;
t342 = -t186 * t253 - t52;
t237 = t266 * t141;
t96 = t140 * t188 + t237;
t341 = -t96 + t328 * (pkin(4) * t192 - pkin(9) * t195);
t340 = -pkin(9) * t154 - t342;
t191 = sin(qJ(5));
t194 = cos(qJ(5));
t258 = t195 * t153;
t113 = -t154 * t194 - t191 * t258;
t252 = qJD(4) * t195;
t339 = t191 * t252 + t113;
t160 = (t188 * t196 + t236) * t189;
t155 = qJD(2) * t160;
t148 = qJD(1) * t155;
t181 = qJD(1) * t190 + qJD(2);
t122 = t154 * t192 + t181 * t195;
t156 = qJD(2) * t201;
t149 = qJD(1) * t156;
t81 = qJD(4) * t122 + t149 * t195;
t123 = -t154 * t195 + t181 * t192;
t151 = -t205 * t256 + qJD(4);
t83 = -t123 * t191 + t151 * t194;
t36 = qJD(5) * t83 + t148 * t191 + t194 * t81;
t317 = t36 / 0.2e1;
t84 = t123 * t194 + t151 * t191;
t37 = -qJD(5) * t84 + t148 * t194 - t191 * t81;
t316 = t37 / 0.2e1;
t82 = qJD(4) * t123 + t149 * t192;
t310 = t82 / 0.2e1;
t176 = qJD(2) * t179;
t200 = (-qJD(2) * t242 + qJD(3) * t196) * t189;
t117 = qJD(1) * t200 + t176;
t263 = t189 * t193;
t126 = -qJD(2) * t326 - qJD(3) * t263;
t198 = t126 * qJD(1);
t66 = t117 * t266 + t188 * t198;
t124 = pkin(2) * t181 + t140;
t75 = t188 * t124 + t237;
t69 = pkin(8) * t181 + t75;
t171 = (-pkin(2) * t196 - pkin(1)) * t189;
t167 = qJD(1) * t171 + qJD(3);
t90 = -t153 * pkin(3) + t154 * pkin(8) + t167;
t254 = qJD(2) * t189;
t241 = qJD(1) * t254;
t229 = t193 * t241;
t210 = pkin(2) * t229;
t91 = pkin(3) * t148 - pkin(8) * t149 + t210;
t13 = t192 * t91 + t195 * t66 + t90 * t252 - t253 * t69;
t11 = pkin(9) * t148 + t13;
t65 = t117 * t188 - t266 * t198;
t27 = pkin(4) * t82 - pkin(9) * t81 + t65;
t42 = t192 * t90 + t195 * t69;
t33 = pkin(9) * t151 + t42;
t74 = t124 * t266 - t132;
t68 = -t181 * pkin(3) - t74;
t40 = -t122 * pkin(4) - t123 * pkin(9) + t68;
t9 = -t191 * t33 + t194 * t40;
t1 = qJD(5) * t9 + t11 * t194 + t191 * t27;
t338 = t1 * mrSges(6,2);
t10 = t191 * t40 + t194 * t33;
t2 = -qJD(5) * t10 - t11 * t191 + t194 * t27;
t337 = t2 * mrSges(6,1);
t62 = mrSges(5,1) * t148 - mrSges(5,3) * t81;
t8 = -mrSges(6,1) * t37 + mrSges(6,2) * t36;
t336 = -t62 + t8;
t249 = t266 * pkin(2);
t187 = -t249 - pkin(3);
t170 = -t195 * pkin(4) - t192 * pkin(9) + t187;
t264 = t186 * t195;
t131 = t170 * t191 + t194 * t264;
t335 = -qJD(5) * t131 + t340 * t191 + t194 * t341;
t130 = t170 * t194 - t191 * t264;
t334 = qJD(5) * t130 + t191 * t341 - t340 * t194;
t50 = -mrSges(6,1) * t83 + mrSges(6,2) * t84;
t89 = mrSges(5,1) * t151 - mrSges(5,3) * t123;
t285 = t50 - t89;
t270 = t154 * mrSges(4,3);
t333 = mrSges(4,1) * t181 + mrSges(5,1) * t122 - mrSges(5,2) * t123 + t270;
t159 = t188 * t263 - t189 * t235;
t112 = t159 * pkin(3) - t160 * pkin(8) + t171;
t137 = pkin(2) * t190 + t184 - t230;
t257 = pkin(7) * t262 + t183;
t152 = qJ(3) * t262 + t257;
t104 = t188 * t137 + t266 * t152;
t93 = pkin(8) * t190 + t104;
t332 = t192 * t112 + t195 * t93;
t251 = qJD(5) * t192;
t331 = t194 * t251 + t339;
t114 = -t154 * t191 + t194 * t258;
t330 = t191 * t251 - t194 * t252 + t114;
t14 = -qJD(4) * t42 - t192 * t66 + t195 * t91;
t329 = t13 * t195 - t14 * t192;
t226 = t1 * t194 - t191 * t2;
t216 = Ifges(6,5) * t194 - Ifges(6,6) * t191;
t279 = Ifges(6,4) * t194;
t219 = -Ifges(6,2) * t191 + t279;
t280 = Ifges(6,4) * t191;
t222 = Ifges(6,1) * t194 - t280;
t224 = mrSges(6,1) * t191 + mrSges(6,2) * t194;
t225 = t10 * t191 + t194 * t9;
t121 = qJD(5) - t122;
t79 = Ifges(6,4) * t83;
t31 = Ifges(6,1) * t84 + Ifges(6,5) * t121 + t79;
t268 = t194 * t31;
t304 = t121 / 0.2e1;
t306 = t84 / 0.2e1;
t308 = t83 / 0.2e1;
t296 = Ifges(6,4) * t84;
t30 = Ifges(6,2) * t83 + Ifges(6,6) * t121 + t296;
t318 = -t30 / 0.2e1;
t41 = -t192 * t69 + t195 * t90;
t32 = -pkin(4) * t151 - t41;
t327 = -t225 * mrSges(6,3) + t191 * t318 + t268 / 0.2e1 + t32 * t224 + t219 * t308 + t222 * t306 + t216 * t304;
t243 = t193 * t254;
t233 = pkin(2) * t243;
t109 = pkin(3) * t155 - pkin(8) * t156 + t233;
t180 = qJD(2) * t184;
t125 = t180 + t200;
t72 = t125 * t266 + t188 * t126;
t22 = -qJD(4) * t332 + t109 * t195 - t192 * t72;
t325 = t9 * mrSges(6,1) - t10 * mrSges(6,2);
t275 = Ifges(6,3) * t121;
t294 = Ifges(6,6) * t83;
t295 = Ifges(6,5) * t84;
t29 = t275 + t294 + t295;
t276 = Ifges(5,6) * t151;
t283 = Ifges(5,4) * t123;
t60 = Ifges(5,2) * t122 + t276 + t283;
t324 = t60 / 0.2e1 - t29 / 0.2e1 - t325;
t323 = -0.2e1 * pkin(1);
t6 = Ifges(6,4) * t36 + Ifges(6,2) * t37 + Ifges(6,6) * t82;
t322 = t6 / 0.2e1;
t321 = Ifges(6,1) * t317 + Ifges(6,4) * t316 + Ifges(6,5) * t310;
t319 = t29 / 0.2e1;
t315 = -t60 / 0.2e1;
t119 = Ifges(5,4) * t122;
t277 = Ifges(5,5) * t151;
t61 = Ifges(5,1) * t123 + t119 + t277;
t313 = -t61 / 0.2e1;
t312 = t81 / 0.2e1;
t311 = -t82 / 0.2e1;
t309 = -t83 / 0.2e1;
t307 = -t84 / 0.2e1;
t305 = -t121 / 0.2e1;
t135 = t160 * t192 - t190 * t195;
t303 = -t135 / 0.2e1;
t136 = t160 * t195 + t190 * t192;
t302 = t136 / 0.2e1;
t300 = -t154 / 0.2e1;
t298 = t190 / 0.2e1;
t297 = -t193 / 0.2e1;
t289 = t41 * mrSges(5,3);
t288 = t42 * mrSges(5,3);
t284 = Ifges(3,4) * t193;
t282 = Ifges(5,4) * t192;
t281 = Ifges(5,4) * t195;
t278 = Ifges(3,5) * t196;
t12 = -pkin(4) * t148 - t14;
t274 = t12 * t192;
t271 = t153 * mrSges(4,3);
t269 = t154 * Ifges(4,4);
t261 = t191 * t192;
t260 = t192 * t153;
t259 = t192 * t194;
t5 = Ifges(6,5) * t36 + Ifges(6,6) * t37 + Ifges(6,3) * t82;
t250 = Ifges(5,5) * t81 - Ifges(5,6) * t82 + Ifges(5,3) * t148;
t71 = t125 * t188 - t266 * t126;
t232 = mrSges(3,3) * t248;
t231 = mrSges(3,3) * t247;
t223 = Ifges(5,1) * t195 - t282;
t221 = Ifges(6,1) * t191 + t279;
t220 = -Ifges(5,2) * t192 + t281;
t218 = Ifges(6,2) * t194 + t280;
t217 = Ifges(5,5) * t195 - Ifges(5,6) * t192;
t215 = Ifges(6,5) * t191 + Ifges(6,6) * t194;
t23 = mrSges(6,1) * t82 - mrSges(6,3) * t36;
t24 = -mrSges(6,2) * t82 + mrSges(6,3) * t37;
t214 = -t191 * t23 + t194 * t24;
t48 = pkin(9) * t159 + t332;
t103 = t137 * t266 - t188 * t152;
t92 = -t190 * pkin(3) - t103;
t53 = t135 * pkin(4) - t136 * pkin(9) + t92;
t18 = t191 * t53 + t194 * t48;
t17 = -t191 * t48 + t194 * t53;
t57 = -mrSges(6,2) * t121 + mrSges(6,3) * t83;
t58 = mrSges(6,1) * t121 - mrSges(6,3) * t84;
t213 = -t191 * t57 - t194 * t58;
t51 = t108 * t195 - t192 * t97;
t54 = t112 * t195 - t192 * t93;
t111 = t136 * t194 + t159 * t191;
t110 = -t136 * t191 + t159 * t194;
t21 = t192 * t109 + t112 * t252 + t195 * t72 - t253 * t93;
t206 = t181 * (-Ifges(3,6) * t193 + t278);
t166 = t257 * qJD(2);
t157 = -pkin(7) * t229 + t176;
t158 = qJD(1) * t166;
t199 = -t158 * mrSges(3,1) - t65 * mrSges(4,1) - t157 * mrSges(3,2) - t66 * mrSges(4,2);
t178 = Ifges(3,4) * t247;
t175 = t241 * t278;
t168 = -pkin(7) * t263 + t184;
t165 = -pkin(7) * t243 + t180;
t164 = t257 * qJD(1);
t163 = -pkin(7) * t248 + t179;
t162 = -t181 * mrSges(3,2) + t231;
t161 = mrSges(3,1) * t181 - t232;
t150 = Ifges(4,4) * t153;
t145 = Ifges(4,5) * t149;
t144 = Ifges(4,6) * t148;
t142 = t149 * mrSges(4,2);
t139 = Ifges(3,1) * t248 + t181 * Ifges(3,5) + t178;
t138 = Ifges(3,6) * t181 + (Ifges(3,2) * t196 + t284) * t256;
t127 = -mrSges(4,2) * t181 + t271;
t115 = -mrSges(4,1) * t153 - mrSges(4,2) * t154;
t107 = -t154 * Ifges(4,1) + t181 * Ifges(4,5) + t150;
t106 = t153 * Ifges(4,2) + t181 * Ifges(4,6) - t269;
t102 = -qJD(4) * t135 + t156 * t195;
t101 = qJD(4) * t136 + t156 * t192;
t88 = -mrSges(5,2) * t151 + mrSges(5,3) * t122;
t73 = pkin(4) * t123 - pkin(9) * t122;
t63 = -mrSges(5,2) * t148 - mrSges(5,3) * t82;
t59 = t123 * Ifges(5,5) + t122 * Ifges(5,6) + t151 * Ifges(5,3);
t49 = mrSges(5,1) * t82 + mrSges(5,2) * t81;
t47 = -pkin(4) * t159 - t54;
t46 = qJD(5) * t110 + t102 * t194 + t155 * t191;
t45 = -qJD(5) * t111 - t102 * t191 + t155 * t194;
t43 = pkin(4) * t154 - t51;
t39 = t81 * Ifges(5,1) - t82 * Ifges(5,4) + t148 * Ifges(5,5);
t38 = t81 * Ifges(5,4) - t82 * Ifges(5,2) + t148 * Ifges(5,6);
t28 = pkin(4) * t101 - pkin(9) * t102 + t71;
t26 = t191 * t73 + t194 * t41;
t25 = -t191 * t41 + t194 * t73;
t16 = -pkin(4) * t155 - t22;
t15 = pkin(9) * t155 + t21;
t4 = -qJD(5) * t18 - t15 * t191 + t194 * t28;
t3 = qJD(5) * t17 + t15 * t194 + t191 * t28;
t7 = [m(5) * (t13 * t332 + t14 * t54 + t21 * t42 + t22 * t41 + t65 * t92 + t68 * t71) + t332 * t63 + m(3) * (t157 * t257 - t158 * t168 - t163 * t166 + t164 * t165) + (t157 * t196 + t158 * t193) * t189 * mrSges(3,3) + (t175 / 0.2e1 + t145 / 0.2e1 - t144 / 0.2e1 + t199) * t190 + m(6) * (t1 * t18 + t10 * t3 + t12 * t47 + t16 * t32 + t17 * t2 + t4 * t9) + (t196 * t139 / 0.2e1 + t138 * t297 + t206 / 0.2e1 + t193 * pkin(2) * t115 + (-t163 * t196 - t164 * t193) * mrSges(3,3)) * t254 + (-t103 * mrSges(4,3) + Ifges(4,1) * t160 - Ifges(4,4) * t159 + Ifges(4,5) * t298) * t149 + (-t104 * mrSges(4,3) + Ifges(5,5) * t302 + Ifges(5,6) * t303 - Ifges(4,4) * t160 + t171 * mrSges(4,1) - Ifges(4,6) * t190 / 0.2e1 + (Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t159) * t148 + t171 * t142 + m(4) * (-t103 * t65 + t104 * t66 + t167 * t233 - t71 * t74 + t72 * t75) + (Ifges(4,1) * t156 - Ifges(4,4) * t155) * t300 + t39 * t302 + t38 * t303 + (Ifges(6,5) * t46 + Ifges(6,6) * t45 + Ifges(6,3) * t101) * t304 + (Ifges(6,1) * t46 + Ifges(6,4) * t45 + Ifges(6,5) * t101) * t306 + (Ifges(6,4) * t46 + Ifges(6,2) * t45 + Ifges(6,6) * t101) * t308 + (Ifges(6,5) * t111 + Ifges(6,6) * t110 + Ifges(6,3) * t135) * t310 + (Ifges(5,4) * t136 - Ifges(5,2) * t135 + Ifges(5,6) * t159) * t311 + (Ifges(5,1) * t136 - Ifges(5,4) * t135 + Ifges(5,5) * t159) * t312 + t101 * t315 + (Ifges(6,4) * t111 + Ifges(6,2) * t110 + Ifges(6,6) * t135) * t316 + (Ifges(6,1) * t111 + Ifges(6,4) * t110 + Ifges(6,5) * t135) * t317 + t101 * t319 + t111 * t321 + t110 * t322 + t159 * t250 / 0.2e1 - t333 * t71 + ((-t168 * mrSges(3,3) + Ifges(3,5) * t298 + (mrSges(3,2) * t323 + 0.3e1 / 0.2e1 * Ifges(3,4) * t196) * t189) * t196 + (-t257 * mrSges(3,3) - Ifges(3,6) * t190 + (mrSges(3,1) * t323 - 0.3e1 / 0.2e1 * t284 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t196) * t189 + (m(4) * t171 + mrSges(4,1) * t159 + mrSges(4,2) * t160) * pkin(2)) * t193) * t241 + (-t155 * t75 - t156 * t74 - t159 * t66 + t160 * t65) * mrSges(4,3) + t17 * t23 + t18 * t24 + t45 * t30 / 0.2e1 + t46 * t31 / 0.2e1 + t32 * (-mrSges(6,1) * t45 + mrSges(6,2) * t46) + t47 * t8 + t16 * t50 + t3 * t57 + t4 * t58 + t54 * t62 + t21 * t88 + t22 * t89 + t92 * t49 + t10 * (-mrSges(6,2) * t101 + mrSges(6,3) * t45) + t9 * (mrSges(6,1) * t101 - mrSges(6,3) * t46) + t102 * t61 / 0.2e1 + t68 * (mrSges(5,1) * t101 + mrSges(5,2) * t102) + t12 * (-mrSges(6,1) * t110 + mrSges(6,2) * t111) + t72 * t127 + t135 * t5 / 0.2e1 + t1 * (-mrSges(6,2) * t135 + mrSges(6,3) * t110) + t2 * (mrSges(6,1) * t135 - mrSges(6,3) * t111) + t65 * (mrSges(5,1) * t135 + mrSges(5,2) * t136) + t155 * t59 / 0.2e1 + t42 * (-mrSges(5,2) * t155 - mrSges(5,3) * t101) + t41 * (mrSges(5,1) * t155 - mrSges(5,3) * t102) + t122 * (Ifges(5,4) * t102 - Ifges(5,2) * t101 + Ifges(5,6) * t155) / 0.2e1 + t123 * (Ifges(5,1) * t102 - Ifges(5,4) * t101 + Ifges(5,5) * t155) / 0.2e1 + t151 * (Ifges(5,5) * t102 - Ifges(5,6) * t101 + Ifges(5,3) * t155) / 0.2e1 - t155 * t106 / 0.2e1 + t156 * t107 / 0.2e1 + t153 * (Ifges(4,4) * t156 - Ifges(4,2) * t155) / 0.2e1 + t13 * (-mrSges(5,2) * t159 - mrSges(5,3) * t135) + t14 * (mrSges(5,1) * t159 - mrSges(5,3) * t136) + t165 * t162 - t166 * t161 + t167 * (mrSges(4,1) * t155 + mrSges(4,2) * t156) + t181 * (Ifges(4,5) * t156 - Ifges(4,6) * t155) / 0.2e1; t63 * t264 + t74 * t271 + t224 * t274 - (t206 + (-Ifges(3,2) * t248 + t139 + t178) * t196) * t256 / 0.2e1 + (t258 * t41 + t260 * t42 + t329) * mrSges(5,3) + (mrSges(6,1) * t331 - mrSges(6,2) * t330) * t32 + (-t1 * t261 - t10 * t331 - t2 * t259 + t330 * t9) * mrSges(6,3) + t175 + (Ifges(6,5) * t114 + Ifges(6,6) * t113) * t305 + ((Ifges(3,1) * t196 - t284) * t297 + pkin(1) * (mrSges(3,1) * t193 + mrSges(3,2) * t196)) * qJD(1) ^ 2 * t189 ^ 2 - t6 * t261 / 0.2e1 - t144 + t145 - t252 * t289 + t339 * t318 + (-t167 * t234 + t74 * t96 - t75 * t97 + (t188 * t66 - t266 * t65) * pkin(2)) * m(4) + t342 * t88 - t151 * (-Ifges(5,3) * t154 + t153 * t217) / 0.2e1 - t122 * (-Ifges(5,6) * t154 + t153 * t220) / 0.2e1 - t123 * (-Ifges(5,5) * t154 + t153 * t223) / 0.2e1 + (t187 * t65 - t41 * t51 - t42 * t52 - t68 * t96) * m(5) + (t1 * t131 + t334 * t10 + t130 * t2 - t32 * t43 + t335 * t9) * m(6) + (((-t192 * t42 - t195 * t41) * qJD(4) + t329) * m(5) + t285 * t252 + t336 * t192 + m(6) * (t252 * t32 + t274)) * t186 + (t231 - t162) * t163 + (Ifges(6,4) * t114 + Ifges(6,2) * t113) * t309 + (-t288 + t315 + t319 + t325) * t253 - t115 * t234 + t106 * t300 + (-t215 * t251 + (Ifges(6,3) * t192 + t195 * t216) * qJD(4)) * t304 + (-t221 * t251 + (Ifges(6,5) * t192 + t195 * t222) * qJD(4)) * t306 + (-t218 * t251 + (Ifges(6,6) * t192 + t195 * t219) * qJD(4)) * t308 + (-Ifges(6,3) * t195 + t192 * t216) * t310 + (Ifges(5,2) * t195 + t282) * t311 - t42 * mrSges(5,2) * t154 + t41 * mrSges(5,1) * t154 + t138 * t248 / 0.2e1 - t75 * t270 + t199 + (Ifges(5,1) * t192 + t281) * t312 + t258 * t313 + (-Ifges(6,6) * t195 + t192 * t219) * t316 + (-Ifges(6,5) * t195 + t192 * t222) * t317 + t259 * t321 + (-t148 * t292 - t149 * t249) * mrSges(4,3) + t333 * t96 + t328 * t68 * (mrSges(5,1) * t192 + mrSges(5,2) * t195) - Ifges(3,6) * t229 + (Ifges(6,5) * t307 + Ifges(6,6) * t309 + Ifges(6,3) * t305 + t324) * t260 - (t191 * t31 + t194 * t30) * t251 / 0.2e1 + (t268 + t61) * t252 / 0.2e1 + (t122 * t220 + t123 * t223 + t151 * t217) * qJD(4) / 0.2e1 + (Ifges(4,1) * t153 + t269 + t59) * t154 / 0.2e1 - (Ifges(4,2) * t154 + t107 + t150) * t153 / 0.2e1 + (Ifges(6,1) * t114 + Ifges(6,4) * t113) * t307 + t195 * t338 + t334 * t57 + t335 * t58 - t195 * t337 - t43 * t50 - t51 * t89 + (t232 + t161) * t164 - t114 * t31 / 0.2e1 - t97 * t127 + t130 * t23 + t131 * t24 - t167 * (-mrSges(4,1) * t154 + mrSges(4,2) * t153) - t181 * (Ifges(4,5) * t153 + Ifges(4,6) * t154) / 0.2e1 + t187 * t49 + t192 * t39 / 0.2e1 - t195 * t5 / 0.2e1 + t195 * t38 / 0.2e1 + t65 * (-mrSges(5,1) * t195 + mrSges(5,2) * t192) + t148 * (Ifges(5,5) * t192 + Ifges(5,6) * t195) / 0.2e1; t148 * mrSges(4,1) - t113 * t58 - t114 * t57 - t153 * t127 + t142 - t333 * t154 + (-t153 * t88 + (-t191 * t58 + t194 * t57 + t88) * qJD(4) - t336) * t195 + (qJD(5) * t213 + t285 * t328 + t214 + t63) * t192 + (-t10 * t114 - t113 * t9 - t32 * t260 + (-t12 + (t10 * t194 - t191 * t9) * qJD(4)) * t195 + (qJD(4) * t32 - qJD(5) * t225 + t226) * t192) * m(6) + (t13 * t192 + t14 * t195 + t154 * t68 + t328 * (-t192 * t41 + t195 * t42)) * m(5) + (-t153 * t75 - t154 * t74 + t210) * m(4); t250 + t218 * t316 + t221 * t317 - t285 * t42 + t215 * t310 + (-t68 * mrSges(5,1) + t276 / 0.2e1 - t294 / 0.2e1 - t295 / 0.2e1 - t275 / 0.2e1 + t288 + t283 / 0.2e1 + t324) * t123 - pkin(4) * t8 - t13 * mrSges(5,2) + t14 * mrSges(5,1) - t26 * t57 - t25 * t58 - t41 * t88 + (t289 + t313 - t68 * mrSges(5,2) - t119 / 0.2e1 - t277 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t123 - t327) * t122 + t327 * qJD(5) + t226 * mrSges(6,3) + t191 * t321 + t194 * t322 + t12 * (-mrSges(6,1) * t194 + mrSges(6,2) * t191) + (-pkin(4) * t12 - t10 * t26 - t25 * t9 - t32 * t42) * m(6) + (t214 + m(6) * t226 + (-m(6) * t225 + t213) * qJD(5)) * pkin(9); -t338 + t337 - t32 * (mrSges(6,1) * t84 + mrSges(6,2) * t83) + (Ifges(6,1) * t83 - t296) * t307 + t30 * t306 + (Ifges(6,5) * t83 - Ifges(6,6) * t84) * t305 - t9 * t57 + t10 * t58 + (t10 * t84 + t83 * t9) * mrSges(6,3) + t5 + (-Ifges(6,2) * t84 + t31 + t79) * t309;];
tauc = t7(:);
