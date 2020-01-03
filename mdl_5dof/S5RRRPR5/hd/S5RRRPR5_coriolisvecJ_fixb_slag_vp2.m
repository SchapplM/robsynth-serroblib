% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:01
% EndTime: 2019-12-31 21:13:13
% DurationCPUTime: 4.95s
% Computational Cost: add. (7261->442), mult. (19271->620), div. (0->0), fcn. (13836->8), ass. (0->222)
t205 = sin(qJ(3));
t206 = sin(qJ(2));
t208 = cos(qJ(3));
t209 = cos(qJ(2));
t179 = t205 * t209 + t208 * t206;
t169 = t179 * qJD(1);
t202 = qJD(2) + qJD(3);
t204 = sin(qJ(5));
t207 = cos(qJ(5));
t252 = t208 * t209;
t178 = -t205 * t206 + t252;
t168 = t178 * qJD(1);
t203 = sin(pkin(9));
t265 = cos(pkin(9));
t217 = t203 * t168 + t169 * t265;
t111 = t202 * t207 - t204 * t217;
t236 = t265 * t168 - t169 * t203;
t119 = qJD(5) - t236;
t231 = mrSges(6,1) * t204 + mrSges(6,2) * t207;
t307 = -pkin(7) - pkin(6);
t192 = t307 * t209;
t184 = qJD(1) * t192;
t170 = t205 * t184;
t191 = t307 * t206;
t183 = qJD(1) * t191;
t176 = qJD(2) * pkin(2) + t183;
t135 = t208 * t176 + t170;
t163 = t169 * qJ(4);
t107 = t135 - t163;
t100 = pkin(3) * t202 + t107;
t173 = t208 * t184;
t136 = t176 * t205 - t173;
t262 = qJ(4) * t168;
t108 = t136 + t262;
t256 = t203 * t108;
t56 = t100 * t265 - t256;
t54 = -t202 * pkin(4) - t56;
t218 = t54 * t231;
t226 = Ifges(6,5) * t207 - Ifges(6,6) * t204;
t283 = Ifges(6,4) * t207;
t228 = -Ifges(6,2) * t204 + t283;
t284 = Ifges(6,4) * t204;
t230 = Ifges(6,1) * t207 - t284;
t295 = t207 / 0.2e1;
t296 = -t204 / 0.2e1;
t112 = t202 * t204 + t207 * t217;
t302 = t112 / 0.2e1;
t280 = t112 * Ifges(6,4);
t41 = t111 * Ifges(6,2) + t119 * Ifges(6,6) + t280;
t109 = Ifges(6,4) * t111;
t42 = Ifges(6,1) * t112 + Ifges(6,5) * t119 + t109;
t318 = t218 + t42 * t295 + t41 * t296 + t119 * t226 / 0.2e1 + t230 * t302 + t111 * t228 / 0.2e1;
t315 = mrSges(5,3) * t236;
t274 = t217 * Ifges(5,4);
t277 = t236 * Ifges(5,4);
t267 = mrSges(5,1) * t202 + mrSges(6,1) * t111 - mrSges(6,2) * t112 - mrSges(5,3) * t217;
t141 = t208 * t183 + t170;
t113 = -t163 + t141;
t140 = -t183 * t205 + t173;
t220 = t140 - t262;
t237 = t265 * t205;
t282 = pkin(2) * qJD(3);
t314 = -t113 * t203 + t220 * t265 + (t203 * t208 + t237) * t282;
t255 = t203 * t205;
t160 = (t208 * t265 - t255) * t282;
t64 = t113 * t265 + t203 * t220;
t313 = t160 - t64;
t180 = t205 * t191;
t146 = -t208 * t192 + t180;
t101 = t265 * t108;
t57 = t203 * t100 + t101;
t55 = pkin(8) * t202 + t57;
t198 = -pkin(2) * t209 - pkin(1);
t190 = qJD(1) * t198;
t144 = -t168 * pkin(3) + qJD(4) + t190;
t61 = -pkin(4) * t236 - pkin(8) * t217 + t144;
t15 = -t204 * t55 + t207 * t61;
t16 = t204 * t61 + t207 * t55;
t312 = -t15 * t204 + t16 * t207;
t212 = (t252 * t307 - t180) * qJD(2) * qJD(1);
t142 = t202 * t178;
t131 = t142 * qJD(1);
t221 = -t131 * qJ(4) - t169 * qJD(4);
t247 = qJD(3) * t208;
t248 = qJD(3) * t205;
t143 = t202 * t179;
t132 = t143 * qJD(1);
t242 = qJD(2) * t307;
t89 = t242 * t169 + t176 * t247 + t184 * t248;
t49 = -qJ(4) * t132 + qJD(4) * t168 + t89;
t14 = t265 * t49 + (-t176 * t248 + t184 * t247 + t212 + t221) * t203;
t251 = qJD(1) * t206;
t200 = pkin(2) * t251;
t116 = pkin(3) * t132 + qJD(2) * t200;
t86 = t131 * t203 + t132 * t265;
t87 = t131 * t265 - t203 * t132;
t25 = pkin(4) * t86 - pkin(8) * t87 + t116;
t2 = qJD(5) * t15 + t14 * t207 + t204 * t25;
t259 = qJD(5) * t16;
t3 = -t14 * t204 + t207 * t25 - t259;
t45 = qJD(5) * t111 + t207 * t87;
t46 = -qJD(5) * t112 - t204 * t87;
t311 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t45 + Ifges(6,6) * t46;
t310 = t45 / 0.2e1;
t309 = t46 / 0.2e1;
t308 = t86 / 0.2e1;
t306 = pkin(1) * mrSges(3,1);
t305 = pkin(1) * mrSges(3,2);
t304 = -t111 / 0.2e1;
t303 = -t112 / 0.2e1;
t301 = -t119 / 0.2e1;
t299 = t168 / 0.2e1;
t298 = t169 / 0.2e1;
t294 = m(4) * t190;
t293 = pkin(3) * t169;
t292 = pkin(3) * t203;
t90 = -qJD(3) * t136 + t212;
t13 = t203 * t49 - t265 * (t221 + t90);
t122 = qJ(4) * t178 + t146;
t145 = t208 * t191 + t192 * t205;
t219 = -qJ(4) * t179 + t145;
t72 = t122 * t203 - t219 * t265;
t291 = t13 * t72;
t290 = t15 * mrSges(6,3);
t289 = t2 * t207;
t288 = t3 * t204;
t287 = mrSges(4,3) * t168;
t286 = Ifges(3,4) * t206;
t285 = Ifges(4,4) * t169;
t281 = t111 * Ifges(6,6);
t279 = t112 * Ifges(6,5);
t278 = t119 * Ifges(6,3);
t276 = t236 * Ifges(5,2);
t275 = t217 * Ifges(5,1);
t271 = t169 * mrSges(4,3);
t270 = t202 * Ifges(5,5);
t269 = t202 * Ifges(5,6);
t268 = t57 * t217;
t264 = Ifges(3,5) * qJD(2);
t263 = Ifges(3,6) * qJD(2);
t261 = qJD(2) * mrSges(3,1);
t260 = qJD(2) * mrSges(3,2);
t258 = t236 * t204;
t257 = t236 * t207;
t197 = pkin(2) * t208 + pkin(3);
t162 = pkin(2) * t237 + t203 * t197;
t250 = qJD(1) * t209;
t249 = qJD(2) * t206;
t246 = qJD(5) * t204;
t245 = qJD(5) * t207;
t241 = t265 * pkin(3);
t240 = t264 / 0.2e1;
t239 = -t263 / 0.2e1;
t238 = t86 * mrSges(5,1) + t87 * mrSges(5,2);
t128 = pkin(2) * t249 + pkin(3) * t143;
t233 = -t2 * t204 - t3 * t207;
t232 = mrSges(6,1) * t207 - mrSges(6,2) * t204;
t229 = Ifges(6,1) * t204 + t283;
t227 = Ifges(6,2) * t207 + t284;
t225 = Ifges(6,5) * t204 + Ifges(6,6) * t207;
t224 = -t15 * t207 - t16 * t204;
t68 = -mrSges(6,2) * t119 + mrSges(6,3) * t111;
t69 = mrSges(6,1) * t119 - mrSges(6,3) * t112;
t222 = -t204 * t69 + t207 * t68;
t138 = -t178 * t265 + t179 * t203;
t139 = t203 * t178 + t179 * t265;
t150 = -t178 * pkin(3) + t198;
t71 = t138 * pkin(4) - t139 * pkin(8) + t150;
t73 = t122 * t265 + t203 * t219;
t29 = -t204 * t73 + t207 * t71;
t30 = t204 * t71 + t207 * t73;
t185 = t206 * t242;
t186 = t209 * t242;
t96 = t208 * t185 + t205 * t186 + t191 * t247 + t192 * t248;
t67 = pkin(4) * t217 - pkin(8) * t236 + t293;
t161 = -pkin(2) * t255 + t197 * t265;
t97 = -t146 * qJD(3) - t185 * t205 + t208 * t186;
t213 = qJD(5) * t224 - t288 + t289;
t211 = -qJ(4) * t142 - qJD(4) * t179 + t97;
t120 = Ifges(4,2) * t168 + Ifges(4,6) * t202 + t285;
t164 = Ifges(4,4) * t168;
t121 = Ifges(4,1) * t169 + Ifges(4,5) * t202 + t164;
t40 = t278 + t279 + t281;
t75 = t269 + t274 + t276;
t76 = t270 + t275 + t277;
t8 = t45 * Ifges(6,4) + t46 * Ifges(6,2) + t86 * Ifges(6,6);
t9 = t45 * Ifges(6,1) + t46 * Ifges(6,4) + t86 * Ifges(6,5);
t210 = -(Ifges(4,5) * t168 + Ifges(5,5) * t236 - Ifges(4,6) * t169 - Ifges(5,6) * t217) * t202 / 0.2e1 - (-Ifges(4,2) * t169 + t121 + t164) * t168 / 0.2e1 - (-Ifges(5,2) * t217 + t277 + t76) * t236 / 0.2e1 - (Ifges(5,1) * t236 - t274 + t40) * t217 / 0.2e1 + (-t232 - mrSges(5,1)) * t13 - t169 * (Ifges(4,1) * t168 - t285) / 0.2e1 - t42 * t257 / 0.2e1 + t41 * t258 / 0.2e1 + t204 * t9 / 0.2e1 + Ifges(4,5) * t131 - Ifges(4,6) * t132 - t144 * (mrSges(5,1) * t217 + mrSges(5,2) * t236) + (Ifges(6,3) * t217 + t226 * t236) * t301 + (Ifges(6,5) * t217 + t230 * t236) * t303 + (Ifges(6,6) * t217 + t228 * t236) * t304 - t190 * (mrSges(4,1) * t169 + mrSges(4,2) * t168) + t217 * t75 / 0.2e1 - t16 * (-mrSges(6,2) * t217 - mrSges(6,3) * t258) - t15 * (mrSges(6,1) * t217 - mrSges(6,3) * t257) + t229 * t310 + t8 * t295 + t120 * t298 + t225 * t308 + t227 * t309 + mrSges(6,3) * t289 + t135 * t287 + t136 * t271 + t56 * t315 - Ifges(5,6) * t86 + Ifges(5,5) * t87 - t89 * mrSges(4,2) + t90 * mrSges(4,1) - t14 * mrSges(5,2) + t318 * qJD(5) - t236 * t218;
t199 = Ifges(3,4) * t250;
t196 = -t241 - pkin(4);
t188 = mrSges(3,3) * t250 - t260;
t187 = -mrSges(3,3) * t251 + t261;
t167 = Ifges(3,1) * t251 + t199 + t264;
t166 = t263 + (t209 * Ifges(3,2) + t286) * qJD(1);
t155 = pkin(8) + t162;
t154 = -pkin(4) - t161;
t149 = mrSges(4,1) * t202 - t271;
t148 = -mrSges(4,2) * t202 + t287;
t147 = t200 + t293;
t134 = -mrSges(4,1) * t168 + mrSges(4,2) * t169;
t114 = -mrSges(5,2) * t202 + t315;
t92 = t142 * t265 - t203 * t143;
t91 = t142 * t203 + t143 * t265;
t83 = Ifges(6,3) * t86;
t81 = -mrSges(5,1) * t236 + mrSges(5,2) * t217;
t65 = t200 + t67;
t62 = -qJ(4) * t143 + qJD(4) * t178 + t96;
t59 = t107 * t265 - t256;
t58 = t107 * t203 + t101;
t33 = pkin(4) * t91 - pkin(8) * t92 + t128;
t24 = t204 * t65 + t207 * t64;
t23 = -t204 * t64 + t207 * t65;
t22 = t204 * t67 + t207 * t59;
t21 = -t204 * t59 + t207 * t67;
t20 = -mrSges(6,2) * t86 + mrSges(6,3) * t46;
t19 = mrSges(6,1) * t86 - mrSges(6,3) * t45;
t18 = t203 * t211 + t265 * t62;
t17 = t203 * t62 - t211 * t265;
t10 = -mrSges(6,1) * t46 + mrSges(6,2) * t45;
t5 = -qJD(5) * t30 - t18 * t204 + t207 * t33;
t4 = qJD(5) * t29 + t18 * t207 + t204 * t33;
t1 = [t202 * (Ifges(4,5) * t142 - Ifges(4,6) * t143) / 0.2e1 + t190 * (mrSges(4,1) * t143 + mrSges(4,2) * t142) + (t226 * t308 + t230 * t310 + t228 * t309 + Ifges(5,1) * t87 - Ifges(5,4) * t86 + t116 * mrSges(5,2) + t8 * t296 + t9 * t295 + (mrSges(5,3) + t231) * t13 + t233 * mrSges(6,3) + (t54 * t232 + t225 * t301 + t227 * t304 + t229 * t303 + t42 * t296 - t207 * t41 / 0.2e1 - t312 * mrSges(6,3)) * qJD(5)) * t139 + (-t131 * t145 - t132 * t146 - t135 * t142 - t136 * t143 + t178 * t89 - t179 * t90) * mrSges(4,3) + (-t132 * t178 - t143 * t299) * Ifges(4,2) + (t131 * t178 - t132 * t179 + t142 * t299 - t143 * t298) * Ifges(4,4) + t198 * (mrSges(4,1) * t132 + mrSges(4,2) * t131) + (-t56 * t92 - t57 * t91 + t72 * t87 - t73 * t86) * mrSges(5,3) + (-t269 / 0.2e1 + t144 * mrSges(5,1) - t274 / 0.2e1 - t276 / 0.2e1 + t278 / 0.2e1 + t281 / 0.2e1 + t279 / 0.2e1 + t15 * mrSges(6,1) - t16 * mrSges(6,2) + t40 / 0.2e1 - t75 / 0.2e1) * t91 - t267 * t17 + t96 * t148 + t97 * t149 + t142 * t121 / 0.2e1 - t143 * t120 / 0.2e1 + t128 * t81 + t18 * t114 + t150 * t238 + (t83 / 0.2e1 - Ifges(5,4) * t87 + t116 * mrSges(5,1) - t14 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t86 + t311) * t138 + (-pkin(6) * t187 + t167 / 0.2e1 + t240 + (-0.2e1 * t305 + 0.3e1 / 0.2e1 * Ifges(3,4) * t209) * qJD(1)) * t209 * qJD(2) + (t131 * t179 + t142 * t298) * Ifges(4,1) + (-pkin(6) * t188 - t166 / 0.2e1 + t239 + (-0.2e1 * t306 - 0.3e1 / 0.2e1 * t286 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t209) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t178 + mrSges(4,2) * t179) + t134 + 0.2e1 * t294) * pkin(2)) * t249 + m(6) * (t15 * t5 + t16 * t4 + t17 * t54 + t2 * t30 + t29 * t3 + t291) + m(5) * (t116 * t150 + t128 * t144 + t14 * t73 - t17 * t56 + t18 * t57 + t291) + m(4) * (t135 * t97 + t136 * t96 + t145 * t90 + t146 * t89) + t4 * t68 + t5 * t69 + t72 * t10 + t29 * t19 + t30 * t20 + (t270 / 0.2e1 + t144 * mrSges(5,2) + t275 / 0.2e1 + t277 / 0.2e1 + t76 / 0.2e1 + t224 * mrSges(6,3) + t318) * t92; -t141 * t148 - t140 * t149 + t154 * t10 - t147 * t81 + (-t161 * t87 - t162 * t86 + t268) * mrSges(5,3) + ((t148 * t208 - t149 * t205) * qJD(3) + (-t131 * t208 - t132 * t205) * mrSges(4,3)) * pkin(2) + t313 * t114 + t210 + (-t160 * t69 + (-qJD(5) * t68 - t19) * t155 + (-t3 - t259) * mrSges(6,3)) * t204 + ((t240 - t167 / 0.2e1 - t199 / 0.2e1 + qJD(1) * t305 + (t187 - t261) * pkin(6)) * t209 + (t239 + t166 / 0.2e1 + (t306 + t286 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t209) * qJD(1) + (t188 + t260) * pkin(6) + (-t134 - t294) * pkin(2)) * t206) * qJD(1) - t24 * t68 - t23 * t69 + (t155 * t20 + t160 * t68 + (-t155 * t69 - t290) * qJD(5)) * t207 - t267 * t314 + (-t13 * t161 + t14 * t162 - t144 * t147 + t313 * t57 - t314 * t56) * m(5) + ((t205 * t89 + t208 * t90 + (-t135 * t205 + t136 * t208) * qJD(3)) * pkin(2) - t135 * t140 - t136 * t141) * m(4) + (t13 * t154 - t15 * t23 + t213 * t155 - t16 * t24 + t312 * t160 + t314 * t54) * m(6); t196 * t10 - t59 * t114 - t135 * t148 + t136 * t149 - t21 * t69 - t22 * t68 - t245 * t290 - t81 * t293 + t210 + t267 * t58 + (-t16 * t246 - t288) * mrSges(6,3) + (t13 * t196 - t15 * t21 - t16 * t22 - t54 * t58) * m(6) + ((-t13 * t265 + t14 * t203) * pkin(3) - t144 * t293 + t56 * t58 - t57 * t59) * m(5) + (-t241 * t87 - t292 * t86 + t268) * mrSges(5,3) + (m(6) * t213 - t204 * t19 + t207 * t20 - t245 * t69 - t246 * t68) * (pkin(8) + t292); t207 * t19 + t204 * t20 + t267 * t217 + t222 * qJD(5) + (-t114 - t222) * t236 + t238 + (t119 * t312 - t54 * t217 - t233) * m(6) + (t217 * t56 - t236 * t57 + t116) * m(5); t83 - t54 * (mrSges(6,1) * t112 + mrSges(6,2) * t111) + (Ifges(6,1) * t111 - t280) * t303 + t41 * t302 + (Ifges(6,5) * t111 - Ifges(6,6) * t112) * t301 - t15 * t68 + t16 * t69 + (t111 * t15 + t112 * t16) * mrSges(6,3) + (-Ifges(6,2) * t112 + t109 + t42) * t304 + t311;];
tauc = t1(:);
