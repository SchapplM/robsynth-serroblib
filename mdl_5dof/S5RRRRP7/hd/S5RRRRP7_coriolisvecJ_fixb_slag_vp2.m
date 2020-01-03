% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:08
% EndTime: 2019-12-31 21:56:23
% DurationCPUTime: 6.05s
% Computational Cost: add. (5569->454), mult. (14048->590), div. (0->0), fcn. (9318->6), ass. (0->215)
t340 = Ifges(5,1) + Ifges(6,1);
t334 = Ifges(5,5) + Ifges(6,4);
t339 = Ifges(5,6) - Ifges(6,6);
t199 = sin(qJ(4));
t202 = cos(qJ(4));
t330 = -t339 * t199 + t334 * t202;
t279 = Ifges(6,5) * t199;
t282 = Ifges(5,4) * t199;
t327 = -t340 * t202 - t279 + t282;
t255 = qJD(2) + qJD(3);
t338 = Ifges(4,6) * t255 / 0.2e1;
t204 = cos(qJ(2));
t196 = -pkin(2) * t204 - pkin(1);
t189 = qJD(1) * t196;
t242 = Ifges(4,5) * t255;
t337 = t189 * mrSges(4,2) + t242 / 0.2e1;
t200 = sin(qJ(3));
t201 = sin(qJ(2));
t203 = cos(qJ(3));
t179 = t200 * t201 - t203 * t204;
t171 = t179 * qJD(1);
t180 = t200 * t204 + t203 * t201;
t172 = t180 * qJD(1);
t284 = Ifges(4,4) * t172;
t336 = t338 - Ifges(4,2) * t171 / 0.2e1 + t284 / 0.2e1;
t152 = t255 * t180;
t139 = t152 * qJD(1);
t151 = t255 * t179;
t138 = t151 * qJD(1);
t156 = t199 * t172 - t202 * t255;
t88 = -qJD(4) * t156 - t202 * t138;
t157 = t202 * t172 + t199 * t255;
t89 = qJD(4) * t157 - t199 * t138;
t333 = (-Ifges(5,4) + Ifges(6,5)) * t89 + t340 * t88 + t334 * t139;
t155 = Ifges(5,4) * t156;
t168 = qJD(4) + t171;
t280 = Ifges(6,5) * t156;
t322 = t340 * t157 + t334 * t168 - t155 + t280;
t313 = -pkin(7) - pkin(6);
t191 = t313 * t204;
t183 = qJD(1) * t191;
t173 = t200 * t183;
t190 = t313 * t201;
t182 = qJD(1) * t190;
t176 = qJD(2) * pkin(2) + t182;
t145 = t203 * t176 + t173;
t124 = -pkin(3) * t255 - t145;
t235 = mrSges(6,1) * t199 - mrSges(6,3) * t202;
t237 = mrSges(5,1) * t199 + mrSges(5,2) * t202;
t39 = t156 * pkin(4) - t157 * qJ(5) + t124;
t331 = t124 * t237 + t39 * t235;
t329 = t334 * t199 + t339 * t202;
t278 = Ifges(6,5) * t202;
t281 = Ifges(5,4) * t202;
t328 = t340 * t199 - t278 + t281;
t113 = t171 * pkin(3) - t172 * pkin(8) + t189;
t174 = t203 * t183;
t146 = t200 * t176 - t174;
t125 = pkin(8) * t255 + t146;
t256 = qJD(4) * t202;
t257 = qJD(4) * t199;
t258 = qJD(2) * t201;
t254 = pkin(2) * t258;
t52 = pkin(3) * t139 + pkin(8) * t138 + qJD(1) * t254;
t247 = qJD(2) * t313;
t240 = qJD(1) * t247;
t177 = t201 * t240;
t216 = t204 * t240;
t78 = qJD(3) * t145 + t203 * t177 + t200 * t216;
t6 = t113 * t256 - t125 * t257 + t199 * t52 + t202 * t78;
t45 = t113 * t199 + t125 * t202;
t7 = -qJD(4) * t45 - t199 * t78 + t202 * t52;
t326 = -t199 * t7 + t202 * t6;
t2 = qJ(5) * t139 + qJD(5) * t168 + t6;
t44 = t113 * t202 - t125 * t199;
t319 = qJD(5) - t44;
t26 = -pkin(4) * t168 + t319;
t4 = -pkin(4) * t139 - t7;
t325 = t199 * t4 + t2 * t202 + t26 * t256;
t224 = Ifges(6,3) * t199 + t278;
t230 = -Ifges(5,2) * t199 + t281;
t299 = t199 / 0.2e1;
t300 = -t199 / 0.2e1;
t305 = t157 / 0.2e1;
t307 = t156 / 0.2e1;
t308 = -t156 / 0.2e1;
t154 = Ifges(6,5) * t157;
t70 = Ifges(6,6) * t168 + Ifges(6,3) * t156 + t154;
t283 = Ifges(5,4) * t157;
t73 = -Ifges(5,2) * t156 + Ifges(5,6) * t168 + t283;
t324 = (-t199 * t45 - t202 * t44) * mrSges(5,3) + t224 * t307 + t230 * t308 + t73 * t300 + t70 * t299 + t331 - t327 * t305 + t330 * t168 / 0.2e1;
t27 = qJ(5) * t168 + t45;
t323 = -t189 * mrSges(4,1) - t44 * mrSges(5,1) + t26 * mrSges(6,1) + t45 * mrSges(5,2) - t27 * mrSges(6,3) + t336;
t24 = mrSges(5,1) * t89 + mrSges(5,2) * t88;
t79 = qJD(3) * t146 + t177 * t200 - t203 * t216;
t321 = m(5) * t79 + t24;
t167 = pkin(4) * t257 - qJ(5) * t256 - qJD(5) * t199;
t268 = t171 * t202;
t269 = t171 * t199;
t239 = -pkin(4) * t269 + qJ(5) * t268;
t320 = -t146 - t239 + t167;
t144 = t179 * pkin(3) - t180 * pkin(8) + t196;
t159 = t190 * t200 - t191 * t203;
t318 = t199 * t144 + t202 * t159;
t317 = t203 * t190 + t191 * t200;
t87 = pkin(3) * t152 + pkin(8) * t151 + t254;
t184 = t201 * t247;
t185 = t204 * t247;
t94 = qJD(3) * t317 + t184 * t203 + t185 * t200;
t13 = -qJD(4) * t318 - t199 * t94 + t202 * t87;
t316 = t88 / 0.2e1;
t315 = -t89 / 0.2e1;
t314 = t89 / 0.2e1;
t311 = pkin(1) * mrSges(3,1);
t310 = pkin(1) * mrSges(3,2);
t309 = t139 / 0.2e1;
t306 = -t157 / 0.2e1;
t304 = -t168 / 0.2e1;
t302 = t171 / 0.2e1;
t298 = -t202 / 0.2e1;
t297 = t202 / 0.2e1;
t296 = m(4) * t189;
t295 = pkin(2) * t200;
t294 = pkin(2) * t203;
t293 = pkin(4) * t172;
t288 = mrSges(5,3) * t156;
t287 = mrSges(5,3) * t157;
t286 = Ifges(4,1) * t172;
t285 = Ifges(3,4) * t201;
t166 = Ifges(4,4) * t171;
t276 = t317 * t79;
t275 = t171 * mrSges(4,3);
t274 = t172 * mrSges(4,3);
t273 = Ifges(3,5) * qJD(2);
t272 = Ifges(3,6) * qJD(2);
t271 = qJD(2) * mrSges(3,1);
t270 = qJD(2) * mrSges(3,2);
t266 = t199 * t203;
t265 = t202 * t203;
t263 = -mrSges(4,1) * t255 + mrSges(5,1) * t156 + mrSges(5,2) * t157 + t274;
t106 = -mrSges(6,2) * t156 + mrSges(6,3) * t168;
t107 = -mrSges(5,2) * t168 - t288;
t262 = t106 + t107;
t108 = mrSges(5,1) * t168 - t287;
t109 = -mrSges(6,1) * t168 + mrSges(6,2) * t157;
t261 = -t108 + t109;
t142 = pkin(3) * t172 + pkin(8) * t171;
t69 = t199 * t142 + t202 * t145;
t260 = qJD(1) * t201;
t119 = pkin(2) * t260 + t142;
t148 = t182 * t203 + t173;
t64 = t199 * t119 + t202 * t148;
t259 = qJD(1) * t204;
t253 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t252 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t251 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t246 = t273 / 0.2e1;
t245 = -t272 / 0.2e1;
t31 = -t139 * mrSges(6,1) + t88 * mrSges(6,2);
t147 = t182 * t200 - t174;
t238 = mrSges(5,1) * t202 - mrSges(5,2) * t199;
t236 = mrSges(6,1) * t202 + mrSges(6,3) * t199;
t229 = Ifges(5,2) * t202 + t282;
t223 = -Ifges(6,3) * t202 + t279;
t222 = pkin(4) * t202 + qJ(5) * t199;
t221 = pkin(4) * t199 - qJ(5) * t202;
t63 = t119 * t202 - t148 * t199;
t68 = t142 * t202 - t145 * t199;
t90 = t144 * t202 - t159 * t199;
t186 = -pkin(3) - t222;
t12 = t144 * t256 - t159 * t257 + t199 * t87 + t202 * t94;
t207 = t7 * mrSges(5,1) - t4 * mrSges(6,1) - t6 * mrSges(5,2) + t2 * mrSges(6,3);
t95 = qJD(3) * t159 + t184 * t200 - t203 * t185;
t29 = -mrSges(6,2) * t89 + mrSges(6,3) * t139;
t30 = mrSges(5,1) * t139 - mrSges(5,3) * t88;
t32 = -mrSges(5,2) * t139 - mrSges(5,3) * t89;
t206 = (t29 + t32) * t202 + (-t30 + t31) * t199 + (-t199 * t262 + t202 * t261) * qJD(4) + m(5) * (-t256 * t44 - t257 * t45 + t326) + m(6) * (-t257 * t27 + t325);
t10 = pkin(4) * t89 - qJ(5) * t88 - qJD(5) * t157 + t79;
t121 = -t166 + t242 + t286;
t19 = t88 * Ifges(6,5) + t139 * Ifges(6,6) + t89 * Ifges(6,3);
t20 = t88 * Ifges(5,4) - t89 * Ifges(5,2) + t139 * Ifges(5,6);
t71 = t157 * Ifges(5,5) - t156 * Ifges(5,6) + t168 * Ifges(5,3);
t72 = t157 * Ifges(6,4) + t168 * Ifges(6,2) + t156 * Ifges(6,6);
t205 = -t78 * mrSges(4,2) - Ifges(4,5) * t138 - Ifges(4,6) * t139 - t10 * t236 - t145 * t275 + t19 * t298 + t20 * t297 + t223 * t314 + t229 * t315 + (-t238 - mrSges(4,1)) * t79 + t328 * t316 + t329 * t309 + (-t166 + t121) * t302 + t333 * t299 + (-t73 / 0.2e1 + t70 / 0.2e1) * t269 + (-t224 * t308 - t230 * t307 - t330 * t304 + t327 * t306 + t331 + t337) * t171 + (t338 + Ifges(5,6) * t307 + Ifges(6,6) * t308 - Ifges(4,2) * t302 + t334 * t306 + (Ifges(5,3) + Ifges(6,2)) * t304 + t323) * t172 + (-t44 * t268 - t45 * t269 + t326) * mrSges(5,3) - (-Ifges(4,1) * t171 - t284 + t71 + t72) * t172 / 0.2e1 + (t26 * t268 + (-t257 - t269) * t27 + t325) * mrSges(6,2) + t324 * qJD(4) + (t268 / 0.2e1 + t256 / 0.2e1) * t322;
t197 = Ifges(3,4) * t259;
t188 = mrSges(3,3) * t259 - t270;
t187 = -mrSges(3,3) * t260 + t271;
t178 = t186 - t294;
t170 = Ifges(3,1) * t260 + t197 + t273;
t169 = t272 + (t204 * Ifges(3,2) + t285) * qJD(1);
t165 = t172 * qJ(5);
t162 = qJD(3) * t295 + t167;
t160 = -mrSges(4,2) * t255 - t275;
t141 = mrSges(4,1) * t171 + mrSges(4,2) * t172;
t136 = Ifges(6,2) * t139;
t134 = Ifges(5,3) * t139;
t99 = mrSges(6,1) * t156 - mrSges(6,3) * t157;
t98 = pkin(4) * t157 + qJ(5) * t156;
t93 = t180 * t221 - t317;
t86 = Ifges(6,4) * t88;
t85 = Ifges(5,5) * t88;
t84 = Ifges(5,6) * t89;
t83 = Ifges(6,6) * t89;
t81 = t147 + t239;
t58 = -pkin(4) * t179 - t90;
t53 = qJ(5) * t179 + t318;
t41 = -t68 - t293;
t40 = t165 + t69;
t34 = -t63 - t293;
t33 = t165 + t64;
t23 = mrSges(6,1) * t89 - mrSges(6,3) * t88;
t18 = -t221 * t151 + (qJD(4) * t222 - qJD(5) * t202) * t180 + t95;
t11 = -pkin(4) * t152 - t13;
t8 = qJ(5) * t152 + qJD(5) * t179 + t12;
t1 = [t53 * t29 + t58 * t31 + t90 * t30 + t318 * t32 + t93 * t23 + t18 * t99 + t8 * t106 + t12 * t107 + t13 * t108 + t11 * t109 - t317 * t24 + t94 * t160 + t196 * (mrSges(4,1) * t139 - mrSges(4,2) * t138) + t263 * t95 + (t138 * t317 - t139 * t159) * mrSges(4,3) + m(6) * (t10 * t93 + t11 * t26 + t18 * t39 + t2 * t53 + t27 * t8 + t4 * t58) + m(4) * (-t145 * t95 + t146 * t94 + t159 * t78 - t276) + m(5) * (t12 * t45 + t124 * t95 + t13 * t44 + t318 * t6 + t7 * t90 - t276) + (t170 / 0.2e1 - pkin(6) * t187 + t246 + (-0.2e1 * t310 + 0.3e1 / 0.2e1 * Ifges(3,4) * t204) * qJD(1)) * t204 * qJD(2) + (-t146 * mrSges(4,3) + t71 / 0.2e1 + t72 / 0.2e1 + t252 * t168 + t253 * t157 + t251 * t156 - t323 - t336) * t152 - (t322 * t297 - t166 / 0.2e1 + t121 / 0.2e1 + t286 / 0.2e1 + (-t199 * t27 + t202 * t26) * mrSges(6,2) + t324 - t145 * mrSges(4,3) + t337) * t151 + (-t78 * mrSges(4,3) + t85 / 0.2e1 - t84 / 0.2e1 + t134 / 0.2e1 + t86 / 0.2e1 + t136 / 0.2e1 + t83 / 0.2e1 + Ifges(4,4) * t138 + t251 * t89 + t253 * t88 + (Ifges(4,2) + t252) * t139 + t207) * t179 + (t230 * t315 + t19 * t299 + t20 * t300 + t10 * t235 + t224 * t314 - Ifges(4,1) * t138 - Ifges(4,4) * t139 + (mrSges(4,3) + t237) * t79 + (-t199 * t6 - t202 * t7) * mrSges(5,3) + (-t199 * t2 + t202 * t4) * mrSges(6,2) + (t39 * t236 + t124 * t238 + t223 * t308 + t229 * t307 + t73 * t298 + (t199 * t44 - t202 * t45) * mrSges(5,3) + (-t199 * t26 - t202 * t27) * mrSges(6,2) + t328 * t306 + t329 * t304 + t322 * t300) * qJD(4) - t327 * t316 + t330 * t309 + (qJD(4) * t70 + t333) * t297) * t180 + (-t169 / 0.2e1 - pkin(6) * t188 + t245 + (-0.2e1 * t311 - 0.3e1 / 0.2e1 * t285 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t204) * qJD(1) + (qJD(1) * (mrSges(4,1) * t179 + mrSges(4,2) * t180) + 0.2e1 * t296 + t141) * pkin(2)) * t258; m(6) * (t10 * t178 + t162 * t39) + t205 + (-t81 + t162) * t99 + ((-t170 / 0.2e1 - t197 / 0.2e1 + t246 + qJD(1) * t310 + (t187 - t271) * pkin(6)) * t204 + (t169 / 0.2e1 + t245 + (t311 + t285 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t204) * qJD(1) + (t188 + t270) * pkin(6) + (-t141 - t296) * pkin(2)) * t201) * qJD(1) + (m(4) * (t200 * t78 - t203 * t79) + (t138 * t203 - t139 * t200) * mrSges(4,3) + (t263 * t200 + (t199 * t261 + t202 * t262 + t160) * t203 + m(4) * (-t145 * t200 + t146 * t203) + m(6) * (t26 * t266 + t265 * t27) + m(5) * (t124 * t200 + t265 * t45 - t266 * t44)) * qJD(3)) * pkin(2) + t206 * (pkin(8) + t295) - t263 * t147 - m(4) * (-t145 * t147 + t146 * t148) - m(5) * (t124 * t147 + t44 * t63 + t45 * t64) - m(6) * (t26 * t34 + t27 * t33 + t39 * t81) + t146 * t274 - t33 * t106 - t64 * t107 - t63 * t108 - t34 * t109 - t148 * t160 + t178 * t23 + t321 * (-pkin(3) - t294); (-t263 + t274) * t146 - m(5) * (t124 * t146 + t44 * t68 + t45 * t69) + t320 * t99 + t205 + t206 * pkin(8) - t40 * t106 - t69 * t107 - t68 * t108 - t41 * t109 - t145 * t160 + t186 * t23 - t321 * pkin(3) + (t10 * t186 - t26 * t41 - t27 * t40 + t320 * t39) * m(6); t207 + t136 + t134 + t86 + t85 - t84 + t83 + (-t261 + t287) * t45 + (-t262 - t288) * t44 + t73 * t305 + (Ifges(6,3) * t157 - t280) * t308 + (t156 * t26 + t157 * t27) * mrSges(6,2) + qJ(5) * t29 - pkin(4) * t31 - t98 * t99 + qJD(5) * t106 - t39 * (mrSges(6,1) * t157 + mrSges(6,3) * t156) - t124 * (mrSges(5,1) * t157 - mrSges(5,2) * t156) + (-t334 * t156 - t339 * t157) * t304 + (-pkin(4) * t4 + qJ(5) * t2 - t26 * t45 + t27 * t319 - t39 * t98) * m(6) + (-Ifges(5,2) * t157 - t155 + t322) * t307 + (-t340 * t156 + t154 - t283 + t70) * t306; -t168 * t106 + t157 * t99 + 0.2e1 * (t4 / 0.2e1 + t39 * t305 + t27 * t304) * m(6) + t31;];
tauc = t1(:);
