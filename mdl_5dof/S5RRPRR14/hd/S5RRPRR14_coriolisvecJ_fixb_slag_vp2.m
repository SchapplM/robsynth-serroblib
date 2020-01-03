% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR14
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR14_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR14_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:48
% EndTime: 2019-12-31 20:36:19
% DurationCPUTime: 14.11s
% Computational Cost: add. (9126->608), mult. (24842->883), div. (0->0), fcn. (19587->10), ass. (0->295)
t221 = sin(pkin(10));
t223 = cos(pkin(10));
t226 = sin(qJ(4));
t229 = cos(qJ(4));
t198 = t221 * t229 + t223 * t226;
t222 = sin(pkin(5));
t230 = cos(qJ(2));
t287 = t222 * t230;
t238 = t198 * t287;
t155 = qJD(1) * t238;
t193 = t198 * qJD(4);
t284 = t155 - t193;
t197 = t221 * t226 - t229 * t223;
t237 = t197 * t287;
t156 = qJD(1) * t237;
t192 = t197 * qJD(4);
t283 = -t156 + t192;
t224 = cos(pkin(5));
t281 = qJD(1) * t224;
t215 = qJD(2) + t281;
t227 = sin(qJ(2));
t282 = qJD(1) * t222;
t274 = t227 * t282;
t172 = t215 * t223 - t221 * t274;
t173 = t215 * t221 + t223 * t274;
t270 = t229 * t172 - t173 * t226;
t370 = t270 / 0.2e1;
t273 = t230 * t282;
t209 = qJD(4) - t273;
t372 = t209 / 0.2e1;
t366 = Ifges(5,4) * t370 + Ifges(5,5) * t372;
t246 = t172 * t226 + t229 * t173;
t371 = t246 / 0.2e1;
t367 = Ifges(5,1) * t371;
t348 = t367 + t366;
t225 = sin(qJ(5));
t228 = cos(qJ(5));
t275 = pkin(1) * t281;
t187 = pkin(7) * t273 + t227 * t275;
t159 = qJ(3) * t215 + t187;
t180 = (-pkin(2) * t230 - qJ(3) * t227 - pkin(1)) * t222;
t165 = qJD(1) * t180;
t103 = -t221 * t159 + t223 * t165;
t69 = -pkin(3) * t273 - t173 * pkin(8) + t103;
t104 = t223 * t159 + t221 * t165;
t79 = pkin(8) * t172 + t104;
t36 = t226 * t69 + t229 * t79;
t32 = pkin(9) * t209 + t36;
t186 = -pkin(7) * t274 + t230 * t275;
t150 = -t215 * pkin(2) + qJD(3) - t186;
t116 = -t172 * pkin(3) + t150;
t45 = -pkin(4) * t270 - pkin(9) * t246 + t116;
t12 = -t225 * t32 + t228 * t45;
t195 = t224 * t227 * pkin(1) + pkin(7) * t287;
t189 = t195 * qJD(2);
t178 = qJD(1) * t189;
t280 = qJD(2) * t222;
t271 = qJD(1) * t280;
t266 = t230 * t271;
t243 = t221 * t266;
t139 = pkin(3) * t243 + t178;
t234 = qJD(2) * t237;
t77 = -qJD(1) * t234 + qJD(4) * t270;
t235 = qJD(2) * t238;
t78 = qJD(1) * t235 + qJD(4) * t246;
t27 = t78 * pkin(4) - t77 * pkin(9) + t139;
t278 = qJD(4) * t229;
t279 = qJD(4) * t226;
t286 = t223 * t230;
t239 = (pkin(3) * t227 - pkin(8) * t286) * t222;
t236 = qJD(2) * t239;
t253 = pkin(2) * t227 - qJ(3) * t230;
t162 = (qJD(2) * t253 - qJD(3) * t227) * t222;
t144 = qJD(1) * t162;
t325 = pkin(1) * t230;
t277 = t224 * t325;
t214 = qJD(2) * t277;
t267 = t227 * t271;
t177 = -pkin(7) * t267 + qJD(1) * t214;
t145 = qJD(3) * t215 + t177;
t96 = t223 * t144 - t221 * t145;
t70 = qJD(1) * t236 + t96;
t97 = t221 * t144 + t223 * t145;
t80 = -pkin(8) * t243 + t97;
t10 = t226 * t70 + t229 * t80 + t69 * t278 - t279 * t79;
t8 = pkin(9) * t267 + t10;
t1 = qJD(5) * t12 + t225 * t27 + t228 * t8;
t377 = t1 * mrSges(6,2);
t13 = t225 * t45 + t228 * t32;
t2 = -qJD(5) * t13 - t225 * t8 + t228 * t27;
t376 = t2 * mrSges(6,1);
t185 = t253 * t282;
t126 = t223 * t185 - t221 * t186;
t102 = qJD(1) * t239 + t126;
t127 = t221 * t185 + t223 * t186;
t269 = t221 * t273;
t111 = -pkin(8) * t269 + t127;
t315 = pkin(8) + qJ(3);
t204 = t315 * t221;
t205 = t315 * t223;
t244 = -t229 * t204 - t205 * t226;
t361 = -qJD(3) * t197 + qJD(4) * t244 - t226 * t102 - t229 * t111;
t115 = qJD(5) - t270;
t296 = Ifges(6,3) * t115;
t93 = t209 * t228 - t225 * t246;
t326 = Ifges(6,6) * t93;
t94 = t209 * t225 + t228 * t246;
t327 = Ifges(6,5) * t94;
t28 = t296 + t326 + t327;
t297 = Ifges(5,6) * t209;
t301 = Ifges(5,2) * t270;
t308 = Ifges(5,4) * t246;
t61 = t297 + t301 + t308;
t375 = -t61 / 0.2e1 + t28 / 0.2e1;
t35 = -t226 * t79 + t229 * t69;
t374 = -t116 * mrSges(5,2) + t35 * mrSges(5,3);
t373 = -t116 * mrSges(5,1) - t12 * mrSges(6,1) + t13 * mrSges(6,2) + t36 * mrSges(5,3);
t369 = -pkin(9) * t274 + t361;
t148 = pkin(3) * t269 + t187;
t368 = -pkin(4) * t284 + t283 * pkin(9) - t148;
t154 = -t204 * t226 + t205 * t229;
t362 = -qJD(3) * t198 - qJD(4) * t154 - t102 * t229 + t111 * t226;
t220 = -pkin(3) * t223 - pkin(2);
t136 = pkin(4) * t197 - pkin(9) * t198 + t220;
t88 = t136 * t225 + t154 * t228;
t365 = -qJD(5) * t88 - t225 * t369 + t228 * t368;
t87 = t136 * t228 - t154 * t225;
t364 = qJD(5) * t87 + t225 * t368 + t228 * t369;
t363 = pkin(4) * t274 - t362;
t252 = t12 * t228 + t13 * t225;
t360 = t252 * mrSges(6,3);
t359 = t1 * t228 - t2 * t225;
t255 = Ifges(6,5) * t228 - Ifges(6,6) * t225;
t306 = Ifges(6,4) * t228;
t257 = -Ifges(6,2) * t225 + t306;
t307 = Ifges(6,4) * t225;
t260 = Ifges(6,1) * t228 - t307;
t262 = mrSges(6,1) * t225 + mrSges(6,2) * t228;
t328 = Ifges(6,4) * t94;
t29 = Ifges(6,2) * t93 + Ifges(6,6) * t115 + t328;
t89 = Ifges(6,4) * t93;
t30 = Ifges(6,1) * t94 + Ifges(6,5) * t115 + t89;
t31 = -pkin(4) * t209 - t35;
t330 = t228 / 0.2e1;
t333 = -t225 / 0.2e1;
t341 = t115 / 0.2e1;
t343 = t94 / 0.2e1;
t345 = t93 / 0.2e1;
t232 = t255 * t341 + t257 * t345 + t260 * t343 + t262 * t31 + t29 * t333 + t30 * t330;
t358 = t348 + t232 + t366;
t293 = t209 * Ifges(5,3);
t294 = t246 * Ifges(5,5);
t295 = t270 * Ifges(5,6);
t60 = t293 + t294 + t295;
t11 = -qJD(4) * t36 - t226 * t80 + t229 * t70;
t179 = qJ(3) * t224 + t195;
t122 = -t221 * t179 + t223 * t180;
t288 = t222 * t227;
t191 = t221 * t224 + t223 * t288;
t85 = -pkin(3) * t287 - t191 * pkin(8) + t122;
t123 = t223 * t179 + t221 * t180;
t190 = -t221 * t288 + t223 * t224;
t99 = pkin(8) * t190 + t123;
t314 = t226 * t85 + t229 * t99;
t272 = t227 * t280;
t188 = -pkin(7) * t272 + t214;
t171 = qJD(3) * t224 + t188;
t109 = t223 * t162 - t221 * t171;
t83 = t109 + t236;
t110 = t221 * t162 + t223 * t171;
t289 = t221 * t230;
t268 = t280 * t289;
t95 = -pkin(8) * t268 + t110;
t18 = -qJD(4) * t314 - t226 * t95 + t229 * t83;
t357 = -0.2e1 * pkin(1);
t42 = -qJD(5) * t94 - t225 * t77 + t228 * t267;
t38 = Ifges(6,6) * t42;
t41 = qJD(5) * t93 + t225 * t267 + t228 * t77;
t39 = Ifges(6,5) * t41;
t5 = Ifges(6,3) * t78 + t38 + t39;
t356 = t5 / 0.2e1;
t7 = t41 * Ifges(6,1) + t42 * Ifges(6,4) + t78 * Ifges(6,5);
t355 = t7 / 0.2e1;
t354 = Ifges(5,2) / 0.2e1;
t352 = -t29 / 0.2e1;
t351 = t41 / 0.2e1;
t350 = t42 / 0.2e1;
t347 = t78 / 0.2e1;
t346 = -t93 / 0.2e1;
t344 = -t94 / 0.2e1;
t342 = -t115 / 0.2e1;
t245 = t229 * t190 - t191 * t226;
t340 = t245 / 0.2e1;
t131 = t190 * t226 + t191 * t229;
t339 = t131 / 0.2e1;
t338 = t190 / 0.2e1;
t337 = t191 / 0.2e1;
t336 = -t221 / 0.2e1;
t335 = t223 / 0.2e1;
t334 = t224 / 0.2e1;
t331 = t227 / 0.2e1;
t323 = t10 * mrSges(5,2);
t322 = t11 * mrSges(5,1);
t320 = t35 * mrSges(5,1);
t319 = t36 * mrSges(5,2);
t318 = t77 * Ifges(5,1);
t317 = t77 * Ifges(5,4);
t316 = t78 * Ifges(5,4);
t313 = mrSges(4,2) * t223;
t311 = Ifges(3,4) * t227;
t310 = Ifges(4,4) * t221;
t309 = Ifges(4,4) * t223;
t305 = Ifges(4,5) * t173;
t304 = Ifges(4,5) * t227;
t302 = Ifges(3,2) * t230;
t300 = Ifges(3,6) * t215;
t299 = Ifges(4,6) * t172;
t298 = Ifges(4,6) * t227;
t292 = t215 * Ifges(3,5);
t101 = mrSges(5,1) * t209 - mrSges(5,3) * t246;
t50 = -mrSges(6,1) * t93 + mrSges(6,2) * t94;
t291 = t101 - t50;
t285 = -mrSges(3,1) * t215 - mrSges(4,1) * t172 + mrSges(4,2) * t173 + mrSges(3,3) * t274;
t160 = mrSges(4,1) * t243 + t266 * t313;
t276 = Ifges(5,5) * t77 - Ifges(5,6) * t78 + Ifges(5,3) * t267;
t40 = t78 * mrSges(5,1) + t77 * mrSges(5,2);
t265 = t376 - t377;
t264 = -t1 * t225 - t2 * t228;
t263 = mrSges(6,1) * t228 - mrSges(6,2) * t225;
t261 = Ifges(4,1) * t223 - t310;
t259 = Ifges(6,1) * t225 + t306;
t258 = -Ifges(4,2) * t221 + t309;
t256 = Ifges(6,2) * t228 + t307;
t254 = Ifges(6,5) * t225 + Ifges(6,6) * t228;
t251 = t12 * t225 - t13 * t228;
t47 = -pkin(9) * t287 + t314;
t216 = pkin(7) * t288;
t182 = t216 + (-pkin(2) - t325) * t224;
t132 = -t190 * pkin(3) + t182;
t59 = -pkin(4) * t245 - t131 * pkin(9) + t132;
t22 = t225 * t59 + t228 * t47;
t21 = -t225 * t47 + t228 * t59;
t57 = -mrSges(6,2) * t115 + mrSges(6,3) * t93;
t58 = mrSges(6,1) * t115 - mrSges(6,3) * t94;
t250 = -t225 * t58 + t228 * t57;
t48 = -t226 * t99 + t229 * t85;
t242 = mrSges(4,1) * t227 - mrSges(4,3) * t286;
t241 = -mrSges(4,2) * t227 - mrSges(4,3) * t289;
t112 = -t225 * t131 - t228 * t287;
t240 = -t228 * t131 + t225 * t287;
t17 = t226 * t83 + t229 * t95 + t85 * t278 - t279 * t99;
t149 = pkin(3) * t268 + t189;
t233 = -t308 / 0.2e1 + t327 / 0.2e1 - t297 / 0.2e1 + t326 / 0.2e1 + t296 / 0.2e1 + t375;
t210 = Ifges(3,4) * t273;
t207 = Ifges(3,5) * t266;
t194 = -t216 + t277;
t184 = -t215 * mrSges(3,2) + mrSges(3,3) * t273;
t167 = t242 * t271;
t166 = t241 * t271;
t158 = Ifges(3,1) * t274 + t210 + t292;
t157 = t300 + (t302 + t311) * t282;
t142 = -mrSges(4,1) * t273 - t173 * mrSges(4,3);
t141 = mrSges(4,2) * t273 + t172 * mrSges(4,3);
t134 = (t230 * t261 + t304) * t271;
t133 = (t230 * t258 + t298) * t271;
t129 = -t156 * t228 + t225 * t274;
t128 = t156 * t225 + t228 * t274;
t108 = Ifges(4,1) * t173 + Ifges(4,4) * t172 - Ifges(4,5) * t273;
t107 = Ifges(4,4) * t173 + Ifges(4,2) * t172 - Ifges(4,6) * t273;
t106 = -Ifges(4,3) * t273 + t299 + t305;
t100 = -mrSges(5,2) * t209 + mrSges(5,3) * t270;
t91 = qJD(4) * t131 + t235;
t90 = qJD(4) * t245 - t234;
t66 = -mrSges(5,2) * t267 - mrSges(5,3) * t78;
t65 = mrSges(5,1) * t267 - mrSges(5,3) * t77;
t64 = pkin(4) * t246 - pkin(9) * t270;
t63 = -mrSges(5,1) * t270 + mrSges(5,2) * t246;
t52 = qJD(5) * t240 - t225 * t90 + t228 * t272;
t51 = qJD(5) * t112 + t225 * t272 + t228 * t90;
t46 = pkin(4) * t287 - t48;
t37 = t91 * pkin(4) - t90 * pkin(9) + t149;
t34 = Ifges(5,5) * t267 - t316 + t318;
t33 = -t78 * Ifges(5,2) + Ifges(5,6) * t267 + t317;
t24 = -mrSges(6,2) * t78 + mrSges(6,3) * t42;
t23 = mrSges(6,1) * t78 - mrSges(6,3) * t41;
t20 = t225 * t64 + t228 * t35;
t19 = -t225 * t35 + t228 * t64;
t16 = -pkin(4) * t272 - t18;
t15 = pkin(9) * t272 + t17;
t14 = -mrSges(6,1) * t42 + mrSges(6,2) * t41;
t9 = -pkin(4) * t267 - t11;
t6 = t41 * Ifges(6,4) + t42 * Ifges(6,2) + t78 * Ifges(6,6);
t4 = -qJD(5) * t22 - t15 * t225 + t228 * t37;
t3 = qJD(5) * t21 + t15 * t228 + t225 * t37;
t25 = [(Ifges(6,4) * t51 + Ifges(6,2) * t52) * t345 + (Ifges(6,1) * t51 + Ifges(6,4) * t52) * t343 + m(5) * (t10 * t314 + t11 * t48 + t116 * t149 + t132 * t139 + t17 * t36 + t18 * t35) + t314 * t66 + (-Ifges(5,4) * t371 + Ifges(6,5) * t343 - Ifges(5,2) * t370 - Ifges(5,6) * t372 + Ifges(6,6) * t345 + Ifges(6,3) * t341 - t373 + t375) * t91 + t285 * t189 + (0.2e1 * t348 - t374) * t90 + (Ifges(6,5) * t51 + Ifges(6,6) * t52) * t341 + t9 * (-mrSges(6,1) * t112 - mrSges(6,2) * t240) - t240 * t355 - t78 * (Ifges(5,4) * t131 + Ifges(5,2) * t245 - Ifges(5,6) * t287) / 0.2e1 + t77 * (Ifges(5,1) * t131 + Ifges(5,4) * t245 - Ifges(5,5) * t287) / 0.2e1 + t139 * (-mrSges(5,1) * t245 + mrSges(5,2) * t131) + (-Ifges(6,1) * t240 + Ifges(6,4) * t112 - Ifges(6,5) * t245) * t351 + (-Ifges(6,5) * t240 + Ifges(6,6) * t112 - Ifges(6,3) * t245) * t347 + (-Ifges(6,4) * t240 + Ifges(6,2) * t112 - Ifges(6,6) * t245) * t350 - t245 * t356 + (t10 * t245 - t11 * t131) * mrSges(5,3) + ((-Ifges(3,6) * t224 + Ifges(5,5) * t339 + Ifges(5,6) * t340 + Ifges(4,5) * t337 + Ifges(4,6) * t338 - t195 * mrSges(3,3) + (mrSges(3,1) * t357 - 0.3e1 / 0.2e1 * t311) * t222) * t227 + (Ifges(3,5) * t334 + (Ifges(4,4) * t191 + Ifges(4,2) * t190) * t336 - t194 * mrSges(3,3) + (Ifges(4,1) * t191 + Ifges(4,4) * t190) * t335 + (mrSges(3,2) * t357 + (0.3e1 / 0.2e1 * Ifges(3,4) + 0.3e1 / 0.2e1 * t221 * Ifges(4,6) - 0.3e1 / 0.2e1 * t223 * Ifges(4,5)) * t230) * t222 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(5,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,3)) * t288) * t230) * t271 + t245 * t377 - t287 * t322 + m(6) * (t1 * t22 + t12 * t4 + t13 * t3 + t16 * t31 + t2 * t21 + t46 * t9) + m(4) * (t103 * t109 + t104 * t110 + t122 * t96 + t123 * t97 + t150 * t189 + t178 * t182) + m(3) * (t177 * t195 - t178 * t194 - t186 * t189 + t187 * t188) - t245 * t376 + (t1 * t112 - t12 * t51 + t13 * t52 + t2 * t240) * mrSges(6,3) - t276 * t287 / 0.2e1 + t96 * (-mrSges(4,1) * t287 - t191 * mrSges(4,3)) + t177 * (-t224 * mrSges(3,2) + mrSges(3,3) * t287) + t97 * (mrSges(4,2) * t287 + t190 * mrSges(4,3)) + (-mrSges(3,1) * t224 - mrSges(4,1) * t190 + mrSges(4,2) * t191 + mrSges(3,3) * t288) * t178 + t21 * t23 + t22 * t24 + t46 * t14 + t16 * t50 + t51 * t30 / 0.2e1 + t52 * t29 / 0.2e1 + t31 * (-mrSges(6,1) * t52 + mrSges(6,2) * t51) + t3 * t57 + t4 * t58 + t48 * t65 + t287 * t323 + t207 * t334 + t134 * t337 + t133 * t338 + t34 * t339 + t33 * t340 + t17 * t100 + t18 * t101 + t112 * t6 / 0.2e1 + ((t107 * t336 + t108 * t335 - t186 * mrSges(3,3) + t158 / 0.2e1 + t150 * (mrSges(4,1) * t221 + t313) + t292 / 0.2e1 + t172 * t258 / 0.2e1 + t173 * t261 / 0.2e1 + (-t103 * t223 - t104 * t221) * mrSges(4,3)) * t230 + (-t187 * mrSges(3,3) - t157 / 0.2e1 + t106 / 0.2e1 + t60 / 0.2e1 + t320 - t319 + t295 / 0.2e1 + t294 / 0.2e1 + t293 / 0.2e1 - t300 / 0.2e1 - t104 * mrSges(4,2) + t103 * mrSges(4,1) + t299 / 0.2e1 + t305 / 0.2e1) * t227) * t280 + t132 * t40 + t110 * t141 + t109 * t142 + t149 * t63 + t123 * t166 + t122 * t167 + t182 * t160 + t188 * t184; (t97 * mrSges(4,3) + qJD(3) * t141 + qJ(3) * t166 + t133 / 0.2e1 - t178 * mrSges(4,1)) * t223 + ((t225 * t192 - t128) * mrSges(6,3) + t284 * mrSges(6,2)) * t13 + ((t228 * t192 + t129) * mrSges(6,3) - t284 * mrSges(6,1)) * t12 + (t283 * t35 + t284 * t36) * mrSges(5,3) + (-mrSges(5,1) * t284 - mrSges(5,2) * t283) * t116 - t285 * t187 + t364 * t57 + (t1 * t88 + t12 * t365 + t13 * t364 + t2 * t87 - t244 * t9 + t31 * t363) * m(6) + t365 * t58 + t361 * t100 + (t10 * t154 + t11 * t244 - t116 * t148 + t139 * t220 + t35 * t362 + t36 * t361) * m(5) + t362 * t101 + t363 * t50 + (-t301 / 0.2e1 + t233) * t193 - (t14 - t65) * t244 - t209 * (-Ifges(5,5) * t156 - Ifges(5,6) * t155) / 0.2e1 - t270 * (-Ifges(5,4) * t156 - Ifges(5,2) * t155) / 0.2e1 - t246 * (-Ifges(5,1) * t156 - Ifges(5,4) * t155) / 0.2e1 + (t356 + t39 / 0.2e1 + t38 / 0.2e1 - t33 / 0.2e1 + t139 * mrSges(5,1) - t10 * mrSges(5,3) - t317 / 0.2e1 + (Ifges(6,3) / 0.2e1 + t354) * t78 + t265) * t197 - (t367 + t358) * t192 + t207 + (-pkin(2) * t178 + (-t103 * t221 + t104 * t223) * qJD(3) + (-t221 * t96 + t223 * t97) * qJ(3) - t103 * t126 - t104 * t127 - t150 * t187) * m(4) + (-t96 * mrSges(4,3) - qJD(3) * t142 - qJ(3) * t167 + t134 / 0.2e1 + t178 * mrSges(4,2)) * t221 + (t107 * t289 / 0.2e1 - t108 * t286 / 0.2e1 - t227 * t320 + t227 * t319 + (t302 * t331 + pkin(1) * (mrSges(3,1) * t227 + mrSges(3,2) * t230) + t230 * (Ifges(4,5) * t286 - Ifges(4,6) * t289 + Ifges(4,3) * t227) / 0.2e1) * t282 + t157 * t331 - t150 * (mrSges(4,1) * t289 + mrSges(4,2) * t286) - t215 * (Ifges(3,5) * t230 - Ifges(3,6) * t227) / 0.2e1 - t104 * t241 - t103 * t242 - t172 * (Ifges(4,4) * t286 - Ifges(4,2) * t289 + t298) / 0.2e1 - t173 * (Ifges(4,1) * t286 - Ifges(4,4) * t289 + t304) / 0.2e1 + (t186 * t230 + t187 * t227) * mrSges(3,3) + ((Ifges(4,5) * t221 / 0.2e1 + Ifges(4,6) * t335 + Ifges(5,5) * t198 / 0.2e1 - Ifges(5,6) * t197 / 0.2e1 - Ifges(3,6)) * t227 + ((Ifges(4,1) * t221 + t309) * t335 + (Ifges(4,2) * t223 + t310) * t336) * t230) * qJD(2) - (t210 + t158) * t230 / 0.2e1 - ((Ifges(3,1) * t230 - t311) * t282 + t106 + 0.2e1 * t60) * t227 / 0.2e1) * t282 + t128 * t352 + t156 * t348 + (Ifges(6,5) * t129 + Ifges(6,6) * t128 + Ifges(6,3) * t155) * t342 + (Ifges(6,1) * t129 + Ifges(6,4) * t128 + Ifges(6,5) * t155) * t344 + (Ifges(6,4) * t129 + Ifges(6,2) * t128 + Ifges(6,6) * t155) * t346 + t87 * t23 + t88 * t24 + (-t11 * mrSges(5,3) + t6 * t333 + t7 * t330 + t34 / 0.2e1 + t318 / 0.2e1 - t316 / 0.2e1 + t139 * mrSges(5,2) + t9 * t262 + t260 * t351 + t257 * t350 + t255 * t347 + t264 * mrSges(6,3) + (mrSges(6,3) * t251 + t228 * t352 + t254 * t342 + t256 * t346 + t259 * t344 + t263 * t31 + t30 * t333) * qJD(5)) * t198 - t129 * t30 / 0.2e1 - t31 * (-mrSges(6,1) * t128 + mrSges(6,2) * t129) - t127 * t141 - t126 * t142 - t148 * t63 + t154 * t66 - t155 * t28 / 0.2e1 + t155 * t61 / 0.2e1 - pkin(2) * t160 - t177 * mrSges(3,2) - t178 * mrSges(3,1) - t186 * t184 + t220 * t40; -t172 * t141 + t173 * t142 + t225 * t24 + t228 * t23 + t291 * t246 + t250 * qJD(5) + (-t100 - t250) * t270 + t40 + t160 + (-t115 * t251 - t246 * t31 - t264) * m(6) + (t246 * t35 - t270 * t36 + t139) * m(5) + (t103 * t173 - t104 * t172 + t178) * m(4); t259 * t351 - t9 * t263 + t276 + (-t233 + t373) * t246 + ((t354 - Ifges(5,1) / 0.2e1) * t246 + t360 - t358 + t374) * t270 + (t232 - t360) * qJD(5) + t254 * t347 + t256 * t350 + t291 * t36 + t359 * mrSges(6,3) - t323 + t322 - pkin(4) * t14 - t20 * t57 - t19 * t58 - t35 * t100 + t225 * t355 + t6 * t330 + (-pkin(4) * t9 - t12 * t19 - t13 * t20 - t31 * t36) * m(6) + (-t225 * t23 + t228 * t24 + (-m(6) * t252 - t225 * t57 - t228 * t58) * qJD(5) + m(6) * t359) * pkin(9); -t31 * (mrSges(6,1) * t94 + mrSges(6,2) * t93) + (Ifges(6,1) * t93 - t328) * t344 + t29 * t343 + (Ifges(6,5) * t93 - Ifges(6,6) * t94) * t342 - t12 * t57 + t13 * t58 + (t12 * t93 + t13 * t94) * mrSges(6,3) + t265 + t5 + (-Ifges(6,2) * t94 + t30 + t89) * t346;];
tauc = t25(:);
