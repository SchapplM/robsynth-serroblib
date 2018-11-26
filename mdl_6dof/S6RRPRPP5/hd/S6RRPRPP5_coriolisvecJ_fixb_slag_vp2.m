% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:59:02
% EndTime: 2018-11-23 16:59:11
% DurationCPUTime: 8.54s
% Computational Cost: add. (3599->546), mult. (8428->666), div. (0->0), fcn. (4365->4), ass. (0->248)
t346 = Ifges(7,4) + Ifges(6,5);
t337 = Ifges(6,4) + Ifges(5,5);
t347 = -Ifges(7,5) + t337;
t343 = Ifges(7,2) + Ifges(6,3);
t345 = Ifges(6,6) - Ifges(7,6);
t182 = cos(qJ(2));
t262 = qJD(1) * t182;
t170 = pkin(7) * t262;
t137 = pkin(3) * t262 + t170;
t178 = qJD(2) * qJ(3);
t119 = t178 + t137;
t179 = sin(qJ(4));
t181 = cos(qJ(4));
t131 = qJD(2) * t179 + t181 * t262;
t240 = t179 * t262;
t260 = qJD(2) * t181;
t132 = -t240 + t260;
t199 = qJ(5) * t132 - t119;
t306 = pkin(4) + pkin(5);
t20 = -t131 * t306 + qJD(6) + t199;
t184 = -pkin(2) - pkin(8);
t180 = sin(qJ(2));
t237 = -qJ(3) * t180 - pkin(1);
t125 = t182 * t184 + t237;
t103 = t125 * qJD(1);
t263 = qJD(1) * t180;
t168 = pkin(3) * t263;
t167 = pkin(7) * t263;
t251 = qJD(3) + t167;
t104 = qJD(2) * t184 + t168 + t251;
t35 = -t179 * t103 + t181 * t104;
t36 = t181 * t103 + t179 * t104;
t200 = t179 * t35 - t181 * t36;
t163 = qJD(4) + t263;
t319 = qJD(5) - t35;
t27 = -pkin(4) * t163 + t319;
t156 = t163 * qJ(5);
t29 = t156 + t36;
t201 = t179 * t27 + t181 * t29;
t18 = qJ(6) * t132 + t35;
t320 = qJD(5) - t18;
t13 = -t163 * t306 + t320;
t19 = qJ(6) * t131 + t36;
t16 = t156 + t19;
t202 = t13 * t179 + t16 * t181;
t222 = mrSges(7,1) * t181 + mrSges(7,2) * t179;
t224 = mrSges(6,1) * t181 + mrSges(6,3) * t179;
t226 = mrSges(5,1) * t181 - mrSges(5,2) * t179;
t294 = t181 / 0.2e1;
t295 = -t181 / 0.2e1;
t297 = -t179 / 0.2e1;
t128 = Ifges(5,4) * t131;
t330 = t346 * t131;
t339 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t316 = t339 * t132 + t347 * t163 - t128 + t330;
t332 = t346 * t132;
t321 = t343 * t131 + t163 * t345 + t332;
t38 = pkin(4) * t131 - t199;
t284 = Ifges(5,4) * t132;
t48 = -Ifges(5,2) * t131 + Ifges(5,6) * t163 + t284;
t344 = -t201 * mrSges(6,2) + t200 * mrSges(5,3) + t202 * mrSges(7,3) + t119 * t226 - t20 * t222 + t38 * t224 + t294 * t321 + t295 * t48 + t297 * t316;
t322 = Ifges(5,6) - Ifges(6,6);
t145 = -pkin(2) * t182 + t237;
t120 = t145 * qJD(1);
t147 = -t170 - t178;
t206 = Ifges(7,5) * t179 - Ifges(7,6) * t181;
t283 = Ifges(5,4) * t179;
t213 = Ifges(5,2) * t181 + t283;
t298 = t163 / 0.2e1;
t299 = -t163 / 0.2e1;
t300 = t132 / 0.2e1;
t302 = t131 / 0.2e1;
t303 = -t131 / 0.2e1;
t282 = Ifges(5,4) * t181;
t341 = t346 * t181;
t314 = t179 * t339 + t282 - t341;
t317 = t179 * t337 + t322 * t181;
t340 = t346 * t179;
t318 = -t181 * t343 + t340;
t324 = qJD(2) / 0.2e1;
t325 = -qJD(2) / 0.2e1;
t326 = -qJD(1) / 0.2e1;
t342 = Ifges(4,5) * t324 + (-Ifges(4,6) * t180 - t182 * Ifges(4,3)) * qJD(1) / 0.2e1 + Ifges(3,6) * t325 + (Ifges(3,4) * t180 + t182 * Ifges(3,2)) * t326 + t206 * t299 + t213 * t303 - t120 * mrSges(4,2) + t147 * mrSges(4,1) + t318 * t302 + t317 * t298 + t314 * t300 - t344;
t307 = pkin(3) + pkin(7);
t338 = -mrSges(3,1) + mrSges(4,2);
t250 = qJD(1) * qJD(2);
t238 = t182 * t250;
t261 = qJD(2) * t180;
t241 = t179 * t261;
t249 = qJD(2) * qJD(4);
t255 = qJD(4) * t182;
t83 = -t179 * t249 + (-t181 * t255 + t241) * qJD(1);
t239 = t180 * t250;
t84 = -qJD(4) * t240 + (-t239 + t249) * t181;
t336 = t238 * t345 + t343 * t84 + t346 * t83;
t272 = qJ(5) * t179;
t197 = t181 * t306 + t272;
t235 = -qJD(5) * t181 + qJD(3);
t264 = t167 + t168;
t335 = -t163 * t197 - t235 - t264;
t265 = qJ(6) + t184;
t236 = qJD(4) * t265;
t252 = t181 * qJD(6);
t169 = pkin(2) * t263;
t203 = pkin(8) * t180 - qJ(3) * t182;
t111 = qJD(1) * t203 + t169;
t52 = -t179 * t111 + t137 * t181;
t334 = t179 * t236 - t252 - (-qJ(6) * t179 * t180 - t182 * t306) * qJD(1) + t52;
t253 = t179 * qJD(6);
t53 = t181 * t111 + t179 * t137;
t40 = qJ(5) * t262 + t53;
t333 = t253 - t40 + (qJ(6) * t263 + t236) * t181;
t331 = t179 * t343 + t341;
t329 = (-Ifges(5,4) + t346) * t84 + t339 * t83 + t347 * t238;
t328 = t181 * t339 - t283 + t340;
t301 = -t132 / 0.2e1;
t271 = qJ(5) * t181;
t313 = t179 * t306 - t271;
t143 = -qJD(2) * pkin(2) + t251;
t166 = Ifges(3,4) * t262;
t243 = Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t233 = Ifges(7,6) / 0.2e1 + t243;
t247 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t234 = Ifges(7,5) / 0.2e1 - t247;
t245 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t277 = Ifges(4,6) * t182;
t311 = -t233 * t131 - t234 * t132 + (Ifges(7,3) / 0.2e1 + t245) * t163 + t143 * mrSges(4,1) + t16 * mrSges(7,2) + t29 * mrSges(6,3) + t35 * mrSges(5,1) + Ifges(3,1) * t263 / 0.2e1 + Ifges(3,5) * t324 + t166 / 0.2e1 + Ifges(4,4) * t325 + (-t180 * Ifges(4,2) - t277) * t326 + Ifges(7,5) * t301 + Ifges(6,6) * t302 - t120 * mrSges(4,3) - t13 * mrSges(7,1) - t27 * mrSges(6,1) - t36 * mrSges(5,2) + (Ifges(7,6) + Ifges(5,6)) * t303 + t337 * t300 + (Ifges(7,3) + Ifges(5,3) + Ifges(6,2)) * t298;
t309 = -t84 / 0.2e1;
t308 = t84 / 0.2e1;
t305 = pkin(1) * mrSges(3,1);
t304 = pkin(1) * mrSges(3,2);
t296 = t179 / 0.2e1;
t293 = t83 * mrSges(7,3);
t59 = mrSges(5,1) * t238 - mrSges(5,3) * t83;
t75 = t83 * mrSges(6,2);
t60 = -mrSges(6,1) * t238 + t75;
t292 = t59 - t60;
t57 = -mrSges(6,2) * t84 + mrSges(6,3) * t238;
t61 = mrSges(7,2) * t238 + mrSges(7,3) * t84;
t291 = t61 + t57;
t69 = mrSges(6,1) * t131 - mrSges(6,3) * t132;
t70 = -mrSges(7,1) * t131 + mrSges(7,2) * t132;
t290 = t69 - t70;
t90 = mrSges(7,2) * t163 + mrSges(7,3) * t131;
t92 = -mrSges(6,2) * t131 + mrSges(6,3) * t163;
t289 = t90 + t92;
t286 = mrSges(5,3) * t131;
t91 = -mrSges(5,2) * t163 - t286;
t288 = t91 + t92;
t285 = mrSges(5,3) * t132;
t94 = mrSges(5,1) * t163 - t285;
t95 = -mrSges(6,1) * t163 + mrSges(6,2) * t132;
t287 = -t94 + t95;
t149 = -mrSges(4,1) * t262 - qJD(2) * mrSges(4,3);
t71 = mrSges(5,1) * t131 + mrSges(5,2) * t132;
t274 = -t149 + t71;
t273 = qJ(5) * t131;
t270 = qJD(2) * mrSges(3,2);
t269 = t179 * t182;
t268 = t179 * t184;
t267 = t181 * t182;
t266 = t181 * t184;
t151 = t307 * t180;
t66 = t181 * t125 + t179 * t151;
t140 = pkin(7) * t239 - qJD(2) * qJD(3);
t152 = t307 * t182;
t259 = qJD(2) * t182;
t258 = qJD(3) * t180;
t257 = qJD(4) * t179;
t256 = qJD(4) * t181;
t254 = qJD(5) * t179;
t248 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t246 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t244 = -0.3e1 / 0.2e1 * Ifges(4,6) - 0.3e1 / 0.2e1 * Ifges(3,4);
t242 = m(4) * pkin(7) + mrSges(4,1);
t55 = t180 * qJ(5) + t66;
t32 = -t84 * mrSges(7,1) + t83 * mrSges(7,2);
t65 = -t179 * t125 + t151 * t181;
t139 = t307 * t259;
t231 = m(4) * t143 + (mrSges(4,1) + mrSges(3,3)) * t263 + t338 * qJD(2);
t230 = m(4) * t147 - mrSges(3,3) * t262 + t149 + t270;
t124 = qJD(1) * t139;
t162 = pkin(2) * t239;
t190 = qJD(2) * t203 - t258;
t73 = qJD(1) * t190 + t162;
t8 = -t103 * t257 + t104 * t256 + t179 * t124 + t181 * t73;
t3 = qJ(5) * t238 + t163 * qJD(5) + t8;
t1 = qJ(6) * t84 + qJD(6) * t131 + t3;
t9 = -t103 * t256 - t104 * t257 + t124 * t181 - t179 * t73;
t2 = -qJ(6) * t83 - qJD(6) * t132 - t238 * t306 - t9;
t229 = t1 * t179 - t181 * t2;
t7 = -pkin(4) * t238 - t9;
t228 = -t179 * t3 + t181 * t7;
t227 = -t179 * t8 - t181 * t9;
t225 = mrSges(5,1) * t179 + mrSges(5,2) * t181;
t223 = mrSges(6,1) * t179 - mrSges(6,3) * t181;
t221 = mrSges(7,1) * t179 - mrSges(7,2) * t181;
t214 = -Ifges(5,2) * t179 + t282;
t205 = pkin(4) * t181 + t272;
t204 = -pkin(4) * t179 + t271;
t172 = pkin(2) * t261;
t97 = t172 + t190;
t15 = -t125 * t256 + t139 * t181 - t151 * t257 - t179 * t97;
t193 = -qJ(3) * t259 - t258;
t192 = pkin(3) * t239 + t140;
t14 = -t125 * t257 + t179 * t139 + t151 * t256 + t181 * t97;
t10 = qJ(5) * t259 + t180 * qJD(5) + t14;
t93 = -mrSges(7,1) * t163 - mrSges(7,3) * t132;
t191 = (t90 + t288) * t181 + (t93 + t287) * t179;
t189 = qJ(5) * t83 + qJD(5) * t132 + t192;
t188 = t9 * mrSges(5,1) - t7 * mrSges(6,1) - t2 * mrSges(7,1) - t8 * mrSges(5,2) + t1 * mrSges(7,2) + t3 * mrSges(6,3);
t160 = Ifges(6,2) * t238;
t159 = Ifges(5,3) * t238;
t144 = qJ(3) - t204;
t142 = t265 * t181;
t141 = t265 * t179;
t138 = t307 * t261;
t134 = (mrSges(4,2) * t182 - mrSges(4,3) * t180) * qJD(1);
t123 = -qJ(3) - t313;
t114 = t172 + t193;
t108 = qJD(4) * t205 + t235;
t100 = qJD(1) * t193 + t162;
t96 = t182 * t205 + t152;
t79 = Ifges(6,4) * t83;
t78 = Ifges(5,5) * t83;
t77 = Ifges(5,6) * t84;
t76 = Ifges(6,6) * t84;
t72 = -t205 * t263 - t264;
t68 = pkin(4) * t132 + t273;
t67 = -t182 * t197 - t152;
t62 = -mrSges(5,2) * t238 - mrSges(5,3) * t84;
t58 = -mrSges(7,1) * t238 - t293;
t56 = -pkin(4) * t180 - t65;
t42 = -pkin(4) * t262 - t52;
t41 = qJ(6) * t267 + t55;
t39 = -t132 * t306 - t273;
t37 = qJ(6) * t269 - t180 * t306 - t65;
t33 = mrSges(5,1) * t84 + mrSges(5,2) * t83;
t31 = mrSges(6,1) * t84 - mrSges(6,3) * t83;
t30 = (qJD(4) * t204 + t254) * t182 + (-t205 - t307) * t261;
t23 = t83 * Ifges(5,4) - t84 * Ifges(5,2) + Ifges(5,6) * t238;
t17 = (qJD(4) * t313 - t254) * t182 + (t197 + t307) * t261;
t12 = -pkin(4) * t259 - t15;
t11 = pkin(4) * t84 - t189;
t6 = t182 * t252 + (-t179 * t255 - t180 * t260) * qJ(6) + t10;
t5 = -qJ(6) * t241 + (qJ(6) * t256 - qJD(2) * t306 + t253) * t182 - t15;
t4 = -t306 * t84 + t189;
t21 = [m(5) * (-t119 * t138 + t14 * t36 + t15 * t35 - t152 * t192 + t65 * t9 + t66 * t8) + t152 * t33 + t114 * t134 - t138 * t71 + t6 * t90 + t14 * t91 + t10 * t92 + t5 * t93 + t15 * t94 + t12 * t95 + t96 * t31 + t17 * t70 + t56 * t60 + t41 * t61 + t65 * t59 + t66 * t62 + t67 * t32 + t30 * t69 + t55 * t57 + t37 * t58 + (t23 * t295 + t213 * t308 - t192 * t226 + t11 * t224 - t4 * t222 + t100 * mrSges(4,2) - t242 * t140 + (t1 * t181 + t179 * t2) * mrSges(7,3) + (t179 * t9 - t181 * t8) * mrSges(5,3) + (-t179 * t7 - t181 * t3) * mrSges(6,2) + (-t119 * t225 - t38 * t223 + t20 * t221 + (-Ifges(7,5) * t181 - Ifges(7,6) * t179) * t299 + t214 * t302 + t48 * t296 + (t13 * t181 - t16 * t179) * mrSges(7,3) + (t179 * t36 + t181 * t35) * mrSges(5,3) + (t179 * t29 - t181 * t27) * mrSges(6,2) + t331 * t303 + (t179 * t322 - t181 * t337) * t298 + t328 * t301 + t316 * t295) * qJD(4) + (t231 * pkin(7) + (-t145 * mrSges(4,3) - 0.2e1 * t304 - t244 * t182 - t233 * t267 + t234 * t269 + (0.3e1 / 0.2e1 * Ifges(3,1) + Ifges(7,3) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + t242 * pkin(7) + t245) * t180) * qJD(1) + t248 * qJD(2) + t311) * qJD(2) + t318 * t309 + t336 * t294 - t314 * t83 / 0.2e1 + (qJD(4) * t321 + t329) * t297) * t182 + (t159 / 0.2e1 + t160 / 0.2e1 + t79 / 0.2e1 - t77 / 0.2e1 + t78 / 0.2e1 + t76 / 0.2e1 - t100 * mrSges(4,3) + (-Ifges(7,6) - t243) * t84 + (-Ifges(7,5) + t247) * t83 + (t230 * pkin(7) + t246 * qJD(2) + (-t145 * mrSges(4,2) + t180 * t244 - 0.2e1 * t305) * qJD(1) + t342) * qJD(2) + t188) * t180 + m(4) * (t100 * t145 + t114 * t120) + m(6) * (t10 * t29 + t11 * t96 + t12 * t27 + t3 * t55 + t30 * t38 + t56 * t7) + m(7) * (t1 * t41 + t13 * t5 + t16 * t6 + t17 * t20 + t2 * t37 + t4 * t67); ((-m(5) * t200 + m(6) * t201 + t179 * t287 + t181 * t288) * t184 + t206 * t298 + t213 * t302 + t318 * t303 + t317 * t299 + t314 * t301 + t344) * qJD(4) + (t1 * t141 + t123 * t4 + t13 * t334 - t142 * t2 + t16 * t333 + t20 * t335) * m(7) + t335 * t70 + t336 * t296 + t333 * t90 + t334 * t93 + t331 * t308 + t328 * t83 / 0.2e1 + t329 * t294 + (-m(4) * t120 - t134) * (-qJ(3) * t262 + t169) - t192 * t225 + m(5) * (-qJ(3) * t192 + t119 * qJD(3) + t266 * t9 + t268 * t8) + t264 * t71 - m(5) * (-t119 * t264 + t35 * t52 + t36 * t53) - t140 * mrSges(4,3) + t141 * t61 - t142 * t58 + t144 * t31 + t123 * t32 - t53 * t91 - t40 * t92 - t52 * t94 - t42 * t95 + qJ(3) * t33 + m(4) * (-qJ(3) * t140 - qJD(3) * t147) + t228 * mrSges(6,2) + t229 * mrSges(7,3) + t11 * t223 + t227 * mrSges(5,3) - t4 * t221 + ((-t166 / 0.2e1 + ((-m(4) * pkin(2) + t338) * qJD(2) - t231) * pkin(7) + (-pkin(2) * mrSges(4,1) - t179 * t233 - t181 * t234 + t248) * qJD(2) + (t304 - t277 / 0.2e1) * qJD(1) - t311) * t182 + ((-t230 + t270) * pkin(7) + (t305 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t180) * qJD(1) + (-qJ(3) * mrSges(4,1) + t246) * qJD(2) + (Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t262 - t342) * t180) * qJD(1) + (t108 - t72) * t69 + m(6) * (t38 * t108 + t11 * t144 - t266 * t7 + t268 * t3) + t274 * qJD(3) - m(6) * (t27 * t42 + t29 * t40 + t38 * t72) + (t292 * t181 + (t57 + t62) * t179) * t184 + t214 * t309 + t23 * t297; (-t58 + t292) * t181 + (t62 + t291) * t179 + (-t274 - t290) * qJD(2) + t191 * qJD(4) + (t242 * t259 + (t134 + t191) * t180) * qJD(1) - m(4) * (-qJD(2) * t147 - t120 * t263) + (qJD(2) * t20 + t163 * t202 + t229) * m(7) + (-qJD(2) * t38 + t163 * t201 - t228) * m(6) + (-qJD(2) * t119 - t163 * t200 - t227) * m(5); -t306 * t58 + (qJ(5) * t1 - t13 * t19 + t16 * t320 - t2 * t306 - t20 * t39) * m(7) + t159 + t160 - t119 * (mrSges(5,1) * t132 - mrSges(5,2) * t131) - t38 * (mrSges(6,1) * t132 + mrSges(6,3) * t131) - t20 * (-mrSges(7,1) * t132 - mrSges(7,2) * t131) + t79 - t77 + t78 + t76 - t18 * t90 - t19 * t93 - Ifges(7,5) * t83 - Ifges(7,6) * t84 - t39 * t70 - pkin(4) * t60 - t68 * t69 + (-t131 * t337 - t132 * t322) * t299 + (-t131 * t339 - t284 + t321 + t332) * t301 + t188 + (-Ifges(5,2) * t132 - t128 + t316) * t302 + (-pkin(4) * t7 + qJ(5) * t3 - t27 * t36 + t29 * t319 - t38 * t68) * m(6) + (-t13 * t131 - t132 * t16) * mrSges(7,3) + (t131 * t27 + t132 * t29) * mrSges(6,2) + (t285 - t287) * t36 + (-t286 - t288) * t35 + t289 * qJD(5) + t291 * qJ(5) + (-Ifges(7,5) * t131 + Ifges(7,6) * t132) * t298 + t48 * t300 + (t132 * t343 - t330) * t303 + Ifges(7,3) * t238; -t293 + t75 - t289 * t163 + t290 * t132 + (-mrSges(6,1) - mrSges(7,1)) * t238 + (-t132 * t20 - t16 * t163 + t2) * m(7) + (t132 * t38 - t163 * t29 + t7) * m(6); -t131 * t90 + t132 * t93 + 0.2e1 * (t4 / 0.2e1 + t13 * t300 + t16 * t303) * m(7) + t32;];
tauc  = t21(:);
