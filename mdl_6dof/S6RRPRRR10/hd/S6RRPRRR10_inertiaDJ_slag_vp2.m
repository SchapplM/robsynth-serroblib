% Calculate time derivative of joint inertia matrix for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:25
% EndTime: 2019-03-09 14:18:43
% DurationCPUTime: 7.57s
% Computational Cost: add. (16832->700), mult. (42850->1038), div. (0->0), fcn. (44514->12), ass. (0->281)
t257 = sin(pkin(6));
t355 = 0.2e1 * t257;
t256 = sin(pkin(12));
t258 = cos(pkin(12));
t262 = sin(qJ(4));
t266 = cos(qJ(4));
t223 = t256 * t266 + t258 * t262;
t261 = sin(qJ(5));
t296 = qJD(5) * t261;
t222 = t256 * t262 - t258 * t266;
t216 = t222 * qJD(4);
t265 = cos(qJ(5));
t301 = t265 * t216;
t268 = t223 * t296 + t301;
t251 = -pkin(3) * t258 - pkin(2);
t167 = pkin(4) * t222 - pkin(10) * t223 + t251;
t316 = pkin(9) + qJ(3);
t235 = t316 * t256;
t236 = t316 * t258;
t188 = -t235 * t262 + t236 * t266;
t175 = t265 * t188;
t118 = t167 * t261 + t175;
t354 = qJD(5) * t118;
t260 = sin(qJ(6));
t264 = cos(qJ(6));
t274 = t260 * t261 - t264 * t265;
t325 = -t274 / 0.2e1;
t225 = t260 * t265 + t261 * t264;
t324 = t225 / 0.2e1;
t353 = t261 / 0.2e1;
t318 = t265 / 0.2e1;
t158 = t274 * t223;
t352 = -t235 * t266 - t236 * t262;
t237 = -mrSges(6,1) * t265 + mrSges(6,2) * t261;
t351 = -m(6) * pkin(4) + t237;
t279 = mrSges(6,1) * t261 + mrSges(6,2) * t265;
t226 = t279 * qJD(5);
t350 = qJD(5) + qJD(6);
t349 = 2 * m(4);
t348 = 2 * m(5);
t347 = 0.2e1 * m(6);
t346 = 2 * m(7);
t345 = -2 * mrSges(3,3);
t344 = -2 * mrSges(5,3);
t150 = qJD(3) * t223 + qJD(4) * t188;
t343 = 0.2e1 * t150;
t342 = -0.2e1 * t352;
t259 = cos(pkin(6));
t263 = sin(qJ(2));
t304 = t257 * t263;
t214 = -t256 * t304 + t258 * t259;
t215 = t256 * t259 + t258 * t304;
t153 = t214 * t262 + t215 * t266;
t267 = cos(qJ(2));
t303 = t257 * t267;
t139 = -t153 * t261 - t265 * t303;
t270 = -t153 * t265 + t261 * t303;
t275 = t214 * t266 - t215 * t262;
t63 = -Ifges(6,4) * t270 + Ifges(6,2) * t139 - Ifges(6,6) * t275;
t340 = -t63 / 0.2e1;
t299 = qJD(2) * t257;
t286 = t267 * t299;
t120 = qJD(4) * t275 - t222 * t286;
t287 = t263 * t299;
t75 = qJD(5) * t270 - t120 * t261 + t265 * t287;
t339 = t75 / 0.2e1;
t84 = t139 * t264 + t260 * t270;
t338 = t84 / 0.2e1;
t85 = t139 * t260 - t264 * t270;
t337 = t85 / 0.2e1;
t217 = t223 * qJD(4);
t295 = qJD(5) * t265;
t269 = -t216 * t261 + t223 * t295;
t92 = -Ifges(6,1) * t268 - Ifges(6,4) * t269 + Ifges(6,5) * t217;
t336 = t92 / 0.2e1;
t335 = -pkin(11) - pkin(10);
t334 = t139 / 0.2e1;
t313 = Ifges(6,4) * t261;
t278 = Ifges(6,1) * t265 - t313;
t143 = Ifges(6,5) * t222 + t223 * t278;
t333 = t143 / 0.2e1;
t157 = t225 * t223;
t332 = -t157 / 0.2e1;
t331 = -t158 / 0.2e1;
t178 = t350 * t274;
t330 = -t178 / 0.2e1;
t179 = t350 * t225;
t329 = -t179 / 0.2e1;
t185 = Ifges(7,4) * t225 - Ifges(7,2) * t274;
t328 = t185 / 0.2e1;
t186 = Ifges(7,1) * t225 - Ifges(7,4) * t274;
t327 = t186 / 0.2e1;
t326 = -t223 / 0.2e1;
t229 = t278 * qJD(5);
t323 = t229 / 0.2e1;
t239 = Ifges(6,2) * t265 + t313;
t322 = -t239 / 0.2e1;
t312 = Ifges(6,4) * t265;
t240 = Ifges(6,1) * t261 + t312;
t321 = t240 / 0.2e1;
t320 = t258 / 0.2e1;
t319 = -t261 / 0.2e1;
t317 = pkin(1) * t267;
t219 = pkin(1) * t259 * t263 + pkin(8) * t303;
t203 = qJ(3) * t259 + t219;
t204 = (-pkin(2) * t267 - qJ(3) * t263 - pkin(1)) * t257;
t147 = -t203 * t256 + t204 * t258;
t109 = -pkin(3) * t303 - pkin(9) * t215 + t147;
t148 = t203 * t258 + t204 * t256;
t125 = pkin(9) * t214 + t148;
t70 = t109 * t262 + t125 * t266;
t61 = -pkin(10) * t303 + t70;
t246 = pkin(8) * t304;
t206 = t246 + (-pkin(2) - t317) * t259;
t161 = -pkin(3) * t214 + t206;
t83 = -pkin(4) * t275 - pkin(10) * t153 + t161;
t36 = t261 * t83 + t265 * t61;
t315 = Ifges(4,4) * t256;
t314 = Ifges(4,4) * t258;
t311 = Ifges(6,6) * t261;
t208 = t219 * qJD(2);
t280 = t256 * t286;
t182 = pkin(3) * t280 + t208;
t310 = t182 * mrSges(5,1);
t309 = t182 * mrSges(5,2);
t308 = t150 * t352;
t307 = t223 * t261;
t306 = t223 * t265;
t302 = t258 * t267;
t127 = -Ifges(7,5) * t178 - Ifges(7,6) * t179;
t192 = (-qJD(3) * t263 + (pkin(2) * t263 - qJ(3) * t267) * qJD(2)) * t257;
t292 = t259 * t317;
t207 = -pkin(8) * t287 + qJD(2) * t292;
t199 = qJD(3) * t259 + t207;
t136 = t192 * t256 + t199 * t258;
t300 = -Ifges(5,5) * t216 - Ifges(5,6) * t217;
t195 = mrSges(4,2) * t258 * t286 + mrSges(4,1) * t280;
t298 = qJD(4) * t262;
t297 = qJD(4) * t266;
t294 = qJD(6) * t260;
t293 = qJD(6) * t264;
t121 = qJD(4) * t153 + t223 * t286;
t74 = qJD(5) * t139 + t120 * t265 + t261 * t287;
t22 = qJD(6) * t84 + t260 * t75 + t264 * t74;
t23 = -qJD(6) * t85 - t260 * t74 + t264 * t75;
t6 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t121;
t25 = Ifges(6,5) * t74 + Ifges(6,6) * t75 + Ifges(6,3) * t121;
t78 = -t179 * t223 + t216 * t274;
t79 = t158 * t350 + t225 * t216;
t37 = Ifges(7,5) * t78 + Ifges(7,6) * t79 + Ifges(7,3) * t217;
t291 = pkin(5) * t296;
t64 = -Ifges(6,1) * t270 + Ifges(6,4) * t139 - Ifges(6,5) * t275;
t290 = t64 * t318;
t289 = Ifges(5,5) * t120 - Ifges(5,6) * t121 + Ifges(5,3) * t287;
t288 = qJD(5) * t335;
t253 = Ifges(6,5) * t295;
t284 = t127 / 0.2e1 - Ifges(6,6) * t296 / 0.2e1 + t253 / 0.2e1;
t283 = Ifges(6,5) * t353 + Ifges(7,5) * t324 + Ifges(6,6) * t318 + Ifges(7,6) * t325;
t35 = -t261 * t61 + t265 * t83;
t71 = mrSges(5,1) * t121 + mrSges(5,2) * t120;
t162 = mrSges(5,1) * t217 - mrSges(5,2) * t216;
t126 = mrSges(7,1) * t179 - t178 * mrSges(7,2);
t69 = t109 * t266 - t125 * t262;
t149 = -qJD(3) * t222 + qJD(4) * t352;
t166 = pkin(4) * t217 + pkin(10) * t216;
t282 = -t149 * t261 + t166 * t265;
t117 = t167 * t265 - t188 * t261;
t135 = t192 * t258 - t199 * t256;
t60 = pkin(4) * t303 - t69;
t277 = -Ifges(6,2) * t261 + t312;
t107 = (pkin(3) * t263 - pkin(9) * t302) * t299 + t135;
t124 = -pkin(9) * t280 + t136;
t32 = t107 * t262 + t109 * t297 + t124 * t266 - t125 * t298;
t30 = pkin(10) * t287 + t32;
t54 = t121 * pkin(4) - t120 * pkin(10) + t182;
t10 = t261 * t54 + t265 * t30 + t295 * t83 - t296 * t61;
t11 = -qJD(5) * t36 - t261 * t30 + t265 * t54;
t276 = t10 * t265 - t11 * t261;
t24 = -pkin(5) * t275 + pkin(11) * t270 + t35;
t28 = pkin(11) * t139 + t36;
t12 = t24 * t264 - t260 * t28;
t13 = t24 * t260 + t264 * t28;
t87 = pkin(5) * t222 - pkin(11) * t306 + t117;
t98 = -pkin(11) * t307 + t118;
t50 = -t260 * t98 + t264 * t87;
t51 = t260 * t87 + t264 * t98;
t241 = t335 * t261;
t242 = t335 * t265;
t193 = t241 * t264 + t242 * t260;
t194 = t241 * t260 - t242 * t264;
t230 = t261 * t288;
t231 = t265 * t288;
t131 = qJD(6) * t193 + t230 * t264 + t231 * t260;
t132 = -qJD(6) * t194 - t230 * t260 + t231 * t264;
t273 = mrSges(7,1) * t132 - t131 * mrSges(7,2) + t127;
t33 = t107 * t266 - t109 * t298 - t124 * t262 - t125 * t297;
t4 = pkin(5) * t121 - pkin(11) * t74 + t11;
t5 = pkin(11) * t75 + t10;
t2 = qJD(6) * t12 + t260 * t4 + t264 * t5;
t3 = -qJD(6) * t13 - t260 * t5 + t264 * t4;
t272 = mrSges(7,1) * t3 - t2 * mrSges(7,2) + t6;
t45 = pkin(11) * t301 + pkin(5) * t217 + (-t175 + (pkin(11) * t223 - t167) * t261) * qJD(5) + t282;
t55 = t149 * t265 + t166 * t261 + t167 * t295 - t188 * t296;
t49 = -pkin(11) * t269 + t55;
t15 = qJD(6) * t50 + t260 * t45 + t264 * t49;
t16 = -qJD(6) * t51 - t260 * t49 + t264 * t45;
t271 = mrSges(7,1) * t16 - mrSges(7,2) * t15 + t37;
t31 = -pkin(4) * t287 - t33;
t90 = -Ifges(6,5) * t268 - Ifges(6,6) * t269 + Ifges(6,3) * t217;
t252 = -pkin(5) * t265 - pkin(4);
t244 = Ifges(3,5) * t286;
t228 = t277 * qJD(5);
t220 = (-mrSges(7,1) * t260 - mrSges(7,2) * t264) * qJD(6) * pkin(5);
t218 = -t246 + t292;
t201 = (mrSges(4,1) * t263 - mrSges(4,3) * t302) * t299;
t200 = (-mrSges(4,3) * t256 * t267 - mrSges(4,2) * t263) * t299;
t191 = -mrSges(4,1) * t303 - mrSges(4,3) * t215;
t190 = mrSges(4,2) * t303 + mrSges(4,3) * t214;
t183 = mrSges(7,1) * t274 + mrSges(7,2) * t225;
t177 = Ifges(5,1) * t223 - Ifges(5,4) * t222;
t176 = Ifges(5,4) * t223 - Ifges(5,2) * t222;
t171 = (t263 * Ifges(4,5) + (Ifges(4,1) * t258 - t315) * t267) * t299;
t170 = (t263 * Ifges(4,6) + (-Ifges(4,2) * t256 + t314) * t267) * t299;
t169 = mrSges(6,1) * t222 - mrSges(6,3) * t306;
t168 = -mrSges(6,2) * t222 - mrSges(6,3) * t307;
t165 = t279 * t223;
t164 = -Ifges(5,1) * t216 - Ifges(5,4) * t217;
t163 = -Ifges(5,4) * t216 - Ifges(5,2) * t217;
t151 = pkin(5) * t307 - t352;
t145 = -mrSges(5,1) * t303 - mrSges(5,3) * t153;
t144 = mrSges(5,2) * t303 + mrSges(5,3) * t275;
t142 = Ifges(6,6) * t222 + t223 * t277;
t141 = Ifges(6,3) * t222 + (Ifges(6,5) * t265 - t311) * t223;
t138 = mrSges(7,1) * t222 + mrSges(7,3) * t158;
t137 = -mrSges(7,2) * t222 - mrSges(7,3) * t157;
t134 = -mrSges(6,2) * t217 - mrSges(6,3) * t269;
t133 = mrSges(6,1) * t217 + mrSges(6,3) * t268;
t129 = -Ifges(7,1) * t178 - Ifges(7,4) * t179;
t128 = -Ifges(7,4) * t178 - Ifges(7,2) * t179;
t104 = mrSges(6,1) * t269 - mrSges(6,2) * t268;
t102 = -mrSges(5,2) * t287 - mrSges(5,3) * t121;
t101 = mrSges(5,1) * t287 - mrSges(5,3) * t120;
t100 = mrSges(7,1) * t157 - mrSges(7,2) * t158;
t99 = pkin(5) * t269 + t150;
t97 = Ifges(5,1) * t153 + Ifges(5,4) * t275 - Ifges(5,5) * t303;
t96 = Ifges(5,4) * t153 + Ifges(5,2) * t275 - Ifges(5,6) * t303;
t95 = -Ifges(7,1) * t158 - Ifges(7,4) * t157 + Ifges(7,5) * t222;
t94 = -Ifges(7,4) * t158 - Ifges(7,2) * t157 + Ifges(7,6) * t222;
t93 = -Ifges(7,5) * t158 - Ifges(7,6) * t157 + Ifges(7,3) * t222;
t91 = -Ifges(6,4) * t268 - Ifges(6,2) * t269 + Ifges(6,6) * t217;
t89 = -mrSges(6,1) * t275 + mrSges(6,3) * t270;
t88 = mrSges(6,2) * t275 + mrSges(6,3) * t139;
t86 = -mrSges(6,1) * t139 - mrSges(6,2) * t270;
t68 = -mrSges(7,2) * t217 + mrSges(7,3) * t79;
t67 = mrSges(7,1) * t217 - mrSges(7,3) * t78;
t66 = Ifges(5,1) * t120 - Ifges(5,4) * t121 + Ifges(5,5) * t287;
t65 = Ifges(5,4) * t120 - Ifges(5,2) * t121 + Ifges(5,6) * t287;
t62 = -Ifges(6,5) * t270 + Ifges(6,6) * t139 - Ifges(6,3) * t275;
t58 = -mrSges(7,1) * t275 - mrSges(7,3) * t85;
t57 = mrSges(7,2) * t275 + mrSges(7,3) * t84;
t56 = t282 - t354;
t48 = -pkin(5) * t139 + t60;
t47 = -mrSges(6,2) * t121 + mrSges(6,3) * t75;
t46 = mrSges(6,1) * t121 - mrSges(6,3) * t74;
t44 = -mrSges(7,1) * t84 + mrSges(7,2) * t85;
t43 = -mrSges(7,1) * t79 + mrSges(7,2) * t78;
t42 = Ifges(7,1) * t85 + Ifges(7,4) * t84 - Ifges(7,5) * t275;
t41 = Ifges(7,4) * t85 + Ifges(7,2) * t84 - Ifges(7,6) * t275;
t40 = Ifges(7,5) * t85 + Ifges(7,6) * t84 - Ifges(7,3) * t275;
t39 = Ifges(7,1) * t78 + Ifges(7,4) * t79 + Ifges(7,5) * t217;
t38 = Ifges(7,4) * t78 + Ifges(7,2) * t79 + Ifges(7,6) * t217;
t34 = -mrSges(6,1) * t75 + mrSges(6,2) * t74;
t27 = Ifges(6,1) * t74 + Ifges(6,4) * t75 + Ifges(6,5) * t121;
t26 = Ifges(6,4) * t74 + Ifges(6,2) * t75 + Ifges(6,6) * t121;
t19 = -pkin(5) * t75 + t31;
t18 = -mrSges(7,2) * t121 + mrSges(7,3) * t23;
t17 = mrSges(7,1) * t121 - mrSges(7,3) * t22;
t9 = -mrSges(7,1) * t23 + mrSges(7,2) * t22;
t8 = Ifges(7,1) * t22 + Ifges(7,4) * t23 + Ifges(7,5) * t121;
t7 = Ifges(7,4) * t22 + Ifges(7,2) * t23 + Ifges(7,6) * t121;
t1 = [(t12 * t3 + t13 * t2 + t19 * t48) * t346 + (t10 * t36 + t11 * t35 + t31 * t60) * t347 + (t161 * t182 + t32 * t70 + t33 * t69) * t348 + (t135 * t147 + t136 * t148 + t206 * t208) * t349 + 0.2e1 * m(3) * (t207 * t219 - t208 * t218) + (t40 + t62 - t96) * t121 + t259 * t244 - t270 * t27 - (t25 + t6 - t65 + 0.2e1 * t310) * t275 - t289 * t303 + 0.2e1 * t207 * (-mrSges(3,2) * t259 + mrSges(3,3) * t303) - 0.2e1 * t208 * (mrSges(3,1) * t259 - mrSges(3,3) * t304) + ((t219 * t345 + Ifges(4,5) * t215 + Ifges(5,5) * t153 - 0.2e1 * Ifges(3,6) * t259 + Ifges(4,6) * t214 + Ifges(5,6) * t275 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t263) * t355) * t263 + (-t256 * (Ifges(4,4) * t215 + Ifges(4,2) * t214) + t218 * t345 + t258 * (Ifges(4,1) * t215 + Ifges(4,4) * t214) + Ifges(3,5) * t259 + (-pkin(1) * mrSges(3,2) + (-Ifges(4,5) * t258 + Ifges(4,6) * t256 + Ifges(3,4)) * t267) * t355 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) - Ifges(5,3)) * t304) * t267) * t299 + (t66 + 0.2e1 * t309) * t153 + 0.2e1 * t12 * t17 + 0.2e1 * t13 * t18 + t23 * t41 + t22 * t42 + 0.2e1 * t19 * t44 + 0.2e1 * t35 * t46 + 0.2e1 * t36 * t47 + 0.2e1 * t48 * t9 + 0.2e1 * t2 * t57 + 0.2e1 * t3 * t58 + 0.2e1 * t60 * t34 + t74 * t64 + t75 * t63 + t84 * t7 + t85 * t8 + 0.2e1 * t31 * t86 + 0.2e1 * t10 * t88 + 0.2e1 * t11 * t89 + 0.2e1 * t69 * t101 + 0.2e1 * t70 * t102 + t120 * t97 + t139 * t26 + 0.2e1 * t32 * t144 + 0.2e1 * t33 * t145 + 0.2e1 * t161 * t71 + 0.2e1 * t136 * t190 + 0.2e1 * t135 * t191 + 0.2e1 * t148 * t200 + 0.2e1 * t147 * t201 + 0.2e1 * t206 * t195 + t214 * t170 + t215 * t171 + 0.2e1 * t208 * (-mrSges(4,1) * t214 + mrSges(4,2) * t215); -(t34 - t101) * t352 + m(6) * (t10 * t118 + t11 * t117 + t150 * t60 - t31 * t352 + t35 * t56 + t36 * t55) + m(5) * (t149 * t70 - t150 * t69 + t182 * t251 + t188 * t32 + t33 * t352) + t142 * t339 + t8 * t331 + t7 * t332 + t74 * t333 + t91 * t334 + t39 * t337 + t38 * t338 + m(7) * (t12 * t16 + t13 * t15 + t151 * t19 + t2 * t51 + t3 * t50 + t48 * t99) + (t170 / 0.2e1 - t208 * mrSges(4,1) + t136 * mrSges(4,3) + qJD(3) * t190 + qJ(3) * t200) * t258 + (t171 / 0.2e1 + t208 * mrSges(4,2) - t135 * mrSges(4,3) - qJD(3) * t191 - qJ(3) * t201) * t256 + (t86 - t145) * t150 + (t40 / 0.2e1 + t62 / 0.2e1 - t96 / 0.2e1 - t70 * mrSges(5,3)) * t217 + (t93 / 0.2e1 + t141 / 0.2e1 - t176 / 0.2e1) * t121 - t270 * t336 - (t37 / 0.2e1 + t90 / 0.2e1 - t163 / 0.2e1) * t275 + t244 + (t6 / 0.2e1 + t25 / 0.2e1 - t65 / 0.2e1 + t310 - t32 * mrSges(5,3)) * t222 - (t290 + t97 / 0.2e1 - t69 * mrSges(5,3) + t63 * t319) * t216 + (-t267 * t300 / 0.2e1 + ((-Ifges(3,6) + Ifges(5,5) * t223 / 0.2e1 - Ifges(5,6) * t222 / 0.2e1 + Ifges(4,5) * t256 / 0.2e1 + Ifges(4,6) * t320) * t263 + (-t256 * (Ifges(4,2) * t258 + t315) / 0.2e1 + (Ifges(4,1) * t256 + t314) * t320) * t267) * qJD(2)) * t257 + (t27 * t318 + t66 / 0.2e1 + t309 - t33 * mrSges(5,3) + t26 * t319 + (t265 * t340 + t319 * t64) * qJD(5)) * t223 + t48 * t43 + t50 * t17 + t51 * t18 + t15 * t57 + t16 * t58 + t12 * t67 + t13 * t68 + t78 * t42 / 0.2e1 + t79 * t41 / 0.2e1 + t55 * t88 + t56 * t89 + t23 * t94 / 0.2e1 + t22 * t95 / 0.2e1 + t99 * t44 + t19 * t100 + t60 * t104 + t117 * t46 + t118 * t47 + t35 * t133 + t36 * t134 + t2 * t137 + t3 * t138 + t149 * t144 + t151 * t9 + t161 * t162 + t153 * t164 / 0.2e1 + t31 * t165 + t10 * t168 + t11 * t169 + t120 * t177 / 0.2e1 + t188 * t102 - pkin(2) * t195 - t207 * mrSges(3,2) - t208 * mrSges(3,1) + t251 * t71 + m(4) * (-pkin(2) * t208 + (-t147 * t256 + t148 * t258) * qJD(3) + (-t135 * t256 + t136 * t258) * qJ(3)); (qJ(3) * t349 + 0.2e1 * mrSges(4,3)) * (t256 ^ 2 + t258 ^ 2) * qJD(3) + t104 * t342 + t165 * t343 + (t15 * t51 + t151 * t99 + t16 * t50) * t346 + (t117 * t56 + t118 * t55 - t308) * t347 + (t149 * t188 - t308) * t348 - (mrSges(5,3) * t342 - t142 * t261 + t143 * t265 + t177) * t216 + (mrSges(5,3) * t343 - t261 * t91 + t265 * t92 + t164 + (-t142 * t265 - t143 * t261) * qJD(5)) * t223 + (t188 * t344 + t141 - t176 + t93) * t217 + (t149 * t344 - t163 + t37 + t90) * t222 + 0.2e1 * t50 * t67 + 0.2e1 * t51 * t68 + t79 * t94 + t78 * t95 + 0.2e1 * t99 * t100 + 0.2e1 * t117 * t133 + 0.2e1 * t118 * t134 + 0.2e1 * t15 * t137 + 0.2e1 * t16 * t138 + 0.2e1 * t151 * t43 - t157 * t38 - t158 * t39 + 0.2e1 * t55 * t168 + 0.2e1 * t56 * t169 + 0.2e1 * t251 * t162; -t274 * t17 - t178 * t57 - t179 * t58 + t225 * t18 + t261 * t47 + t265 * t46 + (-t261 * t89 + t265 * t88) * qJD(5) + m(7) * (-t12 * t179 - t13 * t178 + t2 * t225 - t274 * t3) + m(6) * (t10 * t261 + t11 * t265 + (-t261 * t35 + t265 * t36) * qJD(5)) + m(5) * t182 + m(4) * t208 + t71 + t195; t265 * t133 + t261 * t134 - t178 * t137 - t179 * t138 - t274 * t67 + t225 * t68 + (t168 * t265 - t169 * t261) * qJD(5) + m(7) * (t15 * t225 - t16 * t274 - t178 * t51 - t179 * t50) + m(6) * (t261 * t55 + t265 * t56 + (-t117 * t261 + t118 * t265) * qJD(5)) + t162; (-t178 * t225 + t179 * t274) * t346; t283 * t121 + ((-t261 * t36 - t265 * t35) * qJD(5) + t276) * mrSges(6,3) + t239 * t339 + t74 * t321 + t8 * t324 + t7 * t325 + t22 * t327 + t23 * t328 + t41 * t329 + t42 * t330 + t228 * t334 + t129 * t337 + t128 * t338 + t26 * t318 + t27 * t353 + (t12 * t178 - t13 * t179 - t2 * t274 - t225 * t3) * mrSges(7,3) - t270 * t323 - t284 * t275 + t289 + t351 * t31 + m(7) * (t12 * t132 + t13 * t131 + t19 * t252 + t193 * t3 + t194 * t2 + t291 * t48) + (m(6) * (-t295 * t35 - t296 * t36 + t276) - t89 * t295 - t88 * t296 - t261 * t46 + t265 * t47) * pkin(10) + (t290 + (pkin(5) * t44 + t340) * t261) * qJD(5) - t32 * mrSges(5,2) + t33 * mrSges(5,1) - pkin(4) * t34 + t48 * t126 + t131 * t57 + t132 * t58 + t19 * t183 + t193 * t17 + t194 * t18 + t60 * t226 + t252 * t9; t300 - t352 * t226 + t283 * t217 + t284 * t222 + t39 * t324 + t38 * t325 + t78 * t327 + t79 * t328 + t94 * t329 + t95 * t330 + t129 * t331 + t128 * t332 + (t91 / 0.2e1 + t55 * mrSges(6,3) + t223 * t323 - t216 * t321 + (-t117 * mrSges(6,3) + t223 * t322 + t333) * qJD(5) + (-qJD(5) * t169 + m(6) * (-qJD(5) * t117 + t55) + t134) * pkin(10)) * t265 + (t336 - t56 * mrSges(6,3) + t228 * t326 - t216 * t322 + (-t118 * mrSges(6,3) + t240 * t326 - t142 / 0.2e1 + (m(7) * t151 + t100) * pkin(5)) * qJD(5) + (-qJD(5) * t168 - t133 + m(6) * (-t56 - t354)) * pkin(10)) * t261 + (-t15 * t274 - t16 * t225 + t178 * t50 - t179 * t51) * mrSges(7,3) + m(7) * (t131 * t51 + t132 * t50 + t15 * t194 + t16 * t193 + t252 * t99) + (-mrSges(5,1) + t351) * t150 - pkin(4) * t104 + t131 * t137 + t132 * t138 - t149 * mrSges(5,2) + t151 * t126 + t99 * t183 + t193 * t67 + t194 * t68 + t252 * t43; m(7) * (t131 * t225 - t132 * t274 - t178 * t194 - t179 * t193); -t178 * t186 + t225 * t129 - t179 * t185 - t274 * t128 + 0.2e1 * t183 * t291 + 0.2e1 * t252 * t126 + (t131 * t194 + t132 * t193 + t252 * t291) * t346 - t239 * t296 - 0.2e1 * pkin(4) * t226 + t261 * t229 + (qJD(5) * t240 + t228) * t265 + 0.2e1 * (-t131 * t274 - t132 * t225 + t178 * t193 - t179 * t194) * mrSges(7,3); t11 * mrSges(6,1) - t10 * mrSges(6,2) + (m(7) * (-t12 * t294 + t13 * t293 + t2 * t260 + t264 * t3) + t57 * t293 + t260 * t18 - t58 * t294 + t264 * t17) * pkin(5) + t272 + t25; t56 * mrSges(6,1) - t55 * mrSges(6,2) + (m(7) * (t15 * t260 + t16 * t264 + t293 * t51 - t294 * t50) + t137 * t293 + t260 * t68 - t138 * t294 + t264 * t67) * pkin(5) + t90 + t271; -t226 + m(7) * (-t178 * t260 - t179 * t264 + (t225 * t264 + t260 * t274) * qJD(6)) * pkin(5) - t126; t253 + (pkin(10) * t237 - t311) * qJD(5) + (m(7) * (t131 * t260 + t132 * t264 + (-t193 * t260 + t194 * t264) * qJD(6)) + (t264 * t178 - t260 * t179 + (t225 * t260 - t264 * t274) * qJD(6)) * mrSges(7,3)) * pkin(5) + t273; 0.2e1 * t220; t272; t271; -t126; t273; t220; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
