% Calculate time derivative of joint inertia matrix for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:40
% EndTime: 2019-12-31 22:47:01
% DurationCPUTime: 7.05s
% Computational Cost: add. (11294->686), mult. (34355->1018), div. (0->0), fcn. (33891->12), ass. (0->293)
t239 = sin(pkin(5));
t360 = 0.2e1 * t239;
t249 = cos(qJ(2));
t241 = cos(pkin(5));
t322 = pkin(1) * t241;
t233 = t249 * t322;
t226 = qJD(2) * t233;
t245 = sin(qJ(2));
t240 = cos(pkin(6));
t267 = t239 * (-pkin(9) * t240 - pkin(8));
t256 = t245 * t267;
t147 = pkin(2) * t241 + t233 + t256;
t307 = t147 * t240;
t359 = qJD(2) * t256 + qJD(3) * t307 + t226;
t242 = sin(qJ(5));
t246 = cos(qJ(5));
t243 = sin(qJ(4));
t290 = qJD(5) * t243;
t247 = cos(qJ(4));
t292 = qJD(4) * t247;
t251 = -t242 * t290 + t246 * t292;
t244 = sin(qJ(3));
t303 = t240 * t244;
t238 = sin(pkin(6));
t248 = cos(qJ(3));
t305 = t238 * t248;
t194 = pkin(2) * t303 + pkin(9) * t305;
t175 = pkin(10) * t240 + t194;
t176 = (-pkin(3) * t248 - pkin(10) * t244 - pkin(2)) * t238;
t357 = t247 * t175 + t243 * t176;
t358 = qJD(4) * t357;
t325 = t242 / 0.2e1;
t323 = t246 / 0.2e1;
t356 = -m(5) * pkin(3) - mrSges(5,1) * t247 + mrSges(5,2) * t243;
t306 = t238 * t244;
t228 = pkin(9) * t306;
t302 = t240 * t248;
t192 = pkin(2) * t302 - t228;
t298 = t248 * t249;
t301 = t244 * t245;
t355 = t240 * t298 - t301;
t321 = pkin(9) * t238;
t165 = (-pkin(2) * t249 - t245 * t321 - pkin(1)) * t239;
t100 = -t147 * t238 + t240 * t165;
t140 = -t355 * t239 - t241 * t305;
t299 = t245 * t248;
t300 = t244 * t249;
t252 = t240 * t300 + t299;
t141 = t252 * t239 + t241 * t306;
t63 = pkin(3) * t140 - pkin(10) * t141 + t100;
t304 = t239 * t249;
t189 = -t238 * t304 + t241 * t240;
t232 = t245 * t322;
t195 = pkin(8) * t304 + t232;
t137 = (t238 * t241 + t240 * t304) * pkin(9) + t195;
t74 = t248 * t137 + t147 * t303 + t165 * t306;
t67 = pkin(10) * t189 + t74;
t319 = t243 * t63 + t247 * t67;
t354 = qJD(4) * t319;
t73 = -t244 * t137 + (t165 * t238 + t307) * t248;
t296 = qJD(3) * t238;
t184 = (pkin(3) * t244 - pkin(10) * t248) * t296;
t185 = t192 * qJD(3);
t80 = t184 * t247 - t185 * t243 - t358;
t297 = qJD(2) * t239;
t279 = t245 * t297;
t268 = t238 * t279;
t149 = (t249 * t267 - t232) * qJD(2);
t168 = (pkin(2) * t245 - t249 * t321) * t297;
t294 = qJD(3) * t248;
t276 = t238 * t294;
t295 = qJD(3) * t244;
t37 = -t137 * t295 + t149 * t303 + t165 * t276 + t168 * t306 + t359 * t248;
t35 = pkin(10) * t268 + t37;
t277 = t238 * t295;
t104 = t241 * t277 + (t252 * qJD(3) + (t240 * t299 + t300) * qJD(2)) * t239;
t105 = t241 * t276 + (t355 * qJD(3) + (-t240 * t301 + t298) * qJD(2)) * t239;
t106 = -t149 * t238 + t240 * t168;
t46 = pkin(3) * t104 - pkin(10) * t105 + t106;
t9 = -t243 * t35 + t247 * t46 - t354;
t353 = 2 * m(4);
t352 = 0.2e1 * m(5);
t351 = 2 * m(6);
t350 = 0.2e1 * pkin(10);
t349 = -2 * mrSges(3,3);
t348 = -2 * mrSges(4,3);
t108 = t141 * t247 + t189 * t243;
t59 = t108 * qJD(4) + t105 * t243 - t247 * t268;
t107 = t141 * t243 - t189 * t247;
t60 = -t107 * qJD(4) + t105 * t247 + t243 * t268;
t18 = Ifges(5,1) * t60 - Ifges(5,4) * t59 + Ifges(5,5) * t104;
t346 = t18 / 0.2e1;
t52 = Ifges(5,1) * t108 - Ifges(5,4) * t107 + Ifges(5,5) * t140;
t345 = t52 / 0.2e1;
t190 = -t247 * t240 + t243 * t306;
t150 = -t190 * qJD(4) + t247 * t276;
t191 = t240 * t243 + t247 * t306;
t151 = t191 * qJD(4) + t243 * t276;
t97 = Ifges(5,1) * t150 - Ifges(5,4) * t151 + Ifges(5,5) * t277;
t344 = t97 / 0.2e1;
t122 = Ifges(5,1) * t191 - Ifges(5,4) * t190 - Ifges(5,5) * t305;
t343 = t122 / 0.2e1;
t312 = Ifges(6,4) * t242;
t218 = Ifges(6,2) * t246 + t312;
t311 = Ifges(6,4) * t246;
t263 = -Ifges(6,2) * t242 + t311;
t126 = -t218 * t290 + (Ifges(6,6) * t243 + t263 * t247) * qJD(4);
t342 = t126 / 0.2e1;
t220 = Ifges(6,1) * t242 + t311;
t264 = Ifges(6,1) * t246 - t312;
t127 = -t220 * t290 + (Ifges(6,5) * t243 + t264 * t247) * qJD(4);
t341 = t127 / 0.2e1;
t152 = -t191 * t242 - t246 * t305;
t340 = t152 / 0.2e1;
t253 = -t191 * t246 + t242 * t305;
t339 = -t253 / 0.2e1;
t178 = -Ifges(6,6) * t247 + t263 * t243;
t338 = t178 / 0.2e1;
t179 = -Ifges(6,5) * t247 + t264 * t243;
t337 = t179 / 0.2e1;
t289 = qJD(5) * t246;
t236 = Ifges(6,5) * t289;
t291 = qJD(5) * t242;
t202 = -Ifges(6,6) * t291 + t236;
t336 = t202 / 0.2e1;
t204 = t263 * qJD(5);
t335 = t204 / 0.2e1;
t206 = t264 * qJD(5);
t334 = t206 / 0.2e1;
t314 = Ifges(5,4) * t243;
t207 = (Ifges(5,1) * t247 - t314) * qJD(4);
t333 = t207 / 0.2e1;
t332 = Ifges(6,5) * t325 + Ifges(6,6) * t323;
t331 = Ifges(5,5) * t243 / 0.2e1 + Ifges(5,6) * t247 / 0.2e1;
t330 = t218 / 0.2e1;
t329 = t220 / 0.2e1;
t313 = Ifges(5,4) * t247;
t221 = Ifges(5,1) * t243 + t313;
t328 = t221 / 0.2e1;
t327 = t240 / 0.2e1;
t326 = -t242 / 0.2e1;
t324 = -t246 / 0.2e1;
t320 = pkin(10) * t247;
t318 = m(6) * qJD(4);
t317 = mrSges(6,3) * t243;
t316 = Ifges(4,4) * t244;
t315 = Ifges(4,4) * t248;
t310 = Ifges(6,6) * t242;
t187 = -pkin(8) * t279 + t226;
t309 = t187 * mrSges(3,2);
t188 = t195 * qJD(2);
t308 = t188 * mrSges(3,1);
t293 = qJD(4) * t243;
t288 = qJD(5) * t247;
t78 = t108 * t246 + t140 * t242;
t23 = -t78 * qJD(5) + t104 * t246 - t242 * t60;
t77 = -t108 * t242 + t140 * t246;
t24 = t77 * qJD(5) + t104 * t242 + t246 * t60;
t3 = Ifges(6,5) * t24 + Ifges(6,6) * t23 + Ifges(6,3) * t59;
t17 = Ifges(5,4) * t60 - Ifges(5,2) * t59 + Ifges(5,6) * t104;
t286 = t3 / 0.2e1 - t17 / 0.2e1;
t16 = Ifges(5,5) * t60 - Ifges(5,6) * t59 + Ifges(5,3) * t104;
t91 = t152 * qJD(5) + t150 * t246 + t242 * t277;
t92 = t253 * qJD(5) - t150 * t242 + t246 * t277;
t42 = Ifges(6,5) * t91 + Ifges(6,6) * t92 + Ifges(6,3) * t151;
t30 = Ifges(6,5) * t78 + Ifges(6,6) * t77 + Ifges(6,3) * t107;
t51 = Ifges(5,4) * t108 - Ifges(5,2) * t107 + Ifges(5,6) * t140;
t285 = t30 / 0.2e1 - t51 / 0.2e1;
t96 = Ifges(5,4) * t150 - Ifges(5,2) * t151 + Ifges(5,6) * t277;
t284 = t42 / 0.2e1 - t96 / 0.2e1;
t10 = -mrSges(6,1) * t23 + mrSges(6,2) * t24;
t7 = -pkin(4) * t104 - t9;
t283 = -m(6) * t7 - t10;
t121 = Ifges(5,4) * t191 - Ifges(5,2) * t190 - Ifges(5,6) * t305;
t85 = -Ifges(6,5) * t253 + Ifges(6,6) * t152 + Ifges(6,3) * t190;
t281 = -t121 / 0.2e1 + t85 / 0.2e1;
t56 = Ifges(4,5) * t105 - Ifges(4,6) * t104 + Ifges(4,3) * t268;
t95 = Ifges(5,5) * t150 - Ifges(5,6) * t151 + Ifges(5,3) * t277;
t47 = -mrSges(6,1) * t92 + mrSges(6,2) * t91;
t76 = -pkin(4) * t277 - t80;
t280 = -m(6) * t76 - t47;
t250 = t242 * t292 + t243 * t289;
t125 = t251 * Ifges(6,5) - Ifges(6,6) * t250 + Ifges(6,3) * t293;
t205 = (-Ifges(5,2) * t243 + t313) * qJD(4);
t273 = -t205 / 0.2e1 + t125 / 0.2e1;
t177 = -Ifges(6,3) * t247 + (Ifges(6,5) * t246 - t310) * t243;
t219 = Ifges(5,2) * t247 + t314;
t272 = -t219 / 0.2e1 + t177 / 0.2e1;
t210 = (pkin(4) * t243 - pkin(11) * t247) * qJD(4);
t213 = -pkin(4) * t247 - pkin(11) * t243 - pkin(3);
t116 = t213 * t289 + t210 * t242 + (-t242 * t288 - t246 * t293) * pkin(10);
t169 = t213 * t246 - t242 * t320;
t271 = -qJD(5) * t169 + t116;
t117 = -t213 * t291 + t210 * t246 + (t242 * t293 - t246 * t288) * pkin(10);
t170 = t213 * t242 + t246 * t320;
t270 = -qJD(5) * t170 - t117;
t269 = -t137 * t294 - t165 * t277 - t359 * t244;
t26 = pkin(11) * t140 + t319;
t66 = -pkin(3) * t189 - t73;
t33 = pkin(4) * t107 - pkin(11) * t108 + t66;
t11 = -t242 * t26 + t246 * t33;
t36 = -t149 * t302 + (-pkin(3) * t279 - t168 * t248) * t238 - t269;
t13 = pkin(4) * t59 - pkin(11) * t60 + t36;
t8 = t243 * t46 + t247 * t35 + t63 * t292 - t67 * t293;
t6 = pkin(11) * t104 + t8;
t1 = t11 * qJD(5) + t13 * t242 + t246 * t6;
t12 = t242 * t33 + t246 * t26;
t2 = -t12 * qJD(5) + t13 * t246 - t242 * t6;
t266 = t1 * t246 - t2 * t242;
t214 = -mrSges(6,1) * t246 + mrSges(6,2) * t242;
t265 = mrSges(6,1) * t242 + mrSges(6,2) * t246;
t174 = t228 + (-pkin(2) * t248 - pkin(3)) * t240;
t109 = pkin(4) * t190 - pkin(11) * t191 + t174;
t111 = -pkin(11) * t305 + t357;
t68 = t109 * t246 - t111 * t242;
t79 = -t175 * t293 + t176 * t292 + t243 * t184 + t247 * t185;
t75 = pkin(11) * t277 + t79;
t186 = t194 * qJD(3);
t88 = pkin(4) * t151 - pkin(11) * t150 + t186;
t19 = t68 * qJD(5) + t242 * t88 + t246 * t75;
t69 = t109 * t242 + t111 * t246;
t20 = -t69 * qJD(5) - t242 * t75 + t246 * t88;
t262 = t19 * t246 - t20 * t242;
t28 = -t243 * t67 + t247 * t63;
t118 = -t175 * t243 + t176 * t247;
t31 = Ifges(6,4) * t78 + Ifges(6,2) * t77 + Ifges(6,6) * t107;
t32 = Ifges(6,1) * t78 + Ifges(6,4) * t77 + Ifges(6,5) * t107;
t255 = t31 * t326 + t32 * t323;
t86 = -Ifges(6,4) * t253 + Ifges(6,2) * t152 + Ifges(6,6) * t190;
t87 = -Ifges(6,1) * t253 + Ifges(6,4) * t152 + Ifges(6,5) * t190;
t254 = t87 * t323 + t86 * t326;
t237 = Ifges(5,5) * t292;
t225 = Ifges(3,5) * t249 * t297;
t224 = Ifges(4,5) * t276;
t209 = -mrSges(6,1) * t247 - t246 * t317;
t208 = mrSges(6,2) * t247 - t242 * t317;
t203 = -Ifges(5,6) * t293 + t237;
t201 = (mrSges(5,1) * t243 + mrSges(5,2) * t247) * qJD(4);
t200 = t265 * qJD(5);
t199 = -mrSges(4,2) * t240 + mrSges(4,3) * t305;
t198 = mrSges(4,1) * t240 - mrSges(4,3) * t306;
t196 = t265 * t243;
t193 = -pkin(8) * t239 * t245 + t233;
t183 = (Ifges(4,1) * t248 - t316) * t296;
t182 = (-Ifges(4,2) * t244 + t315) * t296;
t181 = -Ifges(4,6) * t277 + t224;
t180 = (mrSges(4,1) * t244 + mrSges(4,2) * t248) * t296;
t172 = Ifges(4,5) * t240 + (Ifges(4,1) * t244 + t315) * t238;
t171 = Ifges(4,6) * t240 + (Ifges(4,2) * t248 + t316) * t238;
t162 = -mrSges(6,2) * t293 - mrSges(6,3) * t250;
t161 = mrSges(6,1) * t293 - mrSges(6,3) * t251;
t158 = -mrSges(5,1) * t305 - mrSges(5,3) * t191;
t157 = mrSges(5,2) * t305 - mrSges(5,3) * t190;
t136 = t250 * mrSges(6,1) + mrSges(6,2) * t251;
t131 = mrSges(5,1) * t190 + mrSges(5,2) * t191;
t124 = -mrSges(5,2) * t277 - mrSges(5,3) * t151;
t123 = mrSges(5,1) * t277 - mrSges(5,3) * t150;
t120 = Ifges(5,5) * t191 - Ifges(5,6) * t190 - Ifges(5,3) * t305;
t115 = mrSges(6,1) * t190 + mrSges(6,3) * t253;
t114 = -mrSges(6,2) * t190 + mrSges(6,3) * t152;
t113 = mrSges(4,1) * t189 - mrSges(4,3) * t141;
t112 = -mrSges(4,2) * t189 - mrSges(4,3) * t140;
t110 = pkin(4) * t305 - t118;
t99 = -mrSges(6,1) * t152 - mrSges(6,2) * t253;
t98 = mrSges(5,1) * t151 + mrSges(5,2) * t150;
t94 = mrSges(4,1) * t268 - mrSges(4,3) * t105;
t93 = -mrSges(4,2) * t268 - mrSges(4,3) * t104;
t84 = Ifges(4,1) * t141 - Ifges(4,4) * t140 + Ifges(4,5) * t189;
t83 = Ifges(4,4) * t141 - Ifges(4,2) * t140 + Ifges(4,6) * t189;
t82 = mrSges(5,1) * t140 - mrSges(5,3) * t108;
t81 = -mrSges(5,2) * t140 - mrSges(5,3) * t107;
t72 = -mrSges(6,2) * t151 + mrSges(6,3) * t92;
t71 = mrSges(6,1) * t151 - mrSges(6,3) * t91;
t70 = mrSges(5,1) * t107 + mrSges(5,2) * t108;
t65 = mrSges(4,1) * t104 + mrSges(4,2) * t105;
t58 = Ifges(4,1) * t105 - Ifges(4,4) * t104 + Ifges(4,5) * t268;
t57 = Ifges(4,4) * t105 - Ifges(4,2) * t104 + Ifges(4,6) * t268;
t50 = Ifges(5,5) * t108 - Ifges(5,6) * t107 + Ifges(5,3) * t140;
t49 = mrSges(6,1) * t107 - mrSges(6,3) * t78;
t48 = -mrSges(6,2) * t107 + mrSges(6,3) * t77;
t44 = Ifges(6,1) * t91 + Ifges(6,4) * t92 + Ifges(6,5) * t151;
t43 = Ifges(6,4) * t91 + Ifges(6,2) * t92 + Ifges(6,6) * t151;
t41 = -mrSges(6,1) * t77 + mrSges(6,2) * t78;
t40 = mrSges(5,1) * t104 - mrSges(5,3) * t60;
t39 = -mrSges(5,2) * t104 - mrSges(5,3) * t59;
t38 = (t149 * t240 + t168 * t238) * t248 + t269;
t27 = mrSges(5,1) * t59 + mrSges(5,2) * t60;
t25 = -pkin(4) * t140 - t28;
t15 = mrSges(6,1) * t59 - mrSges(6,3) * t24;
t14 = -mrSges(6,2) * t59 + mrSges(6,3) * t23;
t5 = Ifges(6,1) * t24 + Ifges(6,4) * t23 + Ifges(6,5) * t59;
t4 = Ifges(6,4) * t24 + Ifges(6,2) * t23 + Ifges(6,6) * t59;
t21 = [(0.2e1 * (t187 * t249 + t188 * t245) * mrSges(3,3) + ((t193 * t349 + Ifges(3,5) * t241 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t249) * t360) * t249 + (-0.2e1 * Ifges(3,6) * t241 + t238 * (Ifges(4,5) * t141 - Ifges(4,6) * t140 + Ifges(4,3) * t189) + t195 * t349 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t245 + (Ifges(3,1) - Ifges(3,2)) * t249) * t360) * t245) * qJD(2)) * t239 + 0.2e1 * t319 * t39 + (t28 * t9 + t319 * t8 + t36 * t66) * t352 + 0.2e1 * m(3) * (t187 * t195 - t188 * t193) + (t225 - 0.2e1 * t308 - 0.2e1 * t309) * t241 + (t30 - t51) * t59 + t189 * t56 + t140 * t16 - t140 * t57 + 0.2e1 * t106 * (t140 * mrSges(4,1) + t141 * mrSges(4,2)) + t141 * t58 + 0.2e1 * t37 * t112 + 0.2e1 * t38 * t113 + t108 * t18 + t105 * t84 + 0.2e1 * t100 * t65 + 0.2e1 * t74 * t93 + 0.2e1 * t73 * t94 + t78 * t5 + 0.2e1 * t8 * t81 + 0.2e1 * t9 * t82 + t77 * t4 + 0.2e1 * t36 * t70 + t60 * t52 + 0.2e1 * t66 * t27 + 0.2e1 * t1 * t48 + 0.2e1 * t2 * t49 + 0.2e1 * t28 * t40 + 0.2e1 * t7 * t41 + t24 * t32 + t23 * t31 + 0.2e1 * t25 * t10 + 0.2e1 * t12 * t14 + 0.2e1 * t11 * t15 + (t3 - t17) * t107 + (t50 - t83) * t104 + (t100 * t106 + t37 * t74 + t38 * t73) * t353 + (t1 * t12 + t11 * t2 + t25 * t7) * t351; m(5) * (t118 * t9 + t174 * t36 + t186 * t66 + t28 * t80 + t319 * t79 + t357 * t8) + t357 * t39 - t309 - t308 + t319 * t124 + t225 + t194 * t93 + t38 * t198 + t37 * t199 + t189 * t181 / 0.2e1 + t192 * t94 + t185 * t112 + t105 * t172 / 0.2e1 + t174 * t27 + t100 * t180 + t141 * t183 / 0.2e1 + t8 * t157 + t9 * t158 + t28 * t123 + t36 * t131 + t1 * t114 + t2 * t115 + t118 * t40 + t110 * t10 + t66 * t98 + t7 * t99 + t92 * t31 / 0.2e1 + t23 * t86 / 0.2e1 + t24 * t87 / 0.2e1 + t91 * t32 / 0.2e1 + t78 * t44 / 0.2e1 + t79 * t81 + t80 * t82 + t76 * t41 + t77 * t43 / 0.2e1 + t68 * t15 + t69 * t14 + t11 * t71 + t12 * t72 + t25 * t47 + t19 * t48 + t20 * t49 + t5 * t339 + t4 * t340 + t60 * t343 + t108 * t344 + t150 * t345 + t191 * t346 + t56 * t327 + (-t171 / 0.2e1 + t120 / 0.2e1) * t104 + (t70 - t113) * t186 + m(4) * (t185 * t74 - t186 * t73 + t192 * t38 + t194 * t37) + (-t182 / 0.2e1 + t95 / 0.2e1) * t140 - Ifges(3,6) * t279 + t281 * t59 + m(6) * (t1 * t69 + t11 * t20 + t110 * t7 + t12 * t19 + t2 * t68 + t25 * t76) + t284 * t107 + t285 * t151 + t286 * t190 + ((t58 / 0.2e1 + t106 * mrSges(4,2)) * t244 + (t57 / 0.2e1 - t16 / 0.2e1 - t106 * mrSges(4,1)) * t248 + (Ifges(4,3) * t327 + (Ifges(4,5) * t244 + Ifges(4,6) * t248) * t238 / 0.2e1) * t279 + (-m(4) * t106 - t65) * pkin(2) + ((-t73 * mrSges(4,3) + t84 / 0.2e1) * t248 + (-t74 * mrSges(4,3) + t50 / 0.2e1 - t83 / 0.2e1) * t244) * qJD(3)) * t238; 0.2e1 * t357 * t124 + (t118 * t80 + t174 * t186 + t357 * t79) * t352 - t253 * t44 + 0.2e1 * (-t198 + t131) * t186 + t240 * t181 + 0.2e1 * t185 * t199 + t191 * t97 + 0.2e1 * t174 * t98 + 0.2e1 * t79 * t157 + 0.2e1 * t80 * t158 + t150 * t122 + t152 * t43 + 0.2e1 * t118 * t123 + 0.2e1 * t19 * t114 + 0.2e1 * t20 * t115 + 0.2e1 * t110 * t47 + 0.2e1 * t76 * t99 + t92 * t86 + t91 * t87 + 0.2e1 * t68 * t71 + 0.2e1 * t69 * t72 + (t185 * t194 - t186 * t192) * t353 + (t110 * t76 + t19 * t69 + t20 * t68) * t351 + (t42 - t96) * t190 + (t85 - t121) * t151 + (-0.2e1 * pkin(2) * t180 + t244 * t183 + (t182 - t95) * t248 + ((t192 * t348 + t172) * t248 + (t194 * t348 + t120 - t171) * t244) * qJD(3)) * t238; (t5 * t323 + t4 * t326 - t9 * mrSges(5,3) + t346 + (t31 * t324 + t32 * t326) * qJD(5) + (-mrSges(5,3) * t319 + t285) * qJD(4) + (-qJD(4) * t81 - t40 + m(5) * (-t9 - t354) - t283) * pkin(10)) * t243 + t356 * t36 + m(6) * (t1 * t170 + t11 * t117 + t116 * t12 + t169 * t2) + t7 * t196 + t66 * t201 + t140 * t203 / 0.2e1 + t1 * t208 + t2 * t209 + t170 * t14 + t11 * t161 + t12 * t162 + t169 * t15 + t25 * t136 + t116 * t48 + t117 * t49 + t38 * mrSges(4,1) - t37 * mrSges(4,2) - pkin(3) * t27 + t24 * t337 + t23 * t338 + t78 * t341 + t77 * t342 + t104 * t331 + t108 * t333 + t60 * t328 + t56 + t272 * t59 + t273 * t107 + (t8 * mrSges(5,3) + (-t28 * mrSges(5,3) + t255 + t345) * qJD(4) + (t39 + (t41 - t82) * qJD(4) + t25 * t318 + m(5) * (-qJD(4) * t28 + t8)) * pkin(10) - t286) * t247; (-mrSges(4,1) + t356) * t186 + (t44 * t323 + t43 * t326 - t80 * mrSges(5,3) + t344 + (t86 * t324 + t87 * t326) * qJD(5) + (-mrSges(5,3) * t357 + t281) * qJD(4) + (-qJD(4) * t157 - t123 + m(5) * (-t80 - t358) - t280) * pkin(10)) * t243 + t224 + t76 * t196 + t174 * t201 + t19 * t208 + t20 * t209 - t185 * mrSges(4,2) + t170 * t72 + t68 * t161 + t69 * t162 + t169 * t71 + t110 * t136 + t116 * t114 + t117 * t115 - pkin(3) * t98 + t91 * t337 + t92 * t338 + t127 * t339 + t126 * t340 + t191 * t333 + t150 * t328 + m(6) * (t116 * t69 + t117 * t68 + t169 * t20 + t170 * t19) + t272 * t151 + t273 * t190 + (-t248 * t203 / 0.2e1 + (t331 - Ifges(4,6)) * t295) * t238 + (t79 * mrSges(5,3) + (-t118 * mrSges(5,3) + t254 + t343) * qJD(4) + (t124 + (-t158 + t99) * qJD(4) + m(5) * (-qJD(4) * t118 + t79) + t110 * t318) * pkin(10) - t284) * t247; 0.2e1 * t117 * t209 + 0.2e1 * t169 * t161 + (t116 * t170 + t117 * t169) * t351 + 0.2e1 * t116 * t208 + 0.2e1 * t170 * t162 - 0.2e1 * pkin(3) * t201 + (-t125 + t205 + (-t178 * t242 + t179 * t246 + t196 * t350 + t221) * qJD(4)) * t247 + (t136 * t350 - t242 * t126 + t246 * t127 + t207 + (-t178 * t246 - t179 * t242) * qJD(5) + (pkin(10) ^ 2 * t247 * t351 + t177 - t219) * qJD(4)) * t243; t4 * t323 + t5 * t325 + t7 * t214 + t59 * t332 + t23 * t330 + t24 * t329 + t25 * t200 + t107 * t336 + t77 * t335 + t78 * t334 + t9 * mrSges(5,1) - t8 * mrSges(5,2) + t255 * qJD(5) + t283 * pkin(4) + ((-t11 * t246 - t12 * t242) * qJD(5) + t266) * mrSges(6,3) + (m(6) * (-t11 * t289 - t12 * t291 + t266) + t246 * t14 - t242 * t15 - t49 * t289 - t48 * t291) * pkin(11) + t16; t43 * t323 + t44 * t325 + t76 * t214 + t151 * t332 + t92 * t330 + t91 * t329 + t110 * t200 + t190 * t336 + t152 * t335 - t253 * t334 - t79 * mrSges(5,2) + t80 * mrSges(5,1) + t254 * qJD(5) + t280 * pkin(4) + ((-t242 * t69 - t246 * t68) * qJD(5) + t262) * mrSges(6,3) + (m(6) * (-t68 * t289 - t69 * t291 + t262) + t246 * t72 - t242 * t71 - t115 * t289 - t114 * t291) * pkin(11) + t95; -pkin(4) * t136 + t237 + (-t202 / 0.2e1 + (-m(6) * pkin(4) - mrSges(5,1) + t214) * qJD(4) * pkin(10)) * t247 + (qJD(5) * t337 + t342 + t292 * t329 + t271 * mrSges(6,3) + (m(6) * t271 - qJD(5) * t209 + t162) * pkin(11)) * t246 + (-qJD(5) * t178 / 0.2e1 + t341 - t218 * t292 / 0.2e1 + t270 * mrSges(6,3) + (m(6) * t270 - qJD(5) * t208 - t161) * pkin(11)) * t242 + (t206 * t323 + t204 * t326 + pkin(10) * t200 + (t218 * t324 + t220 * t326) * qJD(5) + (pkin(10) * mrSges(5,2) - Ifges(5,6) + t332) * qJD(4)) * t243; -0.2e1 * pkin(4) * t200 + t204 * t246 + t206 * t242 + (-t218 * t242 + t220 * t246) * qJD(5); mrSges(6,1) * t2 - mrSges(6,2) * t1 + t3; mrSges(6,1) * t20 - mrSges(6,2) * t19 + t42; mrSges(6,1) * t117 - mrSges(6,2) * t116 + t125; t236 + (t214 * pkin(11) - t310) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t21(1), t21(2), t21(4), t21(7), t21(11); t21(2), t21(3), t21(5), t21(8), t21(12); t21(4), t21(5), t21(6), t21(9), t21(13); t21(7), t21(8), t21(9), t21(10), t21(14); t21(11), t21(12), t21(13), t21(14), t21(15);];
Mq = res;
