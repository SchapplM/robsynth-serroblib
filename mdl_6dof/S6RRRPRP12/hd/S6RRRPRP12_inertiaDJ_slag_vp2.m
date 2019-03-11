% Calculate time derivative of joint inertia matrix for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:45
% EndTime: 2019-03-09 17:52:02
% DurationCPUTime: 7.43s
% Computational Cost: add. (5616->664), mult. (14738->905), div. (0->0), fcn. (13057->8), ass. (0->265)
t218 = cos(qJ(5));
t326 = Ifges(6,5) + Ifges(7,4);
t329 = t326 * t218;
t328 = Ifges(5,1) + Ifges(4,3);
t327 = Ifges(4,5) - Ifges(5,4);
t325 = -Ifges(4,6) + Ifges(5,5);
t215 = sin(qJ(5));
t216 = sin(qJ(3));
t278 = qJD(3) * t216;
t257 = t218 * t278;
t219 = cos(qJ(3));
t274 = qJD(5) * t219;
t227 = t215 * t274 + t257;
t294 = qJ(4) * t216;
t315 = pkin(3) + pkin(10);
t148 = -t219 * t315 - pkin(2) - t294;
t314 = pkin(4) + pkin(9);
t188 = t314 * t216;
t320 = t218 * t148 + t215 * t188;
t324 = qJD(5) * t320;
t323 = 0.2e1 * qJ(4);
t214 = cos(pkin(6));
t220 = cos(qJ(2));
t213 = sin(pkin(6));
t280 = qJD(2) * t213;
t261 = t220 * t280;
t217 = sin(qJ(2));
t293 = t213 * t217;
t264 = t216 * t293;
t103 = -qJD(3) * t264 + (qJD(3) * t214 + t261) * t219;
t142 = t214 * t216 + t219 * t293;
t102 = qJD(3) * t142 + t216 * t261;
t141 = -t214 * t219 + t264;
t279 = qJD(2) * t217;
t262 = t213 * t279;
t292 = t213 * t220;
t263 = t218 * t292;
t52 = -qJD(5) * t263 - t102 * t218 + (qJD(5) * t141 + t262) * t215;
t104 = t141 * t218 + t215 * t292;
t53 = qJD(5) * t104 + t102 * t215 + t218 * t262;
t7 = Ifges(6,5) * t53 - Ifges(6,6) * t52 + Ifges(6,3) * t103;
t8 = Ifges(7,4) * t53 + Ifges(7,2) * t103 + Ifges(7,6) * t52;
t322 = t7 + t8;
t200 = pkin(8) * t293;
t313 = pkin(1) * t220;
t144 = t214 * t313 - t200;
t247 = pkin(3) * t278 - qJD(4) * t216;
t122 = (pkin(10) * t216 - qJ(4) * t219) * qJD(3) + t247;
t277 = qJD(3) * t219;
t170 = t314 * t277;
t51 = -t122 * t215 + t170 * t218 - t324;
t145 = t214 * t217 * pkin(1) + pkin(8) * t292;
t127 = pkin(9) * t214 + t145;
t128 = (-pkin(2) * t220 - pkin(9) * t217 - pkin(1)) * t213;
t136 = (pkin(2) * t217 - pkin(9) * t220) * t280;
t137 = t144 * qJD(2);
t32 = -t127 * t277 - t128 * t278 + t136 * t219 - t216 * t137;
t19 = pkin(4) * t103 - t262 * t315 - t32;
t138 = t145 * qJD(2);
t223 = -qJ(4) * t103 - qJD(4) * t142 + t138;
t23 = t102 * t315 + t223;
t72 = -t216 * t127 + t128 * t219;
t66 = pkin(3) * t292 - t72;
t43 = pkin(4) * t142 + pkin(10) * t292 + t66;
t126 = t200 + (-pkin(2) - t313) * t214;
t224 = -qJ(4) * t142 + t126;
t54 = t141 * t315 + t224;
t308 = t215 * t43 + t218 * t54;
t4 = -qJD(5) * t308 + t19 * t218 - t215 * t23;
t319 = 2 * m(6);
t318 = 2 * m(7);
t317 = -2 * mrSges(3,3);
t316 = 2 * mrSges(3,3);
t312 = pkin(9) * t216;
t211 = t219 * pkin(9);
t25 = -mrSges(6,2) * t103 - mrSges(6,3) * t52;
t28 = -mrSges(7,2) * t52 + mrSges(7,3) * t103;
t310 = t25 + t28;
t26 = mrSges(6,1) * t103 - mrSges(6,3) * t53;
t27 = -t103 * mrSges(7,1) + t53 * mrSges(7,2);
t309 = t26 - t27;
t307 = Ifges(4,4) * t216;
t306 = Ifges(4,4) * t219;
t305 = Ifges(6,4) * t215;
t304 = Ifges(6,4) * t218;
t303 = Ifges(7,5) * t215;
t302 = Ifges(7,5) * t218;
t301 = Ifges(5,6) * t216;
t300 = Ifges(5,6) * t219;
t299 = Ifges(7,6) * t215;
t298 = t137 * mrSges(3,2);
t297 = t138 * mrSges(3,1);
t296 = t138 * mrSges(4,1);
t295 = t138 * mrSges(4,2);
t291 = t215 * t219;
t290 = t215 * t315;
t289 = t218 * t219;
t288 = t218 * t315;
t114 = -mrSges(6,2) * t277 + mrSges(6,3) * t227;
t117 = mrSges(7,2) * t227 + mrSges(7,3) * t277;
t287 = t114 + t117;
t258 = t215 * t278;
t259 = t218 * t274;
t226 = t258 - t259;
t115 = mrSges(6,1) * t277 - mrSges(6,3) * t226;
t275 = qJD(5) * t218;
t116 = mrSges(7,2) * t258 + (-mrSges(7,1) * qJD(3) - mrSges(7,2) * t275) * t219;
t286 = t115 - t116;
t73 = t219 * t127 + t216 * t128;
t238 = -Ifges(7,3) * t218 + t303;
t129 = Ifges(7,6) * t216 - t219 * t238;
t240 = Ifges(6,2) * t218 + t305;
t132 = Ifges(6,6) * t216 - t219 * t240;
t285 = -t129 + t132;
t241 = Ifges(7,1) * t215 - t302;
t133 = Ifges(7,4) * t216 - t219 * t241;
t242 = Ifges(6,1) * t215 + t304;
t134 = Ifges(6,5) * t216 - t219 * t242;
t284 = t133 + t134;
t165 = mrSges(6,1) * t216 + mrSges(6,3) * t291;
t166 = -mrSges(7,1) * t216 - mrSges(7,2) * t291;
t283 = -t165 + t166;
t167 = -mrSges(6,2) * t216 - mrSges(6,3) * t289;
t168 = -mrSges(7,2) * t289 + mrSges(7,3) * t216;
t282 = t167 + t168;
t281 = Ifges(7,4) * t258 + Ifges(7,2) * t277;
t189 = t219 * pkin(4) + t211;
t276 = qJD(5) * t215;
t273 = qJD(6) * t215;
t272 = 0.2e1 * t213;
t6 = Ifges(7,5) * t53 + Ifges(7,6) * t103 + Ifges(7,3) * t52;
t9 = Ifges(6,4) * t53 - Ifges(6,2) * t52 + Ifges(6,6) * t103;
t271 = t6 / 0.2e1 - t9 / 0.2e1;
t10 = Ifges(7,1) * t53 + Ifges(7,4) * t103 + Ifges(7,5) * t52;
t11 = Ifges(6,1) * t53 - Ifges(6,4) * t52 + Ifges(6,5) * t103;
t269 = -t11 / 0.2e1 - t10 / 0.2e1;
t105 = t141 * t215 - t263;
t37 = Ifges(7,5) * t105 + Ifges(7,6) * t142 - Ifges(7,3) * t104;
t40 = Ifges(6,4) * t105 + Ifges(6,2) * t104 + Ifges(6,6) * t142;
t268 = -t37 / 0.2e1 + t40 / 0.2e1;
t41 = Ifges(7,1) * t105 + Ifges(7,4) * t142 - Ifges(7,5) * t104;
t42 = Ifges(6,1) * t105 + Ifges(6,4) * t104 + Ifges(6,5) * t142;
t267 = t41 / 0.2e1 + t42 / 0.2e1;
t180 = Ifges(7,3) * t215 + t302;
t82 = -t180 * t274 + (Ifges(7,6) * t219 + t216 * t238) * qJD(3);
t183 = -Ifges(6,2) * t215 + t304;
t85 = -t183 * t274 + (Ifges(6,6) * t219 + t216 * t240) * qJD(3);
t266 = t85 / 0.2e1 - t82 / 0.2e1;
t185 = Ifges(7,1) * t218 + t303;
t86 = -t185 * t274 + (Ifges(7,4) * t219 + t216 * t241) * qJD(3);
t186 = Ifges(6,1) * t218 - t305;
t87 = -t186 * t274 + (Ifges(6,5) * t219 + t216 * t242) * qJD(3);
t265 = t86 / 0.2e1 + t87 / 0.2e1;
t256 = t292 / 0.2e1;
t255 = t129 / 0.2e1 - t132 / 0.2e1;
t254 = t133 / 0.2e1 + t134 / 0.2e1;
t157 = t238 * qJD(5);
t160 = t240 * qJD(5);
t253 = -t157 / 0.2e1 + t160 / 0.2e1;
t206 = Ifges(7,6) * t275;
t239 = -Ifges(6,5) * t215 - Ifges(6,6) * t218;
t252 = t239 * qJD(5) / 0.2e1 - Ifges(7,4) * t276 / 0.2e1 + t206 / 0.2e1;
t162 = t241 * qJD(5);
t163 = t242 * qJD(5);
t251 = t162 / 0.2e1 + t163 / 0.2e1;
t250 = t180 / 0.2e1 - t183 / 0.2e1;
t249 = -Ifges(6,6) * t215 / 0.2e1 + t299 / 0.2e1 + t329 / 0.2e1;
t248 = t185 / 0.2e1 + t186 / 0.2e1;
t79 = t103 * mrSges(5,1) + mrSges(5,2) * t262;
t246 = Ifges(6,5) * t258 + Ifges(6,6) * t227 + Ifges(6,3) * t277;
t245 = (-m(7) * t315 - mrSges(7,2)) * t215;
t244 = mrSges(6,1) * t218 - mrSges(6,2) * t215;
t177 = t215 * mrSges(6,1) + t218 * mrSges(6,2);
t243 = mrSges(7,1) * t218 + mrSges(7,3) * t215;
t176 = t215 * mrSges(7,1) - t218 * mrSges(7,3);
t175 = t219 * mrSges(5,2) - t216 * mrSges(5,3);
t237 = -pkin(3) * t219 - t294;
t236 = pkin(5) * t218 + qJ(6) * t215;
t235 = -pkin(5) * t215 + qJ(6) * t218;
t12 = qJ(6) * t142 + t308;
t14 = -t215 * t54 + t218 * t43;
t13 = -pkin(5) * t142 - t14;
t234 = t12 * t218 + t13 * t215;
t233 = -t14 * t215 + t218 * t308;
t94 = -t148 * t215 + t188 * t218;
t65 = qJ(4) * t292 - t73;
t228 = t325 * t102 + t327 * t103 + t262 * t328;
t3 = t215 * t19 + t218 * t23 + t43 * t275 - t276 * t54;
t31 = -t127 * t278 + t128 * t277 + t216 * t136 + t219 * t137;
t50 = t218 * t122 - t148 * t276 + t215 * t170 + t188 * t275;
t68 = mrSges(7,2) * t104 + mrSges(7,3) * t142;
t69 = -mrSges(6,2) * t142 + mrSges(6,3) * t104;
t70 = mrSges(6,1) * t142 - mrSges(6,3) * t105;
t71 = -mrSges(7,1) * t142 + mrSges(7,2) * t105;
t225 = (t68 + t69) * t218 + (-t70 + t71) * t215;
t55 = -pkin(4) * t141 - t65;
t24 = -qJ(4) * t262 + qJD(4) * t292 - t31;
t20 = -pkin(4) * t102 - t24;
t222 = m(7) * t235 - t176 - t177;
t209 = Ifges(4,5) * t277;
t208 = Ifges(5,5) * t278;
t193 = Ifges(3,5) * t261;
t187 = Ifges(4,1) * t216 + t306;
t184 = Ifges(4,2) * t219 + t307;
t179 = -Ifges(5,2) * t216 - t300;
t178 = -Ifges(5,3) * t219 - t301;
t172 = -pkin(2) + t237;
t171 = qJ(4) - t235;
t169 = t314 * t278;
t164 = (Ifges(4,1) * t219 - t307) * qJD(3);
t161 = (-Ifges(4,2) * t216 + t306) * qJD(3);
t156 = (-Ifges(5,2) * t219 + t301) * qJD(3);
t155 = (Ifges(5,3) * t216 - t300) * qJD(3);
t154 = (mrSges(4,1) * t216 + mrSges(4,2) * t219) * qJD(3);
t153 = (-mrSges(5,2) * t216 - mrSges(5,3) * t219) * qJD(3);
t152 = t244 * qJD(5);
t151 = t243 * qJD(5);
t147 = t244 * t219;
t146 = t243 * t219;
t140 = -qJ(4) * t277 + t247;
t135 = qJD(5) * t236 - qJD(6) * t218 + qJD(4);
t131 = Ifges(7,2) * t216 + (-Ifges(7,4) * t215 + Ifges(7,6) * t218) * t219;
t130 = Ifges(6,3) * t216 + t219 * t239;
t118 = t219 * t236 + t189;
t110 = -mrSges(4,1) * t292 - mrSges(4,3) * t142;
t109 = mrSges(4,2) * t292 - mrSges(4,3) * t141;
t108 = mrSges(5,1) * t142 - mrSges(5,2) * t292;
t107 = mrSges(5,1) * t141 + mrSges(5,3) * t292;
t92 = -mrSges(6,1) * t227 + mrSges(6,2) * t226;
t91 = -mrSges(7,1) * t227 - mrSges(7,3) * t226;
t90 = -pkin(5) * t216 - t94;
t89 = qJ(6) * t216 + t320;
t88 = -mrSges(5,2) * t141 - mrSges(5,3) * t142;
t84 = -Ifges(7,4) * t259 - Ifges(7,6) * t227 + t281;
t83 = -Ifges(6,5) * t259 + t246;
t81 = mrSges(4,1) * t262 - mrSges(4,3) * t103;
t80 = -mrSges(4,2) * t262 - mrSges(4,3) * t102;
t78 = mrSges(5,1) * t102 - mrSges(5,3) * t262;
t77 = Ifges(4,1) * t142 - Ifges(4,4) * t141 - Ifges(4,5) * t292;
t76 = Ifges(4,4) * t142 - Ifges(4,2) * t141 - Ifges(4,6) * t292;
t75 = -Ifges(5,4) * t292 - Ifges(5,2) * t142 + Ifges(5,6) * t141;
t74 = -Ifges(5,5) * t292 - Ifges(5,6) * t142 + Ifges(5,3) * t141;
t67 = (qJD(5) * t235 + t273) * t219 + (-t236 - t314) * t278;
t64 = pkin(3) * t141 + t224;
t63 = -mrSges(6,1) * t104 + mrSges(6,2) * t105;
t62 = -mrSges(7,1) * t104 - mrSges(7,3) * t105;
t61 = -mrSges(5,2) * t102 - mrSges(5,3) * t103;
t60 = mrSges(4,1) * t102 + mrSges(4,2) * t103;
t59 = Ifges(4,1) * t103 - Ifges(4,4) * t102 + Ifges(4,5) * t262;
t58 = Ifges(4,4) * t103 - Ifges(4,2) * t102 + Ifges(4,6) * t262;
t57 = Ifges(5,4) * t262 - Ifges(5,2) * t103 + Ifges(5,6) * t102;
t56 = Ifges(5,5) * t262 - Ifges(5,6) * t103 + Ifges(5,3) * t102;
t39 = Ifges(7,4) * t105 + Ifges(7,2) * t142 - Ifges(7,6) * t104;
t38 = Ifges(6,5) * t105 + Ifges(6,6) * t104 + Ifges(6,3) * t142;
t36 = -pkin(5) * t277 - t51;
t33 = qJ(6) * t277 + qJD(6) * t216 + t50;
t30 = pkin(3) * t102 + t223;
t29 = -pkin(3) * t262 - t32;
t21 = -pkin(5) * t104 - qJ(6) * t105 + t55;
t17 = mrSges(6,1) * t52 + mrSges(6,2) * t53;
t16 = mrSges(7,1) * t52 - mrSges(7,3) * t53;
t5 = pkin(5) * t52 - qJ(6) * t53 - qJD(6) * t105 + t20;
t2 = -pkin(5) * t103 - t4;
t1 = qJ(6) * t103 + qJD(6) * t142 + t3;
t15 = [(t74 - t76) * t102 + (t41 + t42) * t53 + (-t40 + t37) * t52 + (t14 * t4 + t20 * t55 + t3 * t308) * t319 + 0.2e1 * t308 * t25 + 0.2e1 * m(3) * (t137 * t145 - t138 * t144) + 0.2e1 * m(4) * (t126 * t138 + t31 * t73 + t32 * t72) + 0.2e1 * m(5) * (t24 * t65 + t29 * t66 + t30 * t64) + (t10 + t11) * t105 + (t9 - t6) * t104 + (t77 + t38 + t39 - t75) * t103 + 0.2e1 * t55 * t17 + 0.2e1 * t5 * t62 + 0.2e1 * t20 * t63 + 0.2e1 * t64 * t61 + 0.2e1 * t1 * t68 + 0.2e1 * t3 * t69 + 0.2e1 * t4 * t70 + 0.2e1 * t2 * t71 + 0.2e1 * t65 * t78 + 0.2e1 * t66 * t79 + 0.2e1 * t73 * t80 + 0.2e1 * t72 * t81 + 0.2e1 * t30 * t88 + 0.2e1 * t24 * t107 + 0.2e1 * t29 * t108 + 0.2e1 * t31 * t109 + 0.2e1 * t32 * t110 + 0.2e1 * t126 * t60 + (t1 * t12 + t13 * t2 + t21 * t5) * t318 + 0.2e1 * t21 * t16 + 0.2e1 * t14 * t26 + 0.2e1 * t13 * t27 + 0.2e1 * t12 * t28 + (t138 * t217 * t316 + (t137 * t316 - t228) * t220 + ((t144 * t317 + Ifges(3,5) * t214 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t220) * t272) * t220 + (t145 * t317 - 0.2e1 * Ifges(3,6) * t214 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t217) * t272 + t327 * t142 + t325 * t141 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - t328) * t292) * t217) * qJD(2)) * t213 + (-t57 + t59 + 0.2e1 * t295 + t322) * t142 + (t56 - t58 + 0.2e1 * t296) * t141 + (t193 - 0.2e1 * t297 - 0.2e1 * t298) * t214; -t297 + (t83 / 0.2e1 + t84 / 0.2e1 - t156 / 0.2e1 + t164 / 0.2e1) * t142 + (t155 / 0.2e1 - t161 / 0.2e1) * t141 - t298 + (t130 / 0.2e1 + t131 / 0.2e1 - t179 / 0.2e1 + t187 / 0.2e1) * t103 + (t178 / 0.2e1 - t184 / 0.2e1) * t102 + t308 * t114 + t193 + m(7) * (t1 * t89 + t118 * t5 + t12 * t33 + t13 * t36 + t2 * t90 + t21 * t67) - pkin(2) * t60 + t67 * t62 + t33 * t68 + t50 * t69 + t51 * t70 + t36 * t71 + t89 * t28 + t90 * t27 + t21 * t91 + t55 * t92 + t94 * t26 + t14 * t115 + t13 * t116 + t12 * t117 + t118 * t16 + t140 * t88 + t5 * t146 + t20 * t147 + t64 * t153 + t126 * t154 + t4 * t165 + t2 * t166 + t3 * t167 + t1 * t168 - t169 * t63 + t172 * t61 + t30 * t175 + t189 * t17 + t254 * t53 + t255 * t52 + t265 * t105 + t266 * t104 + ((-t75 / 0.2e1 + t77 / 0.2e1 + t38 / 0.2e1 + t39 / 0.2e1 + Ifges(5,4) * t256 + t66 * mrSges(5,1) - t72 * mrSges(4,3)) * t219 + (Ifges(4,6) * t256 + t74 / 0.2e1 - t76 / 0.2e1 - t73 * mrSges(4,3) + t65 * mrSges(5,1) + t268 * t218 + t267 * t215) * t216 + ((t108 - t110) * t219 + (t107 - t109) * t216 + m(5) * (t216 * t65 + t219 * t66) + m(4) * (-t216 * t73 - t219 * t72)) * pkin(9)) * qJD(3) + (-Ifges(3,6) * t279 + (-t208 / 0.2e1 - t209 / 0.2e1) * t220) * t213 + (-t32 * mrSges(4,3) + t29 * mrSges(5,1) - t57 / 0.2e1 + t59 / 0.2e1 + t7 / 0.2e1 + t8 / 0.2e1 + t295 + (-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t262 + (t79 - t81) * pkin(9)) * t216 + (t31 * mrSges(4,3) - t24 * mrSges(5,1) - t56 / 0.2e1 + t58 / 0.2e1 - t296 + t271 * t218 + t269 * t215 + (-Ifges(5,5) / 0.2e1 + Ifges(4,6) / 0.2e1) * t262 + (-t78 + t80) * pkin(9) + (t215 * t268 - t218 * t267) * qJD(5)) * t219 + m(4) * (-pkin(2) * t138 + t211 * t31 - t312 * t32) + m(5) * (t140 * t64 + t172 * t30 - t211 * t24 + t29 * t312) + m(6) * (t14 * t51 - t169 * t55 + t189 * t20 + t3 * t320 + t308 * t50 + t4 * t94) + t320 * t25; -0.2e1 * pkin(2) * t154 + 0.2e1 * t320 * t114 + 0.2e1 * t94 * t115 + 0.2e1 * t90 * t116 + 0.2e1 * t89 * t117 + 0.2e1 * t118 * t91 + 0.2e1 * t67 * t146 - 0.2e1 * t169 * t147 + 0.2e1 * t172 * t153 + 0.2e1 * t51 * t165 + 0.2e1 * t36 * t166 + 0.2e1 * t50 * t167 + 0.2e1 * t33 * t168 + 0.2e1 * t189 * t92 + 0.2e1 * (m(5) * t172 + t175) * t140 + (-t169 * t189 + t320 * t50 + t51 * t94) * t319 + (t118 * t67 + t33 * t89 + t36 * t90) * t318 + (-t156 + t164 + t83 + t84 + (t215 * t284 + t218 * t285 + t178 - t184) * qJD(3)) * t216 + (-t155 + t161 + (t82 - t85) * t218 + (-t86 - t87) * t215 + (t215 * t285 - t218 * t284) * qJD(5) + (t130 + t131 - t179 + t187) * qJD(3)) * t219; (t63 - t107) * qJD(4) + (-t78 + t17) * qJ(4) + m(5) * (-pkin(3) * t29 - qJ(4) * t24 - qJD(4) * t65) + m(6) * (qJ(4) * t20 + qJD(4) * t55 - t288 * t4 - t290 * t3) + m(7) * (-t1 * t290 + t135 * t21 + t171 * t5 + t2 * t288) - pkin(3) * t79 + t228 + t135 * t62 + t21 * t151 + t55 * t152 + t171 * t16 + t5 * t176 + t20 * t177 - t24 * mrSges(5,3) + t29 * mrSges(5,2) - t31 * mrSges(4,2) + t32 * mrSges(4,1) + t248 * t53 + t249 * t103 + t250 * t52 - t251 * t105 + t252 * t142 - t253 * t104 + ((-t12 * mrSges(7,2) - mrSges(6,3) * t308 - t268) * t218 + (-t13 * mrSges(7,2) + t14 * mrSges(6,3) - t267) * t215 - (m(6) * t233 + m(7) * t234 + t225) * t315) * qJD(5) + (t2 * mrSges(7,2) - t4 * mrSges(6,3) - t309 * t315 - t269) * t218 + (-t1 * mrSges(7,2) - t3 * mrSges(6,3) - t310 * t315 + t271) * t215; qJ(4) * t92 + t118 * t151 + t135 * t146 + t189 * t152 + t171 * t91 + t67 * t176 + t208 + t209 + t252 * t216 + m(7) * (t135 * t118 + t171 * t67) + ((-mrSges(5,1) * qJ(4) - Ifges(4,6)) * t216 + (-pkin(3) * mrSges(5,1) - Ifges(5,4) + t249) * t219 + (m(5) * t237 - t219 * mrSges(4,1) + t216 * mrSges(4,2) + t175) * pkin(9)) * qJD(3) + (t36 * mrSges(7,2) - t51 * mrSges(6,3) + t253 * t219 - t250 * t278 + (-t89 * mrSges(7,2) - mrSges(6,3) * t320 - t219 * t248 + t255) * qJD(5) - (t282 * qJD(5) + m(7) * (qJD(5) * t89 - t36) + m(6) * (t51 + t324) + t286) * t315 + t265) * t218 + (-t50 * mrSges(6,3) - t33 * mrSges(7,2) + t251 * t219 + t248 * t278 + (-t90 * mrSges(7,2) + t94 * mrSges(6,3) - t219 * t250 - t254) * qJD(5) - (t283 * qJD(5) + m(7) * (qJD(5) * t90 + t33) + m(6) * (-qJD(5) * t94 + t50) + t287) * t315 - t266) * t215 - (m(6) * qJ(4) + t177) * t169 + (m(5) * t211 + m(6) * t189 + t219 * mrSges(5,1) + t147) * qJD(4); t152 * t323 + 0.2e1 * t151 * t171 + (-t162 - t163) * t218 + (-t157 + t160) * t215 + ((t180 - t183) * t218 + (-t185 - t186) * t215) * qJD(5) + 0.2e1 * (m(7) * t171 + t176) * t135 + (0.2e1 * mrSges(5,3) + 0.2e1 * t177 + (m(5) + m(6)) * t323) * qJD(4); t309 * t218 + t310 * t215 + t225 * qJD(5) + m(7) * (qJD(5) * t234 + t1 * t215 - t2 * t218) + m(6) * (qJD(5) * t233 + t215 * t3 + t218 * t4) + m(5) * t29 + t79; t286 * t218 + t287 * t215 + (m(5) * pkin(9) + mrSges(5,1)) * t277 + (t215 * t283 + t218 * t282) * qJD(5) + m(7) * (t215 * t33 - t218 * t36 + (t215 * t90 + t218 * t89) * qJD(5)) + m(6) * (t215 * t50 + t218 * t51 + (-t215 * t94 + t218 * t320) * qJD(5)); 0; 0; -pkin(5) * t27 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t12) + qJD(6) * t68 + qJ(6) * t28 + t1 * mrSges(7,3) - t3 * mrSges(6,2) + t4 * mrSges(6,1) - t2 * mrSges(7,1) + t322; -pkin(5) * t116 + m(7) * (-pkin(5) * t36 + qJ(6) * t33 + qJD(6) * t89) + qJD(6) * t168 + qJ(6) * t117 + t33 * mrSges(7,3) - t36 * mrSges(7,1) - t50 * mrSges(6,2) + t51 * mrSges(6,1) - Ifges(7,6) * t257 + (-t299 - t329) * t274 + t246 + t281; t206 + qJD(6) * t245 + ((-mrSges(7,2) * qJ(6) - Ifges(6,6)) * t218 + (mrSges(7,2) * pkin(5) - t326) * t215 - t222 * t315) * qJD(5); m(7) * t273 + qJD(5) * t222; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t2 + t27; m(7) * t36 + t116; qJD(5) * t245; m(7) * t276; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
