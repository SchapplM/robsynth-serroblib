% Calculate time derivative of joint inertia matrix for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:25:05
% EndTime: 2019-03-09 19:25:22
% DurationCPUTime: 8.20s
% Computational Cost: add. (9820->691), mult. (24919->966), div. (0->0), fcn. (23628->10), ass. (0->267)
t350 = Ifges(5,4) + Ifges(4,5);
t349 = Ifges(5,2) + Ifges(4,3);
t348 = Ifges(5,6) - Ifges(4,6);
t228 = cos(qJ(5));
t223 = sin(qJ(6));
t227 = cos(qJ(6));
t277 = t223 ^ 2 + t227 ^ 2;
t256 = t277 * mrSges(7,3);
t347 = (mrSges(6,2) - t256) * t228;
t224 = sin(qJ(5));
t225 = sin(qJ(3));
t229 = cos(qJ(3));
t163 = t224 * t225 + t228 * t229;
t121 = (qJD(3) - qJD(5)) * t163;
t164 = -t224 * t229 + t225 * t228;
t264 = qJD(6) * t223;
t237 = -t227 * t121 + t164 * t264;
t182 = -t227 * mrSges(7,1) + mrSges(7,2) * t223;
t295 = mrSges(6,1) - t182;
t338 = -m(7) * pkin(5) - t295;
t344 = m(7) * pkin(11);
t343 = t223 / 0.2e1;
t315 = t227 / 0.2e1;
t231 = -pkin(3) - pkin(4);
t179 = -t224 * qJ(4) + t228 * t231;
t148 = t228 * qJD(4) + qJD(5) * t179;
t342 = t148 * mrSges(6,2);
t222 = cos(pkin(6));
t221 = sin(pkin(6));
t226 = sin(qJ(2));
t285 = t221 * t226;
t158 = t222 * t225 + t229 * t285;
t230 = cos(qJ(2));
t276 = qJD(2) * t221;
t258 = t230 * t276;
t126 = qJD(3) * t158 + t225 * t258;
t157 = -t222 * t229 + t225 * t285;
t127 = -qJD(3) * t157 + t229 * t258;
t284 = t221 * t230;
t313 = pkin(1) * t222;
t280 = pkin(8) * t284 + t226 * t313;
t153 = t280 * qJD(2);
t52 = t126 * pkin(3) - t127 * qJ(4) - t158 * qJD(4) + t153;
t246 = -t229 * pkin(3) - t225 * qJ(4);
t180 = t228 * qJ(4) + t224 * t231;
t273 = qJD(3) * t229;
t274 = qJD(3) * t225;
t341 = t273 * t350 + t348 * t274;
t340 = t277 * t228;
t145 = pkin(9) * t222 + t280;
t146 = (-pkin(2) * t230 - pkin(9) * t226 - pkin(1)) * t221;
t92 = -t225 * t145 + t146 * t229;
t84 = pkin(3) * t284 - t92;
t65 = pkin(4) * t284 - pkin(10) * t158 + t84;
t93 = t229 * t145 + t225 * t146;
t83 = -qJ(4) * t284 + t93;
t70 = pkin(10) * t157 + t83;
t307 = t224 * t65 + t228 * t70;
t339 = qJD(5) * t307;
t275 = qJD(2) * t226;
t259 = t221 * t275;
t151 = (pkin(2) * t226 - pkin(9) * t230) * t276;
t160 = -pkin(8) * t285 + t230 * t313;
t152 = t160 * qJD(2);
t56 = -t145 * t273 - t146 * t274 + t151 * t229 - t225 * t152;
t35 = -pkin(10) * t127 + t231 * t259 - t56;
t55 = -t145 * t274 + t146 * t273 + t225 * t151 + t229 * t152;
t47 = qJ(4) * t259 - qJD(4) * t284 + t55;
t36 = pkin(10) * t126 + t47;
t9 = -t224 * t36 + t228 * t35 - t339;
t263 = qJD(6) * t227;
t171 = Ifges(7,4) * t263 - Ifges(7,2) * t264;
t173 = Ifges(7,1) * t263 - Ifges(7,4) * t264;
t303 = Ifges(7,4) * t223;
t186 = Ifges(7,2) * t227 + t303;
t302 = Ifges(7,4) * t227;
t188 = Ifges(7,1) * t223 + t302;
t234 = -(t186 * t223 - t188 * t227) * qJD(6) + t227 * t171 + t223 * t173;
t337 = t348 * t126 + t127 * t350 + t349 * t259;
t336 = 2 * m(6);
t335 = 0.2e1 * m(7);
t334 = -2 * mrSges(3,3);
t333 = -2 * mrSges(6,3);
t324 = pkin(9) - pkin(10);
t178 = t324 * t274;
t191 = t324 * t229;
t252 = qJD(3) * t191;
t261 = t324 * t225;
t133 = t228 * t191 + t224 * t261;
t270 = qJD(5) * t133;
t81 = -t224 * t178 - t228 * t252 + t270;
t332 = 0.2e1 * t81;
t132 = t224 * t191 - t228 * t261;
t331 = 0.2e1 * t132;
t330 = 0.2e1 * t153;
t247 = mrSges(7,1) * t223 + mrSges(7,2) * t227;
t166 = t247 * qJD(6);
t329 = -0.2e1 * t166;
t242 = t157 * t228 - t158 * t224;
t49 = qJD(5) * t242 + t126 * t224 + t127 * t228;
t104 = t157 * t224 + t158 * t228;
t86 = t104 * t227 + t223 * t284;
t23 = -qJD(6) * t86 - t223 * t49 - t227 * t259;
t328 = t23 / 0.2e1;
t266 = qJD(5) * t228;
t268 = qJD(5) * t224;
t122 = t224 * t273 + t225 * t266 - t228 * t274 - t229 * t268;
t238 = t223 * t121 + t164 * t263;
t38 = -Ifges(7,4) * t237 - Ifges(7,2) * t238 + Ifges(7,6) * t122;
t327 = t38 / 0.2e1;
t39 = -Ifges(7,1) * t237 - Ifges(7,4) * t238 + Ifges(7,5) * t122;
t326 = t39 / 0.2e1;
t89 = Ifges(7,5) * t163 + (Ifges(7,1) * t227 - t303) * t164;
t325 = t89 / 0.2e1;
t323 = -t164 / 0.2e1;
t205 = Ifges(7,6) * t264;
t322 = Ifges(7,5) * t263 / 0.2e1 - t205 / 0.2e1;
t321 = t171 / 0.2e1;
t320 = t173 / 0.2e1;
t319 = Ifges(7,5) * t343 + Ifges(7,6) * t315;
t318 = -t186 / 0.2e1;
t317 = t188 / 0.2e1;
t316 = -t223 / 0.2e1;
t312 = -qJD(6) / 0.2e1;
t311 = qJD(6) / 0.2e1;
t85 = -t104 * t223 + t227 * t284;
t24 = qJD(6) * t85 - t223 * t259 + t227 * t49;
t10 = -mrSges(7,1) * t23 + mrSges(7,2) * t24;
t42 = -mrSges(6,1) * t259 - mrSges(6,3) * t49;
t310 = t10 - t42;
t43 = -mrSges(7,1) * t85 + mrSges(7,2) * t86;
t91 = mrSges(6,1) * t284 - mrSges(6,3) * t104;
t309 = t43 - t91;
t48 = qJD(5) * t104 - t126 * t228 + t127 * t224;
t308 = Ifges(6,5) * t49 - Ifges(6,6) * t48;
t306 = mrSges(7,3) * t164;
t305 = Ifges(4,4) * t225;
t304 = Ifges(4,4) * t229;
t301 = Ifges(5,5) * t225;
t300 = Ifges(5,5) * t229;
t299 = Ifges(7,5) * t227;
t298 = t132 * t81;
t297 = t152 * mrSges(3,2);
t296 = t153 * mrSges(3,1);
t26 = pkin(11) * t284 + t307;
t144 = -t222 * pkin(2) - t160;
t82 = t157 * pkin(3) - t158 * qJ(4) + t144;
t69 = -pkin(4) * t157 - t82;
t29 = -pkin(5) * t242 - pkin(11) * t104 + t69;
t12 = -t223 * t26 + t227 * t29;
t294 = qJD(6) * t12;
t13 = t223 * t29 + t227 * t26;
t293 = qJD(6) * t13;
t181 = -pkin(2) + t246;
t162 = t229 * pkin(4) - t181;
t98 = pkin(5) * t163 - pkin(11) * t164 + t162;
t66 = -t133 * t223 + t227 * t98;
t292 = qJD(6) * t66;
t67 = t133 * t227 + t223 * t98;
t291 = qJD(6) * t67;
t149 = t224 * qJD(4) + qJD(5) * t180;
t290 = t132 * t149;
t289 = t148 * t224;
t288 = t149 * t228;
t282 = t228 * t166;
t281 = Ifges(6,5) * t121 - Ifges(6,6) * t122;
t279 = qJ(4) * t273 + t225 * qJD(4);
t272 = qJD(4) * t229;
t271 = qJD(5) * t132;
t269 = qJD(5) * t223;
t267 = qJD(5) * t227;
t177 = -pkin(11) + t180;
t265 = qJD(6) * t177;
t262 = 0.2e1 * t221;
t3 = Ifges(7,5) * t24 + Ifges(7,6) * t23 + Ifges(7,3) * t48;
t140 = t231 * t274 + t279;
t59 = pkin(5) * t122 - pkin(11) * t121 + t140;
t80 = -t228 * t178 + t224 * t252 - t271;
t16 = t223 * t59 + t227 * t80 + t292;
t255 = -t16 + t292;
t17 = -t223 * t80 + t227 * t59 - t291;
t254 = t17 + t291;
t253 = t277 * t148;
t251 = t121 * t317 + t327;
t250 = t121 * t318 + t326;
t40 = -pkin(4) * t126 - t52;
t11 = pkin(5) * t48 - pkin(11) * t49 + t40;
t8 = t224 * t35 + t228 * t36 + t65 * t266 - t268 * t70;
t6 = -pkin(11) * t259 + t8;
t1 = t11 * t223 + t227 * t6 + t294;
t2 = t11 * t227 - t223 * t6 - t293;
t248 = t1 * t227 - t2 * t223;
t183 = -t229 * mrSges(5,1) - t225 * mrSges(5,3);
t14 = -mrSges(7,2) * t48 + mrSges(7,3) * t23;
t15 = mrSges(7,1) * t48 - mrSges(7,3) * t24;
t245 = t227 * t14 - t223 * t15;
t27 = -t224 * t70 + t228 * t65;
t102 = -mrSges(5,1) * t259 + t127 * mrSges(5,2);
t31 = Ifges(7,4) * t86 + Ifges(7,2) * t85 - Ifges(7,6) * t242;
t32 = Ifges(7,1) * t86 + Ifges(7,4) * t85 - Ifges(7,5) * t242;
t239 = t31 * t316 + t315 * t32;
t235 = -t12 * t263 - t13 * t264 + t248;
t233 = -t80 * mrSges(6,2) + t122 * t319 + t132 * t166 + t163 * t322 + t281;
t37 = -Ifges(7,5) * t237 - Ifges(7,6) * t238 + Ifges(7,3) * t122;
t25 = -pkin(5) * t284 - t27;
t7 = pkin(5) * t259 - t9;
t232 = t9 * mrSges(6,1) - t8 * mrSges(6,2) - Ifges(6,3) * t259 + t25 * t166 + t7 * t182 + t186 * t328 + t24 * t317 - t242 * t322 + t319 * t48 + t320 * t86 + t321 * t85 + t308;
t195 = Ifges(3,5) * t258;
t190 = Ifges(4,1) * t225 + t304;
t189 = Ifges(5,1) * t225 - t300;
t187 = Ifges(4,2) * t229 + t305;
t185 = -Ifges(5,3) * t229 + t301;
t176 = pkin(5) - t179;
t175 = (Ifges(4,1) * t229 - t305) * qJD(3);
t174 = (Ifges(5,1) * t229 + t301) * qJD(3);
t172 = (-Ifges(4,2) * t225 + t304) * qJD(3);
t170 = (Ifges(5,3) * t225 + t300) * qJD(3);
t168 = (mrSges(4,1) * t225 + mrSges(4,2) * t229) * qJD(3);
t167 = (mrSges(5,1) * t225 - mrSges(5,3) * t229) * qJD(3);
t156 = pkin(3) * t274 - t279;
t131 = mrSges(5,1) * t284 + mrSges(5,2) * t158;
t130 = -mrSges(4,1) * t284 - mrSges(4,3) * t158;
t129 = mrSges(4,2) * t284 - mrSges(4,3) * t157;
t128 = -mrSges(5,2) * t157 - mrSges(5,3) * t284;
t125 = Ifges(6,1) * t164 - Ifges(6,4) * t163;
t124 = Ifges(6,4) * t164 - Ifges(6,2) * t163;
t123 = mrSges(6,1) * t163 + mrSges(6,2) * t164;
t110 = mrSges(7,1) * t163 - t227 * t306;
t109 = -mrSges(7,2) * t163 - t223 * t306;
t106 = t247 * t164;
t105 = mrSges(5,1) * t157 - mrSges(5,3) * t158;
t101 = mrSges(4,1) * t259 - mrSges(4,3) * t127;
t100 = -mrSges(4,2) * t259 - mrSges(4,3) * t126;
t99 = -mrSges(5,2) * t126 + mrSges(5,3) * t259;
t97 = Ifges(4,1) * t158 - Ifges(4,4) * t157 - Ifges(4,5) * t284;
t96 = Ifges(5,1) * t158 - Ifges(5,4) * t284 + Ifges(5,5) * t157;
t95 = Ifges(4,4) * t158 - Ifges(4,2) * t157 - Ifges(4,6) * t284;
t94 = Ifges(5,5) * t158 - Ifges(5,6) * t284 + Ifges(5,3) * t157;
t90 = -mrSges(6,2) * t284 + mrSges(6,3) * t242;
t88 = Ifges(7,6) * t163 + (-Ifges(7,2) * t223 + t302) * t164;
t87 = Ifges(7,3) * t163 + (-Ifges(7,6) * t223 + t299) * t164;
t79 = mrSges(4,1) * t126 + mrSges(4,2) * t127;
t78 = mrSges(5,1) * t126 - mrSges(5,3) * t127;
t77 = Ifges(6,1) * t121 - Ifges(6,4) * t122;
t76 = Ifges(6,4) * t121 - Ifges(6,2) * t122;
t75 = mrSges(6,1) * t122 + mrSges(6,2) * t121;
t74 = Ifges(4,1) * t127 - Ifges(4,4) * t126 + Ifges(4,5) * t259;
t73 = Ifges(5,1) * t127 + Ifges(5,4) * t259 + Ifges(5,5) * t126;
t72 = Ifges(4,4) * t127 - Ifges(4,2) * t126 + Ifges(4,6) * t259;
t71 = Ifges(5,5) * t127 + Ifges(5,6) * t259 + Ifges(5,3) * t126;
t63 = -mrSges(6,1) * t242 + mrSges(6,2) * t104;
t61 = -mrSges(7,2) * t122 - mrSges(7,3) * t238;
t60 = mrSges(7,1) * t122 + mrSges(7,3) * t237;
t58 = Ifges(6,1) * t104 + Ifges(6,4) * t242 + Ifges(6,5) * t284;
t57 = Ifges(6,4) * t104 + Ifges(6,2) * t242 + Ifges(6,6) * t284;
t54 = -mrSges(7,1) * t242 - mrSges(7,3) * t86;
t53 = mrSges(7,2) * t242 + mrSges(7,3) * t85;
t51 = -pkin(3) * t259 - t56;
t50 = mrSges(7,1) * t238 - mrSges(7,2) * t237;
t41 = mrSges(6,2) * t259 - mrSges(6,3) * t48;
t30 = Ifges(7,5) * t86 + Ifges(7,6) * t85 - Ifges(7,3) * t242;
t20 = mrSges(6,1) * t48 + mrSges(6,2) * t49;
t19 = Ifges(6,1) * t49 - Ifges(6,4) * t48 - Ifges(6,5) * t259;
t18 = Ifges(6,4) * t49 - Ifges(6,2) * t48 - Ifges(6,6) * t259;
t5 = Ifges(7,1) * t24 + Ifges(7,4) * t23 + Ifges(7,5) * t48;
t4 = Ifges(7,4) * t24 + Ifges(7,2) * t23 + Ifges(7,6) * t48;
t21 = [0.2e1 * m(4) * (t144 * t153 + t55 * t93 + t56 * t92) + 0.2e1 * m(5) * (t47 * t83 + t51 * t84 + t52 * t82) + (t1 * t13 + t12 * t2 + t25 * t7) * t335 + (mrSges(4,1) * t157 + mrSges(4,2) * t158) * t330 - (-t18 + t3) * t242 + (t195 - 0.2e1 * t296 - 0.2e1 * t297) * t222 + (-t57 + t30) * t48 + (mrSges(3,3) * t226 * t330 + (0.2e1 * mrSges(3,3) * t152 + t308 - t337) * t230 + ((t160 * t334 + Ifges(3,5) * t222 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t230) * t262) * t230 + (t280 * t334 - Ifges(6,5) * t104 - 0.2e1 * Ifges(3,6) * t222 - Ifges(6,6) * t242 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t226) * t262 + t350 * t158 + t348 * t157 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - 0.2e1 * Ifges(6,3) - t349) * t284) * t226) * qJD(2)) * t221 + (t27 * t9 + t307 * t8 + t40 * t69) * t336 + 0.2e1 * t307 * t41 + 0.2e1 * m(3) * (t152 * t280 - t153 * t160) + (t96 + t97) * t127 + (t94 - t95) * t126 + t157 * t71 - t157 * t72 + t158 * t74 + t158 * t73 + 0.2e1 * t144 * t79 + 0.2e1 * t47 * t128 + 0.2e1 * t55 * t129 + 0.2e1 * t56 * t130 + 0.2e1 * t51 * t131 + t104 * t19 + 0.2e1 * t52 * t105 + 0.2e1 * t83 * t99 + 0.2e1 * t93 * t100 + 0.2e1 * t92 * t101 + 0.2e1 * t84 * t102 + t85 * t4 + t86 * t5 + 0.2e1 * t8 * t90 + 0.2e1 * t9 * t91 + 0.2e1 * t82 * t78 + 0.2e1 * t69 * t20 + t49 * t58 + 0.2e1 * t40 * t63 + 0.2e1 * t1 * t53 + 0.2e1 * t2 * t54 + 0.2e1 * t27 * t42 + 0.2e1 * t7 * t43 + t23 * t31 + t24 * t32 + 0.2e1 * t25 * t10 + 0.2e1 * t13 * t14 + 0.2e1 * t12 * t15; t195 + t24 * t325 + t86 * t326 + t85 * t327 + t88 * t328 - (-t76 / 0.2e1 + t37 / 0.2e1) * t242 + (-m(4) * t153 - t79) * pkin(2) + (t5 * t315 + t4 * t316 - t9 * mrSges(6,3) + t19 / 0.2e1 + (-t227 * t31 / 0.2e1 + t32 * t316) * qJD(6)) * t164 + t309 * t81 + t310 * t132 + (t55 * mrSges(4,3) + t47 * mrSges(5,2) - t153 * mrSges(4,1) - t71 / 0.2e1 + t72 / 0.2e1) * t229 + ((t100 + t99) * t229 + (-t101 + t102) * t225 + ((-t130 + t131) * t229 + (-t128 - t129) * t225) * qJD(3) + m(5) * (t51 * t225 + t47 * t229 + t273 * t84 - t274 * t83) + m(4) * (-t56 * t225 + t55 * t229 - t273 * t92 - t274 * t93)) * pkin(9) + (-t56 * mrSges(4,3) + t51 * mrSges(5,2) + t153 * mrSges(4,2) + t73 / 0.2e1 + t74 / 0.2e1) * t225 + m(7) * (t1 * t67 + t12 * t17 + t13 * t16 + t132 * t7 + t2 * t66 + t25 * t81) + (-t27 * mrSges(6,3) + t58 / 0.2e1 + t239) * t121 + m(6) * (-t132 * t9 + t133 * t8 + t140 * t69 + t162 * t40 - t27 * t81 + t307 * t80) + (-t307 * mrSges(6,3) - t57 / 0.2e1 + t30 / 0.2e1) * t122 + (-t124 / 0.2e1 + t87 / 0.2e1) * t48 + ((t96 / 0.2e1 + t97 / 0.2e1 + t84 * mrSges(5,2) - t92 * mrSges(4,3)) * t229 + (t94 / 0.2e1 - t95 / 0.2e1 - t83 * mrSges(5,2) - t93 * mrSges(4,3)) * t225) * qJD(3) + t181 * t78 + t52 * t183 + t82 * t167 + t144 * t168 + t162 * t20 + t156 * t105 + t133 * t41 + t140 * t63 + t40 * t123 + t49 * t125 / 0.2e1 + t104 * t77 / 0.2e1 + t7 * t106 + t1 * t109 + t2 * t110 + t80 * t90 + t69 * t75 + t66 * t15 + t67 * t14 + t12 * t60 + t13 * t61 + t25 * t50 + t16 * t53 + t17 * t54 + (-t8 * mrSges(6,3) + t3 / 0.2e1 - t18 / 0.2e1) * t163 + (t174 / 0.2e1 + t175 / 0.2e1) * t158 + (-t172 / 0.2e1 + t170 / 0.2e1) * t157 + (t189 / 0.2e1 + t190 / 0.2e1) * t127 + (t185 / 0.2e1 - t187 / 0.2e1) * t126 - t296 - t297 + ((Ifges(6,5) * t323 + Ifges(6,6) * t163 / 0.2e1 - Ifges(3,6) + (-Ifges(5,6) / 0.2e1 + Ifges(4,6) / 0.2e1) * t229 + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t225) * t275 + (t281 / 0.2e1 - t341 / 0.2e1) * t230) * t221 + m(5) * (t156 * t82 + t181 * t52); -0.2e1 * pkin(2) * t168 + t106 * t332 + 0.2e1 * t16 * t109 + 0.2e1 * t17 * t110 + 0.2e1 * t140 * t123 + t50 * t331 + 0.2e1 * t162 * t75 + 0.2e1 * t181 * t167 + 0.2e1 * t66 * t60 + 0.2e1 * t67 * t61 + (t172 - t170) * t229 + (t174 + t175) * t225 + 0.2e1 * (m(5) * t181 + t183) * t156 + (t133 * t80 + t140 * t162 + t298) * t336 + (t16 * t67 + t17 * t66 + t298) * t335 + (t333 * t80 + t37 - t76) * t163 + (t133 * t333 - t124 + t87) * t122 + (mrSges(6,3) * t331 - t223 * t88 + t227 * t89 + t125) * t121 + ((t189 + t190) * t229 + (t185 - t187) * t225) * qJD(3) + (mrSges(6,3) * t332 - t223 * t38 + t227 * t39 + t77 + (-t223 * t89 - t227 * t88) * qJD(6)) * t164; t337 + (m(7) * (-t12 * t148 - t13 * t265 - t177 * t2) - t177 * t15 - t148 * t54 + t31 * t311 - t53 * t265 - t5 / 0.2e1 + (t2 + t293) * mrSges(7,3)) * t223 + (t32 * t312 + m(7) * (t1 * t177 - t12 * t265 + t13 * t148) + t177 * t14 + t148 * t53 - t54 * t265 - t4 / 0.2e1 + (-t1 + t294) * mrSges(7,3)) * t227 + t309 * t149 + m(6) * (t148 * t307 - t149 * t27 + t179 * t9 + t180 * t8) - t232 + t180 * t41 + t176 * t10 + t179 * t42 + t148 * t90 + qJD(4) * t128 + qJ(4) * t99 - pkin(3) * t102 - t55 * mrSges(4,2) + t56 * mrSges(4,1) - t51 * mrSges(5,1) + t47 * mrSges(5,3) + m(7) * (t149 * t25 + t176 * t7) + m(5) * (-pkin(3) * t51 + qJ(4) * t47 + qJD(4) * t83); (qJD(3) * t246 + t272) * mrSges(5,2) + t295 * t81 + m(6) * (t133 * t148 - t179 * t81 + t180 * t80 + t290) + (-t121 * t179 - t122 * t180 - t148 * t163 + t149 * t164) * mrSges(6,3) + m(7) * (t176 * t81 + t290) + (m(5) * t272 + (m(5) * t246 - t229 * mrSges(4,1) + t225 * mrSges(4,2) + t183) * qJD(3)) * pkin(9) + (t89 * t312 + m(7) * (t148 * t67 + t16 * t177 - t265 * t66) + t177 * t61 + t148 * t109 - t110 * t265 + (-t173 / 0.2e1 + t186 * t311) * t164 + t255 * mrSges(7,3) - t251) * t227 + (m(7) * (-t148 * t66 - t17 * t177 - t265 * t67) - t177 * t60 - t148 * t110 + t88 * t311 - t109 * t265 + (t188 * t311 + t321) * t164 + t254 * mrSges(7,3) - t250) * t223 + t176 * t50 + t149 * t106 - t233 + t341; t176 * t329 + 0.2e1 * t342 + (t149 * t176 + t177 * t253) * t335 + (t148 * t180 - t149 * t179) * t336 + t234 + 0.2e1 * t295 * t149 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4) - 0.2e1 * t256 * t148; m(5) * t51 + ((-t223 * t54 + t227 * t53 + t90) * qJD(5) + m(7) * (-t12 * t269 + t13 * t267 - t7) + m(6) * (t9 + t339) - t310) * t228 + (t41 + (-t223 * t53 - t227 * t54) * qJD(6) + t309 * qJD(5) + m(7) * (qJD(5) * t25 + t235) + m(6) * (-qJD(5) * t27 + t8) + t245) * t224 + t102; (m(5) * pkin(9) + mrSges(5,2)) * t273 + (-t121 * mrSges(6,3) - t50 + (-t163 * mrSges(6,3) + t109 * t227 - t110 * t223) * qJD(5) + m(7) * (t267 * t67 - t269 * t66 - t81) + m(6) * (-t81 + t270)) * t228 + (qJD(5) * t106 - t223 * t60 + t227 * t61 + (-t109 * t223 - t110 * t227) * qJD(6) + (qJD(5) * t164 - t122) * mrSges(6,3) + m(7) * (t16 * t227 - t17 * t223 - t263 * t66 - t264 * t67 + t271) + m(6) * (t80 + t271)) * t224; t282 + m(7) * (t277 * t289 - t288) + m(6) * (-t288 + t289) + (t295 * t224 + t347 + m(7) * (t176 * t224 + t177 * t340) + m(6) * (-t179 * t224 + t180 * t228)) * qJD(5); (-0.1e1 + t277) * t224 * t266 * t335; t232 + (-m(7) * t7 - t10) * pkin(5) + (m(7) * t235 - t263 * t54 - t264 * t53 + t245) * pkin(11) + t4 * t315 + t5 * t343 + t239 * qJD(6) + ((-t12 * t227 - t13 * t223) * qJD(6) + t248) * mrSges(7,3); -pkin(5) * t50 + t338 * t81 + (t164 * t320 + t16 * mrSges(7,3) + (-t66 * mrSges(7,3) + t164 * t318 + t325) * qJD(6) + (-m(7) * t255 - qJD(6) * t110 + t61) * pkin(11) + t251) * t227 + (t171 * t323 - t17 * mrSges(7,3) + (-t88 / 0.2e1 + t188 * t323 - t67 * mrSges(7,3)) * qJD(6) + (-m(7) * t254 - qJD(6) * t109 - t60) * pkin(11) + t250) * t223 + t233; -t342 + (pkin(5) + t176) * t166 + (mrSges(7,3) + t344) * t253 + t338 * t149 - t234; -t282 + (t338 * t224 + t340 * t344 - t347) * qJD(5); pkin(5) * t329 + t234; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t3; mrSges(7,1) * t17 - mrSges(7,2) * t16 + t37; t205 - t247 * t148 + (t177 * t182 - t299) * qJD(6); (t224 * t264 - t227 * t266) * mrSges(7,2) + (-t223 * t266 - t224 * t263) * mrSges(7,1); -t205 + (pkin(11) * t182 + t299) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
