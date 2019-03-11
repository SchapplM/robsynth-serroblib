% Calculate time derivative of joint inertia matrix for
% S6RRRPRP11
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:47
% EndTime: 2019-03-09 17:42:01
% DurationCPUTime: 6.71s
% Computational Cost: add. (5611->676), mult. (14765->912), div. (0->0), fcn. (13115->8), ass. (0->257)
t323 = Ifges(5,1) + Ifges(4,3);
t322 = Ifges(4,5) - Ifges(5,4);
t321 = Ifges(6,5) + Ifges(7,5);
t320 = -Ifges(4,6) + Ifges(5,5);
t319 = Ifges(6,6) + Ifges(7,6);
t226 = cos(qJ(5));
t223 = sin(qJ(5));
t227 = cos(qJ(3));
t275 = qJD(5) * t227;
t261 = t223 * t275;
t224 = sin(qJ(3));
t279 = qJD(3) * t224;
t233 = t226 * t279 + t261;
t259 = t223 * t279;
t260 = t226 * t275;
t232 = t259 - t260;
t309 = -2 * mrSges(7,3);
t222 = cos(pkin(6));
t228 = cos(qJ(2));
t221 = sin(pkin(6));
t281 = qJD(2) * t221;
t262 = t228 * t281;
t225 = sin(qJ(2));
t288 = t221 * t225;
t264 = t224 * t288;
t103 = -qJD(3) * t264 + (qJD(3) * t222 + t262) * t227;
t144 = t222 * t224 + t227 * t288;
t102 = qJD(3) * t144 + t224 * t262;
t143 = -t222 * t227 + t264;
t287 = t221 * t228;
t234 = -t143 * t223 + t226 * t287;
t280 = qJD(2) * t225;
t263 = t221 * t280;
t54 = qJD(5) * t234 + t102 * t226 - t223 * t263;
t104 = t143 * t226 + t223 * t287;
t55 = qJD(5) * t104 + t102 * t223 + t226 * t263;
t8 = Ifges(7,5) * t55 + Ifges(7,6) * t54 + Ifges(7,3) * t103;
t9 = Ifges(6,5) * t55 + Ifges(6,6) * t54 + Ifges(6,3) * t103;
t318 = t8 + t9;
t289 = qJ(4) * t224;
t308 = pkin(3) + pkin(10);
t151 = -t227 * t308 - pkin(2) - t289;
t307 = pkin(4) + pkin(9);
t192 = t307 * t224;
t95 = t226 * t151 + t223 * t192;
t207 = pkin(8) * t288;
t306 = pkin(1) * t228;
t147 = t222 * t306 - t207;
t148 = t222 * t225 * pkin(1) + pkin(8) * t287;
t127 = pkin(9) * t222 + t148;
t128 = (-pkin(2) * t228 - pkin(9) * t225 - pkin(1)) * t221;
t137 = (pkin(2) * t225 - pkin(9) * t228) * t281;
t138 = t147 * qJD(2);
t278 = qJD(3) * t227;
t35 = -t127 * t278 - t128 * t279 + t137 * t227 - t224 * t138;
t20 = pkin(4) * t103 - t263 * t308 - t35;
t139 = t148 * qJD(2);
t230 = -qJ(4) * t103 - qJD(4) * t144 + t139;
t23 = t102 * t308 + t230;
t276 = qJD(5) * t226;
t277 = qJD(5) * t223;
t73 = -t224 * t127 + t128 * t227;
t68 = pkin(3) * t287 - t73;
t45 = pkin(4) * t144 + pkin(10) * t287 + t68;
t126 = t207 + (-pkin(2) - t306) * t222;
t231 = -qJ(4) * t144 + t126;
t56 = t143 * t308 + t231;
t3 = t223 * t20 + t226 * t23 + t45 * t276 - t277 * t56;
t15 = t223 * t45 + t226 * t56;
t4 = -qJD(5) * t15 + t226 * t20 - t223 * t23;
t314 = t223 * t3 + t226 * t4;
t313 = 2 * m(6);
t312 = 2 * m(7);
t311 = -2 * mrSges(3,3);
t310 = 2 * mrSges(3,3);
t305 = pkin(9) * t224;
t219 = t227 * pkin(9);
t301 = Ifges(4,4) * t224;
t300 = Ifges(4,4) * t227;
t299 = Ifges(6,4) * t223;
t298 = Ifges(6,4) * t226;
t297 = Ifges(7,4) * t223;
t296 = Ifges(7,4) * t226;
t295 = Ifges(5,6) * t224;
t294 = Ifges(5,6) * t227;
t293 = t138 * mrSges(3,2);
t292 = t139 * mrSges(3,1);
t291 = t139 * mrSges(4,1);
t290 = t139 * mrSges(4,2);
t286 = t223 * t227;
t285 = t226 * t227;
t284 = qJ(6) + t308;
t74 = t227 * t127 + t224 * t128;
t240 = Ifges(7,2) * t226 + t297;
t131 = Ifges(7,6) * t224 - t227 * t240;
t241 = Ifges(6,2) * t226 + t299;
t132 = Ifges(6,6) * t224 - t227 * t241;
t283 = t131 + t132;
t242 = Ifges(7,1) * t223 + t296;
t133 = Ifges(7,5) * t224 - t227 * t242;
t243 = Ifges(6,1) * t223 + t298;
t134 = Ifges(6,5) * t224 - t227 * t243;
t282 = t133 + t134;
t193 = t227 * pkin(4) + t219;
t274 = t226 * qJD(6);
t273 = 0.2e1 * t221;
t10 = Ifges(7,4) * t55 + Ifges(7,2) * t54 + Ifges(7,6) * t103;
t11 = Ifges(6,4) * t55 + Ifges(6,2) * t54 + Ifges(6,6) * t103;
t271 = -t10 / 0.2e1 - t11 / 0.2e1;
t12 = Ifges(7,1) * t55 + Ifges(7,4) * t54 + Ifges(7,5) * t103;
t13 = Ifges(6,1) * t55 + Ifges(6,4) * t54 + Ifges(6,5) * t103;
t270 = -t12 / 0.2e1 - t13 / 0.2e1;
t41 = -Ifges(7,4) * t234 + Ifges(7,2) * t104 + Ifges(7,6) * t144;
t42 = -Ifges(6,4) * t234 + Ifges(6,2) * t104 + Ifges(6,6) * t144;
t269 = t41 / 0.2e1 + t42 / 0.2e1;
t43 = -Ifges(7,1) * t234 + Ifges(7,4) * t104 + Ifges(7,5) * t144;
t44 = -Ifges(6,1) * t234 + Ifges(6,4) * t104 + Ifges(6,5) * t144;
t268 = t43 / 0.2e1 + t44 / 0.2e1;
t186 = -Ifges(7,2) * t223 + t296;
t86 = -t186 * t275 + (Ifges(7,6) * t227 + t224 * t240) * qJD(3);
t187 = -Ifges(6,2) * t223 + t298;
t87 = -t187 * t275 + (Ifges(6,6) * t227 + t224 * t241) * qJD(3);
t267 = t86 / 0.2e1 + t87 / 0.2e1;
t189 = Ifges(7,1) * t226 - t297;
t88 = -t189 * t275 + (Ifges(7,5) * t227 + t224 * t242) * qJD(3);
t190 = Ifges(6,1) * t226 - t299;
t89 = -t190 * t275 + (Ifges(6,5) * t227 + t224 * t243) * qJD(3);
t266 = t88 / 0.2e1 + t89 / 0.2e1;
t265 = m(7) * pkin(5) + mrSges(7,1);
t257 = t287 / 0.2e1;
t256 = t131 / 0.2e1 + t132 / 0.2e1;
t255 = t133 / 0.2e1 + t134 / 0.2e1;
t238 = -Ifges(7,5) * t223 - Ifges(7,6) * t226;
t239 = -Ifges(6,5) * t223 - Ifges(6,6) * t226;
t254 = (t238 + t239) * qJD(5) / 0.2e1;
t163 = t240 * qJD(5);
t164 = t241 * qJD(5);
t253 = -t163 / 0.2e1 - t164 / 0.2e1;
t166 = t242 * qJD(5);
t167 = t243 * qJD(5);
t252 = -t166 / 0.2e1 - t167 / 0.2e1;
t251 = -t319 * t223 / 0.2e1 + t321 * t226 / 0.2e1;
t250 = t186 / 0.2e1 + t187 / 0.2e1;
t249 = t189 / 0.2e1 + t190 / 0.2e1;
t16 = -t54 * mrSges(7,1) + t55 * mrSges(7,2);
t14 = -t223 * t56 + t226 * t45;
t81 = t103 * mrSges(5,1) + mrSges(5,2) * t263;
t176 = t284 * t226;
t248 = qJ(6) * t227 - t151;
t247 = pkin(3) * t279 - qJD(4) * t224;
t246 = Ifges(7,5) * t259 + Ifges(7,6) * t233 + Ifges(7,3) * t278;
t245 = Ifges(6,5) * t259 + Ifges(6,6) * t233 + Ifges(6,3) * t278;
t155 = mrSges(7,1) * t276 - mrSges(7,2) * t277;
t244 = mrSges(6,1) * t226 - mrSges(6,2) * t223;
t179 = t227 * mrSges(5,2) - t224 * mrSges(5,3);
t237 = -pkin(3) * t227 - t289;
t236 = -t14 * t223 + t15 * t226;
t67 = qJ(4) * t287 - t74;
t235 = t320 * t102 + t322 * t103 + t263 * t323;
t34 = -t127 * t279 + t128 * t278 + t224 * t137 + t227 * t138;
t122 = (pkin(10) * t224 - qJ(4) * t227) * qJD(3) + t247;
t174 = t307 * t278;
t52 = t226 * t122 - t151 * t277 + t223 * t174 + t192 * t276;
t57 = -pkin(4) * t143 - t67;
t24 = -qJ(4) * t263 + qJD(4) * t287 - t34;
t21 = -pkin(4) * t102 - t24;
t92 = -mrSges(7,1) * t233 + mrSges(7,2) * t232;
t217 = Ifges(4,5) * t278;
t216 = Ifges(5,5) * t279;
t212 = pkin(5) * t223 + qJ(4);
t206 = pkin(5) * t276 + qJD(4);
t197 = Ifges(3,5) * t262;
t191 = Ifges(4,1) * t224 + t300;
t188 = Ifges(4,2) * t227 + t301;
t183 = -Ifges(5,2) * t224 - t294;
t182 = -Ifges(5,3) * t227 - t295;
t181 = mrSges(6,1) * t223 + mrSges(6,2) * t226;
t180 = mrSges(7,1) * t223 + mrSges(7,2) * t226;
t177 = -pkin(2) + t237;
t175 = t284 * t223;
t173 = t307 * t279;
t172 = -mrSges(6,2) * t224 - mrSges(6,3) * t285;
t171 = -mrSges(7,2) * t224 - mrSges(7,3) * t285;
t170 = mrSges(6,1) * t224 + mrSges(6,3) * t286;
t169 = mrSges(7,1) * t224 + mrSges(7,3) * t286;
t168 = (Ifges(4,1) * t227 - t301) * qJD(3);
t165 = (-Ifges(4,2) * t224 + t300) * qJD(3);
t160 = (-Ifges(5,2) * t227 + t295) * qJD(3);
t159 = (Ifges(5,3) * t224 - t294) * qJD(3);
t158 = (mrSges(4,1) * t224 + mrSges(4,2) * t227) * qJD(3);
t157 = (-mrSges(5,2) * t224 - mrSges(5,3) * t227) * qJD(3);
t156 = t244 * qJD(5);
t154 = t226 * t192;
t150 = t244 * t227;
t149 = (mrSges(7,1) * t226 - mrSges(7,2) * t223) * t227;
t146 = t226 * t174;
t142 = pkin(5) * t285 + t193;
t141 = -qJ(4) * t278 + t247;
t136 = -qJD(5) * t176 - t223 * qJD(6);
t135 = t277 * t284 - t274;
t130 = Ifges(6,3) * t224 + t227 * t239;
t129 = Ifges(7,3) * t224 + t227 * t238;
t118 = mrSges(6,1) * t278 - mrSges(6,3) * t232;
t117 = mrSges(7,1) * t278 - mrSges(7,3) * t232;
t116 = -mrSges(6,2) * t278 + mrSges(6,3) * t233;
t115 = -mrSges(7,2) * t278 + mrSges(7,3) * t233;
t111 = -mrSges(4,1) * t287 - mrSges(4,3) * t144;
t110 = mrSges(4,2) * t287 - mrSges(4,3) * t143;
t109 = mrSges(5,1) * t144 - mrSges(5,2) * t287;
t108 = mrSges(5,1) * t143 + mrSges(5,3) * t287;
t106 = -pkin(5) * t261 + (-pkin(5) * t226 - t307) * t279;
t94 = -t151 * t223 + t154;
t93 = -mrSges(6,1) * t233 + mrSges(6,2) * t232;
t91 = -mrSges(5,2) * t143 - mrSges(5,3) * t144;
t90 = -qJ(6) * t285 + t95;
t85 = -Ifges(6,5) * t260 + t245;
t84 = -Ifges(7,5) * t260 + t246;
t83 = mrSges(4,1) * t263 - mrSges(4,3) * t103;
t82 = -mrSges(4,2) * t263 - mrSges(4,3) * t102;
t80 = mrSges(5,1) * t102 - mrSges(5,3) * t263;
t79 = pkin(5) * t224 + t223 * t248 + t154;
t78 = Ifges(4,1) * t144 - Ifges(4,4) * t143 - Ifges(4,5) * t287;
t77 = Ifges(4,4) * t144 - Ifges(4,2) * t143 - Ifges(4,6) * t287;
t76 = -Ifges(5,4) * t287 - Ifges(5,2) * t144 + Ifges(5,6) * t143;
t75 = -Ifges(5,5) * t287 - Ifges(5,6) * t144 + Ifges(5,3) * t143;
t72 = mrSges(6,1) * t144 + mrSges(6,3) * t234;
t71 = mrSges(7,1) * t144 + mrSges(7,3) * t234;
t70 = -mrSges(6,2) * t144 + mrSges(6,3) * t104;
t69 = -mrSges(7,2) * t144 + mrSges(7,3) * t104;
t66 = pkin(3) * t143 + t231;
t65 = -mrSges(6,1) * t104 - mrSges(6,2) * t234;
t64 = -mrSges(7,1) * t104 - mrSges(7,2) * t234;
t63 = -mrSges(5,2) * t102 - mrSges(5,3) * t103;
t62 = mrSges(4,1) * t102 + mrSges(4,2) * t103;
t61 = Ifges(4,1) * t103 - Ifges(4,4) * t102 + Ifges(4,5) * t263;
t60 = Ifges(4,4) * t103 - Ifges(4,2) * t102 + Ifges(4,6) * t263;
t59 = Ifges(5,4) * t263 - Ifges(5,2) * t103 + Ifges(5,6) * t102;
t58 = Ifges(5,5) * t263 - Ifges(5,6) * t103 + Ifges(5,3) * t102;
t53 = -qJD(5) * t95 - t122 * t223 + t146;
t40 = -Ifges(6,5) * t234 + Ifges(6,6) * t104 + Ifges(6,3) * t144;
t39 = -Ifges(7,5) * t234 + Ifges(7,6) * t104 + Ifges(7,3) * t144;
t33 = pkin(3) * t102 + t230;
t32 = qJ(6) * t233 - t227 * t274 + t52;
t31 = -pkin(3) * t263 - t35;
t30 = -pkin(5) * t104 + t57;
t29 = pkin(5) * t278 + t146 + t248 * t276 + (-qJ(6) * t279 - qJD(5) * t192 + qJD(6) * t227 - t122) * t223;
t28 = mrSges(6,1) * t103 - mrSges(6,3) * t55;
t27 = mrSges(7,1) * t103 - mrSges(7,3) * t55;
t26 = -mrSges(6,2) * t103 + mrSges(6,3) * t54;
t25 = -mrSges(7,2) * t103 + mrSges(7,3) * t54;
t17 = -mrSges(6,1) * t54 + mrSges(6,2) * t55;
t7 = qJ(6) * t104 + t15;
t6 = pkin(5) * t144 + qJ(6) * t234 + t14;
t5 = -pkin(5) * t54 + t21;
t2 = qJ(6) * t54 + qJD(6) * t104 + t3;
t1 = pkin(5) * t103 - qJ(6) * t55 + qJD(6) * t234 + t4;
t18 = [(t75 - t77) * t102 + 0.2e1 * m(3) * (t138 * t148 - t139 * t147) + 0.2e1 * m(4) * (t126 * t139 + t34 * t74 + t35 * t73) + 0.2e1 * m(5) * (t24 * t67 + t31 * t68 + t33 * t66) + (t43 + t44) * t55 + (t41 + t42) * t54 + (t10 + t11) * t104 + 0.2e1 * t14 * t28 + 0.2e1 * t30 * t16 + 0.2e1 * t7 * t25 + 0.2e1 * t15 * t26 + 0.2e1 * t6 * t27 + (t139 * t225 * t310 + (t138 * t310 - t235) * t228 + ((t147 * t311 + Ifges(3,5) * t222 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t228) * t273) * t228 + (t148 * t311 - 0.2e1 * Ifges(3,6) * t222 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t225) * t273 + t322 * t144 + t320 * t143 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - t323) * t287) * t225) * qJD(2)) * t221 + (t1 * t6 + t2 * t7 + t30 * t5) * t312 + (t14 * t4 + t15 * t3 + t21 * t57) * t313 - (t12 + t13) * t234 + 0.2e1 * t33 * t91 + 0.2e1 * t24 * t108 + 0.2e1 * t31 * t109 + 0.2e1 * t34 * t110 + 0.2e1 * t35 * t111 + 0.2e1 * t126 * t62 + 0.2e1 * t5 * t64 + 0.2e1 * t21 * t65 + 0.2e1 * t66 * t63 + 0.2e1 * t2 * t69 + 0.2e1 * t3 * t70 + 0.2e1 * t1 * t71 + 0.2e1 * t4 * t72 + 0.2e1 * t67 * t80 + 0.2e1 * t68 * t81 + 0.2e1 * t74 * t82 + 0.2e1 * t73 * t83 + (-t59 + t61 + 0.2e1 * t290 + t318) * t144 + 0.2e1 * t57 * t17 + (t58 - t60 + 0.2e1 * t291) * t143 + (t197 - 0.2e1 * t292 - 0.2e1 * t293) * t222 + (t39 + t40 + t78 - t76) * t103; m(6) * (t14 * t53 + t15 * t52 - t173 * t57 + t193 * t21 + t3 * t95 + t4 * t94) - t292 - t293 + (t84 / 0.2e1 + t85 / 0.2e1 - t160 / 0.2e1 + t168 / 0.2e1) * t144 + (t159 / 0.2e1 - t165 / 0.2e1) * t143 + (t129 / 0.2e1 + t130 / 0.2e1 - t183 / 0.2e1 + t191 / 0.2e1) * t103 + (t182 / 0.2e1 - t188 / 0.2e1) * t102 + t197 + m(7) * (t1 * t79 + t106 * t30 + t142 * t5 + t2 * t90 + t29 * t6 + t32 * t7) - t266 * t234 + t90 * t25 + t30 * t92 + t57 * t93 + t94 * t28 + t95 * t26 + t106 * t64 + t7 * t115 + t15 * t116 + t6 * t117 + t14 * t118 - pkin(2) * t62 + t32 * t69 + t52 * t70 + t29 * t71 + t53 * t72 + t79 * t27 + t141 * t91 + t142 * t16 + t5 * t149 + t21 * t150 + t66 * t157 + t126 * t158 + t1 * t169 + t4 * t170 + t2 * t171 + t3 * t172 - t173 * t65 + t177 * t63 + t33 * t179 + t193 * t17 + t255 * t55 + t256 * t54 + t267 * t104 + ((Ifges(5,4) * t257 - t76 / 0.2e1 + t78 / 0.2e1 + t39 / 0.2e1 + t40 / 0.2e1 + t68 * mrSges(5,1) - t73 * mrSges(4,3)) * t227 + (Ifges(4,6) * t257 + t75 / 0.2e1 - t77 / 0.2e1 + t67 * mrSges(5,1) - t74 * mrSges(4,3) + t269 * t226 + t268 * t223) * t224 + ((t109 - t111) * t227 + (t108 - t110) * t224 + m(5) * (t224 * t67 + t227 * t68) + m(4) * (-t224 * t74 - t227 * t73)) * pkin(9)) * qJD(3) + (-Ifges(3,6) * t280 + (-t217 / 0.2e1 - t216 / 0.2e1) * t228) * t221 + (-t35 * mrSges(4,3) + t31 * mrSges(5,1) + t9 / 0.2e1 + t8 / 0.2e1 - t59 / 0.2e1 + t61 / 0.2e1 + t290 + (Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1) * t263 + (t81 - t83) * pkin(9)) * t224 + (t34 * mrSges(4,3) - t24 * mrSges(5,1) - t58 / 0.2e1 + t60 / 0.2e1 - t291 + t271 * t226 + t270 * t223 + (Ifges(4,6) / 0.2e1 - Ifges(5,5) / 0.2e1) * t263 + (-t80 + t82) * pkin(9) + (t223 * t269 - t226 * t268) * qJD(5)) * t227 + m(4) * (-pkin(2) * t139 + t34 * t219 - t35 * t305) + m(5) * (t141 * t66 + t177 * t33 - t24 * t219 + t31 * t305); -0.2e1 * pkin(2) * t158 + 0.2e1 * t106 * t149 + 0.2e1 * t90 * t115 + 0.2e1 * t95 * t116 + 0.2e1 * t79 * t117 + 0.2e1 * t94 * t118 + 0.2e1 * t142 * t92 - 0.2e1 * t173 * t150 + 0.2e1 * t177 * t157 + 0.2e1 * t29 * t169 + 0.2e1 * t53 * t170 + 0.2e1 * t32 * t171 + 0.2e1 * t52 * t172 + 0.2e1 * t193 * t93 + 0.2e1 * (m(5) * t177 + t179) * t141 + (-t173 * t193 + t52 * t95 + t53 * t94) * t313 + (t106 * t142 + t29 * t79 + t32 * t90) * t312 + (-t160 + t168 + t84 + t85 + (t223 * t282 + t226 * t283 + t182 - t188) * qJD(3)) * t224 + (-t159 + t165 + (-t86 - t87) * t226 + (-t88 - t89) * t223 + (t223 * t283 - t226 * t282) * qJD(5) + (t129 + t130 - t183 + t191) * qJD(3)) * t227; (-t108 + t65) * qJD(4) + (t17 - t80) * qJ(4) + t235 + t31 * mrSges(5,2) - t34 * mrSges(4,2) + t35 * mrSges(4,1) - t24 * mrSges(5,3) + m(5) * (-pkin(3) * t31 - qJ(4) * t24 - qJD(4) * t67) - t252 * t234 + m(6) * (qJ(4) * t21 + qJD(4) * t57 - t308 * t314) + ((-t15 * mrSges(6,3) - t7 * mrSges(7,3) - t269) * t226 + (t14 * mrSges(6,3) + t6 * mrSges(7,3) - t268) * t223 - (m(6) * t236 - t223 * t72 + t226 * t70) * t308) * qJD(5) + (-t4 * mrSges(6,3) - t1 * mrSges(7,3) - t28 * t308 - t270) * t226 + (-t3 * mrSges(6,3) - t2 * mrSges(7,3) - t26 * t308 + t271) * t223 + m(7) * (-t1 * t176 + t135 * t6 + t136 * t7 - t175 * t2 + t206 * t30 + t212 * t5) - pkin(3) * t81 + t135 * t71 + t136 * t69 + t30 * t155 + t57 * t156 - t175 * t25 - t176 * t27 + t5 * t180 + t21 * t181 + t206 * t64 + t212 * t16 + t249 * t55 + t250 * t54 + t251 * t103 + t253 * t104 + t254 * t144; qJ(4) * t93 + t106 * t180 - t175 * t115 - t176 * t117 + t135 * t169 + t136 * t171 + t142 * t155 + t206 * t149 + t193 * t156 + t212 * t92 + t216 + t217 + t254 * t224 + m(7) * (t106 * t212 + t135 * t79 + t136 * t90 + t142 * t206 - t175 * t32 - t176 * t29) + ((-mrSges(5,1) * qJ(4) - Ifges(4,6)) * t224 + (-pkin(3) * mrSges(5,1) - Ifges(5,4) + t251) * t227 + (m(5) * t237 - t227 * mrSges(4,1) + t224 * mrSges(4,2) + t179) * pkin(9)) * qJD(3) + (-t29 * mrSges(7,3) - t53 * mrSges(6,3) - (m(6) * t53 + t118) * t308 - t253 * t227 + t250 * t279 + (-t95 * mrSges(6,3) - t90 * mrSges(7,3) - (m(6) * t95 + t172) * t308 - t249 * t227 - t256) * qJD(5) + t266) * t226 + (-t52 * mrSges(6,3) - t32 * mrSges(7,3) - (m(6) * t52 + t116) * t308 - t252 * t227 + t249 * t279 + (t94 * mrSges(6,3) + t79 * mrSges(7,3) - (-m(6) * t94 - t170) * t308 + t250 * t227 - t255) * qJD(5) - t267) * t223 - (m(6) * qJ(4) + t181) * t173 + (m(5) * t219 + m(6) * t193 + t227 * mrSges(5,1) + t150) * qJD(4); 0.2e1 * t206 * t180 + 0.2e1 * t212 * t155 + (-t135 * t176 - t136 * t175 + t206 * t212) * t312 + 0.2e1 * qJ(4) * t156 + (t135 * t309 - t166 - t167) * t226 + (t136 * t309 + t163 + t164) * t223 + 0.2e1 * (mrSges(5,3) + t181 + (m(5) + m(6)) * qJ(4)) * qJD(4) + ((-t175 * t309 - t186 - t187) * t226 + (t176 * t309 - t189 - t190) * t223) * qJD(5); (t27 + t28) * t226 + (t25 + t26) * t223 + ((t69 + t70) * t226 + (-t71 - t72) * t223) * qJD(5) + m(7) * (t1 * t226 + t2 * t223 + (-t223 * t6 + t226 * t7) * qJD(5)) + m(6) * (qJD(5) * t236 + t314) + m(5) * t31 + t81; (t117 + t118) * t226 + (t115 + t116) * t223 + (m(5) * pkin(9) + mrSges(5,1)) * t278 + ((t171 + t172) * t226 + (-t169 - t170) * t223) * qJD(5) + m(7) * (t223 * t32 + t226 * t29 + (-t223 * t79 + t226 * t90) * qJD(5)) + m(6) * (t223 * t52 + t226 * t53 + (-t223 * t94 + t226 * t95) * qJD(5)); m(7) * (t135 * t226 + t136 * t223 + (-t175 * t226 + t176 * t223) * qJD(5)); 0; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 + (m(7) * t1 + t27) * pkin(5) + t318; mrSges(6,1) * t53 + mrSges(7,1) * t29 - mrSges(6,2) * t52 - mrSges(7,2) * t32 - t321 * t260 + (m(7) * t29 + t117) * pkin(5) + t245 + t246; -t136 * mrSges(7,2) + t265 * t135 + ((mrSges(6,2) * t308 - t319) * t226 + (mrSges(6,1) * t308 + (mrSges(7,3) * pkin(5)) - t321) * t223) * qJD(5); ((-mrSges(6,2) - mrSges(7,2)) * t226 + (-mrSges(6,1) - t265) * t223) * qJD(5); 0; m(7) * t5 + t16; m(7) * t106 + t92; m(7) * t206 + t155; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t18(1) t18(2) t18(4) t18(7) t18(11) t18(16); t18(2) t18(3) t18(5) t18(8) t18(12) t18(17); t18(4) t18(5) t18(6) t18(9) t18(13) t18(18); t18(7) t18(8) t18(9) t18(10) t18(14) t18(19); t18(11) t18(12) t18(13) t18(14) t18(15) t18(20); t18(16) t18(17) t18(18) t18(19) t18(20) t18(21);];
Mq  = res;
