% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPPR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:49
% EndTime: 2019-03-09 02:58:54
% DurationCPUTime: 2.04s
% Computational Cost: add. (3220->308), mult. (5741->397), div. (0->0), fcn. (4258->4), ass. (0->181)
t178 = sin(qJ(6));
t171 = t178 ^ 2;
t180 = cos(qJ(6));
t173 = t180 ^ 2;
t240 = t171 + t173;
t309 = mrSges(7,3) * t240;
t165 = t180 * mrSges(7,1);
t275 = t178 * mrSges(7,2);
t213 = t165 - t275;
t315 = qJD(6) * t213;
t175 = qJ(4) + pkin(5);
t272 = t180 * mrSges(7,2);
t276 = t178 * mrSges(7,1);
t212 = t272 + t276;
t314 = t175 * t212;
t179 = sin(qJ(3));
t311 = -t173 / 0.2e1;
t227 = t311 - t171 / 0.2e1;
t310 = mrSges(7,3) * t227;
t313 = t179 * t310;
t297 = pkin(3) + pkin(4);
t312 = t297 * t179;
t197 = t272 / 0.2e1 + t276 / 0.2e1;
t284 = m(7) * t179;
t308 = mrSges(6,1) + t213;
t170 = -pkin(8) - t297;
t307 = t170 * t179;
t306 = t175 * t179;
t169 = Ifges(7,6) * t178;
t280 = Ifges(7,5) * t180;
t305 = -t169 + t280;
t181 = cos(qJ(3));
t289 = t181 / 0.2e1;
t304 = t289 * t309;
t183 = -pkin(1) - pkin(7);
t245 = qJ(5) + t183;
t133 = t245 * t179;
t68 = t170 * t181 - t306;
t43 = t133 * t180 + t178 * t68;
t271 = t180 * t43;
t42 = -t133 * t178 + t180 * t68;
t274 = t178 * t42;
t204 = t271 - t274;
t303 = -t275 / 0.2e1 + t165 / 0.2e1;
t166 = t181 * mrSges(5,3);
t241 = t181 * mrSges(6,1) + t179 * mrSges(6,2);
t232 = t166 + t241;
t252 = t179 * t180;
t270 = t181 * mrSges(7,1);
t129 = -mrSges(7,3) * t252 + t270;
t250 = t180 * t129;
t256 = t178 * t179;
t268 = t181 * mrSges(7,2);
t128 = -mrSges(7,3) * t256 - t268;
t258 = t178 * t128;
t163 = t181 * qJ(4);
t242 = t163 - qJ(2);
t107 = -t242 + t312;
t286 = m(6) * t107;
t225 = -pkin(3) * t179 + t163;
t137 = qJ(2) - t225;
t288 = m(5) * t137;
t134 = t245 * t181;
t65 = pkin(5) * t181 + t242 + t307;
t40 = t134 * t178 + t180 * t65;
t41 = -t134 * t180 + t178 * t65;
t302 = -m(7) * (t178 * t41 + t180 * t40) - t232 - t250 - t258 + t286 + t288;
t172 = t179 ^ 2;
t301 = 2 * qJD(3);
t300 = m(6) / 0.2e1;
t299 = m(7) / 0.2e1;
t298 = -mrSges(7,1) / 0.2e1;
t296 = (0.1e1 - t240) * t181 * t284;
t295 = -t170 / 0.2e1;
t294 = -t178 / 0.2e1;
t293 = t178 / 0.2e1;
t292 = -t179 / 0.2e1;
t291 = -t180 / 0.2e1;
t290 = t180 / 0.2e1;
t287 = m(6) * qJ(4);
t285 = m(7) * t175;
t283 = mrSges(5,1) + mrSges(4,1);
t282 = Ifges(7,4) * t178;
t281 = Ifges(7,4) * t180;
t279 = Ifges(7,5) * t181;
t278 = Ifges(7,2) * t178;
t277 = Ifges(7,6) * t181;
t273 = t179 * mrSges(5,1);
t269 = t181 * mrSges(4,2);
t108 = t213 * t179;
t140 = -Ifges(7,2) * t180 - t282;
t109 = t179 * t140;
t142 = -Ifges(7,1) * t178 - t281;
t110 = t179 * t142;
t205 = t178 * t40 - t180 * t41;
t209 = -Ifges(7,5) * t178 - Ifges(7,6) * t180;
t210 = -t278 + t281;
t87 = t210 * t179 + t277;
t211 = Ifges(7,1) * t180 - t282;
t88 = t211 * t179 + t279;
t4 = t133 * t108 + t40 * t128 - t41 * t129 + (t205 * mrSges(7,3) + t110 * t290 + t209 * t289 + t87 * t291 + (t109 + t88) * t294) * t179;
t267 = t4 * qJD(1);
t255 = t178 * t181;
t198 = t179 * mrSges(7,2) - mrSges(7,3) * t255;
t249 = t180 * t181;
t199 = -t179 * mrSges(7,1) - mrSges(7,3) * t249;
t185 = (t178 * t43 + t180 * t42) * t299 + t198 * t293 + t199 * t290;
t224 = -t179 * mrSges(6,1) + t181 * mrSges(6,2);
t226 = m(7) * t240;
t231 = t213 * t292;
t7 = -t231 + 0.2e1 * (t285 / 0.4e1 + t287 / 0.2e1) * t179 + 0.2e1 * ((t297 / 0.4e1 + pkin(3) / 0.4e1 + pkin(4) / 0.4e1) * m(6) - t170 * t226 / 0.4e1) * t181 - t185 - t224 + t304;
t266 = t7 * qJD(1);
t265 = qJ(4) * t179;
t174 = t181 ^ 2;
t239 = t174 + t172;
t251 = t180 * t128;
t257 = t178 * t129;
t261 = t133 * t179;
t10 = (-t251 + t257) * t181 + t212 * t172 + t239 * mrSges(6,3) + m(7) * (t205 * t181 + t261) + m(6) * (t134 * t181 + t261);
t264 = qJD(1) * t10;
t11 = (-t273 - t302) * t181;
t263 = qJD(1) * t11;
t13 = t269 + mrSges(3,3) + t283 * t179 + (m(4) + m(3)) * qJ(2) + t302;
t262 = qJD(1) * t13;
t193 = -t258 / 0.2e1 - t250 / 0.2e1;
t254 = t179 * t108;
t14 = -t254 / 0.2e1 + (t193 + t313) * t181 - t303;
t260 = t14 * qJD(1);
t214 = t268 / 0.2e1 - t128 / 0.2e1;
t215 = t270 / 0.2e1 + t129 / 0.2e1;
t17 = -t214 * t178 + t215 * t180 - t313;
t259 = t17 * qJD(1);
t253 = t212 * t179;
t20 = t215 * t178 + t214 * t180;
t248 = t20 * qJD(1);
t187 = (t240 * t174 + t172) * t299 + t239 * t300;
t202 = t300 + t226 / 0.2e1;
t37 = t187 + t202;
t247 = t37 * qJD(1);
t44 = (-0.2e1 * t227 * m(7) + m(6)) * t181;
t246 = t44 * qJD(1);
t243 = t256 * t298 - mrSges(7,2) * t252 / 0.2e1;
t238 = qJD(3) * t179;
t237 = Ifges(4,4) + Ifges(6,4) - Ifges(5,5);
t230 = t249 / 0.2e1;
t229 = -t142 / 0.4e1 + t210 / 0.4e1;
t228 = -t211 / 0.4e1 - t140 / 0.4e1;
t223 = t240 * t181;
t220 = m(5) * t183 - mrSges(5,2) + mrSges(6,3);
t219 = t240 * t299;
t191 = t181 * t212;
t1 = -(-t297 * t181 - t265) * t241 + t107 * t224 - t133 * t191 - t43 * t128 - t42 * t129 - t41 * t198 - t40 * t199 - m(7) * (t133 * t134 + t40 * t42 + t41 * t43) + (t166 - t273 - t288) * (pkin(3) * t181 + t265) + (-qJ(2) * mrSges(4,1) - t137 * mrSges(5,1) - t297 * t286 + t88 * t291 + t87 * t293 + (-t280 / 0.2e1 + t169 / 0.2e1 + t237) * t181) * t181 + (qJ(2) * mrSges(4,2) - t137 * mrSges(5,3) - t134 * t212 - qJ(4) * t286 + (Ifges(7,1) * t311 + Ifges(4,1) + Ifges(5,1) - Ifges(6,1) - Ifges(4,2) + Ifges(6,2) - Ifges(5,3) + Ifges(7,3) + (t281 - t278 / 0.2e1) * t178) * t181 + (-t237 + t305) * t179) * t179;
t192 = -t257 / 0.2e1 + t251 / 0.2e1;
t6 = m(7) * (t133 - t204) * t289 + (t191 / 0.2e1 + (t134 - t205) * t299 + t192) * t179;
t207 = -t1 * qJD(1) + t6 * qJD(2);
t203 = t170 * t310;
t201 = t42 * t298 + t43 * mrSges(7,2) / 0.2e1;
t200 = t6 * qJD(1) + qJD(2) * t296;
t195 = -t109 / 0.4e1 - t88 / 0.4e1 + t129 * t295;
t194 = -t110 / 0.4e1 + t87 / 0.4e1 + t128 * t295;
t29 = -t314 + (-t142 / 0.2e1 + t210 / 0.2e1) * t180 + (t211 / 0.2e1 + t140 / 0.2e1) * t178;
t188 = -t133 * t212 / 0.2e1 + t175 * t108 / 0.2e1 + t181 * t169 / 0.4e1;
t3 = (Ifges(7,3) / 0.2e1 + t203) * t179 + (t277 / 0.2e1 + t229 * t179 + t194) * t178 + (-0.3e1 / 0.4e1 * t279 + t228 * t179 + t195) * t180 + t188 + t201;
t49 = t253 / 0.2e1 + t243;
t190 = t3 * qJD(1) - t49 * qJD(2) + t29 * qJD(3);
t19 = 0.2e1 * (t133 / 0.4e1 + t274 / 0.4e1 - t271 / 0.4e1) * m(7);
t58 = -m(5) * qJ(4) - mrSges(5,3) - t285 - t287 - t308;
t67 = (-0.1e1 / 0.2e1 - t227) * t284;
t189 = qJD(1) * t19 - qJD(2) * t67 - qJD(3) * t58;
t50 = -t253 / 0.2e1 + t243;
t47 = t284 / 0.2e1 + (m(6) + m(5) + t219) * t179;
t45 = t181 * t219 - t223 * t299;
t36 = t187 - t202;
t21 = t197 * t181 + t192;
t18 = t303 * t181 + t292 * t309 + t193;
t16 = t204 * t299 + t198 * t290 + t199 * t294 + 0.2e1 * (m(7) / 0.4e1 + t300) * t133 + (t220 + t197) * t179;
t15 = t254 / 0.2e1 + t128 * t255 / 0.2e1 + t129 * t230 - t303 + t304 * t179;
t9 = -t231 - t181 * t310 + (-t170 * t223 + t306) * t299 + t185;
t5 = t6 * qJD(3);
t2 = Ifges(7,5) * t230 - Ifges(7,6) * t255 / 0.2e1 + Ifges(7,3) * t292 + (-t279 / 0.4e1 + t195) * t180 + t194 * t178 + (t229 * t178 + t228 * t180 + t203) * t179 + t188 - t201;
t8 = [qJD(2) * t13 - qJD(3) * t1 + qJD(4) * t11 + qJD(5) * t10 + qJD(6) * t4, qJD(5) * t36 + qJD(6) * t15 + t262 + t5, t16 * qJD(4) + t9 * qJD(5) + t2 * qJD(6) + ((t134 * t175 + t170 * t204) * t299 + (qJ(4) * t134 - t133 * t297) * t300) * t301 + (pkin(3) * mrSges(5,2) - t297 * mrSges(6,3) - Ifges(5,4) - Ifges(4,5) - Ifges(6,6) + t212 * t170 + (-m(5) * pkin(3) - t283) * t183 - t209) * t238 + t207 + (t133 * mrSges(6,2) + (t142 * t290 + Ifges(5,6) - Ifges(4,6) - Ifges(6,5) + t314 + t210 * t291 + (-mrSges(4,2) + mrSges(5,3)) * t183 + t220 * qJ(4) + (t140 + t211) * t294) * t181 + t308 * t134 - t204 * mrSges(7,3)) * qJD(3), qJD(3) * t16 + qJD(5) * t45 + qJD(6) * t18 + t263, qJD(2) * t36 + qJD(3) * t9 + qJD(4) * t45 + qJD(6) * t21 + t264, t267 + t15 * qJD(2) + t2 * qJD(3) + t18 * qJD(4) + t21 * qJD(5) + (-mrSges(7,1) * t41 - mrSges(7,2) * t40 + t179 * t209) * qJD(6); qJD(5) * t37 - qJD(6) * t14 - t262 + t5, qJD(3) * t296 (t181 * t213 + t232 - t269) * qJD(3) + t47 * qJD(4) + t50 * qJD(6) + (-t283 - t309) * t238 + ((t175 * t181 + t240 * t307) * t299 + m(5) * t225 / 0.2e1 + (t163 - t312) * t300) * t301 + t200, t47 * qJD(3), t247, t50 * qJD(3) + t181 * t315 - t260; qJD(4) * t19 + qJD(5) * t7 + qJD(6) * t3 - t207, -t67 * qJD(4) - t49 * qJD(6) - t200, -qJD(4) * t58 + qJD(6) * t29, t189, t266 (-t170 * t213 - t305) * qJD(6) + t190; -qJD(3) * t19 - qJD(5) * t44 - qJD(6) * t17 - t263, t67 * qJD(3), -t189, 0, -t246, -t259 - t315; -qJD(2) * t37 - qJD(3) * t7 + qJD(4) * t44 - qJD(6) * t20 - t264, -t247, -t266, t246, 0, -qJD(6) * t212 - t248; qJD(2) * t14 - qJD(3) * t3 + qJD(4) * t17 + qJD(5) * t20 - t267, qJD(3) * t49 + t260, -t190, t259, t248, 0;];
Cq  = t8;
