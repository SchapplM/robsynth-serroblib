% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:30
% EndTime: 2019-12-31 18:53:36
% DurationCPUTime: 3.11s
% Computational Cost: add. (6249->328), mult. (12840->449), div. (0->0), fcn. (13251->6), ass. (0->173)
t332 = Ifges(5,5) + Ifges(6,4);
t335 = Ifges(5,6) - Ifges(6,6);
t191 = sin(qJ(4));
t187 = t191 ^ 2;
t193 = cos(qJ(4));
t188 = t193 ^ 2;
t334 = t187 + t188;
t333 = -mrSges(5,1) - mrSges(6,1);
t331 = Ifges(6,2) + Ifges(5,3);
t189 = sin(pkin(8));
t190 = cos(pkin(8));
t192 = sin(qJ(3));
t300 = cos(qJ(3));
t167 = t189 * t192 - t300 * t190;
t245 = -pkin(2) * t190 - pkin(1);
t169 = t300 * t189 + t192 * t190;
t296 = pkin(7) * t169;
t110 = pkin(3) * t167 + t245 - t296;
t294 = pkin(6) + qJ(2);
t171 = t294 * t190;
t235 = t294 * t189;
t125 = t300 * t171 - t192 * t235;
t49 = t110 * t193 - t191 * t125;
t36 = -t167 * pkin(4) - t49;
t292 = t36 + t49;
t266 = t167 * t191;
t112 = -mrSges(5,2) * t169 + mrSges(5,3) * t266;
t119 = mrSges(6,2) * t266 + mrSges(6,3) * t169;
t330 = t112 + t119;
t264 = t169 * t191;
t255 = mrSges(5,3) * t264;
t113 = -t167 * mrSges(5,2) - t255;
t253 = mrSges(6,2) * t264;
t283 = t167 * mrSges(6,3);
t118 = -t253 + t283;
t329 = t113 + t118;
t263 = t169 * t193;
t251 = mrSges(5,3) * t263;
t116 = t167 * mrSges(5,1) - t251;
t252 = mrSges(6,2) * t263;
t117 = -t167 * mrSges(6,1) + t252;
t328 = t116 - t117;
t183 = Ifges(5,4) * t193;
t179 = Ifges(5,1) * t191 + t183;
t325 = -Ifges(5,2) * t191 + t183;
t327 = t179 + t325;
t182 = Ifges(6,5) * t191;
t326 = Ifges(6,1) * t193 + t182;
t289 = Ifges(5,4) * t191;
t224 = Ifges(5,1) * t193 - t289;
t207 = t224 * t169;
t324 = t332 * t167 + t169 * t326 + t207;
t323 = -t335 * t191 + t332 * t193;
t322 = -t334 * t296 / 0.2e1;
t215 = Ifges(6,3) * t193 - t182;
t321 = t326 - t215;
t177 = Ifges(5,2) * t193 + t289;
t288 = Ifges(6,5) * t193;
t178 = Ifges(6,1) * t191 - t288;
t302 = t193 / 0.2e1;
t306 = -t191 / 0.2e1;
t320 = t177 * t306 + t178 * t302;
t297 = pkin(7) * t167;
t298 = pkin(3) * t169;
t122 = t297 + t298;
t124 = t192 * t171 + t300 * t235;
t57 = t122 * t193 + t191 * t124;
t58 = t191 * t122 - t193 * t124;
t319 = -t191 * t57 + t193 * t58;
t40 = qJ(5) * t169 + t58;
t41 = -t169 * pkin(4) - t57;
t318 = t191 * t41 + t193 * t40;
t247 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t181 = m(6) * qJ(5) + mrSges(6,3);
t309 = t169 / 0.2e1;
t317 = (t332 * t191 + t335 * t193) * t309;
t316 = t169 ^ 2;
t315 = m(5) / 0.2e1;
t314 = -m(6) / 0.2e1;
t313 = m(6) / 0.2e1;
t312 = -mrSges(6,1) / 0.2e1;
t121 = t193 * t125;
t88 = t191 * t110;
t50 = t121 + t88;
t311 = t50 / 0.2e1;
t226 = t193 * mrSges(6,1) + t191 * mrSges(6,3);
t308 = -t226 / 0.2e1;
t307 = -t178 / 0.2e1;
t305 = t191 / 0.2e1;
t303 = -t193 / 0.2e1;
t213 = t191 * pkin(4) - qJ(5) * t193;
t299 = m(6) * t213;
t267 = t167 * qJ(5);
t35 = t267 + t50;
t293 = t35 - t50;
t275 = t193 * mrSges(6,3);
t280 = t191 * mrSges(6,1);
t225 = -t275 + t280;
t108 = t225 * t169;
t276 = t193 * mrSges(5,2);
t281 = t191 * mrSges(5,1);
t227 = t276 + t281;
t109 = t227 * t169;
t268 = t124 * t169;
t61 = t213 * t169 + t124;
t5 = (mrSges(4,3) * t169 + t108 + t109) * t169 + (mrSges(4,3) * t167 + t328 * t191 - t329 * t193) * t167 + m(6) * (t61 * t169 + (-t191 * t36 - t193 * t35) * t167) + m(5) * (t268 + (t191 * t49 - t193 * t50) * t167) + m(4) * (-t125 * t167 + t268) + (m(3) * qJ(2) + mrSges(3,3)) * (t189 ^ 2 + t190 ^ 2);
t285 = qJD(1) * t5;
t106 = t225 * t167;
t107 = t227 * t167;
t265 = t167 * t193;
t114 = t169 * mrSges(5,1) + mrSges(5,3) * t265;
t254 = mrSges(6,2) * t265;
t282 = t169 * mrSges(6,1);
t115 = -t254 - t282;
t164 = t167 * mrSges(4,2);
t248 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t62 = -t213 * t167 + t125;
t216 = Ifges(6,3) * t191 + t288;
t74 = Ifges(6,6) * t169 - t216 * t167;
t75 = Ifges(5,6) * t169 - t167 * t325;
t76 = Ifges(6,4) * t169 - t167 * t326;
t77 = Ifges(5,5) * t169 - t224 * t167;
t1 = -t61 * t106 + t62 * t108 + t50 * t112 + t58 * t113 + t49 * t114 + t36 * t115 + t57 * t116 + t41 * t117 + t40 * t118 + t35 * t119 - t124 * t107 + t125 * t109 - t245 * t164 + m(6) * (t35 * t40 + t36 * t41 + t61 * t62) + m(5) * (t124 * t125 + t49 * t57 + t50 * t58) + (t245 * mrSges(4,1) - Ifges(4,4) * t169 + (t76 / 0.2e1 + t77 / 0.2e1 + t248 * t169) * t193 + (t74 / 0.2e1 - t75 / 0.2e1 + t247 * t169) * t191) * t169 + ((-Ifges(4,1) + Ifges(4,2) + (-Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1) * t188 + ((-Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t191 + (Ifges(5,4) - Ifges(6,5)) * t193) * t191 + t331) * t169 + (Ifges(4,4) - t323) * t167) * t167;
t284 = t1 * qJD(1);
t206 = t216 * t169;
t201 = Ifges(6,6) * t167 + t206;
t202 = Ifges(5,6) * t167 + t169 * t325;
t214 = pkin(4) * t193 + t191 * qJ(5);
t205 = t214 * t169;
t208 = t226 * t169;
t228 = t193 * mrSges(5,1) - t191 * mrSges(5,2);
t209 = t228 * t169;
t241 = -t263 / 0.2e1;
t4 = t202 * t263 / 0.2e1 + t201 * t241 - t124 * t209 + t36 * t253 - t61 * t208 + t35 * t252 + t320 * t316 - (-t193 * t179 + t191 * t215) * t316 / 0.2e1 + t324 * t264 / 0.2e1 + (-m(6) * t61 - t108) * t205 + t167 * t317 + (-m(6) * t36 + t251 + t328) * t50 + (-m(6) * t35 - t255 - t329) * t49;
t272 = t4 * qJD(1);
t197 = (t292 * t191 + t293 * t193) * t314 + t116 * t305 + t117 * t306 + t329 * t303;
t198 = (t276 / 0.2e1 + t281 / 0.2e1 + t280 / 0.2e1 - t275 / 0.2e1 + t213 * t313) * t167;
t6 = t198 + t197;
t271 = t6 * qJD(1);
t261 = t334 * t297;
t170 = -pkin(3) - t214;
t262 = t170 * t169;
t195 = (-t261 - t298) * t315 + (-t261 + t262) * t313 + (mrSges(6,2) + mrSges(5,3)) * t167 * (-t187 / 0.2e1 - t188 / 0.2e1);
t200 = -m(5) * (t191 * t58 + t193 * t57) / 0.2e1 + (t191 * t40 - t193 * t41) * t314;
t237 = t308 - t228 / 0.2e1;
t8 = t164 + (t115 / 0.2e1 - t114 / 0.2e1) * t193 + (-t119 / 0.2e1 - t112 / 0.2e1) * t191 + (-mrSges(4,1) + t237) * t169 + t195 + t200;
t270 = t8 * qJD(1);
t16 = m(6) * (t35 * t167 - t61 * t263) - t108 * t263 + t167 * t118;
t269 = qJD(1) * t16;
t259 = qJD(4) * t191;
t103 = m(6) * t266;
t258 = t103 * qJD(1);
t256 = t41 * t313;
t230 = pkin(4) * mrSges(6,2) - t332;
t229 = -qJ(5) * mrSges(6,2) - t335;
t19 = -t213 * t226 + t216 * t303 - pkin(3) * t227 + (t225 + t299) * t170 + t327 * t302 + (t224 + t321) * t305 + t320;
t194 = t205 * t308 + (t170 * t205 + t213 * t61) * t313 + t61 * t225 / 0.2e1 + t124 * t227 / 0.2e1 - pkin(3) * t209 / 0.2e1 + t170 * t208 / 0.2e1 - t191 * t202 / 0.4e1 + t177 * t241 + t213 * t108 / 0.2e1 + t323 * t167 / 0.4e1 + (t206 + t201) * t191 / 0.4e1 + t322 * mrSges(5,3) + (t207 + t324) * t193 / 0.4e1 + (-t179 / 0.4e1 + t307 - t327 / 0.4e1) * t264 + (-t215 / 0.4e1 + t321 / 0.4e1) * t263 + ((-t293 * t191 + t292 * t193) * t313 + (t117 / 0.2e1 - t116 / 0.2e1) * t193 + t329 * t306) * pkin(7) + ((-t35 / 0.2e1 + t311) * t191 + t292 * t302 + t322) * mrSges(6,2);
t196 = (-pkin(4) * t41 + qJ(5) * t40) * t313 - pkin(4) * t115 / 0.2e1 + qJ(5) * t119 / 0.2e1 + t40 * mrSges(6,3) / 0.2e1 + t41 * t312 + t57 * mrSges(5,1) / 0.2e1 - t58 * mrSges(5,2) / 0.2e1;
t2 = -t194 + t196 + t331 * t309 - t247 * t266 - t332 * t265 / 0.2e1;
t212 = t2 * qJD(1) - t19 * qJD(3);
t123 = (m(6) * t170 - t226) * t191;
t199 = (-t191 * t61 + (-t262 + t297) * t193) * t313 + t108 * t306;
t13 = -t254 + (-t226 * t302 + t312) * t169 + t256 - t199;
t211 = qJD(1) * t13 + qJD(3) * t123;
t18 = t283 + 0.2e1 * (t267 / 0.2e1 + t121 / 0.4e1 + t88 / 0.4e1 - t50 / 0.4e1) * m(6);
t210 = qJD(1) * t18 + qJD(4) * t181;
t172 = (m(6) * pkin(7) + mrSges(6,2)) * t193;
t17 = (0.2e1 * t267 + t50) * t313 + m(6) * t311 + t118;
t14 = -t226 * t241 - t282 / 0.2e1 + t256 + t199;
t9 = t114 * t302 + t115 * t303 + t237 * t169 + t330 * t305 + t195 - t200;
t7 = t198 - t197;
t3 = t194 + (-t247 * t191 - t248 * t193) * t167 + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t169 + t196;
t10 = [qJD(2) * t5 + qJD(3) * t1 - qJD(4) * t4 + qJD(5) * t16, qJD(3) * t9 + qJD(4) * t7 + t285, t9 * qJD(2) + t3 * qJD(4) + t14 * qJD(5) + t284 + (-t125 * mrSges(4,1) + t124 * mrSges(4,2) - Ifges(4,6) * t169 + pkin(3) * t107 - t170 * t106 - t125 * t228 - t62 * t226 + t75 * t302 + t74 * t303 + (t330 * t193 + (-t114 + t115) * t191) * pkin(7) + 0.2e1 * (-pkin(3) * t125 + pkin(7) * t319) * t315 + 0.2e1 * (pkin(7) * t318 + t170 * t62) * t313 + (-Ifges(4,5) + (t307 - t179 / 0.2e1) * t193 + (t215 / 0.2e1 + t177 / 0.2e1) * t191) * t167 + t317 + (t76 + t77) * t305 + t319 * mrSges(5,3) + t318 * mrSges(6,2)) * qJD(3), t7 * qJD(2) + t3 * qJD(3) + t17 * qJD(5) - t272 + ((-m(6) * pkin(4) + t333) * t50 + (-mrSges(5,2) + t181) * t49 + (t191 * t230 + t193 * t229) * t169) * qJD(4), qJD(3) * t14 + qJD(4) * t17 + t269; -qJD(3) * t8 - qJD(4) * t6 + qJD(5) * t103 - t285, 0, -t270, -t271 + (t275 - t276 - t299) * qJD(4) + (m(6) * qJD(5) + t333 * qJD(4)) * t191, m(6) * t259 + t258; qJD(2) * t8 - qJD(4) * t2 - qJD(5) * t13 - t284, t270, qJD(4) * t19 - qJD(5) * t123, t172 * qJD(5) + t229 * t259 - t212 + (-t230 * t193 + (-m(6) * t214 - t226 - t228) * pkin(7)) * qJD(4), qJD(4) * t172 - t211; qJD(2) * t6 + qJD(3) * t2 + qJD(5) * t18 + t272, t271, t212, t181 * qJD(5), t210; -qJD(2) * t103 + qJD(3) * t13 - qJD(4) * t18 - t269, -t258, t211, -t210, 0;];
Cq = t10;
