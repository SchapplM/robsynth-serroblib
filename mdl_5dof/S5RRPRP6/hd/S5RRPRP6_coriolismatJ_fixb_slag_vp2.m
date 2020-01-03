% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:49
% EndTime: 2019-12-31 19:56:57
% DurationCPUTime: 3.98s
% Computational Cost: add. (6727->365), mult. (13571->504), div. (0->0), fcn. (14014->6), ass. (0->186)
t211 = cos(qJ(4));
t294 = Ifges(5,6) + Ifges(6,6);
t338 = t294 * t211;
t329 = Ifges(5,5) + Ifges(6,5);
t301 = sin(qJ(2));
t248 = t301 * pkin(2);
t337 = m(4) * t248;
t210 = sin(qJ(4));
t208 = t210 ^ 2;
t209 = t211 ^ 2;
t253 = t208 + t209;
t269 = sin(pkin(8));
t270 = cos(pkin(8));
t302 = cos(qJ(2));
t176 = t269 * t301 - t270 * t302;
t178 = -t269 * t302 - t270 * t301;
t336 = -t178 * mrSges(4,1) - t176 * mrSges(4,2);
t307 = -t210 / 0.2e1;
t303 = t211 / 0.2e1;
t205 = Ifges(6,4) * t211;
t225 = t210 * Ifges(6,2) - t205;
t74 = Ifges(6,6) * t176 + t178 * t225;
t206 = Ifges(5,4) * t211;
t226 = t210 * Ifges(5,2) - t206;
t76 = Ifges(5,6) * t176 + t178 * t226;
t335 = t74 + t76;
t334 = t329 * t210 + t338;
t236 = t269 * pkin(2);
t200 = t236 + pkin(7);
t256 = qJ(5) + t200;
t173 = t256 * t210;
t174 = t256 * t211;
t262 = t178 * t210;
t202 = -pkin(2) * t302 - pkin(1);
t117 = t176 * pkin(3) + t178 * pkin(7) + t202;
t247 = t301 * pkin(6);
t184 = -qJ(3) * t301 - t247;
t249 = t302 * pkin(6);
t187 = qJ(3) * t302 + t249;
t323 = t269 * t184 + t270 * t187;
t50 = t210 * t117 + t211 * t323;
t40 = qJ(5) * t262 + t50;
t274 = t211 * t40;
t296 = t176 * pkin(4);
t261 = t178 * t211;
t49 = t211 * t117 - t210 * t323;
t39 = qJ(5) * t261 + t49;
t32 = t39 + t296;
t222 = t210 * t32 - t274;
t245 = mrSges(6,3) * t262;
t123 = -mrSges(6,2) * t176 + t245;
t127 = t176 * mrSges(6,1) + mrSges(6,3) * t261;
t255 = t123 * t303 + t127 * t307;
t333 = -m(6) * ((-t173 * t211 + t174 * t210) * t178 - t222) / 0.2e1 - t255;
t332 = -pkin(4) / 0.2e1;
t331 = mrSges(6,3) / 0.2e1;
t330 = -t123 / 0.2e1;
t328 = Ifges(5,3) + Ifges(6,3);
t317 = m(6) * pkin(4);
t239 = mrSges(6,1) + t317;
t327 = t173 * t210 + t174 * t211;
t286 = Ifges(6,4) * t210;
t188 = Ifges(6,2) * t211 + t286;
t287 = Ifges(5,4) * t210;
t190 = Ifges(5,2) * t211 + t287;
t326 = t188 + t190;
t192 = Ifges(6,1) * t210 + t205;
t194 = Ifges(5,1) * t210 + t206;
t325 = t192 + t194;
t203 = Ifges(6,5) * t211;
t204 = Ifges(5,5) * t211;
t324 = t203 + t204;
t118 = -t178 * pkin(3) + t176 * pkin(7) + t248;
t140 = -t270 * t184 + t187 * t269;
t55 = t211 * t118 + t140 * t210;
t56 = t210 * t118 - t140 * t211;
t321 = -t210 * t55 + t211 * t56;
t320 = m(5) / 0.2e1;
t319 = m(6) / 0.2e1;
t318 = m(4) * pkin(2);
t316 = -mrSges(5,1) / 0.2e1;
t315 = mrSges(5,2) / 0.2e1;
t314 = mrSges(6,2) / 0.2e1;
t263 = t176 * t211;
t33 = -t178 * pkin(4) + qJ(5) * t263 + t55;
t313 = -t33 / 0.2e1;
t312 = -t176 / 0.2e1;
t310 = -t178 / 0.2e1;
t182 = -mrSges(6,1) * t211 + t210 * mrSges(6,2);
t309 = -t182 / 0.2e1;
t308 = -t200 / 0.2e1;
t305 = t210 / 0.2e1;
t299 = m(6) * t178;
t237 = t270 * pkin(2);
t201 = -t237 - pkin(3);
t295 = t211 * pkin(4);
t181 = t201 - t295;
t298 = m(6) * t181;
t297 = pkin(4) * t210;
t293 = mrSges(6,2) * t211;
t292 = mrSges(4,3) * t176;
t291 = mrSges(4,3) * t178;
t290 = mrSges(6,3) * t211;
t185 = t210 * mrSges(6,1) + t293;
t115 = t185 * t178;
t275 = t211 * mrSges(5,2);
t186 = t210 * mrSges(5,1) + t275;
t116 = t186 * t178;
t246 = mrSges(5,3) * t262;
t279 = t176 * mrSges(5,2);
t124 = t246 - t279;
t280 = t176 * mrSges(5,1);
t128 = mrSges(5,3) * t261 + t280;
t266 = t140 * t178;
t273 = t211 * t50;
t87 = -pkin(4) * t262 + t140;
t5 = (t115 + t116 + t291) * t178 + (t292 + (-t123 - t124) * t211 + (t127 + t128) * t210) * t176 + m(5) * (-t266 + (t210 * t49 - t273) * t176) + m(6) * (t176 * t222 - t87 * t178) + m(4) * (-t176 * t323 - t266);
t285 = qJD(1) * t5;
t183 = -mrSges(5,1) * t211 + t210 * mrSges(5,2);
t213 = (-t176 * t200 * t253 - t201 * t178) * t320 + (-t327 * t176 - t181 * t178) * t319 + t178 * t309 + t183 * t310 + (-t176 * t269 + t178 * t270) * t318 / 0.2e1 + (mrSges(5,3) + mrSges(6,3)) * t253 * t312;
t264 = t176 * t210;
t121 = mrSges(6,2) * t178 + mrSges(6,3) * t264;
t122 = mrSges(5,2) * t178 + mrSges(5,3) * t264;
t125 = -t178 * mrSges(6,1) + mrSges(6,3) * t263;
t126 = -t178 * mrSges(5,1) + mrSges(5,3) * t263;
t42 = qJ(5) * t264 + t56;
t214 = (t210 * t56 + t211 * t55) * t320 + (t210 * t42 + t211 * t33) * t319 + t337 / 0.2e1 + (t121 + t122) * t305 + (t125 + t126) * t303;
t6 = t213 - t214 - t336;
t284 = qJD(1) * t6;
t113 = t185 * t176;
t114 = t186 * t176;
t241 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t242 = Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t73 = -Ifges(6,6) * t178 + t176 * t225;
t75 = -Ifges(5,6) * t178 + t176 * t226;
t193 = t211 * Ifges(6,1) - t286;
t77 = -Ifges(6,5) * t178 - t176 * t193;
t78 = Ifges(6,5) * t176 - t178 * t193;
t195 = t211 * Ifges(5,1) - t287;
t79 = -Ifges(5,5) * t178 - t176 * t195;
t80 = Ifges(5,5) * t176 - t178 * t195;
t86 = -pkin(4) * t264 + t323;
t1 = -pkin(1) * (mrSges(3,1) * t301 + t302 * mrSges(3,2)) - t87 * t113 - t86 * t115 + t40 * t121 + t50 * t122 + t42 * t123 + t56 * t124 + t32 * t125 + t49 * t126 + t33 * t127 + t55 * t128 - t323 * t116 - t140 * t114 + m(5) * (t140 * t323 + t49 * t55 + t50 * t56) + m(6) * (t32 * t33 + t40 * t42 + t86 * t87) + (mrSges(4,1) * t248 + Ifges(4,4) * t176 + (-t78 / 0.2e1 - t80 / 0.2e1 - t242 * t176) * t211 + (t74 / 0.2e1 + t76 / 0.2e1 + t241 * t176) * t210) * t176 + (-mrSges(4,2) * t248 - Ifges(4,4) * t178 + (-t77 / 0.2e1 - t79 / 0.2e1 + t242 * t178) * t211 + (t73 / 0.2e1 + t75 / 0.2e1 - t241 * t178) * t210 + (-Ifges(4,2) + Ifges(4,1) - t328) * t176) * t178 + (-Ifges(3,2) + Ifges(3,1)) * t302 * t301 + (-t301 ^ 2 + t302 ^ 2) * Ifges(3,4) + (t336 + t337) * t202;
t283 = t1 * qJD(1);
t282 = t174 * t32;
t281 = t174 * t39;
t168 = mrSges(6,1) * t261;
t228 = mrSges(6,2) * t262 - t168;
t238 = m(6) * (-t32 + t39);
t250 = pkin(4) * t261;
t4 = -t115 * t250 - t87 * t228 - t39 * t123 + t32 * t245 - t183 * t266 + t50 * t128 + (m(6) * t87 * t295 + t325 * t303 * t178 - mrSges(5,3) * t273 - mrSges(6,3) * t274 + t334 * t312 + (t326 * t178 + t78 + t80) * t307 - t335 * t211 / 0.2e1) * t178 + (-t124 + t246) * t49 + (t127 - t238) * t40;
t271 = t4 * qJD(1);
t15 = (-t210 * t123 - m(6) * (t210 * t40 + t32 * t211) - t211 * t127) * t178;
t268 = qJD(1) * t15;
t260 = t200 * t211;
t59 = 0.2e1 * (t208 / 0.4e1 + t209 / 0.4e1 + 0.1e1 / 0.4e1) * t299;
t257 = t59 * qJD(1);
t235 = t264 / 0.2e1;
t254 = mrSges(6,1) * t235 + t263 * t314;
t252 = qJD(4) * t210;
t251 = t317 / 0.2e1;
t243 = t297 / 0.2e1;
t240 = t39 / 0.2e1 - t32 / 0.2e1;
t234 = -t263 / 0.2e1;
t233 = t124 * t308;
t227 = -mrSges(6,1) * t264 / 0.2e1 + t319 * t86 + mrSges(6,2) * t234;
t12 = t227 + t333;
t81 = m(6) * t327 + t253 * mrSges(6,3);
t221 = -qJD(1) * t12 + qJD(2) * t81;
t175 = t210 * t239 + t293;
t65 = m(6) * t250 - t228;
t220 = qJD(1) * t65 - qJD(2) * t175;
t8 = (t279 / 0.2e1 - t124 / 0.2e1 + t330) * t211 + (t280 / 0.2e1 + t127 / 0.2e1 + t128 / 0.2e1 + (t296 / 0.2e1 - t240) * m(6)) * t210 + t254;
t219 = t8 * qJD(1);
t17 = t181 * t185 + t201 * t186 + (t192 / 0.2e1 - t225 / 0.2e1 + t194 / 0.2e1 - t226 / 0.2e1) * t211 + (t193 / 0.2e1 - t188 / 0.2e1 + t195 / 0.2e1 - t190 / 0.2e1 + (t182 + t298) * pkin(4)) * t210;
t212 = (t208 / 0.2e1 + t209 / 0.2e1) * t200 * mrSges(5,3) + (t194 / 0.4e1 + t192 / 0.4e1 - t226 / 0.4e1 - t225 / 0.4e1 + t173 * t331 + t181 * t314 + t201 * t315 + (Ifges(5,1) / 0.4e1 + Ifges(6,1) / 0.4e1) * t210) * t210 + (-t195 / 0.4e1 - t193 / 0.4e1 + t190 / 0.4e1 + t188 / 0.4e1 + t174 * t331 + t201 * t316 + (Ifges(5,2) / 0.4e1 + Ifges(6,2) / 0.4e1) * t211 + (Ifges(5,4) / 0.2e1 + Ifges(6,4) / 0.2e1) * t210 + (-t298 / 0.2e1 + t309) * pkin(4)) * t211;
t215 = (t128 * t308 + t80 / 0.4e1 + t78 / 0.4e1 + t240 * mrSges(6,3)) * t211 + t140 * t186 / 0.2e1 + t173 * t330 - t174 * t127 / 0.2e1 - t181 * t168 / 0.2e1 + t87 * t185 / 0.2e1;
t216 = mrSges(6,1) * t313 + t125 * t332 + t42 * t314 + t56 * t315 + t55 * t316;
t3 = (t115 * t332 + t233 - t76 / 0.4e1 - t74 / 0.4e1) * t210 + (t87 * t243 - t282 / 0.2e1 + t281 / 0.2e1 + pkin(4) * t313) * m(6) + (t204 / 0.4e1 + t203 / 0.4e1 + t242 * t211 + (-0.3e1 / 0.4e1 * Ifges(5,6) - 0.3e1 / 0.4e1 * Ifges(6,6)) * t210) * t176 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + t212) * t178 + t215 + t216;
t218 = -t3 * qJD(1) - t17 * qJD(2);
t58 = (-0.1e1 / 0.2e1 + t253 / 0.2e1) * t299;
t13 = t227 - t333;
t9 = t124 * t303 + t238 * t305 + t128 * t307 + (t275 / 0.2e1 + (mrSges(5,1) / 0.2e1 + t251) * t210) * t176 + t254 + t255;
t7 = t213 + t214;
t2 = t215 - t216 + t212 * t178 + (t297 * t87 + t281 - t282) * t319 - t115 * t243 + t210 * t233 + t33 * t251 + (-t294 * t210 + t324) * t176 / 0.4e1 + t328 * t310 - t335 * t210 / 0.4e1 + t294 * t235 + t329 * t234;
t10 = [qJD(2) * t1 + qJD(3) * t5 - qJD(4) * t4 - qJD(5) * t15, t283 + (-mrSges(3,1) * t249 + m(6) * (-t173 * t33 + t174 * t42 + t181 * t86) + (m(5) * t201 - t270 * t318 - mrSges(4,1) + t183) * t323 - (t269 * t318 - mrSges(4,2)) * t140 + t236 * t291 + t237 * t292 + t42 * t290 + mrSges(3,2) * t247 + t122 * t260 + t325 * t234 + t326 * t235 + (m(5) * t321 - t210 * t126) * t200 + t321 * mrSges(5,3) - t33 * t210 * mrSges(6,3) - Ifges(3,6) * t301 + Ifges(3,5) * t302 + (t75 + t73) * t303 + (t79 + t77) * t305 - t173 * t125 + t174 * t121 - Ifges(4,5) * t176 + Ifges(4,6) * t178 - t181 * t113 + t86 * t182 - t201 * t114 + t334 * t310) * qJD(2) + t7 * qJD(3) + t2 * qJD(4) + t13 * qJD(5), qJD(2) * t7 + qJD(4) * t9 + qJD(5) * t58 + t285, t2 * qJD(2) + t9 * qJD(3) - t271 + (-t50 * mrSges(5,1) - t49 * mrSges(5,2) - t39 * mrSges(6,2) + (t338 + (-mrSges(6,3) * pkin(4) + t329) * t210) * t178 - t239 * t40) * qJD(4), qJD(2) * t13 + qJD(3) * t58 - t268; qJD(3) * t6 + qJD(4) * t3 - qJD(5) * t12 - t283, qJD(4) * t17 + qJD(5) * t81, t284, (-mrSges(5,1) * t260 + t173 * mrSges(6,2) - pkin(4) * t290 - t239 * t174 + t324) * qJD(4) + (mrSges(5,2) * t200 - t294) * t252 - t218, t221; -qJD(2) * t6 - qJD(4) * t8 + qJD(5) * t59 - t285, -t284, 0, (-mrSges(5,2) - mrSges(6,2)) * qJD(4) * t211 + (-mrSges(5,1) - t239) * t252 - t219, t257; -qJD(2) * t3 + qJD(3) * t8 + qJD(5) * t65 + t271, -t175 * qJD(5) + t218, t219, 0, t220; qJD(2) * t12 - qJD(3) * t59 - qJD(4) * t65 + t268, qJD(4) * t175 - t221, -t257, -t220, 0;];
Cq = t10;
