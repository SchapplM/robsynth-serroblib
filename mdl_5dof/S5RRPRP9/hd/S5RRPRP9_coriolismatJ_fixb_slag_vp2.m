% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP9
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:45
% EndTime: 2019-12-31 20:05:54
% DurationCPUTime: 4.10s
% Computational Cost: add. (7208->424), mult. (15908->581), div. (0->0), fcn. (15610->6), ass. (0->209)
t354 = Ifges(6,4) + Ifges(5,5);
t361 = Ifges(5,6) - Ifges(6,6);
t340 = m(6) / 0.2e1;
t360 = 0.2e1 * t340;
t245 = cos(pkin(8));
t233 = -pkin(3) * t245 - pkin(2);
t244 = sin(pkin(8));
t321 = sin(qJ(4));
t322 = cos(qJ(4));
t207 = t321 * t244 - t322 * t245;
t253 = t244 * t322 + t245 * t321;
t262 = pkin(4) * t207 - qJ(5) * t253;
t114 = t233 + t262;
t133 = mrSges(6,1) * t207 - mrSges(6,3) * t253;
t359 = m(6) * t114 + t133;
t246 = sin(qJ(2));
t247 = cos(qJ(2));
t221 = t246 * pkin(2) - qJ(3) * t247;
t290 = t244 * t246;
t171 = pkin(6) * t290 + t245 * t221;
t288 = t245 * t246;
t172 = -pkin(6) * t288 + t244 * t221;
t358 = -t171 * t244 + t172 * t245;
t357 = -Ifges(5,6) / 0.2e1;
t185 = t253 * t247;
t332 = t185 / 0.2e1;
t187 = t247 * t207;
t356 = t187 / 0.2e1;
t353 = Ifges(6,2) + Ifges(5,3);
t184 = t253 * t246;
t352 = t184 * mrSges(5,3);
t186 = t207 * t246;
t351 = t186 * mrSges(5,3);
t286 = t247 * qJ(5);
t219 = -pkin(2) * t247 - t246 * qJ(3) - pkin(1);
t205 = t245 * t219;
t126 = -pkin(7) * t288 + t205 + (-pkin(6) * t244 - pkin(3)) * t247;
t273 = t321 * t126;
t287 = t245 * t247;
t159 = pkin(6) * t287 + t244 * t219;
t142 = -pkin(7) * t290 + t159;
t275 = t322 * t142;
t49 = t275 + t273;
t39 = t49 - t286;
t238 = t247 * mrSges(6,3);
t149 = -t184 * mrSges(6,2) - t238;
t150 = mrSges(5,2) * t247 - t352;
t350 = -t149 - t150;
t152 = -mrSges(5,1) * t247 + t351;
t153 = mrSges(6,1) * t247 - t186 * mrSges(6,2);
t349 = -t152 + t153;
t348 = -t354 * t207 - t361 * t253;
t347 = -t354 * t184 + t361 * t186;
t237 = m(6) * qJ(5) + mrSges(6,3);
t346 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t345 = -mrSges(5,2) + t237;
t344 = 0.2e1 * m(6);
t243 = t245 ^ 2;
t343 = m(4) / 0.2e1;
t342 = m(5) / 0.2e1;
t341 = -m(6) / 0.2e1;
t339 = m(4) * pkin(6);
t338 = -mrSges(5,1) / 0.2e1;
t337 = -mrSges(5,3) / 0.2e1;
t336 = Ifges(5,4) * t356 + Ifges(5,2) * t332 + t246 * t357;
t294 = qJ(5) * t207;
t319 = pkin(4) * t253;
t263 = -t294 - t319;
t335 = -t263 / 0.2e1;
t317 = pkin(7) + qJ(3);
t220 = t317 * t245;
t268 = t317 * t244;
t146 = t220 * t321 + t268 * t322;
t334 = -t146 / 0.2e1;
t333 = -t185 / 0.2e1;
t331 = t186 / 0.2e1;
t330 = -t187 / 0.2e1;
t329 = t207 / 0.2e1;
t327 = -t244 / 0.2e1;
t326 = t245 / 0.2e1;
t325 = t246 / 0.2e1;
t320 = pkin(4) * t186;
t240 = t246 * pkin(6);
t241 = t247 * pkin(6);
t316 = Ifges(4,4) * t244;
t315 = Ifges(4,4) * t245;
t314 = Ifges(5,4) * t186;
t313 = Ifges(5,4) * t253;
t312 = Ifges(4,5) * t245;
t311 = Ifges(6,5) * t184;
t310 = Ifges(6,5) * t207;
t309 = Ifges(4,6) * t244;
t289 = t244 * t247;
t158 = -pkin(6) * t289 + t205;
t213 = t247 * mrSges(4,2) - mrSges(4,3) * t290;
t215 = -t247 * mrSges(4,1) - mrSges(4,3) * t288;
t48 = t126 * t322 - t142 * t321;
t42 = t247 * pkin(4) - t48;
t8 = -t349 * t186 + t350 * t184 + m(5) * (-t184 * t49 + t186 * t48) + m(6) * (-t184 * t39 - t186 * t42) + (m(4) * (-t158 * t245 - t159 * t244) - t244 * t213 - t245 * t215) * t246;
t308 = qJD(1) * t8;
t307 = t185 * mrSges(5,1);
t306 = t185 * mrSges(6,1);
t305 = t187 * mrSges(5,2);
t304 = t187 * mrSges(6,2);
t303 = t187 * mrSges(6,3);
t302 = t207 * mrSges(6,2);
t301 = t244 * mrSges(4,1);
t300 = t244 * Ifges(4,2);
t299 = t245 * mrSges(4,2);
t298 = t246 * mrSges(6,1);
t103 = mrSges(6,1) * t184 + mrSges(6,3) * t186;
t104 = t303 + t306;
t105 = -t305 + t307;
t148 = -mrSges(6,2) * t185 + mrSges(6,3) * t246;
t151 = -mrSges(5,2) * t246 - mrSges(5,3) * t185;
t154 = mrSges(5,1) * t246 + mrSges(5,3) * t187;
t155 = -t298 - t304;
t181 = Ifges(4,6) * t246 + (-t300 + t315) * t247;
t182 = Ifges(4,5) * t246 + (t245 * Ifges(4,1) - t316) * t247;
t195 = (t299 + t301) * t247;
t214 = -t246 * mrSges(4,2) - mrSges(4,3) * t289;
t216 = t246 * mrSges(4,1) - mrSges(4,3) * t287;
t217 = pkin(3) * t290 + t240;
t218 = pkin(3) * t289 + t241;
t278 = Ifges(6,6) / 0.2e1 + t357;
t279 = Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t252 = -t185 * t278 + t187 * t279;
t94 = -Ifges(6,1) * t186 - t247 * Ifges(6,4) + t311;
t180 = Ifges(5,4) * t184;
t96 = -Ifges(5,1) * t186 - t247 * Ifges(5,5) - t180;
t276 = t94 / 0.2e1 + t96 / 0.2e1;
t177 = Ifges(6,5) * t186;
t90 = -t247 * Ifges(6,6) + Ifges(6,3) * t184 - t177;
t92 = -Ifges(5,2) * t184 - t247 * Ifges(5,6) - t314;
t277 = t90 / 0.2e1 - t92 / 0.2e1;
t127 = t246 * pkin(3) - pkin(7) * t287 + t171;
t145 = -pkin(7) * t289 + t172;
t51 = t321 * t127 + t322 * t145;
t43 = qJ(5) * t246 + t51;
t50 = t127 * t322 - t145 * t321;
t44 = -t246 * pkin(4) - t50;
t264 = pkin(4) * t184 + qJ(5) * t186;
t69 = t264 + t217;
t70 = pkin(4) * t185 + qJ(5) * t187 + t218;
t91 = -Ifges(6,5) * t187 + Ifges(6,6) * t246 + Ifges(6,3) * t185;
t95 = -Ifges(6,1) * t187 + Ifges(6,4) * t246 + Ifges(6,5) * t185;
t97 = -Ifges(5,1) * t187 - Ifges(5,4) * t185 + Ifges(5,5) * t246;
t3 = t70 * t103 + t69 * t104 + t217 * t105 + t39 * t148 + t43 * t149 + t51 * t150 + t49 * t151 + t50 * t152 + t44 * t153 + t48 * t154 + t42 * t155 + t158 * t216 + t159 * t214 + t171 * t215 + t172 * t213 - t276 * t187 + t277 * t185 - (t95 / 0.2e1 + t97 / 0.2e1 + t218 * mrSges(5,2)) * t186 + (t91 / 0.2e1 + t336 + t218 * mrSges(5,1)) * t184 + m(5) * (t217 * t218 + t48 * t50 + t49 * t51) + m(6) * (t39 * t43 + t42 * t44 + t69 * t70) + m(4) * (t158 * t171 + t159 * t172) + (t182 * t326 + t181 * t327 + pkin(6) * t195 - pkin(1) * mrSges(3,1) + (-Ifges(3,4) + t312 / 0.2e1 - t309 / 0.2e1) * t246 - t279 * t186 + t278 * t184) * t246 + (-pkin(1) * mrSges(3,2) + (t243 * Ifges(4,1) / 0.2e1 - Ifges(3,2) + Ifges(3,1) - Ifges(4,3) + (t299 + t339) * pkin(6) + (pkin(6) * mrSges(4,1) - t315 + t300 / 0.2e1) * t244 - t353) * t246 + t252 + (Ifges(3,4) + t309 - t312) * t247) * t247;
t297 = t3 * qJD(1);
t295 = qJ(5) * t184;
t100 = t295 - t320;
t101 = -t186 * mrSges(6,1) + t184 * mrSges(6,3);
t102 = -t186 * mrSges(5,1) - t184 * mrSges(5,2);
t106 = -Ifges(6,3) * t186 - t311;
t107 = Ifges(5,2) * t186 - t180;
t108 = -Ifges(6,1) * t184 - t177;
t109 = -Ifges(5,1) * t184 + t314;
t4 = t217 * t102 + t69 * t101 - (-t39 * mrSges(6,2) + t109 / 0.2e1 + t108 / 0.2e1 + t277) * t186 + (-t42 * mrSges(6,2) - t107 / 0.2e1 + t106 / 0.2e1 - t276) * t184 + (m(6) * t69 + t103) * t100 + (m(6) * t42 + t349 + t351) * t49 + (m(6) * t39 - t350 + t352) * t48 - t347 * t247 / 0.2e1;
t296 = t4 * qJD(1);
t16 = m(6) * (t69 * t186 - t247 * t39) - t247 * t149 + t186 * t103;
t293 = qJD(1) * t16;
t280 = t337 - mrSges(6,2) / 0.2e1;
t271 = t247 * t329;
t270 = -t149 / 0.2e1 - t150 / 0.2e1;
t269 = -t152 / 0.2e1 + t153 / 0.2e1;
t132 = mrSges(5,1) * t253 - t207 * mrSges(5,2);
t131 = mrSges(6,1) * t253 + t207 * mrSges(6,3);
t147 = t220 * t322 - t268 * t321;
t266 = -t146 * t186 - t147 * t184;
t134 = Ifges(6,3) * t253 - t310;
t200 = Ifges(6,5) * t253;
t135 = Ifges(6,3) * t207 + t200;
t203 = Ifges(5,4) * t207;
t136 = -Ifges(5,2) * t253 - t203;
t137 = -Ifges(5,2) * t207 + t313;
t138 = -Ifges(6,1) * t207 + t200;
t139 = Ifges(6,1) * t253 + t310;
t140 = -Ifges(5,1) * t207 - t313;
t141 = Ifges(5,1) * t253 - t203;
t248 = t270 * t146 + t269 * t147 + (t106 / 0.4e1 - t107 / 0.4e1 - t96 / 0.4e1 - t94 / 0.4e1) * t207 - (t92 / 0.4e1 - t109 / 0.4e1 - t108 / 0.4e1 - t90 / 0.4e1) * t253 + (t147 * t331 + t184 * t334 - (-t49 / 0.2e1 + t39 / 0.2e1) * t253 + (-t48 / 0.2e1 - t42 / 0.2e1) * t207) * mrSges(6,2) + (mrSges(5,3) * t334 + t134 / 0.4e1 - t141 / 0.4e1 - t139 / 0.4e1 - t136 / 0.4e1) * t184 - (t147 * t337 + t140 / 0.4e1 + t138 / 0.4e1 + t135 / 0.4e1 - t137 / 0.4e1) * t186 + (t100 * t114 - t263 * t69 + (-t39 + t49) * t146 + (t42 + t48) * t147) * t340 + t100 * t133 / 0.2e1 + t114 * t101 / 0.2e1 + t103 * t335 + t217 * t132 / 0.2e1 + t233 * t102 / 0.2e1 + t69 * t131 / 0.2e1 - t348 * t247 / 0.4e1;
t250 = (-pkin(4) * t44 + qJ(5) * t43) * t341 + pkin(4) * t155 / 0.2e1 - qJ(5) * t148 / 0.2e1 - t43 * mrSges(6,3) / 0.2e1 + t44 * mrSges(6,1) / 0.2e1 + t50 * t338 + t51 * mrSges(5,2) / 0.2e1;
t1 = (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1) * t246 + t248 + t250 + t252;
t7 = t114 * t131 + t233 * t132 - (-t138 / 0.2e1 - t135 / 0.2e1 - t140 / 0.2e1 + t137 / 0.2e1) * t253 + (-t139 / 0.2e1 + t134 / 0.2e1 - t141 / 0.2e1 - t136 / 0.2e1) * t207 - t359 * t263;
t261 = t1 * qJD(1) + t7 * qJD(2);
t13 = (t207 ^ 2 + t253 ^ 2) * (mrSges(5,3) + mrSges(6,2)) + (m(4) * qJ(3) + mrSges(4,3)) * (t244 ^ 2 + t243) + (m(6) + m(5)) * (t146 * t253 - t147 * t207);
t249 = (-t184 * t280 + t270) * t207 - (-t186 * t280 - t269) * t253 + (-t158 * t244 + t159 * t245) * t343 + (-t207 * t49 - t253 * t48 + t266) * t342 + (-t207 * t39 + t253 * t42 + t266) * t340 + t215 * t327 + t213 * t326;
t255 = -m(5) * t218 / 0.2e1 + t70 * t341;
t6 = -(-mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1) * t187 + (t338 - mrSges(6,1) / 0.2e1) * t185 + (-t339 / 0.2e1 - t299 / 0.2e1 - t301 / 0.2e1) * t247 + t249 + t255;
t260 = qJD(1) * t6 + qJD(2) * t13;
t251 = (t114 * t186 - t147 * t247 - t253 * t69) * t340 + t133 * t331 - t253 * t103 / 0.2e1;
t254 = t44 * t341 + t298 / 0.2e1;
t12 = (t271 + t356) * mrSges(6,2) + t251 + t254;
t24 = t359 * t253;
t259 = -qJD(1) * t12 + qJD(2) * t24;
t19 = (t320 / 0.4e1 - t295 / 0.4e1 - t100 / 0.4e1) * t344 - t101 - t102;
t22 = (-t319 / 0.4e1 - t294 / 0.4e1 + t263 / 0.4e1) * t344 - t131 - t132;
t258 = qJD(1) * t19 + qJD(2) * t22;
t116 = m(6) * t253;
t86 = m(6) * t186;
t257 = -qJD(1) * t86 + qJD(2) * t116;
t18 = -t238 + (t275 / 0.4e1 + t273 / 0.4e1 - t286 / 0.2e1 - t49 / 0.4e1) * t344;
t256 = qJD(1) * t18 + qJD(4) * t237;
t52 = m(6) * t147 - t302;
t47 = m(6) * t335 + t263 * t340;
t17 = t360 * t39 + t149;
t11 = mrSges(6,2) * t271 - t304 / 0.2e1 + t251 - t254;
t5 = -t305 / 0.2e1 + t307 / 0.2e1 + t303 / 0.2e1 + t306 / 0.2e1 + t241 * t343 + mrSges(4,2) * t287 / 0.2e1 + mrSges(4,1) * t289 / 0.2e1 + t249 - t255;
t2 = Ifges(5,6) * t333 + Ifges(6,6) * t332 + t353 * t325 + t354 * t330 + t248 - t250;
t9 = [qJD(2) * t3 + qJD(3) * t8 + qJD(4) * t4 + qJD(5) * t16, t5 * qJD(3) + t2 * qJD(4) + t11 * qJD(5) + t297 + ((t141 + t139) * t330 + t358 * mrSges(4,3) + 0.2e1 * (-pkin(2) * t241 + qJ(3) * t358) * t343 + (t97 + t95) * t253 / 0.2e1 + (-t51 * t207 - t253 * t50) * mrSges(5,3) + t218 * (t207 * mrSges(5,1) + mrSges(5,2) * t253) - t43 * t302 + (t114 * t70 + t146 * t44 + t147 * t43) * t360 + t44 * t253 * mrSges(6,2) + t245 * qJ(3) * t214 - t244 * qJ(3) * t216 + (Ifges(4,5) * t244 + Ifges(4,6) * t245 - t361 * t207 + t354 * t253) * t325 + mrSges(3,2) * t240 + t181 * t326 + t91 * t329 + t135 * t332 + t137 * t333 + t207 * t336 + 0.2e1 * (-t146 * t50 + t147 * t51 + t218 * t233) * t342 + ((Ifges(4,1) * t244 + t315) * t326 + (Ifges(4,2) * t245 + t316) * t327 + Ifges(3,5) + (-mrSges(4,1) * t245 + mrSges(4,2) * t244 - mrSges(3,1)) * pkin(6)) * t247 + t114 * t104 + t70 * t133 + t147 * t148 + t147 * t151 - t146 * t154 + t146 * t155 - pkin(2) * t195 + t233 * t105 + t244 * t182 / 0.2e1 - Ifges(3,6) * t246) * qJD(2), qJD(2) * t5 + t308, t296 + t2 * qJD(2) + (t264 * mrSges(6,2) + t345 * t48 + t346 * t49 + t347) * qJD(4) + t17 * qJD(5), qJD(2) * t11 + qJD(4) * t17 + t293; qJD(3) * t6 + qJD(4) * t1 + qJD(5) * t12 - t297, qJD(3) * t13 + qJD(4) * t7 - qJD(5) * t24, qJD(4) * t47 + t260, t47 * qJD(3) + (t262 * mrSges(6,2) - t345 * t146 + t346 * t147 + t348) * qJD(4) + t52 * qJD(5) + t261, qJD(4) * t52 - t259; -qJD(2) * t6 - qJD(4) * t19 + qJD(5) * t86 - t308, -qJD(4) * t22 - qJD(5) * t116 - t260, 0, -t258, -t257; -qJD(2) * t1 + qJD(3) * t19 + qJD(5) * t18 - t296, qJD(3) * t22 - t261, t258, t237 * qJD(5), t256; -qJD(2) * t12 - qJD(3) * t86 - qJD(4) * t18 - t293, qJD(3) * t116 + t259, t257, -t256, 0;];
Cq = t9;
