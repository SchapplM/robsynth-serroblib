% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:50
% EndTime: 2019-12-05 18:45:03
% DurationCPUTime: 5.27s
% Computational Cost: add. (11344->277), mult. (22135->357), div. (0->0), fcn. (24005->6), ass. (0->161)
t212 = sin(qJ(3));
t213 = sin(qJ(2));
t215 = cos(qJ(3));
t216 = cos(qJ(2));
t192 = -t212 * t213 + t215 * t216;
t193 = -t212 * t216 - t215 * t213;
t208 = -pkin(2) * t216 - pkin(1);
t211 = sin(qJ(4));
t214 = cos(qJ(4));
t152 = t192 * t211 - t193 * t214;
t164 = -t192 * pkin(3) + t208;
t241 = t214 * t192 + t193 * t211;
t340 = t241 * mrSges(6,2);
t346 = Ifges(5,4) + Ifges(6,4);
t347 = t241 / 0.2e1;
t353 = t152 * t346;
t364 = -t152 / 0.2e1;
t375 = Ifges(5,2) + Ifges(6,2);
t376 = Ifges(5,1) + Ifges(6,1);
t82 = -pkin(4) * t241 + t164;
t386 = t164 * (mrSges(5,1) * t152 + mrSges(5,2) * t241) + (t241 * t375 + t353) * t364 + t82 * (t152 * mrSges(6,1) + t340) + (t376 * t241 - t353) * t152 / 0.2e1 + (0.2e1 * t346 * t241 + (-t375 + t376) * t152) * t347;
t410 = t208 * (-mrSges(4,1) * t193 + mrSges(4,2) * t192) + t386;
t362 = Ifges(6,5) + Ifges(5,5);
t360 = Ifges(6,6) * t152;
t361 = Ifges(5,6) * t152;
t366 = -t360 / 0.2e1 - t361 / 0.2e1;
t328 = -pkin(7) - pkin(6);
t197 = t328 * t213;
t198 = t328 * t216;
t161 = t197 * t212 - t215 * t198;
t119 = pkin(8) * t192 + t161;
t338 = t215 * t197 + t212 * t198;
t352 = t193 * pkin(8) + t338;
t382 = -t211 * t119 + t214 * t352;
t389 = t382 * mrSges(5,2);
t387 = -qJ(5) * t152 + t382;
t397 = t387 * mrSges(6,2);
t66 = t214 * t119 + t211 * t352;
t398 = t66 * mrSges(5,1);
t42 = qJ(5) * t241 + t66;
t405 = t42 * mrSges(6,1);
t407 = -t405 / 0.2e1 - t389 / 0.2e1 - t397 / 0.2e1 - t398 / 0.2e1;
t227 = t362 * t347 + t366 + t407;
t341 = mrSges(6,3) * t241;
t246 = -t341 / 0.2e1;
t145 = Ifges(6,5) * t241;
t146 = Ifges(5,5) * t241;
t334 = t145 / 0.2e1 + t146 / 0.2e1 + t366;
t237 = pkin(4) * t246 + t334;
t409 = t227 + t237 + t407;
t367 = -t360 - t361;
t377 = -t161 * mrSges(4,1) - t338 * mrSges(4,2) + Ifges(4,5) * t192 + Ifges(4,6) * t193 + t145 + t146 + t367;
t395 = -t405 - t389 - t397 - t398;
t408 = t377 + t395;
t396 = m(6) * t82 - mrSges(6,1) * t241 + mrSges(6,2) * t152;
t327 = pkin(4) * t42;
t313 = pkin(3) * t214;
t206 = pkin(4) + t313;
t275 = t206 * t42;
t399 = m(5) * t164;
t392 = m(5) * (t211 * t382 - t214 * t66);
t207 = pkin(2) * t215 + pkin(3);
t257 = t212 * t214;
t182 = pkin(2) * t257 + t207 * t211;
t186 = (-t211 * t215 - t257) * pkin(2);
t384 = t182 + t186;
t326 = -t340 / 0.2e1;
t363 = pkin(4) * t152;
t355 = t206 - t313;
t314 = pkin(3) * t211;
t348 = -pkin(3) / 0.2e1;
t354 = -m(6) * (t314 * t387 - t275) / 0.2e1 + t348 * t392 - t407;
t183 = t186 * mrSges(6,1);
t184 = t186 * mrSges(5,1);
t351 = (mrSges(4,1) * t212 + mrSges(4,2) * t215) * pkin(2) - t183 - t184;
t298 = mrSges(5,2) + mrSges(6,2);
t173 = t182 * mrSges(6,1);
t174 = t182 * mrSges(5,1);
t258 = t211 * t212;
t181 = -pkin(2) * t258 + t214 * t207;
t335 = -t298 * t181 - t173 - t174;
t332 = 0.2e1 * m(6);
t331 = m(5) / 0.2e1;
t330 = m(6) / 0.2e1;
t329 = m(6) * pkin(4);
t176 = pkin(4) + t181;
t320 = -t176 / 0.2e1;
t319 = -t181 / 0.2e1;
t318 = -t186 / 0.2e1;
t187 = (t214 * t215 - t258) * pkin(2);
t317 = t187 / 0.2e1;
t316 = m(6) * t182;
t315 = pkin(3) * t193;
t210 = t213 * pkin(2);
t293 = Ifges(4,4) * t193;
t238 = Ifges(4,4) * t192 + (-Ifges(4,1) + Ifges(4,2)) * t193;
t74 = -mrSges(5,1) * t241 + mrSges(5,2) * t152;
t99 = -t315 + t363;
t84 = t210 + t99;
t1 = (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t216 + (Ifges(3,1) - Ifges(3,2)) * t213) * t216 + (-mrSges(4,1) * t210 + t238) * t192 + (-mrSges(4,2) * t210 - t293) * t193 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t213) * t213 + m(4) * t208 * t210 + (t399 + t74) * (t210 - t315) + t396 * t84 + t410;
t290 = t1 * qJD(1);
t278 = t187 * mrSges(5,2);
t277 = t187 * mrSges(6,2);
t2 = t238 * t192 + (-pkin(3) * t74 - t293) * t193 - t315 * t399 + t396 * t99 + t410;
t276 = t2 * qJD(1);
t5 = t396 * t363 + t386;
t273 = t5 * qJD(1);
t15 = m(6) * (-t152 * t387 + t241 * t42) + (t152 ^ 2 + t241 ^ 2) * mrSges(6,3);
t270 = qJD(1) * t15;
t269 = t152 * t182;
t268 = t241 * t182;
t267 = t241 * t176;
t266 = t241 * t181;
t265 = t152 * t176;
t264 = t152 * t206;
t263 = t182 * t187;
t223 = 0.2e1 * t364 * mrSges(6,1) + 0.2e1 * t326;
t19 = (t268 / 0.4e1 - t265 / 0.4e1 - t84 / 0.4e1) * t332 + t223;
t262 = t19 * qJD(1);
t251 = t241 * t314;
t21 = (t251 / 0.4e1 - t264 / 0.4e1 - t99 / 0.4e1) * t332 + t223;
t260 = t21 * qJD(1);
t259 = t211 * t152;
t256 = t214 * t241;
t247 = -mrSges(6,1) - t329;
t45 = t247 * t152 - t340;
t255 = t45 * qJD(1);
t253 = -t176 + t181;
t250 = -t327 / 0.2e1;
t249 = -t313 / 0.2e1;
t248 = -mrSges(5,2) / 0.2e1 - mrSges(6,2) / 0.2e1;
t244 = t253 * t42;
t240 = t259 * t348;
t239 = (-mrSges(6,1) / 0.2e1 - mrSges(5,1) / 0.2e1) * t211;
t23 = -t253 * t316 - t335;
t6 = (Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t152 + (-t244 / 0.2e1 + t250) * m(6) + (-Ifges(5,5) / 0.2e1 - Ifges(6,5) / 0.2e1 + (t176 / 0.2e1 + t319) * mrSges(6,3)) * t241 + t237;
t236 = -t6 * qJD(1) - t23 * qJD(2);
t24 = t298 * t187 - m(6) * (t176 * t186 + t263) - m(5) * (t181 * t186 + t263) + t351;
t218 = (t384 * t382 + (-t181 + t187) * t66) * t331 + (t384 * t387 + (-t176 + t187) * t42) * t330 + t407;
t230 = t152 * t318 + t241 * t317;
t224 = -t269 / 0.2e1 + t230;
t221 = -t266 / 0.2e1 + t224;
t4 = ((t206 / 0.2e1 + t320) * t241 - (-t314 / 0.2e1 + t182 / 0.2e1) * t152 + t230) * mrSges(6,3) + ((t256 / 0.2e1 + t259 / 0.2e1) * pkin(3) + t221) * mrSges(5,3) + t218 + t354;
t235 = t4 * qJD(1) - t24 * qJD(2);
t234 = (t313 / 0.2e1 - t206 / 0.2e1) * t241;
t229 = t340 / 0.2e1 + t326;
t228 = -t183 / 0.2e1 - t184 / 0.2e1 + t318 * t329;
t226 = t42 * t313;
t105 = t298 * t313 + (t355 * m(6) + mrSges(5,1) + mrSges(6,1)) * t314;
t219 = -t173 / 0.2e1 - t174 / 0.2e1 + (-t355 * t182 + t253 * t314) * t330;
t17 = pkin(3) * t239 + t219 + t228 + t298 * (t319 + t249 + t317);
t8 = t227 + (pkin(4) * t347 + t234) * mrSges(6,3) + (-t275 / 0.2e1 + t226 / 0.2e1 + t327 / 0.2e1) * m(6) - t407 - t334;
t225 = t8 * qJD(1) + t17 * qJD(2) - t105 * qJD(3);
t172 = t187 * t314;
t20 = t229 + (t251 - t264 + t99) * t330;
t18 = t229 + (-t265 + t268 + t84) * t330;
t16 = -t278 / 0.2e1 - t277 / 0.2e1 + t248 * t181 + (t248 * t214 + t239) * pkin(3) + t219 - t228;
t9 = mrSges(6,3) * t234 + (t226 - t275 - t327) * t330 + t409;
t7 = t244 * t330 + m(6) * t250 + (t320 + t181 / 0.2e1) * t341 + t409;
t3 = t218 + t206 * t246 + (t240 - t267 / 0.2e1 + t224) * mrSges(6,3) - t354 + t377 + (t249 * t241 + t221 + t240) * mrSges(5,3);
t10 = [qJD(2) * t1 + qJD(3) * t2 + qJD(4) * t5 + qJD(5) * t15, t3 * qJD(3) + t7 * qJD(4) + t18 * qJD(5) + t290 + (Ifges(3,5) * t216 - Ifges(3,6) * t213 + 0.2e1 * (-t176 * t42 + t182 * t387) * t330 + 0.2e1 * (-t181 * t66 + t182 * t382) * t331 + (m(4) * (-t161 * t215 + t212 * t338) + (-t192 * t215 + t193 * t212) * mrSges(4,3)) * pkin(2) + (-mrSges(3,1) * t216 + mrSges(3,2) * t213) * pkin(6) + (-t267 - t269) * mrSges(6,3) + (-t266 - t269) * mrSges(5,3) + t408) * qJD(2), t3 * qJD(2) + t9 * qJD(4) + t20 * qJD(5) + t276 + ((-mrSges(6,3) * t259 + (-t256 - t259) * mrSges(5,3) + m(6) * t211 * t387 + t392) * pkin(3) + (-m(6) * t42 - t341) * t206 + t408) * qJD(3), t7 * qJD(2) + t9 * qJD(3) + t273 + (-t42 * t329 + (-mrSges(6,3) * pkin(4) + t362) * t241 + t367 + t395) * qJD(4), qJD(2) * t18 + qJD(3) * t20 + t270; qJD(3) * t4 - qJD(4) * t6 + qJD(5) * t19 - t290, -qJD(3) * t24 - qJD(4) * t23, (-t277 - t278 - t351) * qJD(3) + t16 * qJD(4) + 0.2e1 * ((t186 * t206 + t172) * t330 + (t186 * t313 + t172) * t331) * qJD(3) + t235, t16 * qJD(3) + (-pkin(4) * t316 + t335) * qJD(4) + t236, t262; -qJD(2) * t4 + qJD(4) * t8 + qJD(5) * t21 - t276, qJD(4) * t17 - t235, -t105 * qJD(4), (-t298 * t214 + (-mrSges(5,1) + t247) * t211) * qJD(4) * pkin(3) + t225, t260; qJD(2) * t6 - qJD(3) * t8 + qJD(5) * t45 - t273, -qJD(3) * t17 - t236, -t225, 0, t255; -qJD(2) * t19 - qJD(3) * t21 - qJD(4) * t45 - t270, -t262, -t260, -t255, 0;];
Cq = t10;
