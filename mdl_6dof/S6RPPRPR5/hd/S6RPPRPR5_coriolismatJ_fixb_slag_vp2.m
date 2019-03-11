% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:33
% EndTime: 2019-03-09 01:48:38
% DurationCPUTime: 3.21s
% Computational Cost: add. (7049->417), mult. (14214->573), div. (0->0), fcn. (13595->6), ass. (0->223)
t229 = sin(pkin(9));
t230 = cos(pkin(9));
t346 = sin(qJ(6));
t347 = cos(qJ(6));
t194 = t346 * t229 - t347 * t230;
t236 = cos(qJ(4));
t171 = t194 * t236;
t255 = -t229 * t347 - t230 * t346;
t310 = t255 * t171;
t169 = t255 * t236;
t311 = t194 * t169;
t61 = -t310 + t311;
t377 = -t61 / 0.2e1;
t228 = t230 ^ 2;
t376 = -t228 / 0.2e1;
t227 = t229 ^ 2;
t297 = t227 + t228;
t280 = t297 * mrSges(6,3);
t235 = sin(qJ(4));
t375 = t235 * t236;
t329 = t171 * mrSges(7,2);
t331 = t169 * mrSges(7,1);
t100 = -t329 - t331;
t323 = t230 * mrSges(6,2);
t326 = t229 * mrSges(6,1);
t275 = t323 + t326;
t184 = t275 * t236;
t374 = t100 + t184;
t344 = pkin(4) * t236;
t209 = qJ(5) * t235 + t344;
t198 = t230 * t209;
t231 = qJ(2) - pkin(7);
t301 = t231 * t236;
t143 = -t229 * t301 + t198;
t144 = t229 * t209 + t230 * t301;
t272 = -t143 * t229 + t144 * t230;
t232 = pkin(1) + qJ(3);
t345 = pkin(4) * t235;
t205 = -qJ(5) * t236 + t232 + t345;
t186 = t230 * t205;
t300 = t235 * t231;
t140 = -t229 * t300 + t186;
t141 = t229 * t205 + t230 * t300;
t305 = t229 * t236;
t201 = -t235 * mrSges(6,2) - mrSges(6,3) * t305;
t302 = t230 * t236;
t203 = t235 * mrSges(6,1) - mrSges(6,3) * t302;
t373 = -m(6) * (t140 * t230 + t141 * t229) - t229 * t201 - t230 * t203;
t251 = (t310 / 0.2e1 - t311 / 0.2e1) * mrSges(7,3);
t137 = -mrSges(7,2) * t235 + t169 * mrSges(7,3);
t139 = mrSges(7,1) * t235 + t171 * mrSges(7,3);
t355 = -t255 / 0.2e1;
t356 = t194 / 0.2e1;
t261 = t137 * t356 + t139 * t355;
t243 = t251 + t261;
t170 = t194 * t235;
t330 = t170 * mrSges(7,2);
t168 = t255 * t235;
t332 = t168 * mrSges(7,1);
t264 = t332 / 0.2e1 + t330 / 0.2e1;
t13 = t243 - t264;
t372 = t13 * qJD(1);
t371 = t235 ^ 2;
t370 = t236 ^ 2;
t369 = -m(6) / 0.2e1;
t368 = m(6) / 0.2e1;
t367 = -m(7) / 0.2e1;
t366 = m(7) / 0.2e1;
t365 = -mrSges(7,3) / 0.2e1;
t50 = t168 * t171 - t169 * t170;
t363 = m(7) * t50;
t362 = -t137 / 0.2e1;
t361 = -t139 / 0.2e1;
t360 = -t168 / 0.2e1;
t359 = t169 / 0.2e1;
t358 = t170 / 0.2e1;
t357 = -t171 / 0.2e1;
t354 = -t229 / 0.2e1;
t353 = t229 / 0.2e1;
t352 = -t230 / 0.2e1;
t351 = t230 / 0.2e1;
t350 = t235 / 0.2e1;
t349 = -t236 / 0.2e1;
t343 = pkin(8) + qJ(5);
t341 = t61 * qJD(4) * t366;
t340 = m(7) * qJD(2);
t339 = m(7) * qJD(3);
t337 = Ifges(6,4) * t229;
t336 = Ifges(6,4) * t230;
t335 = Ifges(7,4) * t171;
t334 = Ifges(7,4) * t255;
t136 = -mrSges(7,2) * t236 - t168 * mrSges(7,3);
t138 = mrSges(7,1) * t236 - t170 * mrSges(7,3);
t325 = t229 * Ifges(6,2);
t164 = Ifges(6,6) * t236 + (t325 - t336) * t235;
t165 = Ifges(6,5) * t236 + (-t230 * Ifges(6,1) + t337) * t235;
t183 = t275 * t235;
t281 = pkin(5) * t229 - t231;
t192 = t281 * t235;
t193 = t281 * t236;
t306 = t229 * t235;
t200 = -mrSges(6,2) * t236 + mrSges(6,3) * t306;
t303 = t230 * t235;
t202 = mrSges(6,1) * t236 + mrSges(6,3) * t303;
t226 = t236 * mrSges(5,1);
t265 = Ifges(7,5) * t358 + Ifges(7,6) * t360;
t322 = t230 * Ifges(6,5);
t324 = t229 * Ifges(6,6);
t282 = -t229 * t231 + pkin(5);
t111 = -pkin(8) * t302 + t235 * t282 + t186;
t125 = -pkin(8) * t305 + t141;
t51 = t111 * t347 - t125 * t346;
t52 = t111 * t346 + t125 * t347;
t123 = pkin(8) * t303 + t236 * t282 + t198;
t132 = pkin(8) * t306 + t144;
t57 = t123 * t347 - t132 * t346;
t58 = t123 * t346 + t132 * t347;
t92 = Ifges(7,4) * t170 - Ifges(7,2) * t168 + Ifges(7,6) * t236;
t93 = Ifges(7,2) * t169 + t235 * Ifges(7,6) - t335;
t94 = Ifges(7,1) * t170 - Ifges(7,4) * t168 + Ifges(7,5) * t236;
t157 = Ifges(7,4) * t169;
t95 = -Ifges(7,1) * t171 + t235 * Ifges(7,5) + t157;
t99 = t330 + t332;
t1 = t52 * t136 + t58 * t137 + t51 * t138 + t57 * t139 + t93 * t360 + t92 * t359 + t95 * t358 + t94 * t357 - t192 * t100 + t193 * t99 + t141 * t200 + t144 * t201 + t140 * t202 + t143 * t203 + t232 * t226 + m(6) * (t140 * t143 + t141 * t144) + m(7) * (-t192 * t193 + t51 * t57 + t52 * t58) + (t231 * t183 + t165 * t351 + t164 * t354 + Ifges(7,5) * t357 + Ifges(7,6) * t359 + (-Ifges(5,4) + t322 / 0.2e1 - t324 / 0.2e1) * t236) * t236 + (t231 * t184 - t232 * mrSges(5,2) + (Ifges(5,4) - t322 + t324) * t235 + (-m(6) * t231 ^ 2 + Ifges(6,1) * t376 - Ifges(5,1) + Ifges(5,2) + Ifges(7,3) + Ifges(6,3) + (t336 - t325 / 0.2e1) * t229) * t236 + t265) * t235;
t333 = t1 * qJD(1);
t328 = t194 * mrSges(7,1);
t327 = t255 * mrSges(7,2);
t101 = Ifges(7,2) * t171 + t157;
t102 = Ifges(7,1) * t169 + t335;
t299 = Ifges(7,5) * t169 + Ifges(7,6) * t171;
t98 = -t171 * mrSges(7,1) + mrSges(7,2) * t169;
t4 = -t52 * t139 + t51 * t137 + t299 * t350 + t193 * t98 - (-t52 * mrSges(7,3) + t102 / 0.2e1 - t93 / 0.2e1) * t171 + (-t51 * mrSges(7,3) + t95 / 0.2e1 + t101 / 0.2e1) * t169;
t321 = t4 * qJD(1);
t207 = -mrSges(6,1) * t230 + mrSges(6,2) * t229;
t320 = t207 - mrSges(5,1);
t277 = m(5) * (t370 + t371);
t242 = t277 / 0.2e1 + (t297 * t371 + t370) * t368 + (t168 ^ 2 + t170 ^ 2 + t370) * t366;
t283 = t194 ^ 2 + t255 ^ 2;
t248 = -m(5) / 0.2e1 + t297 * t369 + t283 * t367;
t33 = -m(4) - t242 + t248;
t319 = qJD(1) * t33;
t273 = -t140 * t229 + t141 * t230;
t307 = t229 * t203;
t15 = -mrSges(4,2) - mrSges(3,3) - t201 * t303 - t168 * t139 + t170 * t137 + (mrSges(5,3) * t235 + t307) * t235 + (mrSges(5,3) * t236 + t374) * t236 - m(7) * (t168 * t51 - t170 * t52 - t193 * t236) - m(6) * (t231 * t370 + t235 * t273) - t231 * t277 + (-m(4) - m(3)) * qJ(2);
t315 = t15 * qJD(1);
t314 = t168 * t169;
t124 = t168 * t255;
t313 = t170 * t171;
t312 = t170 * t194;
t21 = t235 * mrSges(5,1) + t236 * mrSges(5,2) - t255 * t137 - t194 * t139 + mrSges(4,3) + m(7) * (-t194 * t51 - t255 * t52) + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t232 - t373;
t309 = t21 * qJD(1);
t308 = t229 * t202;
t304 = t230 * t200;
t298 = -Ifges(7,5) * t194 + Ifges(7,6) * t255;
t296 = qJD(4) * t235;
t126 = -mrSges(7,1) * t255 - mrSges(7,2) * t194;
t295 = t126 * qJD(6);
t294 = t363 / 0.2e1;
t293 = m(6) * t350;
t286 = -t306 / 0.2e1;
t285 = t126 * t349;
t284 = t376 - t227 / 0.2e1;
t279 = t297 * qJ(5);
t276 = t297 * t368;
t271 = -t168 * t194 + t170 * t255;
t270 = t169 * t255 + t171 * t194;
t268 = qJD(1) * t98 + qJD(4) * t126;
t266 = (t124 + t312) * t366;
t263 = t331 / 0.2e1 + t329 / 0.2e1;
t262 = -t328 / 0.2e1 + t327 / 0.2e1;
t260 = t307 / 0.2e1 + t201 * t352;
t245 = (-t313 / 0.2e1 - t314 / 0.2e1) * mrSges(7,3) - t168 * t362 + t139 * t358 + t98 * t349;
t8 = t245 - t262;
t259 = t8 * qJD(1);
t127 = -t327 + t328;
t241 = (-t143 * t230 - t144 * t229) * t368 + (t194 * t57 + t255 * t58) * t366 + t138 * t356 + t255 * t136 / 0.2e1 + t200 * t354 + t202 * t352;
t206 = t343 * t229;
t208 = t343 * t230;
t133 = -t206 * t347 - t208 * t346;
t134 = -t206 * t346 + t208 * t347;
t220 = -pkin(5) * t230 - pkin(4);
t244 = (t235 * t279 + t344) * t369 + (t133 * t168 - t134 * t170 - t220 * t236) * t367;
t9 = -t226 + (t207 / 0.2e1 + t127 / 0.2e1) * t236 + (-t312 / 0.2e1 - t124 / 0.2e1) * mrSges(7,3) + (mrSges(6,3) * t284 + mrSges(5,2)) * t235 + t241 + t244;
t258 = -t9 * qJD(1) + t339 * t377;
t19 = t171 * t139 + t169 * t137 + m(7) * (t169 * t52 + t171 * t51) + t373 * t236;
t257 = t19 * qJD(1) + qJD(3) * t294;
t253 = t270 * t366;
t40 = t253 + (t366 + (0.1e1 / 0.2e1 - t284) * m(6)) * t236;
t256 = t40 * qJD(1);
t254 = (-t133 * t194 - t134 * t255) * t366;
t190 = Ifges(7,4) * t194;
t128 = Ifges(7,2) * t255 - t190;
t129 = -Ifges(7,2) * t194 - t334;
t130 = -Ifges(7,1) * t194 + t334;
t131 = -Ifges(7,1) * t255 - t190;
t17 = t220 * t126 - (t130 / 0.2e1 - t129 / 0.2e1) * t255 + (-t131 / 0.2e1 - t128 / 0.2e1) * t194;
t23 = t285 - t263;
t239 = -(t102 / 0.4e1 - t93 / 0.4e1) * t255 + (-t95 / 0.4e1 - t101 / 0.4e1) * t194 - (t130 / 0.4e1 - t129 / 0.4e1 + t134 * t365) * t171 + (t131 / 0.4e1 + t128 / 0.4e1 + t133 * t365) * t169 + t133 * t137 / 0.2e1 + t134 * t361 + t193 * t126 / 0.2e1 + t220 * t98 / 0.2e1 + t235 * t298 / 0.4e1;
t246 = Ifges(7,3) * t349 - t57 * mrSges(7,1) / 0.2e1 + t58 * mrSges(7,2) / 0.2e1 - t265;
t3 = t239 + t246;
t252 = t3 * qJD(1) + t23 * qJD(3) + t17 * qJD(4);
t250 = -t192 * t367 - t264;
t38 = m(6) * (-0.1e1 + t297) * t375 + (t313 + t314 - t375) * m(7);
t238 = (t273 * t236 + (t272 - 0.2e1 * t301) * t235) * t369 + (t168 * t57 + t169 * t51 - t170 * t58 - t171 * t52 + t192 * t236 + t193 * t235) * t367 + t138 * t360 + t169 * t361 + t136 * t358 - t171 * t362;
t5 = (t99 / 0.2e1 - t183 / 0.2e1 + t260) * t236 + (-t100 / 0.2e1 - t304 / 0.2e1 + t308 / 0.2e1 - t184 / 0.2e1) * t235 + t254 + t238;
t249 = t5 * qJD(1) - t38 * qJD(3) + t340 * t377;
t240 = t251 + t273 * t368 + (t133 * t171 + t134 * t169 - t194 * t52 + t255 * t51) * t366 - t260 - t261;
t12 = (t231 * t369 + t323 / 0.2e1 + t326 / 0.2e1) * t235 + t240 + t250;
t32 = t283 * mrSges(7,3) + t280 + m(7) * (t133 * t255 - t134 * t194) + m(6) * t279;
t43 = t266 + (t367 + (-0.1e1 / 0.2e1 - t284) * m(6)) * t235;
t247 = t12 * qJD(1) + t43 * qJD(3) + t32 * qJD(4);
t48 = t50 * qJD(5) * t366;
t42 = m(7) * t350 + t235 * t276 + t266 + t293;
t41 = t236 * t276 + t253 + (m(6) + m(7)) * t349;
t34 = t242 + t248;
t24 = t285 + t263;
t14 = t243 + t264;
t11 = t231 * t293 - mrSges(6,2) * t303 / 0.2e1 + mrSges(6,1) * t286 + t240 - t250;
t10 = -t312 * t365 + mrSges(7,3) * t124 / 0.2e1 + t241 - t244 + (t207 + t127) * t349 + t350 * t280;
t7 = t245 + t262;
t6 = t201 * t302 / 0.2e1 + t200 * t303 / 0.2e1 - t203 * t305 / 0.2e1 + t202 * t286 + t254 - t238 + t374 * t350 + (t99 - t183) * t349;
t2 = t239 - t246;
t16 = [-qJD(2) * t15 + qJD(3) * t21 + qJD(4) * t1 + qJD(5) * t19 + qJD(6) * t4, t34 * qJD(3) + t10 * qJD(4) + t41 * qJD(5) + t14 * qJD(6) - t271 * t340 - t315, t34 * qJD(2) + t6 * qJD(4) + t7 * qJD(6) + t271 * t339 + t309 + t48, t333 + t10 * qJD(2) + t6 * qJD(3) + t11 * qJD(5) + t2 * qJD(6) + ((Ifges(6,1) * t229 + t336) * t352 + (Ifges(6,2) * t230 + t337) * t353 - Ifges(5,5) + (-m(6) * pkin(4) + t320) * t231) * t296 + (-mrSges(5,2) * t301 + m(7) * (t133 * t57 + t134 * t58 - t192 * t220) + t134 * t136 + t133 * t138 + t129 * t360 + t131 * t358 + pkin(4) * t183 - t192 * t127 - t194 * t92 / 0.2e1 + t94 * t355 + t220 * t99 + t165 * t353 + t164 * t351 - Ifges(5,6) * t236 + (m(6) * t272 + t304 - t308) * qJ(5) + (Ifges(6,5) * t229 - Ifges(7,5) * t255 + Ifges(6,6) * t230 - Ifges(7,6) * t194) * t236 / 0.2e1 + (-t194 * t58 + t255 * t57) * mrSges(7,3) + t272 * mrSges(6,3)) * qJD(4), t41 * qJD(2) + t11 * qJD(4) + t257, t321 + t14 * qJD(2) + t7 * qJD(3) + t2 * qJD(4) + (-mrSges(7,1) * t52 - mrSges(7,2) * t51 + t299) * qJD(6); qJD(3) * t33 + qJD(4) * t9 + qJD(5) * t40 + qJD(6) * t13 + t315, 0, t319 + t341, -t258, t256, t295 + t372; -qJD(2) * t33 - qJD(4) * t5 + qJD(6) * t8 - t309 + t48, -t319 + t341, t38 * qJD(4), t42 * qJD(5) + t24 * qJD(6) + (t127 + t320) * t296 - t249 + (t270 * mrSges(7,3) + (-mrSges(5,2) + t280) * t236 + 0.2e1 * (t133 * t169 - t134 * t171 + t220 * t235) * t366 + 0.2e1 * (t236 * t279 - t345) * t368) * qJD(4), qJD(1) * t294 + t42 * qJD(4), t24 * qJD(4) + (mrSges(7,1) * t170 - mrSges(7,2) * t168) * qJD(6) + t259; -qJD(2) * t9 + qJD(3) * t5 + qJD(5) * t12 + qJD(6) * t3 - t333, t258, t43 * qJD(5) + t23 * qJD(6) + t249, qJD(5) * t32 + qJD(6) * t17, t247 (-mrSges(7,1) * t134 - mrSges(7,2) * t133 + t298) * qJD(6) + t252; -t40 * qJD(2) - t12 * qJD(4) + t98 * qJD(6) - t257, -t256, -qJD(1) * t363 / 0.2e1 - t43 * qJD(4), -t247 + t295, 0, t268; -qJD(2) * t13 - qJD(3) * t8 - qJD(4) * t3 - qJD(5) * t98 - t321, -t372, -t23 * qJD(4) - t259, -qJD(5) * t126 - t252, -t268, 0;];
Cq  = t16;
