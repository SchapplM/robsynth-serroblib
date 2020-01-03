% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP10
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:58
% EndTime: 2019-12-31 18:51:04
% DurationCPUTime: 3.04s
% Computational Cost: add. (6317->351), mult. (13030->483), div. (0->0), fcn. (13499->6), ass. (0->181)
t335 = -m(6) / 0.2e1;
t210 = cos(qJ(4));
t327 = Ifges(5,6) + Ifges(6,6);
t334 = t327 * t210;
t208 = sin(qJ(4));
t302 = -t208 / 0.2e1;
t298 = t210 / 0.2e1;
t328 = Ifges(5,5) + Ifges(6,5);
t206 = sin(pkin(8));
t207 = cos(pkin(8));
t209 = sin(qJ(3));
t297 = cos(qJ(3));
t182 = t206 * t209 - t297 * t207;
t184 = t297 * t206 + t209 * t207;
t200 = Ifges(6,4) * t210;
t224 = -Ifges(6,2) * t208 + t200;
t75 = Ifges(6,6) * t182 + t224 * t184;
t201 = Ifges(5,4) * t210;
t225 = -Ifges(5,2) * t208 + t201;
t77 = Ifges(5,6) * t182 + t225 * t184;
t333 = t77 + t75;
t286 = -qJ(5) - pkin(7);
t187 = t286 * t208;
t190 = t286 * t210;
t266 = t184 * t208;
t241 = -pkin(2) * t207 - pkin(1);
t113 = pkin(3) * t182 - pkin(7) * t184 + t241;
t287 = pkin(6) + qJ(2);
t186 = t287 * t207;
t234 = t287 * t206;
t131 = t297 * t186 - t209 * t234;
t47 = t208 * t113 + t131 * t210;
t38 = -qJ(5) * t266 + t47;
t273 = t210 * t38;
t289 = t182 * pkin(4);
t262 = t210 * t184;
t46 = t210 * t113 - t208 * t131;
t37 = -qJ(5) * t262 + t46;
t30 = t37 + t289;
t223 = t208 * t30 - t273;
t252 = mrSges(6,3) * t266;
t118 = -mrSges(6,2) * t182 - t252;
t122 = t182 * mrSges(6,1) - mrSges(6,3) * t262;
t260 = t118 * t298 + t122 * t302;
t332 = ((-t187 * t210 + t190 * t208) * t184 - t223) * t335 - t260;
t308 = -t182 / 0.2e1;
t329 = t184 / 0.2e1;
t326 = Ifges(5,3) + Ifges(6,3);
t325 = -t30 + t37;
t318 = m(6) * pkin(4);
t243 = mrSges(6,1) + t318;
t324 = -t187 * t208 - t190 * t210;
t198 = Ifges(6,5) * t210;
t199 = Ifges(5,5) * t210;
t323 = t199 + t198;
t194 = Ifges(6,1) * t208 + t200;
t195 = Ifges(5,1) * t208 + t201;
t253 = mrSges(5,3) * t266;
t275 = t182 * mrSges(5,2);
t119 = -t253 - t275;
t310 = -t119 / 0.2e1;
t130 = t186 * t209 + t297 * t234;
t83 = pkin(4) * t266 + t130;
t313 = t83 / 0.2e1;
t322 = mrSges(6,1) * t313 + pkin(7) * t310;
t246 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t321 = (t224 + t225) * t182 / 0.2e1 + (-t246 - t327 / 0.2e1) * t184;
t320 = m(5) / 0.2e1;
t319 = m(6) / 0.2e1;
t317 = -mrSges(5,1) / 0.2e1;
t316 = mrSges(5,2) / 0.2e1;
t315 = mrSges(6,2) / 0.2e1;
t267 = t182 * t210;
t290 = pkin(7) * t182;
t293 = pkin(3) * t184;
t126 = t290 + t293;
t54 = t210 * t126 + t208 * t130;
t33 = t184 * pkin(4) + qJ(5) * t267 + t54;
t314 = -t33 / 0.2e1;
t312 = m(6) * t83;
t284 = mrSges(6,2) * t210;
t228 = t208 * mrSges(6,1) + t284;
t111 = t228 * t184;
t311 = t111 / 0.2e1;
t120 = t184 * mrSges(6,1) + mrSges(6,3) * t267;
t309 = -t120 / 0.2e1;
t307 = t182 / 0.4e1;
t305 = t187 / 0.2e1;
t188 = -mrSges(6,1) * t210 + t208 * mrSges(6,2);
t304 = t188 / 0.2e1;
t303 = t190 / 0.2e1;
t300 = t208 / 0.2e1;
t295 = m(6) * t184;
t291 = pkin(4) * t210;
t197 = -pkin(3) - t291;
t294 = m(6) * t197;
t292 = pkin(4) * t208;
t285 = mrSges(5,1) * t210;
t281 = Ifges(5,4) * t208;
t280 = Ifges(6,4) * t208;
t274 = t210 * mrSges(5,2);
t191 = t208 * mrSges(5,1) + t274;
t112 = t191 * t184;
t276 = t182 * mrSges(5,1);
t123 = -mrSges(5,3) * t262 + t276;
t269 = t130 * t184;
t272 = t210 * t47;
t5 = (mrSges(4,3) * t184 + t111 + t112) * t184 + (mrSges(4,3) * t182 + (-t118 - t119) * t210 + (t122 + t123) * t208) * t182 + m(5) * (t269 + (t208 * t46 - t272) * t182) + m(6) * (t223 * t182 + t83 * t184) + m(4) * (-t131 * t182 + t269) + (m(3) * qJ(2) + mrSges(3,3)) * (t206 ^ 2 + t207 ^ 2);
t279 = qJD(1) * t5;
t268 = t182 * t208;
t116 = -mrSges(6,2) * t184 + mrSges(6,3) * t268;
t117 = -mrSges(5,2) * t184 + mrSges(5,3) * t268;
t121 = t184 * mrSges(5,1) + mrSges(5,3) * t267;
t169 = t182 * mrSges(4,2);
t204 = t208 ^ 2;
t205 = t210 ^ 2;
t235 = -t205 / 0.2e1 - t204 / 0.2e1;
t231 = t235 * mrSges(5,3);
t257 = t204 + t205;
t212 = (t235 * mrSges(6,3) + t231) * t182 + (-t257 * t290 - t293) * t320 + (-t324 * t182 + t197 * t184) * t319;
t55 = t208 * t126 - t210 * t130;
t41 = qJ(5) * t268 + t55;
t214 = -m(5) * (t208 * t55 + t210 * t54) / 0.2e1 + (t208 * t41 + t210 * t33) * t335;
t229 = -t208 * mrSges(5,2) + t285;
t238 = -t229 / 0.2e1 + t304;
t8 = t169 + (-t121 / 0.2e1 + t309) * t210 + (-t117 / 0.2e1 - t116 / 0.2e1) * t208 + (-mrSges(4,1) + t238) * t184 + t212 + t214;
t278 = qJD(1) * t8;
t109 = t228 * t182;
t110 = t191 * t182;
t226 = Ifges(6,1) * t210 - t280;
t227 = Ifges(5,1) * t210 - t281;
t247 = Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t218 = t247 * t184 + t328 * t329 + (t226 + t227) * t308;
t79 = Ifges(6,5) * t182 + t226 * t184;
t81 = Ifges(5,5) * t182 + t227 * t184;
t84 = -pkin(4) * t268 + t131;
t1 = -t83 * t109 + t84 * t111 + t38 * t116 + t47 * t117 + t41 * t118 + t55 * t119 + t30 * t120 + t46 * t121 + t33 * t122 + t54 * t123 - t130 * t110 + t131 * t112 - t241 * t169 + m(6) * (t30 * t33 + t38 * t41 + t83 * t84) + m(5) * (t130 * t131 + t46 * t54 + t47 * t55) + (t241 * mrSges(4,1) - Ifges(4,4) * t184 + t321 * t208 + t218 * t210) * t184 + (Ifges(4,4) * t182 + (-t79 / 0.2e1 - t81 / 0.2e1 - t247 * t182) * t210 + (t75 / 0.2e1 + t77 / 0.2e1 + t246 * t182) * t208 + (-Ifges(4,1) + Ifges(4,2) + t326) * t184) * t182;
t277 = t1 * qJD(1);
t192 = Ifges(6,2) * t210 + t280;
t193 = Ifges(5,2) * t210 + t281;
t165 = mrSges(6,1) * t262;
t232 = -mrSges(6,2) * t266 + t165;
t242 = m(6) * t325;
t254 = pkin(4) * t262;
t4 = t47 * t123 - t229 * t269 - t111 * t254 - t83 * t232 - t37 * t118 - t30 * t252 + (mrSges(5,3) * t272 + mrSges(6,3) * t273 - t291 * t312 + ((t193 + t192) * t302 - (-t195 - t194) * t210 / 0.2e1) * t184 + (-t328 * t208 - t334) * t308 + (t79 + t81) * t300 + t333 * t298) * t184 + (-t119 - t253) * t46 + (t122 - t242) * t38;
t271 = t4 * qJD(1);
t15 = (t208 * t118 - m(6) * (-t208 * t38 - t30 * t210) + t210 * t122) * t184;
t270 = qJD(1) * t15;
t59 = 0.2e1 * (-t204 / 0.4e1 - t205 / 0.4e1 - 0.1e1 / 0.4e1) * t295;
t261 = t59 * qJD(1);
t240 = t268 / 0.2e1;
t259 = mrSges(6,1) * t240 + t267 * t315;
t256 = qJD(4) * t208;
t255 = t318 / 0.2e1;
t245 = t37 / 0.2e1 - t30 / 0.2e1;
t239 = -t267 / 0.2e1;
t237 = t192 / 0.2e1 + t193 / 0.2e1;
t236 = -t194 / 0.2e1 - t195 / 0.2e1;
t233 = t242 / 0.2e1;
t230 = -mrSges(6,1) * t268 / 0.2e1 + t84 * t319 + mrSges(6,2) * t239;
t12 = t230 + t332;
t93 = m(6) * t324 + t257 * mrSges(6,3);
t222 = -qJD(1) * t12 + qJD(3) * t93;
t171 = t243 * t208 + t284;
t66 = -m(6) * t254 - t232;
t221 = qJD(1) * t66 - qJD(3) * t171;
t220 = Ifges(5,1) / 0.4e1 + Ifges(6,1) / 0.4e1 - Ifges(5,2) / 0.4e1 - Ifges(6,2) / 0.4e1;
t6 = (t275 / 0.2e1 + t310 - t118 / 0.2e1) * t210 + (t276 / 0.2e1 + t123 / 0.2e1 + t122 / 0.2e1 + (t289 / 0.2e1 - t245) * m(6)) * t208 + t259;
t219 = t6 * qJD(1);
t19 = pkin(3) * t191 - t197 * t228 - t188 * t292 + (-t200 / 0.2e1 - t201 / 0.2e1 + t236) * t210 + (-pkin(4) * t294 + (Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t208 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t210 + t237) * t208;
t211 = pkin(7) * t231 + (mrSges(6,3) * t305 - t197 * mrSges(6,2) / 0.2e1 + pkin(3) * t316 - t195 / 0.4e1 - t194 / 0.4e1 - t201 / 0.4e1 - t200 / 0.4e1 - t220 * t208) * t208 + (mrSges(6,3) * t303 - t193 / 0.4e1 - t192 / 0.4e1 + pkin(3) * t317 + t220 * t210 + (-0.3e1 / 0.4e1 * Ifges(5,4) - 0.3e1 / 0.4e1 * Ifges(6,4)) * t208 + (t304 + t294 / 0.2e1) * pkin(4)) * t210;
t213 = (-pkin(7) * t123 / 0.2e1 + t81 / 0.4e1 + t79 / 0.4e1 + mrSges(6,2) * t313 + t245 * mrSges(6,3)) * t210 + t130 * t191 / 0.2e1 + t118 * t305 + t197 * t165 / 0.2e1;
t216 = mrSges(6,1) * t314 + t41 * t315 + t55 * t316 + t54 * t317;
t3 = (m(6) * t314 + t309) * pkin(4) - (t233 - t122 / 0.2e1) * t190 + (-t77 / 0.4e1 - t75 / 0.4e1 + (t312 / 0.2e1 + t311) * pkin(4) + t322) * t208 + (t199 / 0.4e1 + t198 / 0.4e1 + t247 * t210 + (-0.3e1 / 0.4e1 * Ifges(5,6) - 0.3e1 / 0.4e1 * Ifges(6,6)) * t208) * t182 + (-Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + t211) * t184 + t213 + t216;
t217 = -t3 * qJD(1) + t19 * qJD(3);
t58 = (-t257 + 0.1e1) * t295 / 0.2e1;
t13 = t230 - t332;
t9 = t238 * t184 + t212 - t214 + (t116 + t117) * t300 + (t120 + t121) * t298;
t7 = t119 * t298 + t123 * t302 + t208 * t233 + (t274 / 0.2e1 + (mrSges(5,1) / 0.2e1 + t255) * t208) * t182 + t259 + t260;
t2 = t211 * t184 + t213 + (-t325 * t190 + t83 * t292) * t319 + pkin(4) * t120 / 0.2e1 - t216 + t122 * t303 + t33 * t255 + t292 * t311 + t323 * t307 + t326 * t329 + t327 * t240 + t328 * t239 + (-t333 / 0.4e1 - t327 * t307 + t322) * t208;
t10 = [qJD(2) * t5 + qJD(3) * t1 - qJD(4) * t4 - qJD(5) * t15, qJD(3) * t9 + qJD(4) * t7 + qJD(5) * t58 + t279, t9 * qJD(2) + t2 * qJD(4) + t13 * qJD(5) + t277 + (-t131 * mrSges(4,1) + t130 * mrSges(4,2) - Ifges(4,6) * t184 + pkin(3) * t110 - t197 * t109 - t190 * t116 + t187 * t120 - t131 * t229 + t84 * t188 + (t55 * mrSges(5,3) + t41 * mrSges(6,3) + pkin(7) * t117 - t321) * t210 + (-t54 * mrSges(5,3) - t33 * mrSges(6,3) - pkin(7) * t121 + t218) * t208 + 0.2e1 * (-pkin(3) * t131 + (-t208 * t54 + t210 * t55) * pkin(7)) * t320 + 0.2e1 * (t187 * t33 - t190 * t41 + t197 * t84) * t319 + (t237 * t208 + t236 * t210 - Ifges(4,5)) * t182) * qJD(3), t7 * qJD(2) + t2 * qJD(3) - t271 + (-t47 * mrSges(5,1) - t46 * mrSges(5,2) - t37 * mrSges(6,2) + (-t334 + (mrSges(6,3) * pkin(4) - t328) * t208) * t184 - t243 * t38) * qJD(4), qJD(2) * t58 + qJD(3) * t13 - t270; -qJD(3) * t8 - qJD(4) * t6 + qJD(5) * t59 - t279, 0, -t278, (-mrSges(5,2) - mrSges(6,2)) * qJD(4) * t210 + (-mrSges(5,1) - t243) * t256 - t219, t261; qJD(2) * t8 + qJD(4) * t3 - qJD(5) * t12 - t277, t278, -qJD(4) * t19 + qJD(5) * t93, (-t187 * mrSges(6,2) - mrSges(6,3) * t291 - pkin(7) * t285 + t243 * t190 + t323) * qJD(4) + (mrSges(5,2) * pkin(7) - t327) * t256 - t217, t222; qJD(2) * t6 - qJD(3) * t3 + qJD(5) * t66 + t271, t219, -t171 * qJD(5) + t217, 0, t221; -qJD(2) * t59 + qJD(3) * t12 - qJD(4) * t66 + t270, -t261, qJD(4) * t171 - t222, -t221, 0;];
Cq = t10;
