% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:47
% EndTime: 2019-12-05 15:18:53
% DurationCPUTime: 2.33s
% Computational Cost: add. (5071->281), mult. (14899->443), div. (0->0), fcn. (16686->12), ass. (0->165)
t195 = sin(qJ(5));
t189 = t195 ^ 2;
t198 = cos(qJ(5));
t191 = t198 ^ 2;
t308 = mrSges(5,2) + (-t189 - t191) * (m(6) * pkin(9) + mrSges(6,3));
t196 = sin(qJ(4));
t199 = cos(qJ(4));
t265 = t199 * mrSges(5,2);
t173 = t196 * mrSges(5,1) + t265;
t194 = sin(pkin(6));
t280 = cos(qJ(3));
t241 = t194 * t280;
t307 = -t173 * t241 / 0.2e1;
t187 = Ifges(6,5) * t198;
t306 = -Ifges(6,6) * t195 + t187;
t188 = Ifges(6,4) * t198;
t273 = Ifges(6,2) * t195;
t175 = t188 - t273;
t258 = sin(pkin(5));
t228 = cos(pkin(11)) * t258;
t259 = cos(pkin(6));
t260 = cos(pkin(5));
t305 = t260 * t194 + t259 * t228;
t266 = t198 * mrSges(6,2);
t270 = t195 * mrSges(6,1);
t216 = t266 / 0.2e1 + t270 / 0.2e1;
t172 = t266 + t270;
t286 = t172 / 0.2e1;
t304 = t216 + t286;
t181 = t196 * pkin(4) - pkin(9) * t199;
t252 = t195 * t196;
t133 = pkin(8) * t252 + t181 * t198;
t250 = t196 * t198;
t134 = -pkin(8) * t250 + t181 * t195;
t303 = -t133 * t195 + t134 * t198;
t275 = Ifges(6,4) * t195;
t174 = Ifges(6,2) * t198 + t275;
t176 = Ifges(6,1) * t195 + t188;
t177 = Ifges(6,1) * t198 - t275;
t244 = Ifges(6,1) / 0.4e1 - Ifges(6,2) / 0.4e1;
t302 = (-t175 / 0.4e1 - t176 / 0.4e1 - t188 - t244 * t195) * t195 + pkin(8) * t286 + (-t189 / 0.2e1 - t191 / 0.2e1) * pkin(9) * mrSges(6,3) + (-t174 / 0.4e1 + t177 / 0.4e1 + t244 * t198) * t198;
t227 = mrSges(6,1) * t198 - mrSges(6,2) * t195;
t301 = -m(6) * pkin(4) - mrSges(5,1) - t227;
t168 = -pkin(4) * t199 - t196 * pkin(9) - pkin(3);
t251 = t195 * t199;
t130 = -pkin(8) * t251 + t198 * t168;
t249 = t198 * t199;
t131 = pkin(8) * t249 + t195 * t168;
t155 = t172 * t199;
t247 = mrSges(6,3) * t252;
t162 = mrSges(6,2) * t199 - t247;
t246 = mrSges(6,3) * t250;
t164 = -mrSges(6,1) * t199 - t246;
t277 = pkin(8) * t199;
t282 = -t198 / 0.2e1;
t283 = t195 / 0.2e1;
t297 = m(6) / 0.2e1;
t299 = (t130 * t195 - t131 * t198 + t277) * t297 + t162 * t282 + t164 * t283 + t155 / 0.2e1;
t192 = t199 ^ 2;
t298 = m(5) / 0.2e1;
t296 = mrSges(6,1) / 0.2e1;
t295 = -mrSges(6,2) / 0.2e1;
t153 = t227 * t196;
t294 = -t153 / 0.2e1;
t154 = t172 * t196;
t293 = t154 / 0.2e1;
t291 = -t162 / 0.2e1;
t163 = -t196 * mrSges(6,2) - mrSges(6,3) * t251;
t290 = t163 / 0.2e1;
t289 = -t164 / 0.2e1;
t165 = t196 * mrSges(6,1) - mrSges(6,3) * t249;
t288 = t165 / 0.2e1;
t287 = -t227 / 0.2e1;
t285 = t173 / 0.2e1;
t284 = -t195 / 0.2e1;
t281 = t198 / 0.2e1;
t279 = pkin(8) * t192;
t278 = pkin(8) * t196;
t276 = mrSges(6,3) * t196;
t274 = Ifges(6,5) * t199;
t271 = Ifges(6,6) * t199;
t205 = -t194 * t228 + t259 * t260;
t197 = sin(qJ(3));
t236 = sin(pkin(11)) * t258;
t99 = t305 * t197 + t280 * t236;
t70 = t196 * t99 - t199 * t205;
t268 = t196 * t70;
t98 = t197 * t236 - t305 * t280;
t267 = t196 * t98;
t71 = t196 * t205 + t99 * t199;
t47 = -t195 * t71 + t198 * t98;
t264 = t47 * t195;
t48 = t195 * t98 + t198 * t71;
t263 = t48 * t198;
t55 = t99 * t198 + t251 * t98;
t262 = t55 * t195;
t56 = t99 * t195 - t249 * t98;
t261 = t56 * t198;
t253 = t194 * t197;
t151 = t196 * t259 + t199 * t253;
t100 = -t195 * t151 - t198 * t241;
t257 = t100 * t195;
t101 = t151 * t198 - t195 * t241;
t256 = t101 * t198;
t245 = pkin(4) * t294;
t243 = mrSges(6,3) * t284;
t242 = mrSges(6,3) * t281;
t240 = t196 * t280;
t239 = t199 * t280;
t190 = t196 ^ 2;
t232 = t190 * t241;
t231 = t194 * t240;
t150 = t196 * t253 - t199 * t259;
t230 = t150 * t240;
t226 = -Ifges(6,5) * t195 - Ifges(6,6) * t198;
t128 = (-t195 * t239 + t197 * t198) * t194;
t129 = (t195 * t197 + t198 * t239) * t194;
t6 = (t100 * t55 + t101 * t56 + t128 * t47 + t129 * t48 + (-t150 * t98 + t241 * t70) * t196) * t297 + ((-t150 * t196 - t151 * t199) * t98 + (t197 * t98 + t239 * t71 + t240 * t70 - t280 * t99) * t194) * t298;
t8 = m(5) * (-t199 * t71 - t268 + t99) * t98 + (-t268 * t98 + t47 * t55 + t48 * t56) * m(6);
t225 = t8 * qJD(1) + t6 * qJD(2);
t219 = -t263 + t71 + t264;
t11 = m(6) * t219 * t70;
t215 = t151 - t256 + t257;
t7 = (t150 * t219 + t215 * t70) * t297;
t224 = t11 * qJD(1) + t7 * qJD(2);
t32 = m(6) * t215 * t150;
t223 = t7 * qJD(1) + t32 * qJD(2);
t35 = m(5) * (t151 * t239 - t197 * t241 + t230) * t194 + m(6) * (t100 * t128 + t101 * t129 + t194 * t230);
t222 = t6 * qJD(1) + t35 * qJD(2);
t221 = t295 * t56 + t296 * t55;
t220 = -Ifges(5,4) + t306;
t218 = t128 * t296 + t129 * t295;
t217 = -t133 * mrSges(6,1) / 0.2e1 + t134 * mrSges(6,2) / 0.2e1;
t18 = -t154 * t277 - t155 * t278 - t133 * t164 - t130 * t165 - m(6) * (t130 * t133 + t131 * t134) - t134 * t162 - t131 * t163 + pkin(3) * t173 + t220 * t192 + (-t220 * t196 + (-Ifges(6,1) * t191 + Ifges(6,3) - m(6) * pkin(8) ^ 2 - Ifges(5,1) + Ifges(5,2) + (-t273 + 0.2e1 * t188) * t195) * t199) * t196;
t202 = t299 * t70 + (t133 * t47 + t134 * t48 + t278 * t71) * t297 + t47 * t288 + t48 * t290 + t71 * t293;
t206 = m(6) * (pkin(4) * t267 + (t261 - t262) * pkin(9));
t3 = (-t261 / 0.2e1 + t262 / 0.2e1) * mrSges(6,3) + (t285 - t265 / 0.2e1 + (t287 - mrSges(5,1) / 0.2e1) * t196) * t98 - t206 / 0.2e1 + t202;
t200 = (t133 * t100 + t134 * t101 + t151 * t278) * t297 + t100 * t288 + t101 * t290 + t151 * t293 + t307 + t299 * t150;
t201 = (-pkin(4) * t231 + (-t128 * t195 + t129 * t198) * pkin(9)) * t297 + t128 * t243 + t129 * t242 + t231 * t287 + t307;
t9 = -t200 + t201;
t212 = t3 * qJD(1) - t9 * qJD(2) - t18 * qJD(3);
t211 = -t172 / 0.2e1 + t216;
t209 = t100 * t291 + t101 * t164 / 0.2e1 + t150 * t294;
t20 = (t256 / 0.2e1 - t257 / 0.2e1) * t276 + t209 + t218;
t25 = -t153 * t278 + t131 * t164 + ((-Ifges(6,4) * t252 - t274) * t195 + (Ifges(6,4) * t250 + t131 * mrSges(6,3) - t271 + (Ifges(6,1) - Ifges(6,2)) * t252) * t198) * t196 + (-t247 - t162) * t130;
t204 = (-t263 / 0.2e1 + t264 / 0.2e1) * t276 + t47 * t162 / 0.2e1 + t48 * t289 + t70 * t153 / 0.2e1;
t4 = t204 - t221;
t210 = t4 * qJD(1) - t20 * qJD(2) - t25 * qJD(3);
t12 = t211 * t70;
t22 = t245 + (t162 * t284 + t164 * t282) * pkin(9) - t306 * t199 + (-Ifges(6,3) / 0.2e1 + t302) * t196 + t217;
t29 = t211 * t150;
t65 = -pkin(4) * t172 + (t176 / 0.2e1 + t175 / 0.2e1) * t198 + (t177 / 0.2e1 - t174 / 0.2e1) * t195;
t208 = t12 * qJD(1) + t29 * qJD(2) - t22 * qJD(3) - t65 * qJD(4);
t171 = -t199 * mrSges(5,1) + t196 * mrSges(5,2);
t169 = pkin(8) * t232;
t84 = t190 * pkin(8) * t98;
t30 = t304 * t150;
t23 = t245 - t199 * t187 / 0.4e1 - Ifges(6,6) * t251 / 0.2e1 + Ifges(6,5) * t249 / 0.2e1 + (-t274 / 0.4e1 + pkin(9) * t289) * t198 + (t271 / 0.2e1 + pkin(9) * t291) * t195 - t217 + (Ifges(6,3) / 0.2e1 + t302) * t196;
t21 = -t101 * t246 / 0.2e1 + t100 * t247 / 0.2e1 - t209 + t218;
t13 = t304 * t70;
t10 = t200 + t201;
t5 = t204 + t221;
t2 = t206 / 0.2e1 + t56 * t242 + t55 * t243 + t202 + (t285 + t265 / 0.2e1) * t98 + (t227 / 0.2e1 + mrSges(5,1) / 0.2e1) * t267;
t1 = t6 * qJD(3) + t7 * qJD(4);
t14 = [t8 * qJD(3) + t11 * qJD(4), t1, t2 * qJD(4) + t5 * qJD(5) + t225 + (-t99 * mrSges(4,1) + t56 * t162 + t55 * t164 + t99 * t171 + (-t196 * t154 + mrSges(4,2) + (-t190 - t192) * mrSges(5,3)) * t98 + 0.2e1 * (t130 * t55 + t131 * t56 - t84) * t297 + 0.2e1 * (-pkin(3) * t99 - t279 * t98 - t84) * t298) * qJD(3), t2 * qJD(3) + (t301 * t71 + t308 * t70) * qJD(4) + t13 * qJD(5) + t224, t5 * qJD(3) + t13 * qJD(4) + (-mrSges(6,1) * t48 - mrSges(6,2) * t47) * qJD(5); t1, t35 * qJD(3) + t32 * qJD(4), (-mrSges(4,2) * t241 + m(5) * (t169 + (-pkin(3) * t197 + t279 * t280) * t194) + t154 * t231 + m(6) * (t128 * t130 + t129 * t131 + t169) + t129 * t162 + t128 * t164 + (-mrSges(4,1) + t171) * t253 + (t192 * t241 + t232) * mrSges(5,3)) * qJD(3) + t10 * qJD(4) + t21 * qJD(5) + t222, t10 * qJD(3) + (t308 * t150 + t301 * t151) * qJD(4) + t30 * qJD(5) + t223, t21 * qJD(3) + t30 * qJD(4) + (-mrSges(6,1) * t101 - mrSges(6,2) * t100) * qJD(5); qJD(4) * t3 + qJD(5) * t4 - t225, -qJD(4) * t9 - qJD(5) * t20 - t222, -qJD(4) * t18 - qJD(5) * t25, t23 * qJD(5) + t212 + (t303 * mrSges(6,3) - pkin(4) * t155 + (m(6) * t303 + t198 * t163 - t195 * t165) * pkin(9) + (t174 * t284 + t177 * t283 + Ifges(5,5) + t301 * pkin(8) + (t176 + t175) * t281) * t199 + (mrSges(5,2) * pkin(8) - Ifges(5,6) - t226) * t196) * qJD(4), t23 * qJD(4) + (-mrSges(6,1) * t131 - mrSges(6,2) * t130 + t196 * t226) * qJD(5) + t210; -qJD(3) * t3 - qJD(5) * t12 - t224, qJD(3) * t9 - qJD(5) * t29 - t223, qJD(5) * t22 - t212, t65 * qJD(5), (-pkin(9) * t227 + t306) * qJD(5) - t208; -t4 * qJD(3) + t12 * qJD(4), t20 * qJD(3) + t29 * qJD(4), -qJD(4) * t22 - t210, t208, 0;];
Cq = t14;
