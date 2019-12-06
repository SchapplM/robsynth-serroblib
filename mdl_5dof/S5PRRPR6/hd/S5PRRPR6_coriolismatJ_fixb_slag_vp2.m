% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:15
% EndTime: 2019-12-05 16:30:25
% DurationCPUTime: 3.37s
% Computational Cost: add. (5870->374), mult. (14565->568), div. (0->0), fcn. (15234->10), ass. (0->199)
t214 = sin(pkin(10));
t295 = pkin(8) + qJ(4);
t190 = t295 * t214;
t216 = cos(pkin(10));
t192 = t295 * t216;
t217 = sin(qJ(5));
t297 = cos(qJ(5));
t127 = -t297 * t190 - t217 * t192;
t128 = -t217 * t190 + t297 * t192;
t205 = -pkin(4) * t216 - pkin(3);
t218 = sin(qJ(3));
t215 = sin(pkin(5));
t221 = cos(qJ(2));
t263 = t215 * t221;
t255 = t218 * t263;
t219 = sin(qJ(2));
t220 = cos(qJ(3));
t259 = t220 * t221;
t135 = (t214 * t219 + t216 * t259) * t215;
t271 = t135 * t216;
t134 = (-t214 * t259 + t216 * t219) * t215;
t272 = t134 * t214;
t319 = -m(6) / 0.2e1;
t321 = -m(5) / 0.2e1;
t75 = t297 * t134 - t217 * t135;
t76 = t217 * t134 + t297 * t135;
t330 = -(t271 / 0.2e1 - t272 / 0.2e1) * mrSges(5,3) + (-pkin(3) * t255 + (t271 - t272) * qJ(4)) * t321 + (t127 * t75 + t128 * t76 + t205 * t255) * t319;
t318 = m(6) / 0.2e1;
t329 = 0.2e1 * t318;
t320 = m(5) / 0.2e1;
t328 = 0.2e1 * t320;
t315 = mrSges(6,3) / 0.2e1;
t233 = t297 * t214 + t217 * t216;
t156 = t233 * t218;
t327 = t156 * mrSges(6,3);
t211 = t216 ^ 2;
t256 = t214 ^ 2 + t211;
t326 = t256 * mrSges(5,3);
t325 = qJ(4) * t256;
t232 = t217 * t214 - t297 * t216;
t172 = Ifges(6,4) * t232;
t120 = Ifges(6,1) * t233 - t172;
t324 = -Ifges(6,2) * t233 + t120 - t172;
t194 = pkin(3) * t218 - qJ(4) * t220;
t266 = t214 * t218;
t144 = pkin(7) * t266 + t216 * t194;
t262 = t216 * t218;
t145 = -pkin(7) * t262 + t214 * t194;
t323 = -t144 * t214 + t145 * t216;
t281 = t220 * mrSges(4,2);
t195 = t218 * mrSges(4,1) + t281;
t317 = -mrSges(6,1) / 0.2e1;
t316 = mrSges(6,2) / 0.2e1;
t129 = mrSges(6,2) * t220 - t327;
t314 = t129 / 0.2e1;
t155 = t232 * t218;
t131 = -mrSges(6,1) * t220 + t155 * mrSges(6,3);
t313 = t131 / 0.2e1;
t312 = -t155 / 0.2e1;
t311 = -t156 / 0.2e1;
t157 = t233 * t220;
t309 = -t157 / 0.2e1;
t158 = t232 * t220;
t308 = -t158 / 0.2e1;
t264 = t215 * t219;
t275 = cos(pkin(5));
t163 = t218 * t264 - t275 * t220;
t307 = t163 / 0.2e1;
t305 = -t232 / 0.2e1;
t303 = t233 / 0.2e1;
t302 = -t233 / 0.2e1;
t301 = -t214 / 0.2e1;
t300 = t214 / 0.2e1;
t299 = t216 / 0.2e1;
t298 = t218 / 0.2e1;
t296 = t220 * pkin(7);
t294 = Ifges(5,4) * t214;
t293 = Ifges(5,4) * t216;
t292 = Ifges(6,4) * t155;
t291 = Ifges(6,4) * t233;
t290 = Ifges(5,5) * t216;
t289 = Ifges(5,6) * t214;
t288 = t157 * mrSges(6,1);
t287 = t158 * mrSges(6,2);
t286 = t214 * mrSges(5,1);
t285 = t214 * Ifges(5,2);
t284 = t216 * mrSges(5,2);
t282 = t218 * mrSges(4,2);
t164 = t275 * t218 + t220 * t264;
t121 = -t164 * t214 - t216 * t263;
t122 = t164 * t216 - t214 * t263;
t57 = t297 * t121 - t217 * t122;
t280 = t57 * t156;
t58 = t217 * t121 + t297 * t122;
t279 = t58 * t155;
t278 = t75 * t233;
t277 = t76 * t232;
t191 = -mrSges(5,1) * t216 + mrSges(5,2) * t214;
t276 = t191 - mrSges(4,1);
t274 = t121 * t216;
t273 = t122 * t214;
t115 = t163 * t164;
t241 = t214 * t121 - t216 * t122;
t92 = t233 * t163;
t93 = t232 * t163;
t15 = m(5) * (t241 * t163 + t115) + m(6) * (t57 * t92 + t58 * t93 + t115);
t268 = t15 * qJD(1);
t133 = t163 * t255;
t16 = m(5) * (t121 * t134 + t122 * t135 + t133) + m(6) * (t57 * t75 + t58 * t76 + t133) + m(4) * (t163 * t218 + t164 * t220 - t264) * t263;
t267 = t16 * qJD(1);
t265 = t214 * t220;
t261 = t216 * t220;
t258 = -Ifges(6,5) * t156 + Ifges(6,6) * t155;
t257 = -Ifges(6,5) * t232 - Ifges(6,6) * t233;
t189 = -pkin(3) * t220 - qJ(4) * t218 - pkin(2);
t137 = pkin(7) * t261 + t214 * t189;
t254 = pkin(4) * t214 + pkin(7);
t117 = mrSges(6,1) * t233 - mrSges(6,2) * t232;
t252 = t117 * t307;
t250 = t263 / 0.2e1;
t247 = 0.2e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * t164;
t246 = t284 + t286;
t245 = -Ifges(6,1) * t232 - t291;
t118 = mrSges(6,1) * t232 + mrSges(6,2) * t233;
t100 = -t287 + t288;
t130 = -mrSges(6,2) * t218 - mrSges(6,3) * t157;
t132 = mrSges(6,1) * t218 + mrSges(6,3) * t158;
t165 = t246 * t218;
t166 = t246 * t220;
t180 = -mrSges(5,2) * t218 - mrSges(5,3) * t265;
t182 = mrSges(5,1) * t218 - mrSges(5,3) * t261;
t183 = t254 * t218;
t184 = t254 * t220;
t179 = t220 * mrSges(5,2) - mrSges(5,3) * t266;
t181 = -t220 * mrSges(5,1) - mrSges(5,3) * t262;
t234 = t181 * t300 - t216 * t179 / 0.2e1;
t174 = t216 * t189;
t136 = -pkin(7) * t265 + t174;
t240 = -t136 * t214 + t137 * t216;
t114 = -pkin(8) * t262 + t174 + (-pkin(7) * t214 - pkin(4)) * t220;
t123 = -pkin(8) * t266 + t137;
t50 = t297 * t114 - t217 * t123;
t51 = t217 * t114 + t297 * t123;
t116 = pkin(4) * t218 - pkin(8) * t261 + t144;
t124 = -pkin(8) * t265 + t145;
t53 = t297 * t116 - t217 * t124;
t54 = t217 * t116 + t297 * t124;
t99 = mrSges(6,1) * t156 - mrSges(6,2) * t155;
t222 = (t99 / 0.2e1 + t165 / 0.2e1) * t164 + (t100 / 0.2e1 + t166 / 0.2e1 + t234) * t163 + (pkin(7) * t164 * t218 + t121 * t144 + t122 * t145 + (-t240 + t296) * t163) * t320 + (t163 * t184 + t164 * t183 + t50 * t92 + t51 * t93 + t53 * t57 + t54 * t58) * t318 + t121 * t182 / 0.2e1 + t122 * t180 / 0.2e1 + t57 * t132 / 0.2e1 + t58 * t130 / 0.2e1 + t92 * t313 + t93 * t314;
t2 = (t281 / 0.2e1 - t195 / 0.2e1 + (mrSges(4,1) / 0.2e1 - t191 / 0.2e1 - t118 / 0.2e1) * t218) * t263 + (t278 / 0.2e1 + t277 / 0.2e1) * mrSges(6,3) + t222 + t330;
t153 = Ifges(5,6) * t218 + (-t285 + t293) * t220;
t154 = Ifges(5,5) * t218 + (t216 * Ifges(5,1) - t294) * t220;
t236 = Ifges(6,5) * t158 / 0.2e1 + Ifges(6,6) * t157 / 0.2e1;
t88 = -Ifges(6,2) * t156 - t220 * Ifges(6,6) - t292;
t89 = -Ifges(6,4) * t158 - Ifges(6,2) * t157 + Ifges(6,6) * t218;
t152 = Ifges(6,4) * t156;
t90 = -Ifges(6,1) * t155 - t220 * Ifges(6,5) - t152;
t91 = -Ifges(6,1) * t158 - Ifges(6,4) * t157 + Ifges(6,5) * t218;
t3 = t54 * t129 + t51 * t130 + t53 * t131 + t50 * t132 + t91 * t312 + t89 * t311 + t88 * t309 + t90 * t308 + t145 * t179 + t137 * t180 + t144 * t181 + t136 * t182 + t183 * t100 + t184 * t99 - pkin(2) * t195 + m(5) * (t136 * t144 + t137 * t145) + m(6) * (t183 * t184 + t50 * t53 + t51 * t54) + (t154 * t299 + t153 * t301 + pkin(7) * t166 + Ifges(6,5) * t312 + Ifges(6,6) * t311 + (-Ifges(4,4) + t290 / 0.2e1 - t289 / 0.2e1) * t218) * t218 + (pkin(7) * t165 + (m(5) * pkin(7) ^ 2 + t211 * Ifges(5,1) / 0.2e1 + Ifges(4,1) - Ifges(4,2) - Ifges(6,3) - Ifges(5,3) + (-t293 + t285 / 0.2e1) * t214) * t218 + t236 + (Ifges(4,4) + t289 - t290) * t220) * t220;
t244 = t2 * qJD(1) + t3 * qJD(2);
t98 = mrSges(6,1) * t155 + t156 * mrSges(6,2);
t228 = t98 * t307 - t57 * t129 / 0.2e1 + t58 * t313;
t238 = t75 * mrSges(6,1) / 0.2e1 - t76 * mrSges(6,2) / 0.2e1;
t6 = (-t280 / 0.2e1 - t279 / 0.2e1) * mrSges(6,3) + t228 + t238;
t101 = Ifges(6,2) * t155 - t152;
t102 = -Ifges(6,1) * t156 + t292;
t8 = -t51 * t131 + t50 * t129 - t220 * t258 / 0.2e1 - t183 * t98 - (-t50 * mrSges(6,3) + t90 / 0.2e1 + t101 / 0.2e1) * t156 + (t51 * mrSges(6,3) - t102 / 0.2e1 + t88 / 0.2e1) * t155;
t243 = -t6 * qJD(1) + t8 * qJD(2);
t20 = t155 * t131 - t156 * t129 + m(6) * (t155 * t50 - t156 * t51) + (-t214 * t179 - t216 * t181 + m(5) * (-t136 * t216 - t137 * t214)) * t218;
t235 = m(6) * (t155 * t57 - t156 * t58);
t23 = -t235 / 0.2e1 + (m(6) * t250 + (t250 + t274 / 0.2e1 + t273 / 0.2e1) * m(5)) * t218;
t242 = -qJD(1) * t23 + qJD(2) * t20;
t239 = qJD(2) * t98 - qJD(3) * t117;
t237 = t93 * t316 + t92 * t317;
t10 = t252 + t237;
t119 = -Ifges(6,2) * t232 + t291;
t17 = -t205 * t117 + t119 * t303 + t245 * t302 + t324 * t232 / 0.2e1;
t223 = -t128 * t131 / 0.2e1 + t183 * t117 / 0.2e1 - t205 * t98 / 0.2e1 - t220 * t257 / 0.4e1 - t324 * t156 / 0.4e1 - (t101 + t90) * t232 / 0.4e1 - (-t102 / 0.4e1 + t88 / 0.4e1) * t233 + (t327 / 0.2e1 + t314) * t127 + (t128 * t315 + t119 / 0.4e1 - t245 / 0.4e1) * t155;
t227 = -Ifges(6,3) * t218 / 0.2e1 + t53 * t317 + t54 * t316 + t236;
t5 = t223 + t227;
t231 = t10 * qJD(1) + t5 * qJD(2) - t17 * qJD(3);
t224 = (t155 * t302 - t156 * t305) * mrSges(6,3) + t240 * t320 + (t127 * t155 - t128 * t156 - t232 * t51 - t233 * t50) * t318 + t129 * t305 + t131 * t302 - t234;
t229 = t184 * t319 - t288 / 0.2e1 + t287 / 0.2e1;
t12 = (pkin(7) * t321 - t284 / 0.2e1 - t286 / 0.2e1) * t220 + t224 + t229;
t226 = t241 * t320 + (-t232 * t58 - t233 * t57) * t319;
t21 = t247 + t226;
t31 = -m(6) * (-t127 * t233 - t128 * t232) - m(5) * t325 + (-t232 ^ 2 - t233 ^ 2) * mrSges(6,3) - t326;
t230 = -qJD(1) * t21 + qJD(2) * t12 - qJD(3) * t31;
t213 = t220 ^ 2;
t212 = t218 ^ 2;
t193 = t212 * pkin(7) * t263;
t24 = m(5) * (-t273 - t274) * t298 + t235 / 0.2e1 + (t320 + t318) * t255;
t22 = t247 - t226;
t11 = t296 * t320 + mrSges(5,2) * t261 / 0.2e1 + mrSges(5,1) * t265 / 0.2e1 + t224 - t229;
t9 = t252 - t237;
t7 = -t228 + t238 + (t279 + t280) * t315;
t4 = t223 - t227;
t1 = t222 + (-t278 - t277) * t315 + (t191 + t118) * t218 * t250 - t195 * t263 - t330;
t13 = [t16 * qJD(2) + t15 * qJD(3), t1 * qJD(3) + t24 * qJD(4) + t7 * qJD(5) + t267 + (t76 * t129 + t75 * t131 + t134 * t181 + t135 * t179 + ((-t220 * mrSges(4,1) - mrSges(3,1) + t282) * t219 + (-mrSges(3,2) + (t165 + t99) * t218 + (t212 + t213) * mrSges(4,3)) * t221) * t215 + (t134 * t136 + t135 * t137 + t193) * t328 + (t183 * t255 + t50 * t75 + t51 * t76) * t329 + m(4) * (t193 + (pkin(7) * t213 * t221 - pkin(2) * t219) * t215)) * qJD(2), t1 * qJD(2) + t22 * qJD(4) + t9 * qJD(5) + t268 + ((-t232 * t93 - t233 * t92) * mrSges(6,3) + (t118 + t276) * t164 + (mrSges(4,2) - t326) * t163 + (-pkin(3) * t164 - t163 * t325) * t328 + (t127 * t92 + t128 * t93 + t164 * t205) * t329) * qJD(3), qJD(2) * t24 + qJD(3) * t22, t7 * qJD(2) + t9 * qJD(3) + (-mrSges(6,1) * t58 - mrSges(6,2) * t57) * qJD(5); qJD(3) * t2 - qJD(4) * t23 - qJD(5) * t6 - t267, qJD(3) * t3 + qJD(4) * t20 + qJD(5) * t8, t11 * qJD(4) + t4 * qJD(5) + t244 + (pkin(7) * t282 + m(6) * (t127 * t53 + t128 * t54 + t184 * t205) + t128 * t130 + t127 * t132 + t119 * t309 + t120 * t308 - pkin(3) * t166 + t89 * t305 + t91 * t303 + t184 * t118 + t205 * t100 + t154 * t300 + t153 * t299 - Ifges(4,6) * t218 + (m(5) * t323 + t216 * t180 - t214 * t182) * qJ(4) + ((Ifges(5,1) * t214 + t293) * t299 + (Ifges(5,2) * t216 + t294) * t301 + Ifges(4,5) + (-m(5) * pkin(3) + t276) * pkin(7)) * t220 + (Ifges(5,5) * t214 + Ifges(6,5) * t233 + Ifges(5,6) * t216 - Ifges(6,6) * t232) * t298 + (-t232 * t54 - t233 * t53) * mrSges(6,3) + t323 * mrSges(5,3)) * qJD(3), qJD(3) * t11 + t242, t4 * qJD(3) + (-mrSges(6,1) * t51 - mrSges(6,2) * t50 + t258) * qJD(5) + t243; -qJD(2) * t2 - qJD(4) * t21 + qJD(5) * t10 - t268, qJD(4) * t12 + qJD(5) * t5 - t244, -qJD(4) * t31 - qJD(5) * t17, t230, (-mrSges(6,1) * t128 - mrSges(6,2) * t127 + t257) * qJD(5) + t231; qJD(2) * t23 + qJD(3) * t21, -qJD(3) * t12 - qJD(5) * t98 - t242, qJD(5) * t117 - t230, 0, -t239; t6 * qJD(2) - t10 * qJD(3), -qJD(3) * t5 + qJD(4) * t98 - t243, -qJD(4) * t117 - t231, t239, 0;];
Cq = t13;
