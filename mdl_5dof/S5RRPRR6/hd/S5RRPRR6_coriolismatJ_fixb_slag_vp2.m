% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:45
% EndTime: 2019-12-05 18:35:51
% DurationCPUTime: 3.22s
% Computational Cost: add. (10609->336), mult. (21889->433), div. (0->0), fcn. (20619->8), ass. (0->185)
t319 = -mrSges(5,2) / 0.2e1;
t197 = sin(qJ(4));
t198 = sin(qJ(2));
t200 = cos(qJ(4));
t195 = cos(pkin(9));
t201 = cos(qJ(2));
t254 = t195 * t201;
t160 = (-t197 * t254 + t198 * t200) * pkin(1);
t161 = (t197 * t198 + t200 * t254) * pkin(1);
t196 = sin(qJ(5));
t199 = cos(qJ(5));
t295 = -mrSges(5,1) / 0.2e1;
t316 = -mrSges(6,2) / 0.2e1;
t317 = mrSges(6,1) / 0.2e1;
t96 = t160 * t199 - t161 * t196;
t97 = t160 * t196 + t161 * t199;
t307 = t97 * t316 + t96 * t317;
t308 = m(6) * pkin(4);
t318 = t307 - t160 * t295 + t161 * t319 + (t196 * t97 + t199 * t96) * t308 / 0.2e1;
t173 = -t196 * t197 + t199 * t200;
t174 = -t196 * t200 - t199 * t197;
t232 = t174 * mrSges(6,1) - t173 * mrSges(6,2);
t315 = t232 * qJD(5);
t194 = sin(pkin(9));
t258 = t194 * t197;
t170 = mrSges(5,2) * t195 - mrSges(5,3) * t258;
t250 = t200 * t170;
t257 = t194 * t200;
t171 = -mrSges(5,1) * t195 - mrSges(5,3) * t257;
t252 = t197 * t171;
t305 = -t252 / 0.2e1 + t250 / 0.2e1;
t176 = -t195 * pkin(3) - t194 * pkin(7) - pkin(2);
t167 = t200 * t176;
t243 = pkin(8) * t257;
t116 = -t243 + t167 + (-qJ(3) * t197 - pkin(4)) * t195;
t255 = t195 * t200;
t137 = qJ(3) * t255 + t176 * t197;
t183 = pkin(8) * t258;
t124 = t137 - t183;
t264 = t124 * t196;
t66 = t116 * t199 - t264;
t263 = t124 * t199;
t67 = t116 * t196 + t263;
t256 = t195 * t197;
t136 = -qJ(3) * t256 + t167;
t123 = t136 - t243;
t74 = -t123 * t196 - t263;
t75 = t123 * t199 - t264;
t314 = -m(6) * ((t66 - t75) * t174 + (t67 + t74) * t173) / 0.2e1 - t305;
t313 = t173 / 0.2e1;
t153 = t174 * t194;
t274 = t153 * mrSges(6,3);
t235 = -t274 / 0.2e1;
t128 = mrSges(6,2) * t195 + t274;
t152 = t173 * t194;
t275 = t152 * mrSges(6,3);
t129 = -mrSges(6,1) * t195 - t275;
t300 = t128 * t313 + (t129 + t275) * t174 / 0.2e1;
t231 = t173 * t235 + t300;
t154 = t174 * t195;
t155 = t173 * t195;
t249 = t154 * t317 + t155 * t316;
t303 = t231 + t249;
t311 = t303 * qJD(3);
t310 = qJD(5) * t303 + m(6) * (t154 * t173 - t155 * t174) * qJD(3);
t192 = t194 ^ 2;
t193 = t195 ^ 2;
t270 = t200 * mrSges(5,2);
t271 = t197 * mrSges(5,1);
t230 = t270 + t271;
t99 = -mrSges(6,1) * t153 + mrSges(6,2) * t152;
t309 = (t230 * t194 + t99) * t194 + (t192 + t193) * mrSges(4,3);
t289 = pkin(1) * t201;
t165 = t176 - t289;
t156 = t200 * t165;
t287 = t198 * pkin(1);
t186 = qJ(3) + t287;
t117 = -t186 * t256 + t156;
t118 = t165 * t197 + t186 * t255;
t178 = t192 * t186;
t189 = t192 * qJ(3);
t247 = t178 + t189;
t304 = ((t118 + t137) * t200 + (-t117 - t136) * t197) * t195 + t247;
t292 = (t154 * t199 + t155 * t196) * t308;
t301 = t256 * t295 + t255 * t319 + t292 / 0.2e1;
t299 = m(4) / 0.2e1;
t297 = m(5) / 0.2e1;
t296 = m(6) / 0.2e1;
t106 = t117 - t243;
t107 = t118 - t183;
t266 = t107 * t196;
t62 = t106 * t199 - t266;
t294 = t62 / 0.2e1;
t293 = t75 / 0.2e1;
t291 = -t195 / 0.2e1;
t288 = pkin(4) * t196;
t95 = -t243 + t156 + (-t186 * t197 - pkin(4)) * t195;
t57 = t199 * t95 - t266;
t286 = t57 * mrSges(6,2);
t265 = t107 * t199;
t58 = t196 * t95 + t265;
t285 = t58 * mrSges(6,1);
t61 = -t106 * t196 - t265;
t284 = t61 * mrSges(6,1);
t283 = t62 * mrSges(6,2);
t282 = t66 * mrSges(6,2);
t281 = t67 * mrSges(6,1);
t280 = t74 * mrSges(6,1);
t279 = t75 * mrSges(6,2);
t276 = pkin(4) * qJD(4);
t184 = pkin(4) * t258;
t157 = t194 * t186 + t184;
t244 = pkin(4) * t257;
t132 = t157 * t244;
t163 = (mrSges(5,1) * t200 - mrSges(5,2) * t197) * t194;
t248 = Ifges(6,5) * t153 - Ifges(6,6) * t152;
t217 = (Ifges(6,4) * t153 + Ifges(6,5) * t291) * t153 + t248 * t291;
t222 = Ifges(6,6) * t291 + Ifges(6,4) * t152 + (-Ifges(6,1) + Ifges(6,2)) * t153;
t44 = t57 * t274;
t98 = mrSges(6,1) * t152 + mrSges(6,2) * t153;
t207 = t157 * t98 - (t58 * mrSges(6,3) + t222) * t152 - t44 + t217;
t215 = Ifges(5,6) * t195 + pkin(4) * t99 + (-Ifges(5,4) * t200 + (-Ifges(5,1) + Ifges(5,2)) * t197) * t194;
t224 = Ifges(5,4) * t258 + Ifges(5,5) * t195;
t5 = m(6) * (t57 * t61 + t58 * t62 + t132) + t62 * t128 + t61 * t129 - t118 * t171 + t117 * t170 + (t186 * t163 + (t117 * mrSges(5,3) + t224) * t197 + (-t118 * mrSges(5,3) + t215) * t200) * t194 + t207;
t269 = t5 * qJD(1);
t7 = t57 * t128 - t58 * t129 + t207;
t268 = t7 * qJD(1);
t218 = t155 * t128 + t154 * t129 + (t250 - t252) * t195 + t309;
t23 = m(6) * (t57 * t154 + t58 * t155 + t157 * t194) + m(5) * (t178 + (-t117 * t197 + t118 * t200) * t195) + m(4) * (t186 * t193 + t178) + t218;
t267 = qJD(1) * t23;
t246 = t192 * t289;
t168 = t186 * t246;
t209 = ((-mrSges(4,1) * t195 + mrSges(4,2) * t194 - mrSges(3,1)) * t198 + (-mrSges(3,2) + t309) * t201) * pkin(1);
t219 = t97 * t128 + t96 * t129 + t160 * t171 + t161 * t170;
t245 = t194 * t289;
t259 = t193 * t201;
t13 = t209 + m(6) * (t157 * t245 + t57 * t96 + t58 * t97) + m(5) * (t117 * t160 + t118 * t161 + t168) + m(4) * (t186 * pkin(1) * t259 + t168 + (-pkin(2) - t289) * t287) + t219;
t262 = t13 * qJD(1);
t241 = t199 * t274;
t236 = t58 / 0.2e1 + t67 / 0.2e1;
t234 = t136 / 0.2e1 + t117 / 0.2e1;
t233 = -t137 / 0.2e1 - t118 / 0.2e1;
t172 = t194 * qJ(3) + t184;
t149 = t172 * t244;
t56 = t66 * t274;
t204 = (t157 / 0.2e1 + t172 / 0.2e1) * t98 - (t236 * mrSges(6,3) + t222) * t152 - t44 / 0.2e1 - t56 / 0.2e1 + t217;
t202 = t234 * t170 + t233 * t171 + (t74 / 0.2e1 + t61 / 0.2e1) * t129 + (t293 + t294) * t128 + ((qJ(3) / 0.2e1 + t186 / 0.2e1) * t163 + (t234 * mrSges(5,3) + t224) * t197 + (t233 * mrSges(5,3) + t215) * t200) * t194 + (t57 * t74 + t58 * t75 + t61 * t66 + t62 * t67 + t132 + t149) * t296 + t204;
t1 = t202 - t318;
t206 = t172 * t98 - (t67 * mrSges(6,3) + t222) * t152 - t56 + t217;
t6 = t75 * t128 + t74 * t129 + t136 * t170 - t137 * t171 + m(6) * (t66 * t74 + t67 * t75 + t149) + (qJ(3) * t163 + (t136 * mrSges(5,3) + t224) * t197 + (-t137 * mrSges(5,3) + t215) * t200) * t194 + t206;
t229 = t1 * qJD(1) + t6 * qJD(2);
t203 = (t57 / 0.2e1 + t66 / 0.2e1) * t128 - t236 * t129 + t204;
t3 = t203 - t307;
t8 = t66 * t128 - t67 * t129 + t206;
t228 = t3 * qJD(1) + t8 * qJD(2);
t27 = m(6) * (t66 * t154 + t67 * t155 + t172 * t194) + m(5) * (t189 + (-t136 * t197 + t137 * t200) * t195) + m(4) * (qJ(3) * t193 + t189) + t218;
t205 = ((qJ(3) + t186) * t193 + t247) * t299 + ((t157 + t172) * t194 + (t58 + t67) * t155 + (t57 + t66) * t154) * t296 + t218;
t208 = (t200 * t160 + t197 * t161) * t297 + (t173 * t96 - t174 * t97) * t296 + t287 * t299;
t9 = -t205 + t208 - m(5) * t304 / 0.2e1;
t227 = qJD(1) * t9 - qJD(2) * t27;
t211 = ((t57 - t62) * t174 + (t58 + t61) * t173) * t296 + t231 + t305;
t12 = (t271 / 0.2e1 + t270 / 0.2e1) * t195 - t292 / 0.2e1 + t211 - t249;
t25 = t274 * t313 + t249 - t300;
t15 = t25 + t301 + t314;
t226 = -t12 * qJD(1) + t15 * qJD(2);
t20 = t231 - t249;
t225 = -t20 * qJD(1) + t25 * qJD(2);
t133 = t275 * t288;
t220 = -Ifges(5,5) * t258 - Ifges(5,6) * t257 - t133 + t248;
t175 = (mrSges(6,1) * t196 + mrSges(6,2) * t199) * pkin(4);
t213 = (-t196 * t129 / 0.2e1 + (t128 / 0.2e1 + t235) * t199) * pkin(4) - t133 / 0.2e1;
t18 = (-t57 / 0.2e1 + t294) * mrSges(6,2) + (-t58 / 0.2e1 - t61 / 0.2e1) * mrSges(6,1) + t213;
t24 = (-t66 / 0.2e1 + t293) * mrSges(6,2) + (-t67 / 0.2e1 - t74 / 0.2e1) * mrSges(6,1) + t213;
t216 = -t18 * qJD(1) - t24 * qJD(2) + t175 * qJD(4);
t210 = t213 + t248;
t182 = qJ(3) * t246;
t169 = t175 * qJD(5);
t22 = -t282 / 0.2e1 - t281 / 0.2e1 - t279 / 0.2e1 + t280 / 0.2e1 + t210;
t16 = t303 + t301 - t314;
t14 = -t286 / 0.2e1 - t285 / 0.2e1 - t283 / 0.2e1 + t284 / 0.2e1 + t210;
t11 = t211 + t249 + t301;
t10 = t304 * t297 + t205 + t208;
t4 = t203 + t307;
t2 = t202 + t318;
t17 = [qJD(2) * t13 + qJD(3) * t23 + qJD(4) * t5 + qJD(5) * t7, t262 + t219 * qJD(2) + t10 * qJD(3) + t2 * qJD(4) + t4 * qJD(5) + qJD(2) * t209 + 0.2e1 * ((t172 * t245 + t66 * t96 + t67 * t97) * t296 + (t136 * t160 + t137 * t161 + t182) * t297 + (t182 + (-pkin(2) * t198 + qJ(3) * t259) * pkin(1)) * t299) * qJD(2), qJD(2) * t10 + qJD(4) * t11 + t267 + t310, t269 + t2 * qJD(2) + t11 * qJD(3) + (-t118 * mrSges(5,1) - t117 * mrSges(5,2) + t220 - t283 + t284) * qJD(4) + t14 * qJD(5) + (m(6) * (t196 * t62 + t199 * t61) - t241) * t276, t268 + t4 * qJD(2) + t311 + t14 * qJD(4) + (t248 - t285 - t286) * qJD(5); -qJD(3) * t9 + qJD(4) * t1 + qJD(5) * t3 - t262, qJD(3) * t27 + qJD(4) * t6 + qJD(5) * t8, qJD(4) * t16 - t227 + t310, t16 * qJD(3) + (-t137 * mrSges(5,1) - t136 * mrSges(5,2) + t220 - t279 + t280) * qJD(4) + t22 * qJD(5) + (-t241 + m(6) * (t196 * t75 + t199 * t74)) * t276 + t229, t311 + t22 * qJD(4) + (t248 - t281 - t282) * qJD(5) + t228; qJD(2) * t9 + qJD(4) * t12 + qJD(5) * t20 - t267, -qJD(4) * t15 - qJD(5) * t25 + t227, 0, (m(6) * (t199 * pkin(4) * t174 + t173 * t288) - t230 + t232) * qJD(4) + t315 - t226, qJD(4) * t232 - t225 + t315; -qJD(2) * t1 - qJD(3) * t12 + qJD(5) * t18 - t269, qJD(3) * t15 + qJD(5) * t24 - t229, t226, -t169, -t169 - t216; -qJD(2) * t3 - qJD(3) * t20 - qJD(4) * t18 - t268, qJD(3) * t25 - qJD(4) * t24 - t228, t225, t216, 0;];
Cq = t17;
