% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP6
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:12
% EndTime: 2019-12-31 18:42:17
% DurationCPUTime: 2.48s
% Computational Cost: add. (3564->345), mult. (7717->481), div. (0->0), fcn. (6074->6), ass. (0->174)
t191 = cos(qJ(4));
t297 = Ifges(5,6) + Ifges(6,6);
t299 = t297 * t191;
t298 = Ifges(5,5) + Ifges(6,5);
t296 = Ifges(5,3) + Ifges(6,3);
t287 = m(6) * pkin(4);
t295 = mrSges(6,1) + t287;
t183 = Ifges(6,4) * t191;
t189 = sin(qJ(4));
t153 = Ifges(6,1) * t189 + t183;
t184 = Ifges(5,4) * t191;
t154 = Ifges(5,1) * t189 + t184;
t190 = sin(qJ(3));
t262 = t190 * pkin(3);
t192 = cos(qJ(3));
t263 = pkin(7) * t192;
t157 = t262 - t263;
t176 = sin(pkin(8)) * pkin(1) + pkin(6);
t238 = t190 * t176;
t45 = t157 * t191 + t189 * t238;
t237 = t190 * t191;
t46 = t157 * t189 - t176 * t237;
t294 = -t45 * t189 + t46 * t191;
t236 = t191 * t192;
t35 = pkin(4) * t190 - qJ(5) * t236 + t45;
t239 = t189 * t192;
t39 = -qJ(5) * t239 + t46;
t293 = -t189 * t35 + t191 * t39;
t228 = Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t254 = Ifges(6,4) * t189;
t209 = Ifges(6,1) * t191 - t254;
t87 = -Ifges(6,5) * t192 + t190 * t209;
t255 = Ifges(5,4) * t189;
t210 = Ifges(5,1) * t191 - t255;
t89 = -Ifges(5,5) * t192 + t190 * t210;
t292 = -t228 * t192 + t87 / 0.2e1 + t89 / 0.2e1;
t291 = m(5) / 0.2e1;
t290 = -m(6) / 0.2e1;
t289 = m(6) / 0.2e1;
t286 = -mrSges(5,1) / 0.2e1;
t285 = mrSges(6,1) / 0.2e1;
t284 = -mrSges(5,2) / 0.2e1;
t283 = mrSges(5,2) / 0.2e1;
t282 = -mrSges(6,2) / 0.2e1;
t281 = mrSges(6,2) / 0.2e1;
t280 = m(6) * t35;
t265 = pkin(4) * t189;
t214 = t176 + t265;
t99 = t214 * t190;
t279 = m(6) * t99;
t248 = t191 * mrSges(6,2);
t211 = t189 * mrSges(6,1) + t248;
t113 = t211 * t190;
t278 = t113 / 0.2e1;
t114 = t192 * t211;
t277 = -t114 / 0.2e1;
t141 = -mrSges(6,1) * t192 - mrSges(6,3) * t237;
t276 = -t141 / 0.2e1;
t246 = t192 * mrSges(5,1);
t142 = -mrSges(5,3) * t237 - t246;
t275 = -t142 / 0.2e1;
t143 = mrSges(6,1) * t190 - mrSges(6,3) * t236;
t274 = -t143 / 0.2e1;
t259 = -qJ(5) - pkin(7);
t146 = t259 * t189;
t273 = t146 / 0.2e1;
t249 = t191 * mrSges(5,2);
t150 = t189 * mrSges(5,1) + t249;
t272 = t150 / 0.2e1;
t175 = mrSges(6,1) * t237;
t271 = t175 / 0.2e1;
t270 = t189 / 0.2e1;
t269 = t190 / 0.2e1;
t264 = pkin(4) * t191;
t177 = -pkin(3) - t264;
t267 = m(6) * t177;
t266 = m(6) * t190;
t261 = t99 * mrSges(6,2);
t224 = qJ(5) * t237;
t223 = -cos(pkin(8)) * pkin(1) - pkin(2);
t120 = -pkin(3) * t192 - pkin(7) * t190 + t223;
t93 = t191 * t120;
t29 = -t224 + t93 + (-t176 * t189 - pkin(4)) * t192;
t235 = t192 * t176;
t40 = -t189 * t235 + t93;
t33 = t40 - t224;
t258 = -t29 + t33;
t252 = t190 * mrSges(5,1);
t251 = t190 * mrSges(5,2);
t250 = t191 * mrSges(5,1);
t186 = t189 ^ 2;
t187 = t191 ^ 2;
t215 = t187 / 0.2e1 + t186 / 0.2e1;
t219 = t275 + t276;
t240 = t189 * t190;
t137 = mrSges(6,2) * t192 - mrSges(6,3) * t240;
t138 = mrSges(5,2) * t192 - mrSges(5,3) * t240;
t220 = t138 / 0.2e1 + t137 / 0.2e1;
t7 = t192 * t271 + (((t284 + t282) * t192 + t220) * t189 + (t246 / 0.2e1 + (-pkin(4) * t192 + t258) * t290 - t219) * t191 + (mrSges(5,3) + mrSges(6,3)) * t190 * t215) * t190;
t243 = t7 * qJD(1);
t242 = mrSges(5,2) * t189 - mrSges(4,1) - t250;
t41 = t120 * t189 + t191 * t235;
t34 = -qJ(5) * t240 + t41;
t15 = (t189 * t137 - m(6) * (-t189 * t34 - t29 * t191) + t191 * t141) * t190;
t241 = qJD(1) * t15;
t222 = -t239 / 0.2e1;
t234 = mrSges(6,1) * t222 + t236 * t282;
t233 = -t237 * t287 - t175;
t232 = t186 + t187;
t231 = qJD(3) * t192;
t230 = qJD(4) * t190;
t227 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t225 = m(6) * t258;
t221 = t236 / 0.2e1;
t151 = Ifges(6,2) * t191 + t254;
t152 = Ifges(5,2) * t191 + t255;
t218 = -t151 / 0.2e1 - t152 / 0.2e1;
t217 = t153 / 0.2e1 + t154 / 0.2e1;
t181 = Ifges(6,5) * t191;
t182 = Ifges(5,5) * t191;
t216 = -t182 / 0.4e1 - t181 / 0.4e1;
t213 = t232 * mrSges(6,3);
t100 = t214 * t192;
t212 = mrSges(6,2) * t221 + t100 * t289 + t239 * t285;
t208 = -Ifges(5,2) * t189 + t184;
t207 = -Ifges(6,2) * t189 + t183;
t115 = t192 * t150;
t139 = -mrSges(6,2) * t190 - mrSges(6,3) * t239;
t140 = -mrSges(5,3) * t239 - t251;
t144 = -mrSges(5,3) * t236 + t252;
t83 = -Ifges(6,6) * t192 + t190 * t207;
t85 = -Ifges(5,6) * t192 + t190 * t208;
t198 = -t83 / 0.2e1 - t85 / 0.2e1 + t227 * t192;
t84 = Ifges(6,6) * t190 + t192 * t207;
t86 = Ifges(5,6) * t190 + t192 * t208;
t88 = Ifges(6,5) * t190 + t192 * t209;
t90 = Ifges(5,5) * t190 + t192 * t210;
t3 = t100 * t113 + t99 * t114 + t39 * t137 + t46 * t138 + t34 * t139 + t41 * t140 + t35 * t141 + t45 * t142 + t29 * t143 + t40 * t144 + m(6) * (t100 * t99 + t29 * t35 + t34 * t39) + m(5) * (t40 * t45 + t41 * t46) + (t223 * mrSges(4,2) + Ifges(4,4) * t192 + t198 * t189 + t292 * t191) * t192 + (t176 * t115 + t223 * mrSges(4,1) - Ifges(4,4) * t190 + (t88 / 0.2e1 + t90 / 0.2e1 + t228 * t190) * t191 + (-t84 / 0.2e1 - t86 / 0.2e1 - t227 * t190) * t189 + (Ifges(4,1) - Ifges(4,2) + (m(5) * t176 + t150) * t176 - t296) * t192) * t190;
t201 = m(5) * t294;
t6 = (-t115 / 0.2e1 + t277 + t220 * t191 + t219 * t189 + (-t189 * t40 + t191 * t41) * t291 + (-t189 * t29 + t191 * t34 - t100) * t289 - m(5) * t235 / 0.2e1) * t192 + (t278 + (t140 / 0.2e1 + t139 / 0.2e1 + t251 / 0.2e1) * t191 + (-t144 / 0.2e1 + t274 + t252 / 0.2e1) * t189 + t201 / 0.2e1 + (t99 + t293) * t289 + t238 * t291) * t190;
t206 = t3 * qJD(1) + t6 * qJD(2);
t116 = t190 * t151;
t117 = t190 * t152;
t118 = t190 * t153;
t119 = t190 * t154;
t5 = t33 * t137 + t40 * t138 - t41 * t142 + t99 * t175 + (t225 - t141) * t34 + ((t40 * mrSges(5,3) - t261 + t29 * mrSges(6,3) + t116 / 0.2e1 + t117 / 0.2e1 - mrSges(5,2) * t238 - t292) * t189 + (-t41 * mrSges(5,3) - t34 * mrSges(6,3) - t118 / 0.2e1 - t119 / 0.2e1 + mrSges(5,1) * t238 + (t113 + t279) * pkin(4) + t198) * t191) * t190;
t205 = t5 * qJD(1) - t7 * qJD(2);
t22 = 0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * (-0.1e1 + t232) * t192 * t190;
t204 = t6 * qJD(1) + t22 * qJD(2);
t149 = t259 * t191;
t203 = -t146 * t189 - t149 * t191;
t121 = t189 * t295 + t248;
t55 = mrSges(6,2) * t240 + t233;
t202 = qJD(1) * t55 - qJD(3) * t121;
t17 = (t272 + (t284 + t281) * t191 + (t286 + t285) * t189) * t192 + t234;
t147 = -mrSges(6,1) * t191 + t189 * mrSges(6,2);
t193 = t176 * t272 - t215 * pkin(7) * mrSges(5,3) + (mrSges(6,3) * t273 - t154 / 0.4e1 - t153 / 0.4e1 - t184 / 0.4e1 - t183 / 0.4e1 + t177 * t282 + pkin(3) * t283 + (Ifges(5,2) / 0.4e1 + Ifges(6,2) / 0.4e1) * t189) * t189 + (t149 * mrSges(6,3) / 0.2e1 - t152 / 0.4e1 - t151 / 0.4e1 + pkin(3) * t286 + (Ifges(5,1) / 0.4e1 + Ifges(6,1) / 0.4e1) * t191 + (-Ifges(5,4) / 0.4e1 - Ifges(6,4) / 0.4e1) * t189 + (t147 / 0.2e1 + t267 / 0.2e1) * pkin(4)) * t191;
t194 = -(t225 / 0.2e1 + t276) * t149 + (pkin(7) * t275 - t117 / 0.4e1 - t116 / 0.4e1 + t89 / 0.4e1 + t87 / 0.4e1 + t261 / 0.2e1 + (t33 / 0.2e1 - t29 / 0.2e1) * mrSges(6,3)) * t191 + t137 * t273 + t177 * t271;
t196 = -pkin(7) * t138 / 0.2e1 - t119 / 0.4e1 - t118 / 0.4e1 - t85 / 0.4e1 - t83 / 0.4e1 + t99 * t285 + (t279 / 0.2e1 + t278) * pkin(4);
t197 = -t35 * mrSges(6,1) / 0.2e1 + t39 * t281 + t45 * t286 + t46 * t283;
t2 = (-t280 / 0.2e1 + t274) * pkin(4) + t196 * t189 + (-t228 * t191 + (0.3e1 / 0.4e1 * Ifges(5,6) + 0.3e1 / 0.4e1 * Ifges(6,6)) * t189 + t216) * t192 + (-Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + t193) * t190 + t194 + t197;
t9 = pkin(3) * t150 - t177 * t211 - t147 * t265 + (-t183 / 0.2e1 - t184 / 0.2e1 - t217) * t191 + (-pkin(4) * t267 + (Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t189 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t191 - t218) * t189;
t200 = t2 * qJD(1) - t17 * qJD(2) - t9 * qJD(3);
t195 = ((-t146 * t190 + t34) * t191 + (t149 * t190 - t29) * t189) * t290 + t141 * t270 - t191 * t137 / 0.2e1;
t12 = t195 + t212;
t31 = m(6) * t203 + t213;
t57 = (-0.1e1 / 0.2e1 + t215) * t266;
t199 = -qJD(1) * t12 + qJD(2) * t57 + qJD(3) * t31;
t56 = (t232 + 0.1e1) * t266 / 0.2e1;
t18 = t222 * t287 + t234 + t277 + (-t150 / 0.2e1 - t249 / 0.2e1 + (t286 - t287 / 0.2e1) * t189) * t192;
t13 = -t195 + t212;
t4 = qJD(3) * t6 - qJD(4) * t7;
t1 = t216 * t192 + ((Ifges(5,6) / 0.4e1 + Ifges(6,6) / 0.4e1) * t192 + t196) * t189 + t193 * t190 + t194 - t197 + (t143 + t280) * pkin(4) / 0.2e1 + t296 * t269 + t297 * t222 + t298 * t221;
t8 = [qJD(3) * t3 + qJD(4) * t5 - qJD(5) * t15, t4, t1 * qJD(4) + t13 * qJD(5) + (Ifges(4,5) + t217 * t191 + t218 * t189 + (-m(5) * pkin(3) + t242) * t176) * t231 + t206 + (m(6) * (t100 * t177 + t146 * t35 - t149 * t39) + mrSges(4,2) * t238 - pkin(3) * t115 + t146 * t143 + t100 * t147 - t149 * t139 + t177 * t114 - Ifges(4,6) * t190 + (t191 * t140 - t189 * t144 + t201) * pkin(7) + (t88 + t90) * t270 + (t189 * t298 + t299) * t269 + (t84 + t86) * t191 / 0.2e1 + t293 * mrSges(6,3) + t294 * mrSges(5,3)) * qJD(3), t1 * qJD(3) + (-mrSges(5,1) * t41 - mrSges(5,2) * t40 - mrSges(6,2) * t33 - t295 * t34) * qJD(4) + (-t299 + (mrSges(6,3) * pkin(4) - t298) * t189) * t230 + t205, qJD(3) * t13 - t241; t4, t22 * qJD(3), t18 * qJD(4) + t56 * qJD(5) + (t147 + t242) * qJD(3) * t190 + (mrSges(5,3) * t232 - mrSges(4,2) + t213) * t231 + 0.2e1 * ((t232 * t263 - t262) * t291 + (t177 * t190 + t192 * t203) * t289) * qJD(3) + t204, -t243 + t18 * qJD(3) + t233 * qJD(4) + (-t250 + (mrSges(5,2) + mrSges(6,2)) * t189) * t230, t56 * qJD(3); qJD(4) * t2 - qJD(5) * t12 - t206, -qJD(4) * t17 + qJD(5) * t57 - t204, -qJD(4) * t9 + qJD(5) * t31, t200 + (-mrSges(6,2) * t146 - mrSges(6,3) * t264 - pkin(7) * t250 + t181 + t182 + (mrSges(5,2) * pkin(7) - t297) * t189 + t295 * t149) * qJD(4), t199; -qJD(3) * t2 + qJD(5) * t55 - t205, qJD(3) * t17 + t243, -qJD(5) * t121 - t200, 0, t202; qJD(3) * t12 - qJD(4) * t55 + t241, -t57 * qJD(3), qJD(4) * t121 - t199, -t202, 0;];
Cq = t8;
