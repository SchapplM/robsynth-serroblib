% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP13_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP13_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:33
% EndTime: 2019-12-31 18:58:39
% DurationCPUTime: 2.59s
% Computational Cost: add. (3416->368), mult. (7167->499), div. (0->0), fcn. (5421->4), ass. (0->186)
t317 = Ifges(6,4) + Ifges(5,5);
t194 = cos(qJ(4));
t256 = t194 * mrSges(5,2);
t192 = sin(qJ(4));
t265 = t192 * mrSges(5,1);
t148 = t256 + t265;
t319 = -t148 / 0.2e1;
t318 = mrSges(6,2) + mrSges(5,3);
t190 = t192 ^ 2;
t191 = t194 ^ 2;
t308 = t190 + t191;
t316 = Ifges(6,2) + Ifges(5,3);
t193 = sin(qJ(3));
t230 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t315 = t230 * t193;
t195 = cos(qJ(3));
t243 = t192 * t195;
t130 = -t193 * mrSges(5,2) - mrSges(5,3) * t243;
t259 = t193 * mrSges(6,3);
t135 = -mrSges(6,2) * t243 + t259;
t314 = t130 + t135;
t218 = pkin(4) * t194 + qJ(5) * t192;
t142 = -pkin(3) - t218;
t257 = t194 * mrSges(6,1);
t262 = t192 * mrSges(6,3);
t221 = t257 + t262;
t313 = m(6) * t142 - t221;
t255 = t194 * mrSges(6,3);
t264 = t192 * mrSges(6,1);
t147 = -t255 + t264;
t247 = qJ(5) * t194;
t278 = pkin(4) * t192;
t146 = -t247 + t278;
t281 = m(6) * t146;
t312 = t147 + t281;
t186 = Ifges(6,5) * t192;
t311 = Ifges(6,1) * t194 + t186;
t310 = Ifges(6,6) * t192 + t317 * t194;
t189 = Ifges(5,4) * t194;
t309 = -Ifges(5,2) * t192 + t189;
t155 = Ifges(5,1) * t192 + t189;
t196 = -pkin(1) - pkin(6);
t236 = t195 * t196;
t276 = pkin(7) * t193;
t161 = pkin(3) * t195 + t276;
t239 = t194 * t161;
t58 = -t192 * t236 + t239;
t59 = t192 * t161 + t194 * t236;
t307 = -t192 * t58 + t194 * t59;
t45 = qJ(5) * t195 + t59;
t223 = t192 * t196 - pkin(4);
t46 = t195 * t223 - t239;
t215 = t192 * t46 + t194 * t45;
t306 = Ifges(6,3) * t194 - t186;
t305 = t195 ^ 2 / 0.2e1 + t193 ^ 2 / 0.2e1;
t184 = m(6) * qJ(5) + mrSges(6,3);
t302 = m(5) / 0.2e1;
t301 = -m(6) / 0.2e1;
t300 = m(6) / 0.2e1;
t299 = -mrSges(5,1) / 0.2e1;
t298 = -mrSges(6,1) / 0.2e1;
t297 = Ifges(6,6) / 0.2e1;
t296 = t46 / 0.2e1;
t237 = t194 * t196;
t275 = pkin(7) * t195;
t279 = pkin(3) * t193;
t141 = qJ(2) - t275 + t279;
t245 = t192 * t141;
t55 = t193 * t237 + t245;
t295 = t55 / 0.2e1;
t258 = t194 * mrSges(5,1);
t263 = t192 * mrSges(5,2);
t222 = t258 - t263;
t115 = t222 * t195;
t294 = -t115 / 0.2e1;
t242 = t193 * t194;
t171 = mrSges(6,2) * t242;
t252 = t195 * mrSges(6,1);
t132 = -t171 - t252;
t293 = t132 / 0.2e1;
t292 = t146 / 0.2e1;
t291 = t147 / 0.2e1;
t290 = -t192 / 0.2e1;
t289 = t192 / 0.2e1;
t288 = t193 / 0.2e1;
t286 = -t194 / 0.2e1;
t285 = t194 / 0.2e1;
t283 = t195 / 0.2e1;
t280 = m(6) * t194;
t274 = m(6) * qJD(5);
t271 = Ifges(5,4) * t192;
t270 = Ifges(6,5) * t194;
t268 = Ifges(5,6) * t193;
t238 = t194 * t195;
t133 = t193 * mrSges(5,1) - mrSges(5,3) * t238;
t134 = -t193 * mrSges(6,1) + mrSges(6,2) * t238;
t240 = t194 * t141;
t241 = t193 * t196;
t54 = -t192 * t241 + t240;
t214 = t55 * t192 + t54 * t194;
t40 = t245 + (qJ(5) + t237) * t193;
t42 = t193 * t223 - t240;
t216 = t40 * t192 - t42 * t194;
t9 = t193 * mrSges(4,1) + t195 * mrSges(4,2) + mrSges(3,3) + (t133 - t134) * t194 + t314 * t192 + (m(4) + m(3)) * qJ(2) + m(6) * t216 + m(5) * t214;
t266 = qJD(1) * t9;
t113 = t218 * t195;
t114 = t221 * t195;
t118 = t195 * t147;
t119 = t306 * t195;
t151 = Ifges(5,2) * t194 + t271;
t120 = t195 * t151;
t173 = Ifges(6,5) * t238;
t121 = -Ifges(6,1) * t243 + t173;
t122 = t195 * t155;
t172 = Ifges(6,6) * t238;
t91 = Ifges(6,4) * t193 + t195 * t311;
t156 = Ifges(5,1) * t194 - t271;
t93 = Ifges(5,5) * t193 + t195 * t156;
t207 = -t91 / 0.2e1 - t93 / 0.2e1 - t315;
t87 = Ifges(6,6) * t193 + Ifges(6,3) * t243 + t173;
t89 = t195 * t309 + t268;
t228 = -t87 / 0.2e1 + t89 / 0.2e1;
t210 = t146 - t196;
t61 = t210 * t195;
t4 = m(6) * (t113 * t61 + t40 * t54 + t42 * t55) + t113 * t118 + t61 * t114 + t54 * t135 + t172 * t288 + t54 * t130 - t55 * t133 + t55 * t134 + (-t196 * t115 + (-t40 * mrSges(6,2) - t268 / 0.2e1 - t55 * mrSges(5,3) + t121 / 0.2e1 - t122 / 0.2e1 - t228) * t194 + (t54 * mrSges(5,3) - t42 * mrSges(6,2) + t119 / 0.2e1 + t120 / 0.2e1 + t207) * t192) * t195;
t251 = t4 * qJD(1);
t206 = t318 * (-t191 / 0.2e1 - t190 / 0.2e1);
t197 = (t113 * t301 - t114 / 0.2e1 + t294) * t195 + ((t214 - t216) * t300 + t134 * t285 + t133 * t286 + t206 * t195 + t314 * t290) * t193;
t203 = t218 * t301 + t263 / 0.2e1 - t262 / 0.2e1 - t258 / 0.2e1 - t257 / 0.2e1;
t8 = t197 + t203;
t250 = t8 * qJD(1);
t249 = -t222 - mrSges(4,1);
t17 = m(6) * (t193 * t40 - t238 * t61) - t118 * t238 + t193 * t135;
t246 = qJD(1) * t17;
t244 = t192 * t193;
t56 = (0.1e1 / 0.2e1 + t305) * t280;
t235 = t56 * qJD(1);
t234 = t308 * t275;
t232 = qJD(3) * t193;
t231 = m(6) * t296;
t229 = t297 - Ifges(5,6) / 0.2e1;
t226 = t306 / 0.2e1 + t151 / 0.2e1;
t153 = Ifges(6,1) * t192 - t270;
t225 = -t153 / 0.2e1 - t155 / 0.2e1;
t150 = Ifges(6,3) * t192 + t270;
t116 = t147 * t193;
t117 = t193 * t148;
t129 = -t195 * mrSges(5,2) + mrSges(5,3) * t244;
t131 = t195 * mrSges(5,1) + mrSges(5,3) * t242;
t136 = mrSges(6,2) * t244 + t195 * mrSges(6,3);
t60 = t210 * t193;
t86 = Ifges(6,6) * t195 - t150 * t193;
t88 = Ifges(5,6) * t195 - t193 * t309;
t90 = Ifges(6,4) * t195 - t193 * t311;
t92 = Ifges(5,5) * t195 - t156 * t193;
t3 = -t61 * t116 - t60 * t118 + t55 * t129 + t59 * t130 + t54 * t131 + t42 * t132 + t58 * t133 + t46 * t134 + t45 * t135 + t40 * t136 + m(6) * (t40 * t45 + t42 * t46 - t60 * t61) + m(5) * (t54 * t58 + t55 * t59) + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t193 + t207 * t194 + (-t193 * t229 + t228) * t192) * t193 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t195 + t196 * t117 + (t90 / 0.2e1 + t92 / 0.2e1 + t230 * t195) * t194 + (t86 / 0.2e1 - t88 / 0.2e1 + t229 * t195) * t192 + (-Ifges(4,1) + Ifges(4,2) + (-m(5) * t196 + t148) * t196 + t316) * t193) * t195;
t209 = m(5) * t307;
t6 = (t118 / 0.2e1 + (t136 / 0.2e1 + t129 / 0.2e1) * t194 + (t293 - t131 / 0.2e1) * t192 + (t215 + t61) * t300 + t209 / 0.2e1) * t193 + (t116 / 0.2e1 + t148 * t288 + t117 / 0.2e1 + t134 * t289 + t133 * t290 + (t192 * t42 + t194 * t40 + t60) * t300 + (-t192 * t54 + t194 * t55 - 0.2e1 * t241) * t302 + t314 * t285) * t195;
t217 = t3 * qJD(1) + t6 * qJD(2);
t25 = 0.4e1 * (m(6) / 0.4e1 + m(5) / 0.4e1) * (-0.1e1 + t308) * t195 * t193;
t213 = t6 * qJD(1) + t25 * qJD(2);
t205 = (-t192 * t61 + (-t142 * t195 + t276) * t194) * t301 + t118 * t289;
t14 = -t171 + (-t221 * t285 + t298) * t195 + t231 + t205;
t38 = t313 * t192;
t212 = qJD(1) * t14 + qJD(3) * t38;
t21 = t259 + 0.2e1 * (t245 / 0.4e1 - t55 / 0.4e1 + (t237 / 0.4e1 + qJ(5) / 0.2e1) * t193) * m(6);
t211 = qJD(1) * t21 + qJD(4) * t184;
t11 = -pkin(3) * t148 - t146 * t221 + t312 * t142 + (-t150 / 0.2e1 + t309 / 0.2e1 - t225) * t194 + (t311 / 0.2e1 + t156 / 0.2e1 - t226) * t192;
t18 = (t148 / 0.2e1 + t291 + (-mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1) * t194 + (t299 + t298) * t192 + (-t278 / 0.2e1 + t247 / 0.2e1 + t292) * m(6)) * t195;
t198 = t196 * t319 + (t156 / 0.4e1 + t311 / 0.4e1 - t151 / 0.4e1 - t306 / 0.4e1) * t194 + (-t155 / 0.4e1 - t153 / 0.4e1 - t309 / 0.4e1 + t150 / 0.4e1) * t192 + t206 * pkin(7);
t199 = (t113 * t142 + t146 * t61) * t300 + pkin(3) * t294 - t113 * t221 / 0.2e1 + t142 * t114 / 0.2e1 + t118 * t292 + t61 * t291 + t310 * t193 / 0.4e1;
t200 = (-pkin(4) * t46 + qJ(5) * t45) * t301 + pkin(4) * t293 - qJ(5) * t136 / 0.2e1 - t45 * mrSges(6,3) / 0.2e1 + mrSges(6,1) * t296 + t58 * t299 + t59 * mrSges(5,2) / 0.2e1;
t201 = (t295 - t40 / 0.2e1) * mrSges(6,2) + ((-t40 + t55) * t300 - t135 / 0.2e1 - t130 / 0.2e1) * pkin(7) + t121 / 0.4e1 - t122 / 0.4e1 + t87 / 0.4e1 - t89 / 0.4e1;
t202 = -t120 / 0.4e1 - t119 / 0.4e1 + t93 / 0.4e1 + t91 / 0.4e1 + (t54 / 0.2e1 + t42 / 0.2e1) * mrSges(6,2) + ((t42 + t54) * t300 + t134 / 0.2e1 - t133 / 0.2e1) * pkin(7);
t2 = (t202 + t315) * t194 + ((-0.3e1 / 0.4e1 * Ifges(5,6) + t297) * t193 + t201) * t192 + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 + t198) * t195 + t199 + t200;
t208 = t2 * qJD(1) - t18 * qJD(2) + t11 * qJD(3);
t204 = (-m(6) * t218 - t221 - t222) * qJD(4);
t143 = (m(6) * pkin(7) + mrSges(6,2)) * t194;
t57 = (-0.1e1 / 0.2e1 + t305) * t280;
t20 = (t245 + (0.2e1 * qJ(5) + t237) * t193) * t300 + m(6) * t295 + t135;
t19 = (-t256 / 0.2e1 - t265 / 0.2e1 - t264 / 0.2e1 + t255 / 0.2e1 - t281 / 0.2e1 + t319 - t312 / 0.2e1) * t195;
t16 = t221 * t238 / 0.2e1 + t231 - t252 / 0.2e1 - t205;
t7 = t197 - t203;
t5 = t6 * qJD(3);
t1 = (-t268 / 0.4e1 + t201) * t192 + t202 * t194 + t198 * t195 - t200 + t199 + t316 * t283 + (-Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t244 - t317 * t242 / 0.2e1;
t10 = [qJD(2) * t9 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t17, qJD(4) * t7 + qJD(5) * t57 + t266 + t5, t1 * qJD(4) + t16 * qJD(5) + (-Ifges(4,5) + t225 * t194 + t226 * t192 + (-m(5) * pkin(3) + t249) * t196) * t232 + t217 + (pkin(3) * t117 - mrSges(4,2) * t236 - t142 * t116 + t86 * t286 + t88 * t285 - Ifges(4,6) * t195 + ((t129 + t136) * t194 + (-t131 + t132) * t192 + m(6) * t215 + t209) * pkin(7) - t313 * t60 + (t90 + t92) * t289 + ((Ifges(5,6) - Ifges(6,6)) * t194 + t317 * t192) * t283 + t307 * mrSges(5,3) + t215 * mrSges(6,2)) * qJD(3), t7 * qJD(2) + t1 * qJD(3) + t20 * qJD(5) + t251 + (t172 + (-m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1)) * t55 + (-mrSges(5,2) + t184) * t54 + ((-qJ(5) * mrSges(6,2) - Ifges(5,6)) * t194 + (pkin(4) * mrSges(6,2) - t317) * t192) * t195) * qJD(4), qJD(2) * t57 + qJD(3) * t16 + qJD(4) * t20 + t246; qJD(4) * t8 + qJD(5) * t56 - t266 + t5, t25 * qJD(3), t19 * qJD(4) + (-t221 + t249) * t232 + 0.2e1 * ((t234 - t279) * t302 + (t142 * t193 + t234) * t300) * qJD(3) + ((t318 * t308 - mrSges(4,2)) * qJD(3) + t192 * t274) * t195 + t213, t250 + t19 * qJD(3) + (t194 * t274 + t204) * t193, t235 + (qJD(3) * t243 + qJD(4) * t242) * m(6); qJD(4) * t2 - qJD(5) * t14 - t217, -qJD(4) * t18 - t213, qJD(4) * t11 - qJD(5) * t38, (-t218 * mrSges(6,2) - Ifges(5,6) * t192 + t310) * qJD(4) + t143 * qJD(5) + pkin(7) * t204 + t208, qJD(4) * t143 - t212; -qJD(2) * t8 - qJD(3) * t2 + qJD(5) * t21 - t251, qJD(3) * t18 - t250, -t208, t184 * qJD(5), t211; -qJD(2) * t56 + qJD(3) * t14 - qJD(4) * t21 - t246, -t235, t212, -t211, 0;];
Cq = t10;
