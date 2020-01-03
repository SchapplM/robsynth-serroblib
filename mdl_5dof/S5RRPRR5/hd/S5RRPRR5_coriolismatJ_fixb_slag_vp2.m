% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR5
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:35
% EndTime: 2020-01-03 12:03:40
% DurationCPUTime: 2.76s
% Computational Cost: add. (9988->240), mult. (19245->304), div. (0->0), fcn. (20855->8), ass. (0->148)
t194 = sin(qJ(2));
t274 = t194 * pkin(1);
t183 = qJ(3) + t274;
t190 = sin(pkin(9));
t162 = (-pkin(7) - t183) * t190;
t191 = cos(pkin(9));
t187 = t191 * pkin(7);
t163 = t183 * t191 + t187;
t193 = sin(qJ(4));
t196 = cos(qJ(4));
t129 = t162 * t193 + t163 * t196;
t172 = -t190 * t193 + t191 * t196;
t276 = pkin(8) * t172;
t106 = t129 + t276;
t192 = sin(qJ(5));
t195 = cos(qJ(5));
t128 = t196 * t162 - t163 * t193;
t173 = t190 * t196 + t191 * t193;
t169 = t173 * pkin(8);
t294 = t128 - t169;
t306 = -t106 * t192 + t195 * t294;
t65 = t106 * t195 + t192 * t294;
t330 = -t65 * mrSges(6,1) - t306 * mrSges(6,2);
t148 = t172 * t192 + t173 * t195;
t222 = t195 * t172 - t173 * t192;
t333 = Ifges(6,5) * t222 - Ifges(6,6) * t148;
t12 = t333 + t330;
t335 = t12 * qJD(5);
t177 = (-pkin(7) - qJ(3)) * t190;
t179 = qJ(3) * t191 + t187;
t152 = t177 * t193 + t179 * t196;
t121 = t152 + t276;
t151 = t196 * t177 - t179 * t193;
t293 = t151 - t169;
t307 = -t121 * t192 + t195 * t293;
t93 = t121 * t195 + t192 * t293;
t331 = -t93 * mrSges(6,1) - t307 * mrSges(6,2);
t16 = t333 + t331;
t334 = t16 * qJD(5);
t282 = t93 / 0.2e1;
t284 = t65 / 0.2e1;
t332 = t282 + t284;
t316 = -Ifges(6,4) * t222 ^ 2 + (Ifges(6,4) * t148 + (-Ifges(6,1) + Ifges(6,2)) * t222) * t148;
t184 = -t191 * pkin(3) - pkin(2);
t156 = -t172 * pkin(4) + t184;
t300 = t222 * mrSges(6,2);
t311 = t148 * mrSges(6,1);
t317 = t311 + t300;
t324 = t156 * t317;
t197 = cos(qJ(2));
t278 = pkin(1) * t197;
t153 = t156 - t278;
t325 = t153 * t317;
t329 = -t316 + t324 / 0.2e1 + t325 / 0.2e1;
t275 = t173 * pkin(4);
t328 = m(6) * t275;
t286 = m(6) / 0.2e1;
t323 = (t153 + t156) * t173 * t286;
t227 = t311 / 0.2e1;
t319 = (Ifges(5,1) - Ifges(5,2)) * t173;
t157 = t173 * t278;
t158 = t172 * t278;
t109 = -t157 * t195 - t158 * t192;
t110 = -t157 * t192 + t158 * t195;
t235 = t109 * mrSges(6,1) / 0.2e1 - t110 * mrSges(6,2) / 0.2e1;
t318 = t157 * mrSges(5,1) / 0.2e1 + t158 * mrSges(5,2) / 0.2e1 - t235;
t313 = t65 + t93;
t310 = t148 * mrSges(6,3);
t247 = t173 * mrSges(5,3);
t233 = t190 ^ 2 + t191 ^ 2;
t295 = t233 * mrSges(4,3);
t309 = -mrSges(3,2) + t295;
t308 = t129 + t152;
t299 = t222 * mrSges(6,3);
t232 = qJD(1) + qJD(2);
t149 = t173 * mrSges(5,1) + t172 * mrSges(5,2);
t236 = t222 * t192;
t237 = t148 * t195;
t279 = t173 / 0.2e1;
t285 = m(6) * pkin(4);
t29 = (t279 + t237 / 0.2e1 - t236 / 0.2e1) * t285 + t317 + t149;
t297 = t232 * t29;
t46 = 0.2e1 * t227 + t300;
t296 = t232 * t46;
t224 = t233 * qJ(3);
t176 = t184 - t278;
t122 = t176 * t149;
t127 = t184 * t149;
t168 = Ifges(5,4) * t172;
t262 = Ifges(5,4) * t173;
t281 = -t307 / 0.2e1;
t283 = -t306 / 0.2e1;
t43 = t306 * t299;
t44 = t65 * t310;
t54 = t307 * t299;
t55 = t93 * t310;
t291 = (Ifges(5,1) * t172 - t262) * t279 - t173 * (Ifges(5,2) * t172 + t262) / 0.2e1 + (t332 * t148 + (t281 + t283) * t222) * mrSges(6,3) + t122 / 0.2e1 + t127 / 0.2e1 + t43 / 0.2e1 - t44 / 0.2e1 + t54 / 0.2e1 - t55 / 0.2e1 + (0.2e1 * t168 + t319) * t172 / 0.2e1 + (t152 / 0.2e1 + t129 / 0.2e1 - t308 / 0.2e1) * t247 + t329;
t97 = -mrSges(6,1) * t222 + mrSges(6,2) * t148;
t290 = -mrSges(4,1) * t191 - mrSges(5,1) * t172 + mrSges(4,2) * t190 + mrSges(5,2) * t173 - mrSges(3,1) + t97;
t288 = t173 ^ 2;
t287 = m(5) / 0.2e1;
t277 = pkin(2) * t194;
t265 = t29 * qJD(4) + t46 * qJD(5);
t47 = t227 - t311 / 0.2e1;
t79 = (t173 + t236 - t237) * t285 / 0.2e1;
t264 = t79 * qJD(4) + t47 * qJD(5);
t260 = pkin(4) * qJD(4);
t231 = t97 * t275;
t200 = t288 * Ifges(5,4) + (-t168 - t319) * t172 - t231 + t316;
t3 = -t43 + t44 - t153 * t328 - t325 - t122 + (-t65 * t148 + t222 * t306) * mrSges(6,3) + t200;
t246 = t3 * qJD(1);
t7 = t316 - t325;
t245 = t7 * qJD(1);
t204 = t148 * t310 + t222 * t299 + (t172 ^ 2 + t288) * mrSges(5,3) + t295;
t223 = t233 * t183;
t18 = m(6) * (-t148 * t306 + t222 * t65) + m(5) * (-t128 * t173 + t129 * t172) + m(4) * t223 + t204;
t243 = qJD(1) * t18;
t203 = t158 * t172 * mrSges(5,3) - t109 * t310 + t110 * t299 + t157 * t247;
t11 = m(6) * (t109 * t306 + t110 * t65) + m(5) * (-t128 * t157 + t129 * t158) + m(4) * (t223 * t197 - t277) * pkin(1) + t203 + (m(5) * t176 + m(6) * t153 + t290) * t274 + (-m(4) * t274 + t309) * t278;
t240 = t11 * qJD(1);
t228 = t195 * t299;
t210 = m(6) * (t109 * t195 + t110 * t192);
t1 = (-t210 / 0.2e1 + t323 + t173 * t97) * pkin(4) + t291 + t318;
t6 = -t54 + t55 - t324 - t127 - t156 * t328 + (-t93 * t148 + t222 * t307) * mrSges(6,3) + t200;
t220 = t1 * qJD(1) - t6 * qJD(2);
t202 = (-t313 / 0.2e1 + t332) * t310 + t329;
t4 = t202 - t235;
t8 = t316 - t324;
t219 = t4 * qJD(1) - t8 * qJD(2);
t217 = -pkin(4) * t192 * t310 + Ifges(5,5) * t172 - Ifges(5,6) * t173 + t333;
t19 = m(6) * (-t148 * t307 + t222 * t93) + m(5) * (-t151 * t173 + t152 * t172) + m(4) * t224 + t204;
t198 = -m(4) * (t223 + t224) / 0.2e1 - m(5) * ((-t128 - t151) * t173 + t308 * t172) / 0.2e1 - m(6) * (t313 * t222 - (t306 + t307) * t148) / 0.2e1 - t204;
t212 = (t286 + t287 + m(4) / 0.2e1) * t274;
t9 = t212 + t198;
t216 = qJD(1) * t9 - qJD(2) * t19;
t13 = (t283 + t306 / 0.2e1) * mrSges(6,2) + (-t65 / 0.2e1 + t284) * mrSges(6,1);
t17 = (t281 + t307 / 0.2e1) * mrSges(6,2) + (-t93 / 0.2e1 + t282) * mrSges(6,1);
t175 = (mrSges(6,1) * t192 + mrSges(6,2) * t195) * pkin(4);
t206 = -qJD(1) * t13 - qJD(2) * t17 + qJD(4) * t175;
t174 = t175 * qJD(5);
t76 = t79 * qJD(3);
t42 = t46 * qJD(3);
t41 = t47 * qJD(3);
t28 = t29 * qJD(3);
t10 = t212 - t198;
t5 = t202 + t235;
t2 = t231 + t291 - t318 + (t323 + t210 / 0.2e1) * pkin(4);
t14 = [qJD(2) * t11 + qJD(3) * t18 - qJD(4) * t3 - qJD(5) * t7, t10 * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t240 + (t203 + 0.2e1 * (t109 * t307 + t110 * t93) * t286 + 0.2e1 * (-t151 * t157 + t152 * t158) * t287 + (t309 * t197 + m(4) * (t197 * t224 - t277) + (m(5) * t184 + m(6) * t156 + t290) * t194) * pkin(1)) * qJD(2), qJD(2) * t10 + t243 + t264, -t246 + t2 * qJD(2) + t76 + (-t129 * mrSges(5,1) - t128 * mrSges(5,2) + t217 + t330) * qJD(4) + t335 + (-t228 + m(6) * (t192 * t306 - t195 * t65)) * t260, t5 * qJD(2) + t12 * qJD(4) - t245 + t335 + t41; -qJD(3) * t9 + qJD(4) * t1 + qJD(5) * t4 - t240, qJD(3) * t19 - qJD(4) * t6 - qJD(5) * t8, -t216 + t264, t76 + (-t152 * mrSges(5,1) - t151 * mrSges(5,2) + t217 + t331) * qJD(4) + t334 + (-t228 + m(6) * (t192 * t307 - t195 * t93)) * t260 + t220, t16 * qJD(4) + t219 + t334 + t41; qJD(2) * t9 - t243 + t265, t216 + t265, 0, t297, t296; -qJD(2) * t1 + qJD(5) * t13 + t246 - t28, qJD(5) * t17 - t220 - t28, -t297, -t174, -t174 - t206; -qJD(2) * t4 - qJD(4) * t13 + t245 - t42, -qJD(4) * t17 - t219 - t42, -t296, t206, 0;];
Cq = t14;
