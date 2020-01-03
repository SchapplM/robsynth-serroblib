% Calculate vector of inverse dynamics joint torques for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:55
% EndTime: 2019-12-31 19:51:03
% DurationCPUTime: 3.70s
% Computational Cost: add. (3108->324), mult. (4827->400), div. (0->0), fcn. (2963->12), ass. (0->150)
t282 = mrSges(5,1) + mrSges(6,1);
t281 = mrSges(5,2) - mrSges(6,3);
t159 = sin(pkin(8));
t160 = cos(pkin(8));
t214 = t159 ^ 2 + t160 ^ 2;
t283 = mrSges(4,3) * t214;
t276 = Ifges(5,1) + Ifges(6,1);
t269 = Ifges(5,5) + Ifges(6,4);
t165 = cos(qJ(2));
t228 = pkin(1) * qJD(1);
t209 = t165 * t228;
t189 = qJD(3) - t209;
t275 = Ifges(6,5) - Ifges(5,4);
t272 = t160 * mrSges(4,1) - t159 * mrSges(4,2);
t280 = -mrSges(3,1) - t272;
t156 = pkin(8) + qJ(4);
t146 = sin(t156);
t147 = cos(t156);
t279 = t281 * t146 - t282 * t147;
t277 = -mrSges(6,2) - mrSges(5,3) - mrSges(4,3) + mrSges(3,2);
t162 = sin(qJ(4));
t246 = cos(qJ(4));
t109 = t159 * t246 + t162 * t160;
t247 = t109 / 0.2e1;
t268 = Ifges(6,6) - Ifges(5,6);
t157 = qJD(1) + qJD(2);
t172 = -t162 * t159 + t246 * t160;
t93 = t172 * t157;
t244 = Ifges(6,5) * t93;
t90 = Ifges(5,4) * t93;
t94 = t109 * t157;
t274 = t269 * qJD(4) + t276 * t94 - t244 + t90;
t161 = -pkin(7) - qJ(3);
t120 = t161 * t159;
t150 = t160 * pkin(7);
t121 = qJ(3) * t160 + t150;
t77 = t162 * t120 + t121 * t246;
t263 = -qJD(4) * t77 - t189 * t109;
t173 = t120 * t246 - t162 * t121;
t262 = qJD(4) * t173 + t189 * t172;
t100 = t172 * qJD(4);
t153 = qJDD(1) + qJDD(2);
t57 = t100 * t157 + t109 * t153;
t42 = -qJDD(4) * mrSges(6,1) + t57 * mrSges(6,2);
t267 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t57 + t42;
t101 = t109 * qJD(4);
t58 = t101 * t157 - t153 * t172;
t44 = -mrSges(6,2) * t58 + qJDD(4) * mrSges(6,3);
t266 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t58 + t44;
t237 = t93 * mrSges(5,3);
t238 = t93 * mrSges(6,2);
t81 = qJD(4) * mrSges(6,3) + t238;
t265 = -qJD(4) * mrSges(5,2) + t237 + t81;
t235 = t94 * mrSges(5,3);
t236 = t94 * mrSges(6,2);
t264 = t282 * qJD(4) - t235 - t236;
t260 = t283 * t153;
t158 = qJ(1) + qJ(2);
t148 = sin(t158);
t149 = cos(t158);
t259 = g(1) * t149 + g(2) * t148;
t163 = sin(qJ(2));
t210 = t163 * t228;
t111 = qJ(3) * t157 + t210;
t196 = pkin(7) * t157 + t111;
t83 = t196 * t159;
t84 = t196 * t160;
t39 = -t162 * t83 + t246 * t84;
t227 = pkin(1) * qJD(2);
t202 = qJD(1) * t227;
t222 = pkin(1) * qJDD(1);
t106 = t163 * t222 + t165 * t202;
t213 = qJD(3) * t157;
t82 = qJ(3) * t153 + t106 + t213;
t198 = pkin(7) * t153 + t82;
t64 = t198 * t159;
t65 = t198 * t160;
t8 = -qJD(4) * t39 - t162 * t65 - t246 * t64;
t141 = pkin(3) * t160 + pkin(2);
t91 = -t141 * t157 + t189;
t258 = m(5) * t91 - mrSges(5,1) * t93 + mrSges(5,2) * t94 - t157 * t272;
t203 = t246 * t83;
t224 = t162 * t84;
t38 = -t203 - t224;
t30 = -qJD(4) * pkin(4) + qJD(5) - t38;
t257 = -m(6) * t30 + t264;
t31 = qJD(4) * qJ(5) + t39;
t256 = -m(6) * t31 - t265;
t183 = pkin(4) * t147 + qJ(5) * t146;
t255 = t277 * t149 + (-m(6) * (-t141 - t183) - t279 - t280) * t148;
t220 = t147 * t149;
t221 = t146 * t149;
t254 = t277 * t148 + t280 * t149 - t282 * t220 + t281 * t221;
t252 = t93 / 0.2e1;
t251 = -t93 / 0.2e1;
t249 = t94 / 0.2e1;
t245 = Ifges(5,4) * t94;
t243 = pkin(1) * t163;
t164 = sin(qJ(1));
t242 = pkin(1) * t164;
t241 = pkin(1) * t165;
t166 = cos(qJ(1));
t151 = t166 * pkin(1);
t233 = qJD(4) / 0.2e1;
t219 = t149 * t161;
t218 = t153 * t159;
t217 = t153 * t160;
t215 = t149 * pkin(2) + t148 * qJ(3);
t211 = -qJD(4) * t203 - t162 * t64 + t246 * t65;
t208 = t163 * t227;
t207 = t165 * t227;
t15 = t58 * mrSges(5,1) + t57 * mrSges(5,2);
t14 = t58 * mrSges(6,1) - t57 * mrSges(6,3);
t197 = -pkin(2) * t148 + t149 * qJ(3);
t195 = t214 * t82;
t193 = t214 * t111;
t192 = mrSges(3,1) * t157 * t243;
t191 = t149 * t141 - t148 * t161;
t190 = t157 * t209;
t98 = -mrSges(4,1) * t217 + mrSges(4,2) * t218;
t105 = -t163 * t202 + t165 * t222;
t182 = -t141 * t148 - t219;
t177 = pkin(4) * t220 + qJ(5) * t221 + t191;
t176 = qJDD(3) - t105;
t140 = qJ(3) + t243;
t103 = (-pkin(7) - t140) * t159;
t104 = t140 * t160 + t150;
t174 = t103 * t246 - t162 * t104;
t67 = t162 * t103 + t104 * t246;
t40 = pkin(4) * t101 - qJ(5) * t100 - qJD(5) * t109;
t68 = -pkin(4) * t172 - qJ(5) * t109 - t141;
t74 = -t141 * t153 + t176;
t28 = -pkin(4) * t93 - qJ(5) * t94 + t91;
t4 = qJDD(4) * qJ(5) + (qJD(5) - t224) * qJD(4) + t211;
t89 = Ifges(6,5) * t94;
t45 = Ifges(6,6) * qJD(4) - t93 * Ifges(6,3) + t89;
t46 = t93 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t245;
t5 = -qJDD(4) * pkin(4) + qJDD(5) - t8;
t6 = pkin(4) * t58 - qJ(5) * t57 - qJD(5) * t94 + t74;
t7 = -qJD(4) * t224 + t211;
t92 = -pkin(2) * t153 + t176;
t167 = t172 * (Ifges(5,4) * t57 + Ifges(5,6) * qJDD(4)) / 0.2e1 - t92 * t272 + Ifges(3,3) * t153 + t105 * mrSges(3,1) - t106 * mrSges(3,2) + t74 * (-mrSges(5,1) * t172 + mrSges(5,2) * t109) + t6 * (-mrSges(6,1) * t172 - mrSges(6,3) * t109) + (Ifges(4,4) * t159 + Ifges(4,2) * t160) * t217 + (Ifges(4,1) * t159 + Ifges(4,4) * t160) * t218 - t172 * (Ifges(6,5) * t57 + Ifges(6,6) * qJDD(4)) / 0.2e1 + t283 * t82 + (t276 * t109 - t275 * t172) * t57 / 0.2e1 + (t269 * qJDD(4) + t276 * t57) * t247 + (t269 * t109 - t268 * t172) * qJDD(4) / 0.2e1 + (t45 / 0.2e1 + t28 * mrSges(6,1) + t91 * mrSges(5,1) - t46 / 0.2e1 + Ifges(6,3) * t251 - Ifges(5,2) * t252 + t275 * t249 + t268 * t233) * t101 + (-t39 * t101 - t109 * t8 + t172 * t7) * mrSges(5,3) + (-t101 * t31 + t5 * t109 + t172 * t4) * mrSges(6,2) + ((-Ifges(6,3) - Ifges(5,2)) * t172 + 0.2e1 * t275 * t247) * t58 + (t274 / 0.2e1 + t91 * mrSges(5,2) - t28 * mrSges(6,3) + Ifges(5,4) * t252 + Ifges(6,5) * t251 + t269 * t233 + t276 * t249 - t38 * mrSges(5,3) + t30 * mrSges(6,2)) * t100;
t143 = -pkin(2) - t241;
t127 = qJD(3) + t207;
t110 = -pkin(2) * t157 + t189;
t61 = t68 - t241;
t55 = -mrSges(6,1) * t93 - mrSges(6,3) * t94;
t54 = pkin(4) * t94 - qJ(5) * t93;
t37 = t40 + t208;
t1 = [t153 * mrSges(3,1) * t241 + m(6) * (t28 * t37 + t6 * t61) + m(3) * (t105 * t165 + t106 * t163) * pkin(1) + t143 * t98 + t61 * t14 + t167 + t37 * t55 + m(4) * (t127 * t193 + t143 * t92) - qJD(2) * t192 + Ifges(2,3) * qJDD(1) + (m(5) * t74 + t15) * (-t141 - t241) + (m(5) * t7 + m(6) * t4 + t266) * t67 - (-m(5) * t8 + m(6) * t5 + t267) * t174 + (-m(5) * t38 - t257) * (qJD(4) * t67 + t109 * t127) + (m(5) * t39 - t256) * (qJD(4) * t174 + t127 * t172) + (m(4) * t110 + t258) * t208 + t127 * t157 * t283 + (m(4) * t195 + t260) * t140 + (-t153 * t243 - t157 * t207) * mrSges(3,2) + (-mrSges(2,1) * t166 + t164 * mrSges(2,2) - m(3) * t151 - m(4) * (t151 + t215) - m(6) * (t151 + t177) - m(5) * (t151 + t191) + t254) * g(2) + (t164 * mrSges(2,1) + mrSges(2,2) * t166 + m(3) * t242 - m(5) * (t182 - t242) - m(6) * (-t219 - t242) - m(4) * (t197 - t242) + t255) * g(1); mrSges(3,2) * t190 - pkin(2) * t98 + qJD(1) * t192 + t68 * t14 - t141 * t15 + t40 * t55 + t167 + t266 * t77 - t267 * t173 + t260 * qJ(3) + (-t173 * t5 + t262 * t31 - t263 * t30 + t28 * t40 + t4 * t77 + t6 * t68) * m(6) + (-t141 * t74 + t173 * t8 + t262 * t39 + t263 * t38 + t7 * t77) * m(5) + (-pkin(2) * t92 + qJ(3) * t195 + qJD(3) * t193 - (t110 * t163 + t165 * t193) * t228) * m(4) + (-m(6) * t28 - t258 - t55) * t210 + (-m(4) * t215 - m(5) * t191 - m(6) * t177 + t254) * g(2) + (-m(4) * t197 - m(5) * t182 + m(6) * t219 + t255) * g(1) + t262 * t265 + t263 * t264 + (-t190 + t213) * t283; t264 * t94 - t265 * t93 - t157 ^ 2 * t283 + t98 + t14 + t15 + (-g(1) * t148 + g(2) * t149) * (m(4) + m(5) + m(6)) + (-t30 * t94 - t31 * t93 + t6) * m(6) + (t38 * t94 - t39 * t93 + t74) * m(5) + (-t157 * t193 + t92) * m(4); (-pkin(4) * t5 - g(3) * t183 + qJ(5) * t4 + qJD(5) * t31 - t28 * t54) * m(6) - (t276 * t93 - t245 + t45 + t89) * t94 / 0.2e1 - (t268 * t94 + t269 * t93) * qJD(4) / 0.2e1 - t91 * (mrSges(5,1) * t94 + mrSges(5,2) * t93) - t28 * (mrSges(6,1) * t94 - mrSges(6,3) * t93) + (t237 + t256) * t38 + (Ifges(6,3) * t94 + t244) * t252 + ((-qJ(5) * m(6) + t281) * t147 + (pkin(4) * m(6) + t282) * t146) * t259 + t279 * g(3) + t268 * t58 + t269 * t57 + (t235 + t257) * t39 + (Ifges(5,3) + Ifges(6,2)) * qJDD(4) + qJD(5) * t81 - t54 * t55 - pkin(4) * t42 + qJ(5) * t44 + t8 * mrSges(5,1) + t4 * mrSges(6,3) - t5 * mrSges(6,1) - t7 * mrSges(5,2) + (-Ifges(5,2) * t94 + t274 + t90) * t251 + t31 * t236 - t30 * t238 + t46 * t249; -qJD(4) * t81 + t94 * t55 + (g(3) * t147 - t31 * qJD(4) - t146 * t259 + t28 * t94 + t5) * m(6) + t42;];
tau = t1;
