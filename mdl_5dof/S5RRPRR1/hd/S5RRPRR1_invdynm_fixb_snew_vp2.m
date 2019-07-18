% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynm_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:50
% EndTime: 2019-07-18 17:20:59
% DurationCPUTime: 2.88s
% Computational Cost: add. (32142->261), mult. (69302->329), div. (0->0), fcn. (41816->8), ass. (0->103)
t245 = sin(qJ(2));
t249 = cos(qJ(2));
t267 = qJD(1) * qJD(2);
t265 = t249 * t267;
t220 = t245 * qJDD(1) + t265;
t266 = qJD(1) * qJD(3);
t251 = qJD(1) ^ 2;
t272 = t245 * t251;
t246 = sin(qJ(1));
t250 = cos(qJ(1));
t228 = -t250 * g(1) - t246 * g(2);
t273 = t245 * t228;
t258 = qJ(3) * t265 - 0.2e1 * t245 * t266 - t273 + (t272 * t249 + qJDD(2)) * pkin(1);
t276 = -pkin(3) - qJ(3);
t163 = qJDD(2) * pkin(2) + t276 * t220 + (pkin(2) * t272 + pkin(3) * t267 - g(3)) * t249 + t258;
t221 = t249 * qJDD(1) - t245 * t267;
t269 = qJD(1) * t245;
t223 = qJD(2) * pkin(1) - qJ(3) * t269;
t226 = qJD(2) * pkin(2) - pkin(3) * t269;
t199 = -t245 * g(3) + t249 * t228;
t261 = t221 * qJ(3) + 0.2e1 * t249 * t266 + t199;
t275 = t249 ^ 2 * t251;
t278 = -pkin(1) - pkin(2);
t164 = t221 * pkin(3) + t278 * t275 + (-t223 - t226) * qJD(2) + t261;
t244 = sin(qJ(4));
t248 = cos(qJ(4));
t157 = t244 * t163 + t248 * t164;
t208 = (t244 * t249 + t245 * t248) * qJD(1);
t180 = -t208 * qJD(4) - t244 * t220 + t248 * t221;
t207 = (t244 * t245 - t248 * t249) * qJD(1);
t190 = t207 * mrSges(5,1) + t208 * mrSges(5,2);
t236 = qJD(2) + qJD(4);
t194 = t236 * mrSges(5,1) - t208 * mrSges(5,3);
t235 = qJDD(2) + qJDD(4);
t153 = (t207 * t208 + t235) * pkin(4) + t157;
t227 = t246 * g(1) - t250 * g(2);
t260 = t223 * t269 + qJDD(3) - t227;
t166 = t278 * t221 + t226 * t269 + t276 * t275 + t260;
t181 = -t207 * qJD(4) + t248 * t220 + t244 * t221;
t158 = (t207 * t236 - t181) * pkin(4) + t166;
t243 = sin(qJ(5));
t247 = cos(qJ(5));
t151 = -t243 * t153 + t247 * t158;
t191 = -t243 * t208 + t247 * t236;
t162 = t191 * qJD(5) + t247 * t181 + t243 * t235;
t192 = t247 * t208 + t243 * t236;
t171 = -t191 * mrSges(6,1) + t192 * mrSges(6,2);
t179 = qJDD(5) - t180;
t200 = qJD(5) + t207;
t182 = -t200 * mrSges(6,2) + t191 * mrSges(6,3);
t149 = m(6) * t151 + t179 * mrSges(6,1) - t162 * mrSges(6,3) - t192 * t171 + t200 * t182;
t152 = t247 * t153 + t243 * t158;
t161 = -t192 * qJD(5) - t243 * t181 + t247 * t235;
t183 = t200 * mrSges(6,1) - t192 * mrSges(6,3);
t150 = m(6) * t152 - t179 * mrSges(6,2) + t161 * mrSges(6,3) + t191 * t171 - t200 * t183;
t263 = -t243 * t149 + t247 * t150;
t135 = m(5) * t157 - t235 * mrSges(5,2) + t180 * mrSges(5,3) - t207 * t190 - t236 * t194 + t263;
t156 = t248 * t163 - t244 * t164;
t154 = (-t208 ^ 2 - t236 ^ 2) * pkin(4) - t156;
t193 = -t236 * mrSges(5,2) - t207 * mrSges(5,3);
t145 = m(5) * t156 - m(6) * t154 + t235 * mrSges(5,1) + t161 * mrSges(6,1) - t162 * mrSges(6,2) - t181 * mrSges(5,3) + t191 * t182 - t192 * t183 - t208 * t190 + t236 * t193;
t131 = t244 * t135 + t248 * t145;
t277 = t249 * g(3);
t174 = -t220 * qJ(3) + t258 - t277;
t219 = (-mrSges(4,1) * t249 + mrSges(4,2) * t245) * qJD(1);
t268 = qJD(1) * t249;
t225 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t268;
t128 = m(4) * t174 + qJDD(2) * mrSges(4,1) - t220 * mrSges(4,3) + qJD(2) * t225 - t219 * t269 + t131;
t198 = -t273 - t277;
t204 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t245 + Ifges(3,2) * t249) * qJD(1);
t205 = Ifges(4,5) * qJD(2) + (Ifges(4,1) * t245 + Ifges(4,4) * t249) * qJD(1);
t206 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t245 + Ifges(3,4) * t249) * qJD(1);
t175 = -pkin(1) * t275 - qJD(2) * t223 + t261;
t203 = Ifges(4,6) * qJD(2) + (Ifges(4,4) * t245 + Ifges(4,2) * t249) * qJD(1);
t167 = Ifges(6,5) * t192 + Ifges(6,6) * t191 + Ifges(6,3) * t200;
t169 = Ifges(6,1) * t192 + Ifges(6,4) * t191 + Ifges(6,5) * t200;
t142 = -mrSges(6,1) * t154 + mrSges(6,3) * t152 + Ifges(6,4) * t162 + Ifges(6,2) * t161 + Ifges(6,6) * t179 - t192 * t167 + t200 * t169;
t168 = Ifges(6,4) * t192 + Ifges(6,2) * t191 + Ifges(6,6) * t200;
t143 = mrSges(6,2) * t154 - mrSges(6,3) * t151 + Ifges(6,1) * t162 + Ifges(6,4) * t161 + Ifges(6,5) * t179 + t191 * t167 - t200 * t168;
t187 = Ifges(5,4) * t208 - Ifges(5,2) * t207 + Ifges(5,6) * t236;
t188 = Ifges(5,1) * t208 - Ifges(5,4) * t207 + Ifges(5,5) * t236;
t256 = -mrSges(5,1) * t156 + mrSges(5,2) * t157 - Ifges(5,5) * t181 - Ifges(5,6) * t180 - Ifges(5,3) * t235 - pkin(4) * t263 - t247 * t142 - t243 * t143 - t208 * t187 - t207 * t188;
t253 = -mrSges(4,1) * t174 + mrSges(4,2) * t175 - Ifges(4,5) * t220 - Ifges(4,6) * t221 - Ifges(4,3) * qJDD(2) - pkin(2) * t131 - t203 * t269 + t256;
t280 = mrSges(3,1) * t198 - mrSges(3,2) * t199 + Ifges(3,5) * t220 + Ifges(3,6) * t221 + Ifges(3,3) * qJDD(2) + pkin(1) * t128 - (-t245 * t204 + (t205 + t206) * t249) * qJD(1) - t253;
t271 = t247 * t149 + t243 * t150;
t264 = t248 * t135 - t244 * t145;
t137 = m(5) * t166 - t180 * mrSges(5,1) + t181 * mrSges(5,2) + t207 * t193 + t208 * t194 + t271;
t185 = -t221 * pkin(1) - qJ(3) * t275 + t260;
t201 = Ifges(4,3) * qJD(2) + (Ifges(4,5) * t245 + Ifges(4,6) * t249) * qJD(1);
t202 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t245 + Ifges(3,6) * t249) * qJD(1);
t224 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t269;
t186 = Ifges(5,5) * t208 - Ifges(5,6) * t207 + Ifges(5,3) * t236;
t127 = mrSges(5,2) * t166 - mrSges(5,3) * t156 + Ifges(5,1) * t181 + Ifges(5,4) * t180 + Ifges(5,5) * t235 - pkin(4) * t271 - t243 * t142 + t247 * t143 - t207 * t186 - t236 * t187;
t254 = mrSges(6,1) * t151 - mrSges(6,2) * t152 + Ifges(6,5) * t162 + Ifges(6,6) * t161 + Ifges(6,3) * t179 + t192 * t168 - t191 * t169;
t136 = -mrSges(5,1) * t166 + mrSges(5,3) * t157 + Ifges(5,4) * t181 + Ifges(5,2) * t180 + Ifges(5,6) * t235 - t208 * t186 + t236 * t188 - t254;
t257 = -mrSges(4,1) * t185 + mrSges(4,3) * t175 + Ifges(4,4) * t220 + Ifges(4,2) * t221 + Ifges(4,6) * qJDD(2) - pkin(2) * t137 + pkin(3) * t264 + qJD(2) * t205 + t244 * t127 + t248 * t136;
t121 = t257 - pkin(1) * (m(4) * t185 - t221 * mrSges(4,1) + t220 * mrSges(4,2) + t137) + qJ(3) * (m(4) * t175 - qJDD(2) * mrSges(4,2) + t221 * mrSges(4,3) - qJD(2) * t224 + t264) + Ifges(3,4) * t220 + Ifges(3,2) * t221 + mrSges(3,1) * t227 + qJD(2) * t206 + mrSges(3,3) * t199 + Ifges(3,6) * qJDD(2) + ((pkin(1) * t225 + qJ(3) * t219) * t249 + (-pkin(1) * t224 - t201 - t202) * t245) * qJD(1);
t255 = mrSges(4,2) * t185 - mrSges(4,3) * t174 + Ifges(4,1) * t220 + Ifges(4,4) * t221 + Ifges(4,5) * qJDD(2) - pkin(3) * t131 + t248 * t127 - t244 * t136 + t201 * t268;
t123 = t255 + Ifges(3,1) * t220 + Ifges(3,4) * t221 - mrSges(3,2) * t227 - mrSges(3,3) * t198 - qJ(3) * t128 + Ifges(3,5) * qJDD(2) + t202 * t268 + (-t203 - t204) * qJD(2);
t259 = mrSges(2,1) * t227 - mrSges(2,2) * t228 + Ifges(2,3) * qJDD(1) + t249 * t121 + t245 * t123;
t124 = mrSges(2,1) * g(3) + mrSges(2,3) * t228 + t251 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - t280;
t119 = -mrSges(2,2) * g(3) - mrSges(2,3) * t227 + Ifges(2,5) * qJDD(1) - t251 * Ifges(2,6) - t245 * t121 + t249 * t123;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t250 * t119 - t246 * t124, t119, t123, -qJD(2) * t203 + t255, t127, t143; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t246 * t119 + t250 * t124, t124, t121, -t201 * t269 + t257, t136, t142; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t259, t259, t280, -t205 * t268 - t253, -t256, t254;];
m_new  = t1;
