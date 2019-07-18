% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:18
% EndTime: 2019-07-18 17:17:33
% DurationCPUTime: 6.45s
% Computational Cost: add. (86951->287), mult. (176266->367), div. (0->0), fcn. (132333->10), ass. (0->111)
t242 = sin(qJ(1));
t247 = cos(qJ(1));
t229 = -t247 * g(1) - t242 * g(2);
t241 = sin(qJ(2));
t246 = cos(qJ(2));
t212 = -t246 * g(3) - t241 * t229;
t248 = qJD(1) ^ 2;
t204 = (t241 * t246 * t248 + qJDD(2)) * pkin(1) + t212;
t213 = -t241 * g(3) + t246 * t229;
t207 = (-t246 ^ 2 * t248 - qJD(2) ^ 2) * pkin(1) + t213;
t240 = sin(qJ(3));
t245 = cos(qJ(3));
t180 = t240 * t204 + t245 * t207;
t219 = (t240 * t246 + t241 * t245) * qJD(1);
t260 = qJD(1) * qJD(2);
t223 = t241 * qJDD(1) + t246 * t260;
t259 = t241 * t260;
t224 = t246 * qJDD(1) - t259;
t191 = -t219 * qJD(3) - t240 * t223 + t245 * t224;
t261 = qJD(1) * t246;
t262 = qJD(1) * t241;
t218 = -t240 * t262 + t245 * t261;
t198 = -t218 * mrSges(4,1) + t219 * mrSges(4,2);
t236 = qJD(2) + qJD(3);
t209 = t236 * mrSges(4,1) - t219 * mrSges(4,3);
t235 = qJDD(2) + qJDD(3);
t192 = t218 * qJD(3) + t245 * t223 + t240 * t224;
t228 = t242 * g(1) - t247 * g(2);
t201 = -t228 + (-t224 + t259) * pkin(1);
t160 = (-t218 * t236 - t192) * pkin(5) + (t219 * t236 - t191) * pkin(2) + t201;
t199 = -t218 * pkin(2) - t219 * pkin(5);
t234 = t236 ^ 2;
t168 = -t234 * pkin(2) + t235 * pkin(5) + t218 * t199 + t180;
t239 = sin(qJ(4));
t244 = cos(qJ(4));
t149 = t244 * t160 - t239 * t168;
t190 = qJDD(4) - t191;
t205 = -t239 * t219 + t244 * t236;
t206 = t244 * t219 + t239 * t236;
t147 = (t205 * t206 + t190) * pkin(3) + t149;
t150 = t239 * t160 + t244 * t168;
t214 = qJD(4) - t218;
t148 = (-t205 ^ 2 - t214 ^ 2) * pkin(3) + t150;
t238 = sin(qJ(5));
t243 = cos(qJ(5));
t145 = t243 * t147 - t238 * t148;
t170 = -t206 * qJD(4) - t239 * t192 + t244 * t235;
t171 = t205 * qJD(4) + t244 * t192 + t239 * t235;
t181 = t243 * t205 - t238 * t206;
t156 = t181 * qJD(5) + t238 * t170 + t243 * t171;
t182 = t238 * t205 + t243 * t206;
t166 = -t181 * mrSges(6,1) + t182 * mrSges(6,2);
t210 = qJD(5) + t214;
t172 = -t210 * mrSges(6,2) + t181 * mrSges(6,3);
t186 = qJDD(5) + t190;
t142 = m(6) * t145 + t186 * mrSges(6,1) - t156 * mrSges(6,3) - t182 * t166 + t210 * t172;
t146 = t238 * t147 + t243 * t148;
t155 = -t182 * qJD(5) + t243 * t170 - t238 * t171;
t173 = t210 * mrSges(6,1) - t182 * mrSges(6,3);
t143 = m(6) * t146 - t186 * mrSges(6,2) + t155 * mrSges(6,3) + t181 * t166 - t210 * t173;
t133 = t243 * t142 + t238 * t143;
t183 = -t205 * mrSges(5,1) + t206 * mrSges(5,2);
t193 = -t214 * mrSges(5,2) + t205 * mrSges(5,3);
t131 = m(5) * t149 + t190 * mrSges(5,1) - t171 * mrSges(5,3) - t206 * t183 + t214 * t193 + t133;
t194 = t214 * mrSges(5,1) - t206 * mrSges(5,3);
t132 = m(5) * t150 - t190 * mrSges(5,2) + t170 * mrSges(5,3) - t238 * t142 + t243 * t143 + t205 * t183 - t214 * t194;
t258 = -t239 * t131 + t244 * t132;
t120 = m(4) * t180 - t235 * mrSges(4,2) + t191 * mrSges(4,3) + t218 * t198 - t236 * t209 + t258;
t179 = t245 * t204 - t240 * t207;
t208 = -t236 * mrSges(4,2) + t218 * mrSges(4,3);
t167 = -t235 * pkin(2) - t234 * pkin(5) + t219 * t199 - t179;
t151 = (t206 * t214 - t170) * pkin(3) + t167;
t255 = m(6) * t151 - t155 * mrSges(6,1) + t156 * mrSges(6,2) - t181 * t172 + t182 * t173;
t251 = -m(5) * t167 + t170 * mrSges(5,1) - t171 * mrSges(5,2) + t205 * t193 - t206 * t194 - t255;
t137 = m(4) * t179 + t235 * mrSges(4,1) - t192 * mrSges(4,3) - t219 * t198 + t236 * t208 + t251;
t117 = t240 * t120 + t245 * t137;
t216 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t241 + Ifges(3,2) * t246) * qJD(1);
t217 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t241 + Ifges(3,4) * t246) * qJD(1);
t161 = Ifges(6,5) * t182 + Ifges(6,6) * t181 + Ifges(6,3) * t210;
t163 = Ifges(6,1) * t182 + Ifges(6,4) * t181 + Ifges(6,5) * t210;
t134 = -mrSges(6,1) * t151 + mrSges(6,3) * t146 + Ifges(6,4) * t156 + Ifges(6,2) * t155 + Ifges(6,6) * t186 - t182 * t161 + t210 * t163;
t162 = Ifges(6,4) * t182 + Ifges(6,2) * t181 + Ifges(6,6) * t210;
t135 = mrSges(6,2) * t151 - mrSges(6,3) * t145 + Ifges(6,1) * t156 + Ifges(6,4) * t155 + Ifges(6,5) * t186 + t181 * t161 - t210 * t162;
t174 = Ifges(5,5) * t206 + Ifges(5,6) * t205 + Ifges(5,3) * t214;
t176 = Ifges(5,1) * t206 + Ifges(5,4) * t205 + Ifges(5,5) * t214;
t125 = -mrSges(5,1) * t167 + mrSges(5,3) * t150 + Ifges(5,4) * t171 + Ifges(5,2) * t170 + Ifges(5,6) * t190 - pkin(3) * t255 + t243 * t134 + t238 * t135 - t206 * t174 + t214 * t176;
t175 = Ifges(5,4) * t206 + Ifges(5,2) * t205 + Ifges(5,6) * t214;
t127 = mrSges(5,2) * t167 - mrSges(5,3) * t149 + Ifges(5,1) * t171 + Ifges(5,4) * t170 + Ifges(5,5) * t190 - t238 * t134 + t243 * t135 + t205 * t174 - t214 * t175;
t196 = Ifges(4,4) * t219 + Ifges(4,2) * t218 + Ifges(4,6) * t236;
t197 = Ifges(4,1) * t219 + Ifges(4,4) * t218 + Ifges(4,5) * t236;
t253 = -mrSges(4,1) * t179 + mrSges(4,2) * t180 - Ifges(4,5) * t192 - Ifges(4,6) * t191 - Ifges(4,3) * t235 - pkin(2) * t251 - pkin(5) * t258 - t244 * t125 - t239 * t127 - t219 * t196 + t218 * t197;
t263 = mrSges(3,1) * t212 - mrSges(3,2) * t213 + Ifges(3,5) * t223 + Ifges(3,6) * t224 + Ifges(3,3) * qJDD(2) + pkin(1) * t117 + (t241 * t216 - t246 * t217) * qJD(1) - t253;
t122 = t244 * t131 + t239 * t132;
t195 = Ifges(4,5) * t219 + Ifges(4,6) * t218 + Ifges(4,3) * t236;
t114 = mrSges(4,2) * t201 - mrSges(4,3) * t179 + Ifges(4,1) * t192 + Ifges(4,4) * t191 + Ifges(4,5) * t235 - pkin(5) * t122 - t239 * t125 + t244 * t127 + t218 * t195 - t236 * t196;
t254 = -mrSges(6,1) * t145 + mrSges(6,2) * t146 - Ifges(6,5) * t156 - Ifges(6,6) * t155 - Ifges(6,3) * t186 - t182 * t162 + t181 * t163;
t249 = mrSges(5,1) * t149 - mrSges(5,2) * t150 + Ifges(5,5) * t171 + Ifges(5,6) * t170 + Ifges(5,3) * t190 + pkin(3) * t133 + t206 * t175 - t205 * t176 - t254;
t116 = -mrSges(4,1) * t201 + mrSges(4,3) * t180 + Ifges(4,4) * t192 + Ifges(4,2) * t191 + Ifges(4,6) * t235 - pkin(2) * t122 - t219 * t195 + t236 * t197 - t249;
t215 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t241 + Ifges(3,6) * t246) * qJD(1);
t252 = m(4) * t201 - t191 * mrSges(4,1) + t192 * mrSges(4,2) - t218 * t208 + t219 * t209 + t122;
t110 = mrSges(3,1) * t228 + mrSges(3,3) * t213 + Ifges(3,4) * t223 + Ifges(3,2) * t224 + Ifges(3,6) * qJDD(2) - pkin(1) * t252 + qJD(2) * t217 + t240 * t114 + t245 * t116 - t215 * t262;
t112 = -mrSges(3,2) * t228 - mrSges(3,3) * t212 + Ifges(3,1) * t223 + Ifges(3,4) * t224 + Ifges(3,5) * qJDD(2) - qJD(2) * t216 + t245 * t114 - t240 * t116 + t215 * t261;
t256 = mrSges(2,1) * t228 - mrSges(2,2) * t229 + Ifges(2,3) * qJDD(1) + t246 * t110 + t241 * t112;
t227 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t261;
t226 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t262;
t222 = (-mrSges(3,1) * t246 + mrSges(3,2) * t241) * qJD(1);
t118 = qJDD(1) * mrSges(2,1) + t224 * mrSges(3,1) - t248 * mrSges(2,2) - t223 * mrSges(3,2) + (m(2) + m(3)) * t228 + (-t226 * t241 + t227 * t246) * qJD(1) - t252;
t115 = m(2) * t229 - qJDD(1) * mrSges(2,2) - t248 * mrSges(2,1) + t246 * (m(3) * t213 - qJDD(2) * mrSges(3,2) + t224 * mrSges(3,3) - qJD(2) * t226 + t245 * t120 - t240 * t137 + t222 * t261) - t241 * (m(3) * t212 + qJDD(2) * mrSges(3,1) - t223 * mrSges(3,3) + qJD(2) * t227 - t222 * t262 + t117);
t113 = mrSges(2,1) * g(3) + mrSges(2,3) * t229 + t248 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - t263;
t108 = -mrSges(2,2) * g(3) - mrSges(2,3) * t228 + Ifges(2,5) * qJDD(1) - t248 * Ifges(2,6) - t241 * t110 + t246 * t112;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t247 * t108 - t242 * t113 - pkin(4) * (t242 * t115 + t247 * t118), t108, t112, t114, t127, t135; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t242 * t108 + t247 * t113 + pkin(4) * (t247 * t115 - t242 * t118), t113, t110, t116, t125, t134; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t256, t256, t263, -t253, t249, -t254;];
m_new  = t1;
