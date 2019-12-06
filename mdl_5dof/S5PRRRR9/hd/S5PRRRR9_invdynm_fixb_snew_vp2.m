% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:20
% EndTime: 2019-12-05 17:19:39
% DurationCPUTime: 7.92s
% Computational Cost: add. (134409->265), mult. (258510->344), div. (0->0), fcn. (179407->12), ass. (0->117)
t227 = sin(pkin(10));
t229 = cos(pkin(10));
t218 = t227 * g(1) - t229 * g(2);
t219 = -t229 * g(1) - t227 * g(2);
t226 = -g(3) + qJDD(1);
t238 = cos(qJ(2));
t230 = cos(pkin(5));
t234 = sin(qJ(2));
t255 = t230 * t234;
t228 = sin(pkin(5));
t256 = t228 * t234;
t182 = t218 * t255 + t238 * t219 + t226 * t256;
t240 = qJD(2) ^ 2;
t178 = -t240 * pkin(2) + qJDD(2) * pkin(7) + t182;
t197 = -t228 * t218 + t230 * t226;
t233 = sin(qJ(3));
t237 = cos(qJ(3));
t173 = t237 * t178 + t233 * t197;
t214 = (-pkin(3) * t237 - pkin(8) * t233) * qJD(2);
t239 = qJD(3) ^ 2;
t253 = t237 * qJD(2);
t164 = -t239 * pkin(3) + qJDD(3) * pkin(8) + t214 * t253 + t173;
t181 = -t234 * t219 + (t218 * t230 + t226 * t228) * t238;
t177 = -qJDD(2) * pkin(2) - t240 * pkin(7) - t181;
t252 = qJD(2) * qJD(3);
t251 = t237 * t252;
t215 = t233 * qJDD(2) + t251;
t225 = t233 * t252;
t216 = t237 * qJDD(2) - t225;
t167 = (-t215 - t251) * pkin(8) + (-t216 + t225) * pkin(3) + t177;
t232 = sin(qJ(4));
t236 = cos(qJ(4));
t153 = -t232 * t164 + t236 * t167;
t254 = qJD(2) * t233;
t211 = t236 * qJD(3) - t232 * t254;
t189 = t211 * qJD(4) + t232 * qJDD(3) + t236 * t215;
t208 = qJDD(4) - t216;
t212 = t232 * qJD(3) + t236 * t254;
t224 = qJD(4) - t253;
t151 = (t211 * t224 - t189) * pkin(9) + (t211 * t212 + t208) * pkin(4) + t153;
t154 = t236 * t164 + t232 * t167;
t188 = -t212 * qJD(4) + t236 * qJDD(3) - t232 * t215;
t196 = t224 * pkin(4) - t212 * pkin(9);
t207 = t211 ^ 2;
t152 = -t207 * pkin(4) + t188 * pkin(9) - t224 * t196 + t154;
t231 = sin(qJ(5));
t235 = cos(qJ(5));
t149 = t235 * t151 - t231 * t152;
t190 = t235 * t211 - t231 * t212;
t161 = t190 * qJD(5) + t231 * t188 + t235 * t189;
t191 = t231 * t211 + t235 * t212;
t174 = -t190 * mrSges(6,1) + t191 * mrSges(6,2);
t223 = qJD(5) + t224;
t179 = -t223 * mrSges(6,2) + t190 * mrSges(6,3);
t204 = qJDD(5) + t208;
t145 = m(6) * t149 + t204 * mrSges(6,1) - t161 * mrSges(6,3) - t191 * t174 + t223 * t179;
t150 = t231 * t151 + t235 * t152;
t160 = -t191 * qJD(5) + t235 * t188 - t231 * t189;
t180 = t223 * mrSges(6,1) - t191 * mrSges(6,3);
t146 = m(6) * t150 - t204 * mrSges(6,2) + t160 * mrSges(6,3) + t190 * t174 - t223 * t180;
t137 = t235 * t145 + t231 * t146;
t192 = -t211 * mrSges(5,1) + t212 * mrSges(5,2);
t194 = -t224 * mrSges(5,2) + t211 * mrSges(5,3);
t135 = m(5) * t153 + t208 * mrSges(5,1) - t189 * mrSges(5,3) - t212 * t192 + t224 * t194 + t137;
t195 = t224 * mrSges(5,1) - t212 * mrSges(5,3);
t249 = -t231 * t145 + t235 * t146;
t136 = m(5) * t154 - t208 * mrSges(5,2) + t188 * mrSges(5,3) + t211 * t192 - t224 * t195 + t249;
t133 = -t232 * t135 + t236 * t136;
t213 = (-mrSges(4,1) * t237 + mrSges(4,2) * t233) * qJD(2);
t220 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t254;
t131 = m(4) * t173 - qJDD(3) * mrSges(4,2) + t216 * mrSges(4,3) - qJD(3) * t220 + t213 * t253 + t133;
t172 = -t233 * t178 + t237 * t197;
t163 = -qJDD(3) * pkin(3) - t239 * pkin(8) + t214 * t254 - t172;
t155 = -t188 * pkin(4) - t207 * pkin(9) + t212 * t196 + t163;
t245 = m(6) * t155 - t160 * mrSges(6,1) + t161 * mrSges(6,2) - t190 * t179 + t191 * t180;
t147 = -m(5) * t163 + t188 * mrSges(5,1) - t189 * mrSges(5,2) + t211 * t194 - t212 * t195 - t245;
t221 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t253;
t141 = m(4) * t172 + qJDD(3) * mrSges(4,1) - t215 * mrSges(4,3) + qJD(3) * t221 - t213 * t254 + t147;
t125 = t233 * t131 + t237 * t141;
t168 = Ifges(6,5) * t191 + Ifges(6,6) * t190 + Ifges(6,3) * t223;
t170 = Ifges(6,1) * t191 + Ifges(6,4) * t190 + Ifges(6,5) * t223;
t138 = -mrSges(6,1) * t155 + mrSges(6,3) * t150 + Ifges(6,4) * t161 + Ifges(6,2) * t160 + Ifges(6,6) * t204 - t191 * t168 + t223 * t170;
t169 = Ifges(6,4) * t191 + Ifges(6,2) * t190 + Ifges(6,6) * t223;
t139 = mrSges(6,2) * t155 - mrSges(6,3) * t149 + Ifges(6,1) * t161 + Ifges(6,4) * t160 + Ifges(6,5) * t204 + t190 * t168 - t223 * t169;
t183 = Ifges(5,5) * t212 + Ifges(5,6) * t211 + Ifges(5,3) * t224;
t185 = Ifges(5,1) * t212 + Ifges(5,4) * t211 + Ifges(5,5) * t224;
t123 = -mrSges(5,1) * t163 + mrSges(5,3) * t154 + Ifges(5,4) * t189 + Ifges(5,2) * t188 + Ifges(5,6) * t208 - pkin(4) * t245 + pkin(9) * t249 + t235 * t138 + t231 * t139 - t212 * t183 + t224 * t185;
t184 = Ifges(5,4) * t212 + Ifges(5,2) * t211 + Ifges(5,6) * t224;
t126 = mrSges(5,2) * t163 - mrSges(5,3) * t153 + Ifges(5,1) * t189 + Ifges(5,4) * t188 + Ifges(5,5) * t208 - pkin(9) * t137 - t231 * t138 + t235 * t139 + t211 * t183 - t224 * t184;
t202 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t233 + Ifges(4,2) * t237) * qJD(2);
t203 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t233 + Ifges(4,4) * t237) * qJD(2);
t260 = mrSges(4,1) * t172 - mrSges(4,2) * t173 + Ifges(4,5) * t215 + Ifges(4,6) * t216 + Ifges(4,3) * qJDD(3) + pkin(3) * t147 + pkin(8) * t133 + t236 * t123 + t232 * t126 + (t233 * t202 - t237 * t203) * qJD(2);
t110 = -mrSges(3,1) * t197 + mrSges(3,3) * t182 + t240 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t125 - t260;
t250 = t237 * t131 - t233 * t141;
t122 = m(3) * t182 - t240 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t250;
t132 = t236 * t135 + t232 * t136;
t243 = -m(4) * t177 + t216 * mrSges(4,1) - t215 * mrSges(4,2) - t220 * t254 + t221 * t253 - t132;
t128 = m(3) * t181 + qJDD(2) * mrSges(3,1) - t240 * mrSges(3,2) + t243;
t118 = t238 * t122 - t234 * t128;
t261 = pkin(6) * t118 + t110 * t238;
t257 = t128 * t238;
t124 = m(3) * t197 + t125;
t114 = t122 * t255 - t228 * t124 + t230 * t257;
t201 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t233 + Ifges(4,6) * t237) * qJD(2);
t115 = mrSges(4,2) * t177 - mrSges(4,3) * t172 + Ifges(4,1) * t215 + Ifges(4,4) * t216 + Ifges(4,5) * qJDD(3) - pkin(8) * t132 - qJD(3) * t202 - t232 * t123 + t236 * t126 + t201 * t253;
t244 = -mrSges(6,1) * t149 + mrSges(6,2) * t150 - Ifges(6,5) * t161 - Ifges(6,6) * t160 - Ifges(6,3) * t204 - t191 * t169 + t190 * t170;
t241 = mrSges(5,1) * t153 - mrSges(5,2) * t154 + Ifges(5,5) * t189 + Ifges(5,6) * t188 + Ifges(5,3) * t208 + pkin(4) * t137 + t212 * t184 - t211 * t185 - t244;
t119 = -mrSges(4,1) * t177 + mrSges(4,3) * t173 + Ifges(4,4) * t215 + Ifges(4,2) * t216 + Ifges(4,6) * qJDD(3) - pkin(3) * t132 + qJD(3) * t203 - t201 * t254 - t241;
t106 = mrSges(3,1) * t181 - mrSges(3,2) * t182 + Ifges(3,3) * qJDD(2) + pkin(2) * t243 + pkin(7) * t250 + t233 * t115 + t237 * t119;
t108 = mrSges(3,2) * t197 - mrSges(3,3) * t181 + Ifges(3,5) * qJDD(2) - t240 * Ifges(3,6) - pkin(7) * t125 + t237 * t115 - t233 * t119;
t246 = mrSges(2,1) * t218 - mrSges(2,2) * t219 + pkin(1) * t114 + t230 * t106 + t108 * t256 + t261 * t228;
t116 = m(2) * t219 + t118;
t113 = t230 * t124 + (t122 * t234 + t257) * t228;
t111 = m(2) * t218 + t114;
t104 = mrSges(2,2) * t226 - mrSges(2,3) * t218 + t238 * t108 - t234 * t110 + (-t113 * t228 - t114 * t230) * pkin(6);
t103 = -mrSges(2,1) * t226 + mrSges(2,3) * t219 - pkin(1) * t113 - t228 * t106 + (t108 * t234 + t261) * t230;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t229 * t104 - t227 * t103 - qJ(1) * (t229 * t111 + t227 * t116), t104, t108, t115, t126, t139; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t227 * t104 + t229 * t103 + qJ(1) * (-t227 * t111 + t229 * t116), t103, t110, t119, t123, t138; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t246, t246, t106, t260, t241, -t244;];
m_new = t1;
