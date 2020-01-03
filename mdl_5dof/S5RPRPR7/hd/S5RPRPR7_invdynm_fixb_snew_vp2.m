% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:00
% EndTime: 2019-12-31 18:19:09
% DurationCPUTime: 5.44s
% Computational Cost: add. (69743->267), mult. (147827->343), div. (0->0), fcn. (91949->10), ass. (0->110)
t226 = sin(pkin(9));
t228 = cos(pkin(9));
t231 = sin(qJ(3));
t234 = cos(qJ(3));
t196 = (t226 * t231 - t228 * t234) * qJD(1);
t232 = sin(qJ(1));
t235 = cos(qJ(1));
t216 = t232 * g(1) - t235 * g(2);
t207 = qJDD(1) * pkin(1) + t216;
t217 = -t235 * g(1) - t232 * g(2);
t237 = qJD(1) ^ 2;
t209 = -t237 * pkin(1) + t217;
t227 = sin(pkin(8));
t229 = cos(pkin(8));
t187 = t227 * t207 + t229 * t209;
t177 = -t237 * pkin(2) + qJDD(1) * pkin(6) + t187;
t225 = -g(3) + qJDD(2);
t168 = -t231 * t177 + t234 * t225;
t255 = qJD(1) * qJD(3);
t254 = t234 * t255;
t210 = t231 * qJDD(1) + t254;
t156 = (-t210 + t254) * qJ(4) + (t231 * t234 * t237 + qJDD(3)) * pkin(3) + t168;
t169 = t234 * t177 + t231 * t225;
t211 = t234 * qJDD(1) - t231 * t255;
t257 = qJD(1) * t231;
t213 = qJD(3) * pkin(3) - qJ(4) * t257;
t224 = t234 ^ 2;
t157 = -t224 * t237 * pkin(3) + t211 * qJ(4) - qJD(3) * t213 + t169;
t258 = 2 * qJD(4);
t152 = t226 * t156 + t228 * t157 - t196 * t258;
t197 = (t226 * t234 + t228 * t231) * qJD(1);
t179 = t196 * mrSges(5,1) + t197 * mrSges(5,2);
t188 = -t226 * t210 + t228 * t211;
t193 = qJD(3) * mrSges(5,1) - t197 * mrSges(5,3);
t180 = t196 * pkin(4) - t197 * pkin(7);
t236 = qJD(3) ^ 2;
t149 = -t236 * pkin(4) + qJDD(3) * pkin(7) - t196 * t180 + t152;
t186 = t229 * t207 - t227 * t209;
t246 = -qJDD(1) * pkin(2) - t186;
t158 = -t211 * pkin(3) + qJDD(4) + t213 * t257 + (-qJ(4) * t224 - pkin(6)) * t237 + t246;
t189 = t228 * t210 + t226 * t211;
t153 = (qJD(3) * t196 - t189) * pkin(7) + (qJD(3) * t197 - t188) * pkin(4) + t158;
t230 = sin(qJ(5));
t233 = cos(qJ(5));
t146 = -t230 * t149 + t233 * t153;
t190 = t233 * qJD(3) - t230 * t197;
t165 = t190 * qJD(5) + t230 * qJDD(3) + t233 * t189;
t191 = t230 * qJD(3) + t233 * t197;
t167 = -t190 * mrSges(6,1) + t191 * mrSges(6,2);
t195 = qJD(5) + t196;
t170 = -t195 * mrSges(6,2) + t190 * mrSges(6,3);
t185 = qJDD(5) - t188;
t142 = m(6) * t146 + t185 * mrSges(6,1) - t165 * mrSges(6,3) - t191 * t167 + t195 * t170;
t147 = t233 * t149 + t230 * t153;
t164 = -t191 * qJD(5) + t233 * qJDD(3) - t230 * t189;
t171 = t195 * mrSges(6,1) - t191 * mrSges(6,3);
t143 = m(6) * t147 - t185 * mrSges(6,2) + t164 * mrSges(6,3) + t190 * t167 - t195 * t171;
t250 = -t230 * t142 + t233 * t143;
t129 = m(5) * t152 - qJDD(3) * mrSges(5,2) + t188 * mrSges(5,3) - qJD(3) * t193 - t196 * t179 + t250;
t249 = -t228 * t156 + t226 * t157;
t151 = -0.2e1 * qJD(4) * t197 - t249;
t192 = -qJD(3) * mrSges(5,2) - t196 * mrSges(5,3);
t148 = -qJDD(3) * pkin(4) - t236 * pkin(7) + (t258 + t180) * t197 + t249;
t244 = -m(6) * t148 + t164 * mrSges(6,1) - t165 * mrSges(6,2) + t190 * t170 - t191 * t171;
t138 = m(5) * t151 + qJDD(3) * mrSges(5,1) - t189 * mrSges(5,3) + qJD(3) * t192 - t197 * t179 + t244;
t123 = t226 * t129 + t228 * t138;
t202 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t231 + Ifges(4,2) * t234) * qJD(1);
t203 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t231 + Ifges(4,4) * t234) * qJD(1);
t159 = Ifges(6,5) * t191 + Ifges(6,6) * t190 + Ifges(6,3) * t195;
t161 = Ifges(6,1) * t191 + Ifges(6,4) * t190 + Ifges(6,5) * t195;
t135 = -mrSges(6,1) * t148 + mrSges(6,3) * t147 + Ifges(6,4) * t165 + Ifges(6,2) * t164 + Ifges(6,6) * t185 - t191 * t159 + t195 * t161;
t160 = Ifges(6,4) * t191 + Ifges(6,2) * t190 + Ifges(6,6) * t195;
t136 = mrSges(6,2) * t148 - mrSges(6,3) * t146 + Ifges(6,1) * t165 + Ifges(6,4) * t164 + Ifges(6,5) * t185 + t190 * t159 - t195 * t160;
t174 = Ifges(5,4) * t197 - Ifges(5,2) * t196 + Ifges(5,6) * qJD(3);
t175 = Ifges(5,1) * t197 - Ifges(5,4) * t196 + Ifges(5,5) * qJD(3);
t241 = -mrSges(5,1) * t151 + mrSges(5,2) * t152 - Ifges(5,5) * t189 - Ifges(5,6) * t188 - Ifges(5,3) * qJDD(3) - pkin(4) * t244 - pkin(7) * t250 - t233 * t135 - t230 * t136 - t197 * t174 - t196 * t175;
t259 = mrSges(4,1) * t168 - mrSges(4,2) * t169 + Ifges(4,5) * t210 + Ifges(4,6) * t211 + Ifges(4,3) * qJDD(3) + pkin(3) * t123 + (t231 * t202 - t234 * t203) * qJD(1) - t241;
t208 = (-mrSges(4,1) * t234 + mrSges(4,2) * t231) * qJD(1);
t256 = qJD(1) * t234;
t215 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t256;
t121 = m(4) * t168 + qJDD(3) * mrSges(4,1) - t210 * mrSges(4,3) + qJD(3) * t215 - t208 * t257 + t123;
t214 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t257;
t251 = t228 * t129 - t226 * t138;
t122 = m(4) * t169 - qJDD(3) * mrSges(4,2) + t211 * mrSges(4,3) - qJD(3) * t214 + t208 * t256 + t251;
t252 = -t231 * t121 + t234 * t122;
t113 = m(3) * t187 - t237 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t252;
t176 = -t237 * pkin(6) + t246;
t131 = t233 * t142 + t230 * t143;
t243 = m(5) * t158 - t188 * mrSges(5,1) + t189 * mrSges(5,2) + t196 * t192 + t197 * t193 + t131;
t239 = -m(4) * t176 + t211 * mrSges(4,1) - t210 * mrSges(4,2) - t214 * t257 + t215 * t256 - t243;
t125 = m(3) * t186 + qJDD(1) * mrSges(3,1) - t237 * mrSges(3,2) + t239;
t110 = t227 * t113 + t229 * t125;
t115 = t234 * t121 + t231 * t122;
t253 = t229 * t113 - t227 * t125;
t173 = Ifges(5,5) * t197 - Ifges(5,6) * t196 + Ifges(5,3) * qJD(3);
t116 = mrSges(5,2) * t158 - mrSges(5,3) * t151 + Ifges(5,1) * t189 + Ifges(5,4) * t188 + Ifges(5,5) * qJDD(3) - pkin(7) * t131 - qJD(3) * t174 - t230 * t135 + t233 * t136 - t196 * t173;
t240 = mrSges(6,1) * t146 - mrSges(6,2) * t147 + Ifges(6,5) * t165 + Ifges(6,6) * t164 + Ifges(6,3) * t185 + t191 * t160 - t190 * t161;
t117 = -mrSges(5,1) * t158 + mrSges(5,3) * t152 + Ifges(5,4) * t189 + Ifges(5,2) * t188 + Ifges(5,6) * qJDD(3) - pkin(4) * t131 + qJD(3) * t175 - t197 * t173 - t240;
t201 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t231 + Ifges(4,6) * t234) * qJD(1);
t104 = -mrSges(4,1) * t176 + mrSges(4,3) * t169 + Ifges(4,4) * t210 + Ifges(4,2) * t211 + Ifges(4,6) * qJDD(3) - pkin(3) * t243 + qJ(4) * t251 + qJD(3) * t203 + t226 * t116 + t228 * t117 - t201 * t257;
t106 = mrSges(4,2) * t176 - mrSges(4,3) * t168 + Ifges(4,1) * t210 + Ifges(4,4) * t211 + Ifges(4,5) * qJDD(3) - qJ(4) * t123 - qJD(3) * t202 + t228 * t116 - t226 * t117 + t201 * t256;
t245 = mrSges(3,1) * t186 - mrSges(3,2) * t187 + Ifges(3,3) * qJDD(1) + pkin(2) * t239 + pkin(6) * t252 + t234 * t104 + t231 * t106;
t242 = mrSges(2,1) * t216 - mrSges(2,2) * t217 + Ifges(2,3) * qJDD(1) + pkin(1) * t110 + t245;
t108 = m(2) * t217 - t237 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t253;
t107 = m(2) * t216 + qJDD(1) * mrSges(2,1) - t237 * mrSges(2,2) + t110;
t102 = -mrSges(3,1) * t225 + mrSges(3,3) * t187 + t237 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t115 - t259;
t101 = mrSges(3,2) * t225 - mrSges(3,3) * t186 + Ifges(3,5) * qJDD(1) - t237 * Ifges(3,6) - pkin(6) * t115 - t231 * t104 + t234 * t106;
t100 = -mrSges(2,2) * g(3) - mrSges(2,3) * t216 + Ifges(2,5) * qJDD(1) - t237 * Ifges(2,6) - qJ(2) * t110 + t229 * t101 - t227 * t102;
t99 = Ifges(2,6) * qJDD(1) + t237 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t217 + t227 * t101 + t229 * t102 - pkin(1) * (m(3) * t225 + t115) + qJ(2) * t253;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t235 * t100 - t232 * t99 - pkin(5) * (t235 * t107 + t232 * t108), t100, t101, t106, t116, t136; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t232 * t100 + t235 * t99 + pkin(5) * (-t232 * t107 + t235 * t108), t99, t102, t104, t117, t135; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t242, t242, t245, t259, -t241, t240;];
m_new = t1;
