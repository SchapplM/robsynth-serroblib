% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRP8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:22
% EndTime: 2019-12-05 16:58:34
% DurationCPUTime: 4.58s
% Computational Cost: add. (55603->261), mult. (104453->326), div. (0->0), fcn. (67674->10), ass. (0->106)
t224 = sin(pkin(9));
t226 = cos(pkin(9));
t214 = t224 * g(1) - t226 * g(2);
t215 = -t226 * g(1) - t224 * g(2);
t223 = -g(3) + qJDD(1);
t232 = cos(qJ(2));
t227 = cos(pkin(5));
t230 = sin(qJ(2));
t255 = t227 * t230;
t225 = sin(pkin(5));
t256 = t225 * t230;
t165 = t214 * t255 + t232 * t215 + t223 * t256;
t234 = qJD(2) ^ 2;
t163 = -t234 * pkin(2) + qJDD(2) * pkin(7) + t165;
t192 = -t225 * t214 + t227 * t223;
t229 = sin(qJ(3));
t231 = cos(qJ(3));
t157 = t231 * t163 + t229 * t192;
t210 = (-pkin(3) * t231 - pkin(8) * t229) * qJD(2);
t233 = qJD(3) ^ 2;
t250 = t231 * qJD(2);
t153 = -t233 * pkin(3) + qJDD(3) * pkin(8) + t210 * t250 + t157;
t164 = -t230 * t215 + (t214 * t227 + t223 * t225) * t232;
t162 = -qJDD(2) * pkin(2) - t234 * pkin(7) - t164;
t249 = qJD(2) * qJD(3);
t246 = t231 * t249;
t211 = t229 * qJDD(2) + t246;
t247 = t229 * t249;
t212 = t231 * qJDD(2) - t247;
t155 = (-t211 - t246) * pkin(8) + (-t212 + t247) * pkin(3) + t162;
t228 = sin(qJ(4));
t261 = cos(qJ(4));
t150 = t261 * t153 + t228 * t155;
t251 = qJD(2) * t229;
t208 = t228 * qJD(3) + t261 * t251;
t180 = t208 * qJD(4) - t261 * qJDD(3) + t228 * t211;
t220 = qJD(4) - t250;
t189 = t220 * mrSges(5,1) - t208 * mrSges(5,3);
t204 = qJDD(4) - t212;
t207 = -t261 * qJD(3) + t228 * t251;
t184 = t207 * pkin(4) - t208 * qJ(5);
t219 = t220 ^ 2;
t145 = -t219 * pkin(4) + t204 * qJ(5) + 0.2e1 * qJD(5) * t220 - t207 * t184 + t150;
t190 = -t220 * mrSges(6,1) + t208 * mrSges(6,2);
t248 = m(6) * t145 + t204 * mrSges(6,3) + t220 * t190;
t185 = t207 * mrSges(6,1) - t208 * mrSges(6,3);
t252 = -t207 * mrSges(5,1) - t208 * mrSges(5,2) - t185;
t259 = -mrSges(5,3) - mrSges(6,2);
t137 = m(5) * t150 - t204 * mrSges(5,2) + t259 * t180 - t220 * t189 + t252 * t207 + t248;
t149 = -t228 * t153 + t261 * t155;
t181 = -t207 * qJD(4) + t228 * qJDD(3) + t261 * t211;
t188 = -t220 * mrSges(5,2) - t207 * mrSges(5,3);
t147 = -t204 * pkin(4) - t219 * qJ(5) + t208 * t184 + qJDD(5) - t149;
t191 = -t207 * mrSges(6,2) + t220 * mrSges(6,3);
t244 = -m(6) * t147 + t204 * mrSges(6,1) + t220 * t191;
t138 = m(5) * t149 + t204 * mrSges(5,1) + t259 * t181 + t220 * t188 + t252 * t208 + t244;
t133 = t261 * t137 - t228 * t138;
t209 = (-mrSges(4,1) * t231 + mrSges(4,2) * t229) * qJD(2);
t216 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t251;
t129 = m(4) * t157 - qJDD(3) * mrSges(4,2) + t212 * mrSges(4,3) - qJD(3) * t216 + t209 * t250 + t133;
t156 = -t229 * t163 + t231 * t192;
t152 = -qJDD(3) * pkin(3) - t233 * pkin(8) + t210 * t251 - t156;
t148 = -0.2e1 * qJD(5) * t208 + (t207 * t220 - t181) * qJ(5) + (t208 * t220 + t180) * pkin(4) + t152;
t142 = m(6) * t148 + t180 * mrSges(6,1) - t181 * mrSges(6,3) - t208 * t190 + t207 * t191;
t139 = -m(5) * t152 - t180 * mrSges(5,1) - t181 * mrSges(5,2) - t207 * t188 - t208 * t189 - t142;
t217 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t250;
t135 = m(4) * t156 + qJDD(3) * mrSges(4,1) - t211 * mrSges(4,3) + qJD(3) * t217 - t209 * t251 + t139;
t124 = t229 * t129 + t231 * t135;
t170 = Ifges(6,1) * t208 + Ifges(6,4) * t220 + Ifges(6,5) * t207;
t171 = Ifges(5,1) * t208 - Ifges(5,4) * t207 + Ifges(5,5) * t220;
t243 = -mrSges(6,1) * t148 + mrSges(6,2) * t145;
t168 = Ifges(6,4) * t208 + Ifges(6,2) * t220 + Ifges(6,6) * t207;
t254 = -Ifges(5,5) * t208 + Ifges(5,6) * t207 - Ifges(5,3) * t220 - t168;
t130 = -mrSges(5,1) * t152 + mrSges(5,3) * t150 - pkin(4) * t142 + (t170 + t171) * t220 + t254 * t208 + (Ifges(5,6) - Ifges(6,6)) * t204 + (Ifges(5,4) - Ifges(6,5)) * t181 + (-Ifges(5,2) - Ifges(6,3)) * t180 + t243;
t169 = Ifges(5,4) * t208 - Ifges(5,2) * t207 + Ifges(5,6) * t220;
t166 = Ifges(6,5) * t208 + Ifges(6,6) * t220 + Ifges(6,3) * t207;
t239 = mrSges(6,2) * t147 - mrSges(6,3) * t148 + Ifges(6,1) * t181 + Ifges(6,4) * t204 + Ifges(6,5) * t180 + t220 * t166;
t131 = mrSges(5,2) * t152 - mrSges(5,3) * t149 + Ifges(5,1) * t181 - Ifges(5,4) * t180 + Ifges(5,5) * t204 - qJ(5) * t142 - t220 * t169 + t254 * t207 + t239;
t196 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t229 + Ifges(4,2) * t231) * qJD(2);
t197 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t229 + Ifges(4,4) * t231) * qJD(2);
t262 = mrSges(4,1) * t156 - mrSges(4,2) * t157 + Ifges(4,5) * t211 + Ifges(4,6) * t212 + Ifges(4,3) * qJDD(3) + pkin(3) * t139 + pkin(8) * t133 + (t196 * t229 - t197 * t231) * qJD(2) + t261 * t130 + t228 * t131;
t110 = -mrSges(3,1) * t192 + mrSges(3,3) * t165 + t234 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t124 - t262;
t245 = t231 * t129 - t229 * t135;
t122 = m(3) * t165 - t234 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t245;
t132 = t228 * t137 + t261 * t138;
t237 = -m(4) * t162 + t212 * mrSges(4,1) - t211 * mrSges(4,2) - t216 * t251 + t217 * t250 - t132;
t126 = m(3) * t164 + qJDD(2) * mrSges(3,1) - t234 * mrSges(3,2) + t237;
t118 = t232 * t122 - t230 * t126;
t264 = pkin(6) * t118 + t110 * t232;
t240 = mrSges(6,1) * t147 - mrSges(6,3) * t145 - Ifges(6,4) * t181 - Ifges(6,2) * t204 - Ifges(6,6) * t180 - t207 * t170;
t263 = -(-t169 + t166) * t208 + mrSges(5,1) * t149 - mrSges(5,2) * t150 + Ifges(5,5) * t181 - Ifges(5,6) * t180 + Ifges(5,3) * t204 + pkin(4) * (-t181 * mrSges(6,2) - t208 * t185 + t244) + qJ(5) * (-t180 * mrSges(6,2) - t207 * t185 + t248) + t207 * t171 - t240;
t257 = t126 * t232;
t123 = m(3) * t192 + t124;
t114 = t122 * t255 - t225 * t123 + t227 * t257;
t195 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t229 + Ifges(4,6) * t231) * qJD(2);
t115 = mrSges(4,2) * t162 - mrSges(4,3) * t156 + Ifges(4,1) * t211 + Ifges(4,4) * t212 + Ifges(4,5) * qJDD(3) - pkin(8) * t132 - qJD(3) * t196 - t228 * t130 + t261 * t131 + t195 * t250;
t119 = -mrSges(4,1) * t162 + mrSges(4,3) * t157 + Ifges(4,4) * t211 + Ifges(4,2) * t212 + Ifges(4,6) * qJDD(3) - pkin(3) * t132 + qJD(3) * t197 - t195 * t251 - t263;
t106 = mrSges(3,1) * t164 - mrSges(3,2) * t165 + Ifges(3,3) * qJDD(2) + pkin(2) * t237 + pkin(7) * t245 + t229 * t115 + t231 * t119;
t108 = mrSges(3,2) * t192 - mrSges(3,3) * t164 + Ifges(3,5) * qJDD(2) - t234 * Ifges(3,6) - pkin(7) * t124 + t231 * t115 - t229 * t119;
t238 = mrSges(2,1) * t214 - mrSges(2,2) * t215 + pkin(1) * t114 + t227 * t106 + t108 * t256 + t264 * t225;
t116 = m(2) * t215 + t118;
t113 = t227 * t123 + (t122 * t230 + t257) * t225;
t111 = m(2) * t214 + t114;
t104 = mrSges(2,2) * t223 - mrSges(2,3) * t214 + t232 * t108 - t230 * t110 + (-t113 * t225 - t114 * t227) * pkin(6);
t103 = -mrSges(2,1) * t223 + mrSges(2,3) * t215 - pkin(1) * t113 - t225 * t106 + (t108 * t230 + t264) * t227;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t226 * t104 - t224 * t103 - qJ(1) * (t226 * t111 + t224 * t116), t104, t108, t115, t131, -t207 * t168 + t239; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t224 * t104 + t226 * t103 + qJ(1) * (-t224 * t111 + t226 * t116), t103, t110, t119, t130, -t208 * t166 - t240; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t238, t238, t106, t262, t263, Ifges(6,5) * t181 + Ifges(6,6) * t204 + Ifges(6,3) * t180 + t208 * t168 - t220 * t170 - t243;];
m_new = t1;
