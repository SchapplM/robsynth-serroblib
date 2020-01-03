% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:31
% EndTime: 2019-12-31 20:15:35
% DurationCPUTime: 2.28s
% Computational Cost: add. (38496->225), mult. (46694->277), div. (0->0), fcn. (25346->8), ass. (0->96)
t202 = sin(qJ(1));
t206 = cos(qJ(1));
t179 = t202 * g(1) - g(2) * t206;
t173 = qJDD(1) * pkin(1) + t179;
t180 = -g(1) * t206 - g(2) * t202;
t207 = qJD(1) ^ 2;
t174 = -pkin(1) * t207 + t180;
t201 = sin(qJ(2));
t205 = cos(qJ(2));
t152 = t201 * t173 + t205 * t174;
t193 = qJDD(1) + qJDD(2);
t195 = (qJD(1) + qJD(2));
t219 = t193 * qJ(3) + (2 * qJD(3) * t195) + t152;
t231 = (-pkin(2) - pkin(7));
t230 = mrSges(3,1) - mrSges(4,2);
t229 = Ifges(3,5) - Ifges(4,4);
t228 = (-Ifges(3,6) + Ifges(4,5));
t191 = t195 ^ 2;
t138 = (t191 * t231) + t219;
t200 = sin(qJ(4));
t204 = cos(qJ(4));
t224 = qJD(4) * t195;
t168 = -t193 * t200 - t204 * t224;
t222 = t200 * t224;
t169 = t193 * t204 - t222;
t227 = t195 * t200;
t175 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t227;
t226 = t195 * t204;
t176 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t226;
t177 = (qJD(4) * pkin(4)) - pkin(8) * t226;
t198 = t200 ^ 2;
t126 = t177 * t226 - t168 * pkin(4) + (-pkin(8) * t198 + t231) * t191 + t219;
t199 = sin(qJ(5));
t203 = cos(qJ(5));
t161 = (-t199 * t200 + t203 * t204) * t195;
t132 = -t161 * qJD(5) + t168 * t203 - t169 * t199;
t160 = (-t199 * t204 - t200 * t203) * t195;
t133 = t160 * qJD(5) + t168 * t199 + t169 * t203;
t194 = qJD(4) + qJD(5);
t153 = -mrSges(6,2) * t194 + t160 * mrSges(6,3);
t154 = mrSges(6,1) * t194 - t161 * mrSges(6,3);
t217 = m(6) * t126 - t132 * mrSges(6,1) + t133 * mrSges(6,2) - t160 * t153 + t161 * t154;
t115 = -m(5) * t138 + mrSges(5,1) * t168 - t169 * mrSges(5,2) - t175 * t227 - t176 * t226 - t217;
t145 = (pkin(2) * t191) - t219;
t211 = -m(4) * t145 + (t191 * mrSges(4,2)) + t193 * mrSges(4,3) - t115;
t113 = m(3) * t152 - (mrSges(3,1) * t191) - mrSges(3,2) * t193 + t211;
t151 = t205 * t173 - t201 * t174;
t218 = -t191 * qJ(3) + qJDD(3) - t151;
t141 = t193 * t231 + t218;
t135 = t200 * g(3) + t204 * t141;
t125 = (-t169 - t222) * pkin(8) + (-t191 * t200 * t204 + qJDD(4)) * pkin(4) + t135;
t136 = -g(3) * t204 + t200 * t141;
t127 = -pkin(4) * t191 * t198 + pkin(8) * t168 - qJD(4) * t177 + t136;
t122 = t125 * t203 - t127 * t199;
t149 = -mrSges(6,1) * t160 + mrSges(6,2) * t161;
t192 = qJDD(4) + qJDD(5);
t119 = m(6) * t122 + mrSges(6,1) * t192 - t133 * mrSges(6,3) - t161 * t149 + t153 * t194;
t123 = t125 * t199 + t127 * t203;
t120 = m(6) * t123 - mrSges(6,2) * t192 + t132 * mrSges(6,3) + t160 * t149 - t154 * t194;
t108 = t119 * t203 + t120 * t199;
t167 = (mrSges(5,1) * t200 + mrSges(5,2) * t204) * t195;
t105 = m(5) * t135 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t169 + qJD(4) * t175 - t167 * t226 + t108;
t220 = -t119 * t199 + t120 * t203;
t106 = m(5) * t136 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t168 - qJD(4) * t176 - t167 * t227 + t220;
t102 = t204 * t105 + t200 * t106;
t148 = -pkin(2) * t193 + t218;
t216 = -m(4) * t148 + t191 * mrSges(4,3) - t102;
t99 = m(3) * t151 - t191 * mrSges(3,2) + t193 * t230 + t216;
t94 = t113 * t201 + t205 * t99;
t221 = t113 * t205 - t201 * t99;
t103 = -t200 * t105 + t106 * t204;
t142 = Ifges(6,5) * t161 + Ifges(6,6) * t160 + Ifges(6,3) * t194;
t144 = Ifges(6,1) * t161 + Ifges(6,4) * t160 + Ifges(6,5) * t194;
t109 = -mrSges(6,1) * t126 + mrSges(6,3) * t123 + Ifges(6,4) * t133 + Ifges(6,2) * t132 + Ifges(6,6) * t192 - t161 * t142 + t144 * t194;
t143 = Ifges(6,4) * t161 + Ifges(6,2) * t160 + Ifges(6,6) * t194;
t110 = mrSges(6,2) * t126 - mrSges(6,3) * t122 + Ifges(6,1) * t133 + Ifges(6,4) * t132 + Ifges(6,5) * t192 + t160 * t142 - t143 * t194;
t157 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t204 - Ifges(5,6) * t200) * t195;
t159 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t204 - Ifges(5,4) * t200) * t195;
t95 = -mrSges(5,1) * t138 + mrSges(5,3) * t136 + Ifges(5,4) * t169 + Ifges(5,2) * t168 + Ifges(5,6) * qJDD(4) - pkin(4) * t217 + pkin(8) * t220 + qJD(4) * t159 + t203 * t109 + t199 * t110 - t157 * t226;
t158 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t204 - Ifges(5,2) * t200) * t195;
t97 = mrSges(5,2) * t138 - mrSges(5,3) * t135 + Ifges(5,1) * t169 + Ifges(5,4) * t168 + Ifges(5,5) * qJDD(4) - pkin(8) * t108 - qJD(4) * t158 - t109 * t199 + t110 * t203 - t157 * t227;
t215 = mrSges(4,2) * t148 - mrSges(4,3) * t145 + Ifges(4,1) * t193 - pkin(7) * t102 - t200 * t95 + t204 * t97;
t214 = mrSges(6,1) * t122 - mrSges(6,2) * t123 + Ifges(6,5) * t133 + Ifges(6,6) * t132 + Ifges(6,3) * t192 + t161 * t143 - t160 * t144;
t213 = -mrSges(4,1) * t145 - pkin(3) * t115 - pkin(7) * t103 - t200 * t97 - t204 * t95;
t212 = -mrSges(3,2) * t152 + pkin(2) * (-mrSges(4,2) * t193 + t216) + qJ(3) * t211 + mrSges(3,1) * t151 + Ifges(3,3) * t193 + t215;
t210 = mrSges(5,1) * t135 - mrSges(5,2) * t136 + Ifges(5,5) * t169 + Ifges(5,6) * t168 + Ifges(5,3) * qJDD(4) + pkin(4) * t108 + t158 * t226 + t159 * t227 + t214;
t209 = mrSges(2,1) * t179 - mrSges(2,2) * t180 + Ifges(2,3) * qJDD(1) + pkin(1) * t94 + t212;
t208 = mrSges(4,1) * t148 + pkin(3) * t102 + t210;
t101 = -m(4) * g(3) + t103;
t92 = m(2) * t180 - mrSges(2,1) * t207 - qJDD(1) * mrSges(2,2) + t221;
t91 = m(2) * t179 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t207 + t94;
t90 = t208 + (t228 * t191) + t229 * t193 + (-mrSges(3,2) + mrSges(4,3)) * g(3) - mrSges(3,3) * t151 - qJ(3) * t101;
t89 = mrSges(3,3) * t152 - pkin(2) * t101 + g(3) * t230 + t191 * t229 - t193 * t228 + t213;
t88 = -mrSges(2,2) * g(3) - mrSges(2,3) * t179 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t207 - pkin(6) * t94 - t201 * t89 + t205 * t90;
t87 = Ifges(2,6) * qJDD(1) + t207 * Ifges(2,5) + mrSges(2,3) * t180 + t201 * t90 + t205 * t89 - pkin(1) * t103 + pkin(6) * t221 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t206 * t88 - t202 * t87 - pkin(5) * (t202 * t92 + t206 * t91), t88, t90, t215, t97, t110; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t202 * t88 + t206 * t87 + pkin(5) * (-t202 * t91 + t206 * t92), t87, t89, -mrSges(4,3) * g(3) + Ifges(4,4) * t193 - (t191 * Ifges(4,5)) - t208, t95, t109; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t209, t209, t212, mrSges(4,2) * g(3) + t191 * Ifges(4,4) + Ifges(4,5) * t193 - t213, t210, t214;];
m_new = t1;
