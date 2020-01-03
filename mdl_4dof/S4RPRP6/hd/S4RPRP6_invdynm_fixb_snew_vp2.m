% Calculate vector of cutting torques with Newton-Euler for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP6_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:00
% EndTime: 2019-12-31 16:46:01
% DurationCPUTime: 0.61s
% Computational Cost: add. (3865->188), mult. (7299->226), div. (0->0), fcn. (2818->4), ass. (0->70)
t166 = sin(qJ(1));
t168 = cos(qJ(1));
t149 = -t168 * g(1) - t166 * g(2);
t181 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t149;
t169 = qJD(1) ^ 2;
t199 = (-pkin(1) - pkin(5));
t106 = (t199 * t169) + t181;
t165 = sin(qJ(3));
t167 = cos(qJ(3));
t190 = qJD(1) * qJD(3);
t139 = -t165 * qJDD(1) - t167 * t190;
t186 = t165 * t190;
t140 = t167 * qJDD(1) - t186;
t192 = qJD(1) * t165;
t144 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t192;
t191 = qJD(1) * t167;
t147 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t191;
t143 = -(qJD(3) * mrSges(5,2)) - mrSges(5,3) * t192;
t146 = (qJD(3) * mrSges(5,1)) - mrSges(5,3) * t191;
t145 = (qJD(3) * pkin(3)) - qJ(4) * t191;
t162 = t165 ^ 2;
t99 = t145 * t191 - t139 * pkin(3) + qJDD(4) + (-qJ(4) * t162 + t199) * t169 + t181;
t184 = m(5) * t99 + t140 * mrSges(5,2) + t143 * t192 + t146 * t191;
t200 = -m(4) * t106 - t140 * mrSges(4,2) + (mrSges(4,1) + mrSges(5,1)) * t139 - t144 * t192 - t147 * t191 - t184;
t198 = mrSges(2,1) - mrSges(3,2);
t196 = Ifges(2,5) - Ifges(3,4);
t195 = (-Ifges(2,6) + Ifges(3,5));
t148 = t166 * g(1) - t168 * g(2);
t180 = -t169 * qJ(2) + qJDD(2) - t148;
t107 = t199 * qJDD(1) + t180;
t102 = -t167 * g(3) + t165 * t107;
t188 = -2 * qJD(1) * qJD(4);
t96 = -t162 * t169 * pkin(3) + t139 * qJ(4) - qJD(3) * t145 + t165 * t188 + t102;
t194 = m(5) * t96 + t139 * mrSges(5,3);
t116 = (Ifges(5,3) * qJD(3)) + (Ifges(5,5) * t167 - Ifges(5,6) * t165) * qJD(1);
t193 = -t116 - (Ifges(4,3) * qJD(3)) - (Ifges(4,5) * t167 - Ifges(4,6) * t165) * qJD(1);
t101 = t165 * g(3) + t167 * t107;
t95 = t167 * t188 + (-t140 - t186) * qJ(4) + (-t165 * t167 * t169 + qJDD(3)) * pkin(3) + t101;
t187 = m(5) * t95 + qJDD(3) * mrSges(5,1) + qJD(3) * t143;
t137 = (mrSges(5,1) * t165 + mrSges(5,2) * t167) * qJD(1);
t185 = qJD(1) * (-t137 - (mrSges(4,1) * t165 + mrSges(4,2) * t167) * qJD(1));
t87 = m(4) * t101 + qJDD(3) * mrSges(4,1) + qJD(3) * t144 + (-mrSges(4,3) - mrSges(5,3)) * t140 + t167 * t185 + t187;
t88 = m(4) * t102 + t139 * mrSges(4,3) + (-mrSges(4,2) - mrSges(5,2)) * qJDD(3) + (-t146 - t147) * qJD(3) + t165 * t185 + t194;
t82 = -t165 * t87 + t167 * t88;
t81 = t165 * t88 + t167 * t87;
t182 = mrSges(5,2) * t99 - mrSges(5,3) * t95 + Ifges(5,1) * t140 + Ifges(5,4) * t139 + Ifges(5,5) * qJDD(3);
t120 = Ifges(5,5) * qJD(3) + (Ifges(5,1) * t167 - Ifges(5,4) * t165) * qJD(1);
t179 = -mrSges(5,1) * t99 + mrSges(5,3) * t96 + Ifges(5,4) * t140 + Ifges(5,2) * t139 + Ifges(5,6) * qJDD(3) + qJD(3) * t120;
t115 = -qJDD(1) * pkin(1) + t180;
t177 = -m(3) * t115 + (t169 * mrSges(3,3)) - t81;
t118 = Ifges(5,6) * qJD(3) + (Ifges(5,4) * t167 - Ifges(5,2) * t165) * qJD(1);
t176 = mrSges(5,1) * t95 - mrSges(5,2) * t96 + Ifges(5,5) * t140 + Ifges(5,6) * t139 + Ifges(5,3) * qJDD(3) + t118 * t191 + t120 * t192;
t112 = t169 * pkin(1) - t181;
t121 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t167 - Ifges(4,4) * t165) * qJD(1);
t75 = Ifges(4,4) * t140 + Ifges(4,2) * t139 + Ifges(4,6) * qJDD(3) + qJD(3) * t121 - mrSges(4,1) * t106 + mrSges(4,3) * t102 - pkin(3) * (-t139 * mrSges(5,1) + t184) + qJ(4) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t146 + t194) + (-qJ(4) * t137 * t165 + t193 * t167) * qJD(1) + t179;
t119 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t167 - Ifges(4,2) * t165) * qJD(1);
t90 = -t140 * mrSges(5,3) - t137 * t191 + t187;
t77 = mrSges(4,2) * t106 - mrSges(4,3) * t101 + Ifges(4,1) * t140 + Ifges(4,4) * t139 + Ifges(4,5) * qJDD(3) - qJ(4) * t90 + (-t118 - t119) * qJD(3) + t193 * t192 + t182;
t175 = mrSges(3,2) * t115 - mrSges(3,3) * t112 + Ifges(3,1) * qJDD(1) - pkin(5) * t81 - t165 * t75 + t167 * t77;
t174 = -mrSges(3,1) * t112 - pkin(2) * t200 - pkin(5) * t82 - t165 * t77 - t167 * t75;
t173 = -m(3) * t112 + t169 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t200;
t172 = -mrSges(2,2) * t149 + mrSges(2,1) * t148 + Ifges(2,3) * qJDD(1) + t175 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t177) + qJ(2) * t173;
t171 = mrSges(4,1) * t101 - mrSges(4,2) * t102 + Ifges(4,5) * t140 + Ifges(4,6) * t139 + Ifges(4,3) * qJDD(3) + pkin(3) * t90 + t119 * t191 + t121 * t192 + t176;
t170 = mrSges(3,1) * t115 + pkin(2) * t81 + t171;
t83 = m(2) * t149 - t169 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t173;
t80 = -m(3) * g(3) + t82;
t78 = m(2) * t148 - t169 * mrSges(2,2) + t198 * qJDD(1) + t177;
t74 = t170 + (t195 * t169) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t196 * qJDD(1) - qJ(2) * t80 - mrSges(2,3) * t148;
t73 = mrSges(2,3) * t149 - pkin(1) * t80 + t198 * g(3) - t195 * qJDD(1) + t196 * t169 + t174;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t168 * t74 - t166 * t73 - pkin(4) * (t166 * t83 + t168 * t78), t74, t175, t77, -qJD(3) * t118 - t116 * t192 + t182; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t166 * t74 + t168 * t73 + pkin(4) * (-t166 * t78 + t168 * t83), t73, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t169 * Ifges(3,5)) - t170, t75, -t116 * t191 + t179; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t172, t172, mrSges(3,2) * g(3) + t169 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t174, t171, t176;];
m_new = t1;
