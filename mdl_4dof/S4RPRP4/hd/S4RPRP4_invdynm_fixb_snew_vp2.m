% Calculate vector of cutting torques with Newton-Euler for
% S4RPRP4
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP4_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP4_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:39
% EndTime: 2019-12-31 16:43:40
% DurationCPUTime: 0.87s
% Computational Cost: add. (6602->189), mult. (12249->238), div. (0->0), fcn. (5457->6), ass. (0->69)
t160 = sin(qJ(1));
t162 = cos(qJ(1));
t145 = t160 * g(1) - t162 * g(2);
t131 = qJDD(1) * pkin(1) + t145;
t146 = -t162 * g(1) - t160 * g(2);
t164 = qJD(1) ^ 2;
t135 = -t164 * pkin(1) + t146;
t157 = sin(pkin(6));
t158 = cos(pkin(6));
t109 = t157 * t131 + t158 * t135;
t106 = -t164 * pkin(2) + qJDD(1) * pkin(5) + t109;
t159 = sin(qJ(3));
t156 = -g(3) + qJDD(2);
t161 = cos(qJ(3));
t181 = t161 * t156;
t102 = -t159 * t106 + t181;
t103 = t161 * t106 + t159 * t156;
t116 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t159 - Ifges(5,3) * t161) * qJD(1);
t119 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t159 + Ifges(4,2) * t161) * qJD(1);
t133 = (-mrSges(5,1) * t161 - mrSges(5,3) * t159) * qJD(1);
t176 = qJD(1) * qJD(3);
t136 = t159 * qJDD(1) + t161 * t176;
t137 = t161 * qJDD(1) - t159 * t176;
t132 = (-pkin(3) * t161 - qJ(4) * t159) * qJD(1);
t163 = qJD(3) ^ 2;
t101 = -qJDD(3) * pkin(3) - t163 * qJ(4) - t181 + qJDD(4) + (qJD(1) * t132 + t106) * t159;
t177 = qJD(1) * t161;
t99 = -t163 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t132 * t177 + t103;
t170 = -mrSges(5,1) * t101 + mrSges(5,3) * t99 + Ifges(5,4) * t136 + Ifges(5,2) * qJDD(3) - Ifges(5,6) * t137;
t144 = mrSges(5,2) * t177 + qJD(3) * mrSges(5,3);
t172 = -m(5) * t101 + qJDD(3) * mrSges(5,1) + qJD(3) * t144;
t178 = qJD(1) * t159;
t142 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t178;
t173 = m(5) * t99 + qJDD(3) * mrSges(5,3) + qJD(3) * t142 + t133 * t177;
t120 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t159 - Ifges(5,5) * t161) * qJD(1);
t179 = t120 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t159 + Ifges(4,4) * t161) * qJD(1);
t184 = -qJD(1) * ((t116 - t119) * t159 + t179 * t161) + mrSges(4,1) * t102 - mrSges(4,2) * t103 + Ifges(4,5) * t136 + Ifges(4,6) * t137 + Ifges(4,3) * qJDD(3) + pkin(3) * (-t136 * mrSges(5,2) - t133 * t178 + t172) + qJ(4) * (t137 * mrSges(5,2) + t173) + t170;
t182 = mrSges(4,3) + mrSges(5,2);
t134 = (-mrSges(4,1) * t161 + mrSges(4,2) * t159) * qJD(1);
t141 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t178;
t91 = m(4) * t103 - qJDD(3) * mrSges(4,2) - qJD(3) * t141 + t134 * t177 + t182 * t137 + t173;
t143 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t177;
t92 = m(4) * t102 + qJDD(3) * mrSges(4,1) + qJD(3) * t143 - t182 * t136 + (-t133 - t134) * t178 + t172;
t174 = -t159 * t92 + t161 * t91;
t82 = m(3) * t109 - t164 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t174;
t108 = t158 * t131 - t157 * t135;
t105 = -qJDD(1) * pkin(2) - t164 * pkin(5) - t108;
t96 = -t137 * pkin(3) - t136 * qJ(4) + (-0.2e1 * qJD(4) * t159 + (pkin(3) * t159 - qJ(4) * t161) * qJD(3)) * qJD(1) + t105;
t93 = m(5) * t96 - t137 * mrSges(5,1) - t136 * mrSges(5,3) - t142 * t178 - t144 * t177;
t166 = -m(4) * t105 + t137 * mrSges(4,1) - t136 * mrSges(4,2) - t141 * t178 + t143 * t177 - t93;
t86 = m(3) * t108 + qJDD(1) * mrSges(3,1) - t164 * mrSges(3,2) + t166;
t75 = t157 * t82 + t158 * t86;
t84 = t159 * t91 + t161 * t92;
t175 = -t157 * t86 + t158 * t82;
t171 = -mrSges(5,1) * t96 + mrSges(5,2) * t99;
t117 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t159 + Ifges(4,6) * t161) * qJD(1);
t118 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t159 - Ifges(5,6) * t161) * qJD(1);
t78 = -mrSges(4,1) * t105 + mrSges(4,3) * t103 - pkin(3) * t93 + (Ifges(4,2) + Ifges(5,3)) * t137 + (Ifges(4,4) - Ifges(5,5)) * t136 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + t179 * qJD(3) + (-t117 - t118) * t178 + t171;
t168 = mrSges(5,2) * t101 - mrSges(5,3) * t96 + Ifges(5,1) * t136 + Ifges(5,4) * qJDD(3) - Ifges(5,5) * t137 + qJD(3) * t116 + t118 * t177;
t79 = mrSges(4,2) * t105 - mrSges(4,3) * t102 + Ifges(4,1) * t136 + Ifges(4,4) * t137 + Ifges(4,5) * qJDD(3) - qJ(4) * t93 - qJD(3) * t119 + t117 * t177 + t168;
t169 = mrSges(3,1) * t108 - mrSges(3,2) * t109 + Ifges(3,3) * qJDD(1) + pkin(2) * t166 + pkin(5) * t174 + t159 * t79 + t161 * t78;
t167 = mrSges(2,1) * t145 - mrSges(2,2) * t146 + Ifges(2,3) * qJDD(1) + pkin(1) * t75 + t169;
t73 = m(2) * t146 - t164 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t175;
t72 = m(2) * t145 + qJDD(1) * mrSges(2,1) - t164 * mrSges(2,2) + t75;
t71 = -mrSges(3,1) * t156 + mrSges(3,3) * t109 + t164 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t84 - t184;
t70 = mrSges(3,2) * t156 - mrSges(3,3) * t108 + Ifges(3,5) * qJDD(1) - t164 * Ifges(3,6) - pkin(5) * t84 - t159 * t78 + t161 * t79;
t69 = -mrSges(2,2) * g(3) - mrSges(2,3) * t145 + Ifges(2,5) * qJDD(1) - t164 * Ifges(2,6) - qJ(2) * t75 - t157 * t71 + t158 * t70;
t68 = Ifges(2,6) * qJDD(1) + t164 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t146 + t157 * t70 + t158 * t71 - pkin(1) * (m(3) * t156 + t84) + qJ(2) * t175;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t162 * t69 - t160 * t68 - pkin(4) * (t160 * t73 + t162 * t72), t69, t70, t79, t168; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t160 * t69 + t162 * t68 + pkin(4) * (-t160 * t72 + t162 * t73), t68, t71, t78, (-t159 * t116 - t161 * t120) * qJD(1) + t170; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t167, t167, t169, t184, Ifges(5,5) * t136 + Ifges(5,6) * qJDD(3) - Ifges(5,3) * t137 - qJD(3) * t120 + t118 * t178 - t171;];
m_new = t1;
