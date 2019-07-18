% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_invdynm_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:24:57
% EndTime: 2019-07-18 13:24:59
% DurationCPUTime: 0.86s
% Computational Cost: add. (9054->206), mult. (15935->258), div. (0->0), fcn. (11119->8), ass. (0->72)
t143 = sin(qJ(3));
t160 = qJD(1) * t143;
t147 = cos(qJ(3));
t159 = t147 * qJD(1);
t158 = qJD(1) * qJD(3);
t144 = sin(qJ(1));
t148 = cos(qJ(1));
t134 = -t148 * g(1) - t144 * g(2);
t133 = t144 * g(1) - t148 * g(2);
t117 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t134;
t150 = qJD(1) ^ 2;
t123 = -t150 * qJ(2) + qJDD(2) - t133;
t110 = t147 * g(3) + t143 * t117;
t119 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t143 + Ifges(4,6) * t147) * qJD(1);
t120 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t143 + Ifges(4,2) * t147) * qJD(1);
t130 = t143 * qJDD(1) + t147 * t158;
t131 = t147 * qJDD(1) - t143 * t158;
t142 = sin(qJ(4));
t146 = cos(qJ(4));
t127 = t146 * qJD(3) - t142 * t160;
t128 = t142 * qJD(3) + t146 * t160;
t135 = qJD(4) - t159;
t100 = Ifges(5,5) * t128 + Ifges(5,6) * t127 + Ifges(5,3) * t135;
t101 = Ifges(5,4) * t128 + Ifges(5,2) * t127 + Ifges(5,6) * t135;
t104 = -t128 * qJD(4) + t146 * qJDD(3) - t142 * t130;
t105 = t127 * qJD(4) + t142 * qJDD(3) + t146 * t130;
t126 = qJDD(4) - t131;
t141 = sin(qJ(5));
t145 = cos(qJ(5));
t103 = qJDD(5) - t104;
t108 = t145 * t128 + t141 * t135;
t122 = qJD(5) - t127;
t88 = -t108 * qJD(5) - t141 * t105 + t145 * t126;
t107 = -t141 * t128 + t145 * t135;
t89 = t107 * qJD(5) + t145 * t105 + t141 * t126;
t111 = -t143 * g(3) + t147 * t117;
t97 = t146 * t111 + t142 * t123;
t91 = t141 * t110 + t145 * t97;
t92 = Ifges(6,5) * t108 + Ifges(6,6) * t107 + Ifges(6,3) * t122;
t94 = Ifges(6,1) * t108 + Ifges(6,4) * t107 + Ifges(6,5) * t122;
t96 = t142 * t111 - t146 * t123;
t84 = -mrSges(6,1) * t96 + mrSges(6,3) * t91 + Ifges(6,4) * t89 + Ifges(6,2) * t88 + Ifges(6,6) * t103 - t108 * t92 + t122 * t94;
t90 = t145 * t110 - t141 * t97;
t93 = Ifges(6,4) * t108 + Ifges(6,2) * t107 + Ifges(6,6) * t122;
t85 = mrSges(6,2) * t96 - mrSges(6,3) * t90 + Ifges(6,1) * t89 + Ifges(6,4) * t88 + Ifges(6,5) * t103 + t107 * t92 - t122 * t93;
t80 = mrSges(5,2) * t110 + mrSges(5,3) * t96 + Ifges(5,1) * t105 + Ifges(5,4) * t104 + Ifges(5,5) * t126 + t127 * t100 - t135 * t101 - t141 * t84 + t145 * t85;
t102 = Ifges(5,1) * t128 + Ifges(5,4) * t127 + Ifges(5,5) * t135;
t152 = mrSges(6,1) * t90 - mrSges(6,2) * t91 + Ifges(6,5) * t89 + Ifges(6,6) * t88 + Ifges(6,3) * t103 - t107 * t94 + t108 * t93;
t83 = -mrSges(5,1) * t110 + mrSges(5,3) * t97 + Ifges(5,4) * t105 + Ifges(5,2) * t104 + Ifges(5,6) * t126 - t128 * t100 + t135 * t102 - t152;
t75 = mrSges(4,2) * t123 + mrSges(4,3) * t110 + Ifges(4,1) * t130 + Ifges(4,4) * t131 + Ifges(4,5) * qJDD(3) - qJD(3) * t120 + t119 * t159 - t142 * t83 + t146 * t80;
t121 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t143 + Ifges(4,4) * t147) * qJD(1);
t151 = mrSges(5,1) * t96 + mrSges(5,2) * t97 - Ifges(5,5) * t105 - Ifges(5,6) * t104 - Ifges(5,3) * t126 - t128 * t101 + t127 * t102 - t141 * t85 - t145 * t84;
t78 = -mrSges(4,1) * t123 + mrSges(4,3) * t111 + Ifges(4,4) * t130 + Ifges(4,2) * t131 + Ifges(4,6) * qJDD(3) + qJD(3) * t121 - t119 * t160 + t151;
t157 = -mrSges(3,1) * t123 + mrSges(3,3) * t117 + Ifges(3,2) * qJDD(1) + t143 * t75 + t147 * t78;
t156 = mrSges(3,2) * t123 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t150 * Ifges(3,6) - t143 * t78 + t147 * t75;
t106 = -t127 * mrSges(5,1) + t128 * mrSges(5,2);
t112 = -t135 * mrSges(5,2) + t127 * mrSges(5,3);
t113 = t135 * mrSges(5,1) - t128 * mrSges(5,3);
t129 = (-mrSges(4,1) * t147 + mrSges(4,2) * t143) * qJD(1);
t95 = -t107 * mrSges(6,1) + t108 * mrSges(6,2);
t98 = -t122 * mrSges(6,2) + t107 * mrSges(6,3);
t86 = m(6) * t90 + t103 * mrSges(6,1) - t89 * mrSges(6,3) - t108 * t95 + t122 * t98;
t99 = t122 * mrSges(6,1) - t108 * mrSges(6,3);
t87 = m(6) * t91 - t103 * mrSges(6,2) + t88 * mrSges(6,3) + t107 * t95 - t122 * t99;
t76 = m(4) * t111 + t131 * mrSges(4,3) - qJDD(3) * mrSges(4,2) + t129 * t159 - qJD(3) * (qJD(3) * mrSges(4,1) - mrSges(4,3) * t160) + t146 * (m(5) * t97 - t126 * mrSges(5,2) + t104 * mrSges(5,3) + t127 * t106 - t135 * t113 - t141 * t86 + t145 * t87) - t142 * (t126 * mrSges(5,1) + t88 * mrSges(6,1) - t89 * mrSges(6,2) - t105 * mrSges(5,3) - t128 * t106 + t107 * t98 - t108 * t99 + t135 * t112 + (-m(5) - m(6)) * t96);
t81 = -t130 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t129 * t160 + qJD(3) * (-qJD(3) * mrSges(4,2) + mrSges(4,3) * t159) - t105 * mrSges(5,2) + t104 * mrSges(5,1) - t128 * t113 + t127 * t112 - t141 * t87 - t145 * t86 + (-m(4) - m(5)) * t110;
t155 = -mrSges(2,2) * t134 + mrSges(2,1) * t133 + Ifges(2,3) * qJDD(1) + t157 + qJ(2) * (m(3) * t117 - t150 * mrSges(3,1) + qJDD(1) * mrSges(3,3) - t143 * t81 + t147 * t76);
t154 = mrSges(4,1) * t110 + mrSges(4,2) * t111 - Ifges(4,5) * t130 - Ifges(4,6) * t131 - Ifges(4,3) * qJDD(3) - t120 * t160 + t121 * t159 - t142 * t80 - t146 * t83;
t153 = mrSges(3,2) * t117 + t154;
t72 = t153 + (Ifges(3,4) + Ifges(2,5)) * t150 + (mrSges(2,1) + mrSges(3,1)) * g(3) + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + mrSges(2,3) * t134;
t70 = Ifges(2,5) * qJDD(1) - t150 * Ifges(2,6) - mrSges(2,2) * g(3) - mrSges(2,3) * t133 - qJ(2) * (-m(3) * g(3) + t143 * t76 + t147 * t81) + t156;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - t144 * t72 + t148 * t70, t70, t156, t75, t80, t85; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t144 * t70 + t148 * t72, t72, t157, t78, t83, t84; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t155, t155, -mrSges(3,1) * g(3) - t150 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t153, -t154, -t151, t152;];
m_new  = t1;
