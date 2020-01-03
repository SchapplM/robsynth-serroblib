% Calculate vector of cutting torques with Newton-Euler for
% S4RRRR1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:10
% EndTime: 2019-12-31 17:22:11
% DurationCPUTime: 0.98s
% Computational Cost: add. (16532->151), mult. (18278->195), div. (0->0), fcn. (9296->8), ass. (0->68)
t132 = qJD(1) + qJD(2);
t128 = qJD(3) + t132;
t134 = sin(qJ(4));
t138 = cos(qJ(4));
t107 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t134 + Ifges(5,2) * t138) * t128;
t108 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t134 + Ifges(5,4) * t138) * t128;
t131 = qJDD(1) + qJDD(2);
t127 = qJDD(3) + t131;
t152 = qJD(4) * t128;
t112 = t134 * t127 + t138 * t152;
t113 = t138 * t127 - t134 * t152;
t126 = t128 ^ 2;
t137 = sin(qJ(1));
t141 = cos(qJ(1));
t122 = t137 * g(1) - t141 * g(2);
t119 = qJDD(1) * pkin(1) + t122;
t123 = -t141 * g(1) - t137 * g(2);
t142 = qJD(1) ^ 2;
t120 = -t142 * pkin(1) + t123;
t136 = sin(qJ(2));
t140 = cos(qJ(2));
t104 = t140 * t119 - t136 * t120;
t101 = t131 * pkin(2) + t104;
t105 = t136 * t119 + t140 * t120;
t130 = t132 ^ 2;
t102 = -t130 * pkin(2) + t105;
t135 = sin(qJ(3));
t139 = cos(qJ(3));
t98 = t135 * t101 + t139 * t102;
t95 = -t126 * pkin(3) + t127 * pkin(7) + t98;
t92 = -t138 * g(3) - t134 * t95;
t93 = -t134 * g(3) + t138 * t95;
t155 = mrSges(5,1) * t92 - mrSges(5,2) * t93 + Ifges(5,5) * t112 + Ifges(5,6) * t113 + Ifges(5,3) * qJDD(4) + (t107 * t134 - t108 * t138) * t128;
t111 = (-mrSges(5,1) * t138 + mrSges(5,2) * t134) * t128;
t153 = t128 * t138;
t118 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t153;
t154 = t128 * t134;
t90 = m(5) * t92 + qJDD(4) * mrSges(5,1) - t112 * mrSges(5,3) + qJD(4) * t118 - t111 * t154;
t117 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t154;
t91 = m(5) * t93 - qJDD(4) * mrSges(5,2) + t113 * mrSges(5,3) - qJD(4) * t117 + t111 * t153;
t151 = -t134 * t90 + t138 * t91;
t77 = m(4) * t98 - t126 * mrSges(4,1) - t127 * mrSges(4,2) + t151;
t97 = t139 * t101 - t135 * t102;
t94 = -t127 * pkin(3) - t126 * pkin(7) - t97;
t146 = -m(5) * t94 + t113 * mrSges(5,1) - t112 * mrSges(5,2) - t117 * t154 + t118 * t153;
t85 = m(4) * t97 + t127 * mrSges(4,1) - t126 * mrSges(4,2) + t146;
t74 = t135 * t77 + t139 * t85;
t70 = m(3) * t104 + t131 * mrSges(3,1) - t130 * mrSges(3,2) + t74;
t150 = -t135 * t85 + t139 * t77;
t71 = m(3) * t105 - t130 * mrSges(3,1) - t131 * mrSges(3,2) + t150;
t65 = t136 * t71 + t140 * t70;
t79 = t134 * t91 + t138 * t90;
t149 = -t136 * t70 + t140 * t71;
t106 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t134 + Ifges(5,6) * t138) * t128;
t82 = -mrSges(5,1) * t94 + mrSges(5,3) * t93 + Ifges(5,4) * t112 + Ifges(5,2) * t113 + Ifges(5,6) * qJDD(4) + qJD(4) * t108 - t106 * t154;
t83 = mrSges(5,2) * t94 - mrSges(5,3) * t92 + Ifges(5,1) * t112 + Ifges(5,4) * t113 + Ifges(5,5) * qJDD(4) - qJD(4) * t107 + t106 * t153;
t147 = mrSges(4,1) * t97 - mrSges(4,2) * t98 + Ifges(4,3) * t127 + pkin(3) * t146 + pkin(7) * t151 + t134 * t83 + t138 * t82;
t144 = mrSges(3,1) * t104 - mrSges(3,2) * t105 + Ifges(3,3) * t131 + pkin(2) * t74 + t147;
t143 = mrSges(2,1) * t122 - mrSges(2,2) * t123 + Ifges(2,3) * qJDD(1) + pkin(1) * t65 + t144;
t72 = mrSges(4,1) * g(3) + mrSges(4,3) * t98 + t126 * Ifges(4,5) + Ifges(4,6) * t127 - pkin(3) * t79 - t155;
t66 = -mrSges(4,2) * g(3) - mrSges(4,3) * t97 + Ifges(4,5) * t127 - t126 * Ifges(4,6) - pkin(7) * t79 - t134 * t82 + t138 * t83;
t63 = m(2) * t123 - t142 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t149;
t62 = m(2) * t122 + qJDD(1) * mrSges(2,1) - t142 * mrSges(2,2) + t65;
t61 = -mrSges(3,2) * g(3) - mrSges(3,3) * t104 + Ifges(3,5) * t131 - t130 * Ifges(3,6) - pkin(6) * t74 - t135 * t72 + t139 * t66;
t60 = Ifges(3,6) * t131 + t130 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t105 + t135 * t66 + t139 * t72 - pkin(2) * (-m(4) * g(3) + t79) + pkin(6) * t150;
t59 = -mrSges(2,2) * g(3) - mrSges(2,3) * t122 + Ifges(2,5) * qJDD(1) - t142 * Ifges(2,6) - pkin(5) * t65 - t136 * t60 + t140 * t61;
t58 = Ifges(2,6) * qJDD(1) + t142 * Ifges(2,5) + mrSges(2,3) * t123 + t136 * t61 + t140 * t60 - pkin(1) * t79 + pkin(5) * t149 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t141 * t59 - t137 * t58 - pkin(4) * (t137 * t63 + t141 * t62), t59, t61, t66, t83; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t137 * t59 + t141 * t58 + pkin(4) * (-t137 * t62 + t141 * t63), t58, t60, t72, t82; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t143, t143, t144, t147, t155;];
m_new = t1;
