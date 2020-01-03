% Calculate vector of cutting torques with Newton-Euler for
% S4PRRR5
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR5_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:37
% EndTime: 2019-12-31 16:33:38
% DurationCPUTime: 0.78s
% Computational Cost: add. (10183->140), mult. (12884->183), div. (0->0), fcn. (7300->8), ass. (0->65)
t126 = qJD(2) + qJD(3);
t129 = sin(qJ(4));
t132 = cos(qJ(4));
t107 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t129 + Ifges(5,2) * t132) * t126;
t108 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t129 + Ifges(5,4) * t132) * t126;
t125 = qJDD(2) + qJDD(3);
t145 = qJD(4) * t126;
t112 = t129 * t125 + t132 * t145;
t113 = t132 * t125 - t129 * t145;
t128 = sin(pkin(7));
t148 = cos(pkin(7));
t118 = t128 * g(1) - t148 * g(2);
t124 = t126 ^ 2;
t119 = -t148 * g(1) - t128 * g(2);
t127 = -g(3) + qJDD(1);
t131 = sin(qJ(2));
t134 = cos(qJ(2));
t104 = -t131 * t119 + t134 * t127;
t102 = qJDD(2) * pkin(2) + t104;
t105 = t134 * t119 + t131 * t127;
t135 = qJD(2) ^ 2;
t103 = -t135 * pkin(2) + t105;
t130 = sin(qJ(3));
t133 = cos(qJ(3));
t99 = t130 * t102 + t133 * t103;
t96 = -t124 * pkin(3) + t125 * pkin(6) + t99;
t93 = -t132 * t118 - t129 * t96;
t94 = -t129 * t118 + t132 * t96;
t150 = mrSges(5,1) * t93 - mrSges(5,2) * t94 + Ifges(5,5) * t112 + Ifges(5,6) * t113 + Ifges(5,3) * qJDD(4) + (t107 * t129 - t108 * t132) * t126;
t149 = m(3) + m(4);
t111 = (-mrSges(5,1) * t132 + mrSges(5,2) * t129) * t126;
t146 = t126 * t132;
t116 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t146;
t147 = t126 * t129;
t89 = m(5) * t93 + qJDD(4) * mrSges(5,1) - t112 * mrSges(5,3) + qJD(4) * t116 - t111 * t147;
t115 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t147;
t90 = m(5) * t94 - qJDD(4) * mrSges(5,2) + t113 * mrSges(5,3) - qJD(4) * t115 + t111 * t146;
t144 = -t129 * t89 + t132 * t90;
t74 = m(4) * t99 - t124 * mrSges(4,1) - t125 * mrSges(4,2) + t144;
t98 = t133 * t102 - t130 * t103;
t95 = -t125 * pkin(3) - t124 * pkin(6) - t98;
t138 = -m(5) * t95 + t113 * mrSges(5,1) - t112 * mrSges(5,2) - t115 * t147 + t116 * t146;
t85 = m(4) * t98 + t125 * mrSges(4,1) - t124 * mrSges(4,2) + t138;
t71 = t130 * t74 + t133 * t85;
t78 = t129 * t90 + t132 * t89;
t143 = -t130 * t85 + t133 * t74;
t69 = m(3) * t104 + qJDD(2) * mrSges(3,1) - t135 * mrSges(3,2) + t71;
t70 = m(3) * t105 - t135 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t143;
t142 = -t131 * t69 + t134 * t70;
t106 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t129 + Ifges(5,6) * t132) * t126;
t82 = -mrSges(5,1) * t95 + mrSges(5,3) * t94 + Ifges(5,4) * t112 + Ifges(5,2) * t113 + Ifges(5,6) * qJDD(4) + qJD(4) * t108 - t106 * t147;
t83 = mrSges(5,2) * t95 - mrSges(5,3) * t93 + Ifges(5,1) * t112 + Ifges(5,4) * t113 + Ifges(5,5) * qJDD(4) - qJD(4) * t107 + t106 * t146;
t66 = -mrSges(4,2) * t118 - mrSges(4,3) * t98 + Ifges(4,5) * t125 - t124 * Ifges(4,6) - pkin(6) * t78 - t129 * t82 + t132 * t83;
t67 = mrSges(4,1) * t118 + mrSges(4,3) * t99 + t124 * Ifges(4,5) + Ifges(4,6) * t125 - pkin(3) * t78 - t150;
t60 = Ifges(3,6) * qJDD(2) + t135 * Ifges(3,5) + mrSges(3,1) * t118 + mrSges(3,3) * t105 + t130 * t66 + t133 * t67 - pkin(2) * (-m(4) * t118 + t78) + pkin(5) * t143;
t62 = -mrSges(3,2) * t118 - mrSges(3,3) * t104 + Ifges(3,5) * qJDD(2) - t135 * Ifges(3,6) - pkin(5) * t71 - t130 * t67 + t133 * t66;
t140 = -mrSges(2,2) * t119 + pkin(4) * t142 + t131 * t62 + t134 * t60 + mrSges(2,1) * t118 + pkin(1) * (t149 * t118 - t78);
t139 = mrSges(4,1) * t98 - mrSges(4,2) * t99 + Ifges(4,3) * t125 + pkin(3) * t138 + pkin(6) * t144 + t129 * t83 + t132 * t82;
t136 = mrSges(3,1) * t104 - mrSges(3,2) * t105 + Ifges(3,3) * qJDD(2) + pkin(2) * t71 + t139;
t75 = (m(2) + t149) * t118 - t78;
t65 = t131 * t70 + t134 * t69;
t63 = m(2) * t119 + t142;
t58 = -mrSges(2,1) * t127 + mrSges(2,3) * t119 - pkin(1) * t65 - t136;
t57 = mrSges(2,2) * t127 - mrSges(2,3) * t118 - pkin(4) * t65 - t131 * t60 + t134 * t62;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t148 * t57 - t128 * t58 - qJ(1) * (t128 * t63 + t148 * t75), t57, t62, t66, t83; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t128 * t57 + t148 * t58 + qJ(1) * (-t128 * t75 + t148 * t63), t58, t60, t67, t82; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t140, t140, t136, t139, t150;];
m_new = t1;
