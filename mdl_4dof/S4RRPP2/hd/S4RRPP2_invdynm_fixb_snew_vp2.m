% Calculate vector of cutting torques with Newton-Euler for
% S4RRPP2
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
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-05-04 19:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:20:14
% EndTime: 2019-05-04 19:20:14
% DurationCPUTime: 0.32s
% Computational Cost: add. (2751->120), mult. (3151->128), div. (0->0), fcn. (1128->4), ass. (0->51)
t104 = qJDD(1) + qJDD(2);
t105 = (qJD(1) + qJD(2));
t109 = sin(qJ(2));
t111 = cos(qJ(2));
t110 = sin(qJ(1));
t112 = cos(qJ(1));
t84 = t110 * g(1) - t112 * g(2);
t81 = qJDD(1) * pkin(1) + t84;
t114 = qJD(1) ^ 2;
t85 = -t112 * g(1) - t110 * g(2);
t82 = -t114 * pkin(1) + t85;
t77 = t109 * t81 + t111 * t82;
t124 = t104 * qJ(3) + (2 * qJD(3) * t105) + t77;
t134 = t105 ^ 2;
t72 = -(pkin(2) * t134) + t124;
t136 = m(4) * t72 + t104 * mrSges(4,3);
t108 = g(3) + qJDD(4);
t133 = -pkin(2) - pkin(3);
t68 = (t133 * t134) + t124;
t129 = -mrSges(5,3) * t68 - (Ifges(5,5) * t134) + Ifges(5,6) * t104;
t126 = m(5) * t68 - (mrSges(5,1) * t134);
t63 = t104 * mrSges(5,2) + t126;
t135 = mrSges(4,2) * t72 - qJ(4) * t63 + (m(5) * pkin(3) + mrSges(5,1)) * t108 + t129;
t132 = m(5) * t108;
t131 = (mrSges(3,1) + mrSges(4,1));
t58 = m(3) * t77 + (-mrSges(3,2) + mrSges(5,2)) * t104 - (t131 * t134) + t126 + t136;
t76 = -t109 * t82 + t111 * t81;
t120 = -qJ(3) * t134 + qJDD(3) - t76;
t71 = t133 * t104 + t120;
t127 = -m(5) * t71 + mrSges(5,2) * t134;
t74 = -t104 * pkin(2) + t120;
t121 = -m(4) * t74 + t104 * mrSges(4,1) + mrSges(4,3) * t134 + t127;
t60 = m(3) * t76 - t134 * mrSges(3,2) + (mrSges(3,1) + mrSges(5,1)) * t104 + t121;
t53 = t109 * t58 + t111 * t60;
t130 = t104 * mrSges(5,1);
t125 = -t109 * t60 + t111 * t58;
t123 = mrSges(5,1) * t71 - mrSges(5,2) * t68 - Ifges(5,3) * t104;
t122 = mrSges(5,2) * t108 - mrSges(5,3) * t71 - Ifges(5,5) * t104 - Ifges(5,6) * t134;
t64 = -t127 - t130;
t118 = -mrSges(4,1) * t74 + mrSges(4,3) * t72 + Ifges(4,2) * t104 - pkin(3) * t64 - t123;
t117 = mrSges(4,2) * t74 + mrSges(4,3) * g(3) + Ifges(4,4) * t104 + (Ifges(4,6) * t134) - qJ(4) * t64 + t122;
t116 = -mrSges(3,2) * t77 + t118 + qJ(3) * (-mrSges(4,1) * t134 + t136 + t63) + pkin(2) * (t121 + t130) + mrSges(3,1) * t76 + Ifges(3,3) * t104;
t115 = mrSges(2,1) * t84 - mrSges(2,2) * t85 + Ifges(2,3) * qJDD(1) + pkin(1) * t53 + t116;
t86 = -m(4) * g(3) - t132;
t55 = -mrSges(3,2) * g(3) - mrSges(3,3) * t76 + Ifges(3,5) * t104 - (Ifges(3,6) * t134) - qJ(3) * t86 + t117;
t54 = mrSges(3,3) * t77 - pkin(2) * t86 + (Ifges(3,6) - Ifges(4,6)) * t104 + (Ifges(4,4) + Ifges(3,5)) * t134 + t131 * g(3) + t135;
t51 = m(2) * t85 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t125;
t50 = m(2) * t84 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) + t53;
t49 = -mrSges(2,2) * g(3) - mrSges(2,3) * t84 + Ifges(2,5) * qJDD(1) - t114 * Ifges(2,6) - pkin(5) * t53 - t109 * t54 + t111 * t55;
t48 = Ifges(2,6) * qJDD(1) + t114 * Ifges(2,5) + mrSges(2,3) * t85 + t109 * t55 + t111 * t54 + pkin(1) * t132 + pkin(5) * t125 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t112 * t49 - t110 * t48 - pkin(4) * (t110 * t51 + t112 * t50), t49, t55, t117, t122; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t110 * t49 + t112 * t48 + pkin(4) * (-t110 * t50 + t112 * t51), t48, t54, t118, -mrSges(5,1) * t108 - t129; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t115, t115, t116, -mrSges(4,1) * g(3) - Ifges(4,4) * t134 + Ifges(4,6) * t104 - t135, t123;];
m_new  = t1;
