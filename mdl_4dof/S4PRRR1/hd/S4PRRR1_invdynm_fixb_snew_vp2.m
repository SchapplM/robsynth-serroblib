% Calculate vector of cutting torques with Newton-Euler for
% S4PRRR1
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
% Datum: 2019-05-04 19:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:04:23
% EndTime: 2019-05-04 19:04:24
% DurationCPUTime: 0.76s
% Computational Cost: add. (11305->109), mult. (16702->136), div. (0->0), fcn. (11538->8), ass. (0->56)
t117 = sin(qJ(2));
t120 = cos(qJ(2));
t121 = qJD(2) ^ 2;
t116 = sin(qJ(3));
t119 = cos(qJ(3));
t110 = qJD(2) + qJD(3);
t108 = t110 ^ 2;
t109 = qJDD(2) + qJDD(3);
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t105 = qJD(4) + t110;
t103 = t105 ^ 2;
t104 = qJDD(4) + t109;
t113 = sin(pkin(7));
t114 = cos(pkin(7));
t100 = -t114 * g(1) - t113 * g(2);
t99 = t113 * g(1) - t114 * g(2);
t93 = -t117 * t100 + t120 * t99;
t90 = qJDD(2) * pkin(2) + t93;
t94 = t120 * t100 + t117 * t99;
t91 = -t121 * pkin(2) + t94;
t85 = -t116 * t91 + t119 * t90;
t82 = t109 * pkin(3) + t85;
t86 = t116 * t90 + t119 * t91;
t83 = -t108 * pkin(3) + t86;
t80 = -t115 * t83 + t118 * t82;
t77 = m(5) * t80 + t104 * mrSges(5,1) - t103 * mrSges(5,2);
t81 = t115 * t82 + t118 * t83;
t78 = m(5) * t81 - t103 * mrSges(5,1) - t104 * mrSges(5,2);
t71 = t115 * t78 + t118 * t77;
t68 = m(4) * t85 + t109 * mrSges(4,1) - t108 * mrSges(4,2) + t71;
t128 = -t115 * t77 + t118 * t78;
t69 = m(4) * t86 - t108 * mrSges(4,1) - t109 * mrSges(4,2) + t128;
t62 = t116 * t69 + t119 * t68;
t59 = m(3) * t93 + qJDD(2) * mrSges(3,1) - t121 * mrSges(3,2) + t62;
t127 = -t116 * t68 + t119 * t69;
t60 = m(3) * t94 - t121 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t127;
t55 = t117 * t60 + t120 * t59;
t112 = -g(3) + qJDD(1);
t129 = (m(4) + m(5)) * t112;
t126 = -t117 * t59 + t120 * t60;
t125 = mrSges(5,1) * t80 - mrSges(5,2) * t81 + Ifges(5,3) * t104;
t124 = mrSges(4,1) * t85 - mrSges(4,2) * t86 + Ifges(4,3) * t109 + pkin(3) * t71 + t125;
t123 = mrSges(3,1) * t93 - mrSges(3,2) * t94 + Ifges(3,3) * qJDD(2) + pkin(2) * t62 + t124;
t122 = mrSges(2,1) * t99 - mrSges(2,2) * t100 + pkin(1) * t55 + t123;
t73 = mrSges(5,2) * t112 - mrSges(5,3) * t80 + Ifges(5,5) * t104 - t103 * Ifges(5,6);
t72 = -mrSges(5,1) * t112 + mrSges(5,3) * t81 + t103 * Ifges(5,5) + Ifges(5,6) * t104;
t64 = mrSges(4,2) * t112 - mrSges(4,3) * t85 + Ifges(4,5) * t109 - t108 * Ifges(4,6) - pkin(6) * t71 - t115 * t72 + t118 * t73;
t63 = Ifges(4,6) * t109 + t108 * Ifges(4,5) + mrSges(4,3) * t86 + t115 * t73 + t118 * t72 + pkin(6) * t128 + (-pkin(3) * m(5) - mrSges(4,1)) * t112;
t53 = m(2) * t100 + t126;
t52 = m(2) * t99 + t55;
t51 = mrSges(3,2) * t112 - mrSges(3,3) * t93 + Ifges(3,5) * qJDD(2) - t121 * Ifges(3,6) - pkin(5) * t62 - t116 * t63 + t119 * t64;
t50 = -mrSges(3,1) * t112 + mrSges(3,3) * t94 + t121 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t129 + pkin(5) * t127 + t116 * t64 + t119 * t63;
t49 = mrSges(2,2) * t112 - mrSges(2,3) * t99 - pkin(4) * t55 - t117 * t50 + t120 * t51;
t48 = -mrSges(2,1) * t112 + mrSges(2,3) * t100 + t117 * t51 + t120 * t50 - pkin(1) * (m(3) * t112 + t129) + pkin(4) * t126;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t114 * t49 - t113 * t48 - qJ(1) * (t113 * t53 + t114 * t52), t49, t51, t64, t73; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t113 * t49 + t114 * t48 + qJ(1) * (-t113 * t52 + t114 * t53), t48, t50, t63, t72; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t122, t122, t123, t124, t125;];
m_new  = t1;
