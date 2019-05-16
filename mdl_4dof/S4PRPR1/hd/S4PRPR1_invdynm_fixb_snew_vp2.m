% Calculate vector of cutting torques with Newton-Euler for
% S4PRPR1
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-05-04 18:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:56:58
% EndTime: 2019-05-04 18:56:58
% DurationCPUTime: 0.35s
% Computational Cost: add. (3828->112), mult. (5712->128), div. (0->0), fcn. (2808->6), ass. (0->52)
t122 = -pkin(2) - pkin(3);
t106 = sin(qJ(2));
t108 = cos(qJ(2));
t109 = qJD(2) ^ 2;
t105 = sin(qJ(4));
t107 = cos(qJ(4));
t103 = sin(pkin(6));
t104 = cos(pkin(6));
t87 = t103 * g(1) - t104 * g(2);
t88 = -t104 * g(1) - t103 * g(2);
t81 = t106 * t87 + t108 * t88;
t119 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t81;
t72 = t122 * t109 + t119;
t80 = -t106 * t88 + t108 * t87;
t116 = -t109 * qJ(3) + qJDD(3) - t80;
t75 = t122 * qJDD(2) + t116;
t70 = -t105 * t72 + t107 * t75;
t93 = -qJD(2) + qJD(4);
t91 = t93 ^ 2;
t92 = -qJDD(2) + qJDD(4);
t67 = m(5) * t70 + t92 * mrSges(5,1) - (t91 * mrSges(5,2));
t71 = t105 * t75 + t107 * t72;
t68 = m(5) * t71 - t91 * mrSges(5,1) - t92 * mrSges(5,2);
t62 = -t105 * t67 + t107 * t68;
t76 = -t109 * pkin(2) + t119;
t118 = m(4) * t76 + qJDD(2) * mrSges(4,3) + t62;
t57 = m(3) * t81 - qJDD(2) * mrSges(3,2) + (-mrSges(3,1) - mrSges(4,1)) * t109 + t118;
t61 = t105 * t68 + t107 * t67;
t78 = -qJDD(2) * pkin(2) + t116;
t114 = -m(4) * t78 + qJDD(2) * mrSges(4,1) + t109 * mrSges(4,3) - t61;
t58 = m(3) * t80 + qJDD(2) * mrSges(3,1) - t109 * mrSges(3,2) + t114;
t51 = t106 * t57 + t108 * t58;
t121 = -m(5) * pkin(3) - mrSges(4,1);
t120 = -t106 * t58 + t108 * t57;
t117 = mrSges(5,1) * t70 - mrSges(5,2) * t71 + Ifges(5,3) * t92;
t100 = g(3) - qJDD(1);
t64 = -mrSges(5,1) * t100 + mrSges(5,3) * t71 + (t91 * Ifges(5,5)) + Ifges(5,6) * t92;
t65 = mrSges(5,2) * t100 - mrSges(5,3) * t70 + Ifges(5,5) * t92 - t91 * Ifges(5,6);
t115 = mrSges(4,2) * t78 + Ifges(4,4) * qJDD(2) + t109 * Ifges(4,6) - pkin(5) * t61 - t105 * t64 + t107 * t65;
t113 = -mrSges(4,2) * t76 + pkin(5) * t62 + t105 * t65 + t107 * t64;
t112 = -mrSges(4,1) * t78 + mrSges(4,3) * t76 + Ifges(4,2) * qJDD(2) - pkin(3) * t61 - t117;
t111 = -mrSges(3,2) * t81 + t112 + qJ(3) * (-t109 * mrSges(4,1) + t118) + pkin(2) * t114 + mrSges(3,1) * t80 + Ifges(3,3) * qJDD(2);
t110 = mrSges(2,1) * t87 - mrSges(2,2) * t88 + pkin(1) * t51 + t111;
t89 = m(4) * t100;
t86 = -m(5) * t100 - t89;
t53 = -mrSges(3,3) * t80 + Ifges(3,5) * qJDD(2) - t109 * Ifges(3,6) - qJ(3) * t86 + (-mrSges(3,2) + mrSges(4,3)) * t100 + t115;
t52 = mrSges(3,3) * t81 - pkin(2) * t86 + (Ifges(4,4) + Ifges(3,5)) * t109 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + (mrSges(3,1) - t121) * t100 - t113;
t49 = m(2) * t88 + t120;
t48 = m(2) * t87 + t51;
t47 = -mrSges(2,2) * t100 - mrSges(2,3) * t87 - pkin(4) * t51 - t106 * t52 + t108 * t53;
t46 = mrSges(2,3) * t88 + t106 * t53 + t108 * t52 + pkin(1) * t89 + pkin(4) * t120 + (mrSges(2,1) - pkin(1) * (-m(3) - m(5))) * t100;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t104 * t47 - t103 * t46 - qJ(1) * (t103 * t49 + t104 * t48), t47, t53, mrSges(4,3) * t100 + t115, t65; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t103 * t47 + t104 * t46 + qJ(1) * (-t103 * t48 + t104 * t49), t46, t52, t112, t64; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t110, t110, t111, -t109 * Ifges(4,4) + Ifges(4,6) * qJDD(2) + t121 * t100 + t113, t117;];
m_new  = t1;
