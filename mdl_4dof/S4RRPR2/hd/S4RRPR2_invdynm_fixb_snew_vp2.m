% Calculate vector of cutting torques with Newton-Euler for
% S4RRPR2
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
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:32
% EndTime: 2019-07-18 18:16:33
% DurationCPUTime: 0.39s
% Computational Cost: add. (5042->120), mult. (5734->136), div. (0->0), fcn. (2328->6), ass. (0->55)
t129 = -m(4) - m(5);
t128 = -pkin(2) - pkin(3);
t109 = sin(qJ(2));
t112 = cos(qJ(2));
t106 = (qJD(1) + qJD(2));
t104 = t106 ^ 2;
t105 = qJDD(1) + qJDD(2);
t108 = sin(qJ(4));
t111 = cos(qJ(4));
t110 = sin(qJ(1));
t113 = cos(qJ(1));
t88 = t110 * g(1) - t113 * g(2);
t85 = qJDD(1) * pkin(1) + t88;
t115 = qJD(1) ^ 2;
t89 = -t113 * g(1) - t110 * g(2);
t86 = -t115 * pkin(1) + t89;
t81 = t109 * t85 + t112 * t86;
t125 = t105 * qJ(3) + (2 * qJD(3) * t106) + t81;
t72 = t128 * t104 + t125;
t80 = -t109 * t86 + t112 * t85;
t120 = -t104 * qJ(3) + qJDD(3) - t80;
t75 = t128 * t105 + t120;
t70 = -t108 * t72 + t111 * t75;
t102 = qJD(4) - t106;
t98 = t102 ^ 2;
t99 = qJDD(4) - t105;
t67 = m(5) * t70 + t99 * mrSges(5,1) - (t98 * mrSges(5,2));
t71 = t108 * t75 + t111 * t72;
t68 = m(5) * t71 - t98 * mrSges(5,1) - t99 * mrSges(5,2);
t76 = -t104 * pkin(2) + t125;
t124 = m(4) * t76 + t105 * mrSges(4,3) - t108 * t67 + t111 * t68;
t57 = m(3) * t81 - t105 * mrSges(3,2) + (-mrSges(3,1) - mrSges(4,1)) * t104 + t124;
t62 = t108 * t68 + t111 * t67;
t78 = -t105 * pkin(2) + t120;
t119 = -m(4) * t78 + t105 * mrSges(4,1) + t104 * mrSges(4,3) - t62;
t59 = m(3) * t80 + t105 * mrSges(3,1) - t104 * mrSges(3,2) + t119;
t52 = t109 * t57 + t112 * t59;
t127 = -m(5) * pkin(3) - mrSges(4,1);
t126 = -t109 * t59 + t112 * t57;
t123 = mrSges(5,1) * t70 - mrSges(5,2) * t71 + Ifges(5,3) * t99;
t65 = -mrSges(5,1) * g(3) + mrSges(5,3) * t71 + (t98 * Ifges(5,5)) + Ifges(5,6) * t99;
t66 = mrSges(5,2) * g(3) - mrSges(5,3) * t70 + Ifges(5,5) * t99 - t98 * Ifges(5,6);
t122 = mrSges(4,2) * t78 + mrSges(4,3) * g(3) + Ifges(4,4) * t105 + t104 * Ifges(4,6) - t108 * t65 + t111 * t66;
t121 = -mrSges(4,2) * t76 + t108 * t66 + t111 * t65;
t118 = -mrSges(4,1) * t78 + mrSges(4,3) * t76 + Ifges(4,2) * t105 - pkin(3) * t62 - t123;
t117 = -mrSges(3,2) * t81 + t118 + qJ(3) * (-t104 * mrSges(4,1) + t124) + pkin(2) * t119 + mrSges(3,1) * t80 + Ifges(3,3) * t105;
t116 = mrSges(2,1) * t88 - mrSges(2,2) * t89 + Ifges(2,3) * qJDD(1) + pkin(1) * t52 + t117;
t92 = t129 * g(3);
t54 = -mrSges(3,2) * g(3) - mrSges(3,3) * t80 + Ifges(3,5) * t105 - t104 * Ifges(3,6) - qJ(3) * t92 + t122;
t53 = mrSges(3,3) * t81 - pkin(2) * t92 + (Ifges(3,6) - Ifges(4,6)) * t105 + (Ifges(4,4) + Ifges(3,5)) * t104 + (mrSges(3,1) - t127) * g(3) - t121;
t50 = m(2) * t89 - t115 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t126;
t49 = m(2) * t88 + qJDD(1) * mrSges(2,1) - t115 * mrSges(2,2) + t52;
t48 = -mrSges(2,2) * g(3) - mrSges(2,3) * t88 + Ifges(2,5) * qJDD(1) - t115 * Ifges(2,6) - pkin(5) * t52 - t109 * t53 + t112 * t54;
t47 = Ifges(2,6) * qJDD(1) + t115 * Ifges(2,5) + mrSges(2,3) * t89 + t109 * t54 + t112 * t53 + pkin(5) * t126 + (mrSges(2,1) - pkin(1) * (-m(3) + t129)) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t113 * t48 - t110 * t47 - pkin(4) * (t110 * t50 + t113 * t49), t48, t54, t122, t66; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t110 * t48 + t113 * t47 + pkin(4) * (-t110 * t49 + t113 * t50), t47, t53, t118, t65; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t116, t116, t117, -t104 * Ifges(4,4) + Ifges(4,6) * t105 + t127 * g(3) + t121, t123;];
m_new  = t1;
