% Calculate vector of cutting torques with Newton-Euler for
% S4PPRR5
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PPRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR5_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR5_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:42
% EndTime: 2019-12-31 16:19:43
% DurationCPUTime: 0.43s
% Computational Cost: add. (3630->124), mult. (5943->159), div. (0->0), fcn. (3116->6), ass. (0->56)
t104 = sin(qJ(4));
t106 = cos(qJ(4));
t108 = qJD(3) ^ 2;
t101 = -g(3) + qJDD(1);
t105 = sin(qJ(3));
t107 = cos(qJ(3));
t102 = sin(pkin(6));
t103 = cos(pkin(6));
t95 = t102 * g(1) - t103 * g(2);
t93 = qJDD(2) - t95;
t81 = t107 * t101 + t105 * t93;
t78 = -t108 * pkin(3) + qJDD(3) * pkin(5) + t81;
t96 = t103 * g(1) + t102 * g(2);
t75 = -t104 * t78 - t106 * t96;
t76 = -t104 * t96 + t106 * t78;
t83 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t104 + Ifges(5,2) * t106) * qJD(3);
t84 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t104 + Ifges(5,4) * t106) * qJD(3);
t119 = qJD(3) * qJD(4);
t90 = t104 * qJDD(3) + t106 * t119;
t91 = t106 * qJDD(3) - t104 * t119;
t123 = mrSges(5,1) * t75 - mrSges(5,2) * t76 + Ifges(5,5) * t90 + Ifges(5,6) * t91 + Ifges(5,3) * qJDD(4) + (t104 * t83 - t106 * t84) * qJD(3);
t122 = mrSges(2,2) - mrSges(3,3);
t121 = qJD(3) * t104;
t89 = (-mrSges(5,1) * t106 + mrSges(5,2) * t104) * qJD(3);
t120 = qJD(3) * t106;
t98 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t120;
t71 = m(5) * t75 + qJDD(4) * mrSges(5,1) - t90 * mrSges(5,3) + qJD(4) * t98 - t89 * t121;
t97 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t121;
t72 = m(5) * t76 - qJDD(4) * mrSges(5,2) + t91 * mrSges(5,3) - qJD(4) * t97 + t89 * t120;
t61 = t104 * t72 + t106 * t71;
t59 = m(4) * t96 - t61;
t118 = -t104 * t71 + t106 * t72;
t57 = m(4) * t81 - t108 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t118;
t80 = -t105 * t101 + t107 * t93;
t77 = -qJDD(3) * pkin(3) - t108 * pkin(5) - t80;
t112 = -m(5) * t77 + t91 * mrSges(5,1) - t90 * mrSges(5,2) + t98 * t120 - t97 * t121;
t67 = m(4) * t80 + qJDD(3) * mrSges(4,1) - t108 * mrSges(4,2) + t112;
t54 = -t105 * t67 + t107 * t57;
t53 = t105 * t57 + t107 * t67;
t82 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t104 + Ifges(5,6) * t106) * qJD(3);
t65 = -mrSges(5,1) * t77 + mrSges(5,3) * t76 + Ifges(5,4) * t90 + Ifges(5,2) * t91 + Ifges(5,6) * qJDD(4) + qJD(4) * t84 - t82 * t121;
t66 = mrSges(5,2) * t77 - mrSges(5,3) * t75 + Ifges(5,1) * t90 + Ifges(5,4) * t91 + Ifges(5,5) * qJDD(4) - qJD(4) * t83 + t82 * t120;
t48 = -mrSges(4,2) * t96 - mrSges(4,3) * t80 + Ifges(4,5) * qJDD(3) - t108 * Ifges(4,6) - pkin(5) * t61 - t104 * t65 + t106 * t66;
t49 = mrSges(4,1) * t96 + mrSges(4,3) * t81 + t108 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t61 - t123;
t116 = mrSges(3,2) * t93 - pkin(4) * t53 - t105 * t49 + t107 * t48;
t115 = -m(3) * t93 - t53;
t114 = mrSges(4,1) * t80 - mrSges(4,2) * t81 + Ifges(4,3) * qJDD(3) + pkin(3) * t112 + pkin(5) * t118 + t104 * t66 + t106 * t65;
t113 = -pkin(2) * t59 - pkin(4) * t54 - t105 * t48 - t107 * t49;
t110 = mrSges(3,1) * t93 + pkin(2) * t53 + t114;
t109 = pkin(1) * t115 + t122 * t96 + t116 + qJ(2) * (-m(3) * t96 - t59) + mrSges(2,1) * t95;
t55 = (-m(2) - m(3)) * t96 - t59;
t52 = m(3) * t101 + t54;
t50 = m(2) * t95 + t115;
t46 = -mrSges(2,3) * t95 - qJ(2) * t52 + t122 * t101 + t110;
t45 = -pkin(1) * t52 + (-mrSges(3,1) - mrSges(2,3)) * t96 + (-mrSges(2,1) + mrSges(3,2)) * t101 + t113;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t103 * t46 - t102 * t45 - qJ(1) * (t102 * t55 + t103 * t50), t46, -mrSges(3,3) * t96 + t116, t48, t66; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t102 * t46 + t103 * t45 + qJ(1) * (-t102 * t50 + t103 * t55), t45, mrSges(3,3) * t101 - t110, t49, t65; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t109, t109, mrSges(3,1) * t96 - mrSges(3,2) * t101 - t113, t114, t123;];
m_new = t1;
