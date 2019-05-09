% Calculate vector of cutting torques with Newton-Euler for
% S4PPRR2
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
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2019-05-04 18:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:50:36
% EndTime: 2019-05-04 18:50:37
% DurationCPUTime: 0.37s
% Computational Cost: add. (4833->91), mult. (5914->111), div. (0->0), fcn. (3888->6), ass. (0->48)
t101 = qJD(3) ^ 2;
t92 = -g(2) + qJDD(1);
t93 = sin(pkin(6));
t94 = cos(pkin(6));
t81 = t93 * g(1) + t94 * t92;
t82 = -t94 * g(1) + t93 * t92;
t96 = sin(qJ(3));
t98 = cos(qJ(3));
t74 = t98 * t81 - t96 * t82;
t71 = qJDD(3) * pkin(3) + t74;
t75 = t96 * t81 + t98 * t82;
t72 = -t101 * pkin(3) + t75;
t95 = sin(qJ(4));
t97 = cos(qJ(4));
t69 = t97 * t71 - t95 * t72;
t89 = qJD(3) + qJD(4);
t87 = t89 ^ 2;
t88 = qJDD(3) + qJDD(4);
t66 = m(5) * t69 + t88 * mrSges(5,1) - t87 * mrSges(5,2);
t70 = t95 * t71 + t97 * t72;
t67 = m(5) * t70 - t87 * mrSges(5,1) - t88 * mrSges(5,2);
t60 = t97 * t66 + t95 * t67;
t57 = m(4) * t74 + qJDD(3) * mrSges(4,1) - t101 * mrSges(4,2) + t60;
t110 = -t95 * t66 + t97 * t67;
t58 = m(4) * t75 - t101 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t110;
t51 = t98 * t57 + t96 * t58;
t91 = -g(3) + qJDD(2);
t113 = (m(4) + m(5)) * t91;
t48 = m(3) * t81 + t51;
t111 = -t96 * t57 + t98 * t58;
t49 = m(3) * t82 + t111;
t112 = -t93 * t48 + t94 * t49;
t109 = m(3) * t91 + t113;
t108 = qJ(1) * m(2) + mrSges(1,2) + mrSges(2,3);
t107 = mrSges(5,1) * t69 - mrSges(5,2) * t70 + Ifges(5,3) * t88;
t61 = -mrSges(5,1) * t91 + mrSges(5,3) * t70 + t87 * Ifges(5,5) + Ifges(5,6) * t88;
t62 = mrSges(5,2) * t91 - mrSges(5,3) * t69 + Ifges(5,5) * t88 - t87 * Ifges(5,6);
t52 = Ifges(4,6) * qJDD(3) + t101 * Ifges(4,5) + mrSges(4,3) * t75 + t95 * t62 + t97 * t61 + pkin(5) * t110 + (-pkin(3) * m(5) - mrSges(4,1)) * t91;
t53 = mrSges(4,2) * t91 - mrSges(4,3) * t74 + Ifges(4,5) * qJDD(3) - t101 * Ifges(4,6) - pkin(5) * t60 - t95 * t61 + t97 * t62;
t41 = -mrSges(3,1) * t91 + mrSges(3,3) * t82 - pkin(2) * t113 + pkin(4) * t111 + t98 * t52 + t96 * t53;
t44 = mrSges(3,2) * t91 - mrSges(3,3) * t81 - pkin(4) * t51 - t96 * t52 + t98 * t53;
t106 = mrSges(2,1) * g(3) + mrSges(2,2) * g(1) - pkin(1) * t109 + qJ(2) * t112 + t94 * t41 + t93 * t44;
t46 = t94 * t48 + t93 * t49;
t105 = mrSges(2,2) * t92 - qJ(2) * t46 - t93 * t41 + t94 * t44;
t104 = mrSges(4,1) * t74 - mrSges(4,2) * t75 + Ifges(4,3) * qJDD(3) + pkin(3) * t60 + t107;
t103 = mrSges(3,1) * t81 - mrSges(3,2) * t82 + pkin(2) * t51 + t104;
t102 = mrSges(2,1) * t92 + pkin(1) * t46 + t103;
t1 = [mrSges(1,3) * g(2) - t108 * g(3) + qJ(1) * t109 + t105, -mrSges(2,3) * g(3) + t105, t44, t53, t62; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t106, -mrSges(2,3) * g(1) - t102, t41, t52, t61; -mrSges(1,1) * g(2) + t108 * g(1) - qJ(1) * t112 + t102, t106, t103, t104, t107;];
m_new  = t1;
