% Calculate vector of cutting torques with Newton-Euler for
% S4PRPP1
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
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2019-05-04 18:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:53:21
% EndTime: 2019-05-04 18:53:21
% DurationCPUTime: 0.26s
% Computational Cost: add. (1789->110), mult. (2673->116), div. (0->0), fcn. (1128->4), ass. (0->46)
t96 = qJD(2) ^ 2;
t116 = t96 * mrSges(5,3);
t115 = mrSges(5,2) + mrSges(4,3);
t114 = Ifges(4,4) - Ifges(5,5);
t113 = Ifges(4,5) - Ifges(3,6);
t112 = -pkin(2) - qJ(4);
t92 = sin(pkin(5));
t93 = cos(pkin(5));
t74 = t92 * g(1) - t93 * g(2);
t75 = -t93 * g(1) - t92 * g(2);
t94 = sin(qJ(2));
t95 = cos(qJ(2));
t68 = t94 * t74 + t95 * t75;
t101 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t68;
t62 = t112 * t96 + qJDD(4) + t101;
t111 = m(5) * t62 + qJDD(2) * mrSges(5,2);
t63 = t96 * pkin(2) - t101;
t105 = -m(4) * t63 + t96 * mrSges(4,2) + qJDD(2) * mrSges(4,3) + t111;
t52 = m(3) * t68 - qJDD(2) * mrSges(3,2) + (-mrSges(3,1) - mrSges(5,3)) * t96 + t105;
t67 = t95 * t74 - t94 * t75;
t104 = -t96 * qJ(3) + qJDD(3) - t67;
t61 = -(2 * qJD(4) * qJD(2)) + t112 * qJDD(2) + t104;
t56 = m(5) * t61 - t96 * mrSges(5,2) - qJDD(2) * mrSges(5,3);
t65 = -qJDD(2) * pkin(2) + t104;
t103 = -m(4) * t65 + t96 * mrSges(4,3) - t56;
t53 = m(3) * t67 - t96 * mrSges(3,2) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + t103;
t46 = t94 * t52 + t95 * t53;
t87 = -g(3) + qJDD(1);
t73 = (m(4) + m(5)) * t87;
t110 = mrSges(5,1) * t61 + t96 * Ifges(5,4) + Ifges(5,5) * qJDD(2);
t109 = m(5) * qJ(4) - mrSges(4,2);
t108 = t95 * t52 - t94 * t53;
t107 = mrSges(5,1) * t62 - mrSges(5,3) * t87 - Ifges(5,4) * qJDD(2);
t106 = mrSges(5,2) * t62 - mrSges(5,3) * t61 + Ifges(5,1) * qJDD(2);
t102 = -mrSges(4,1) * t65 - pkin(3) * t56 - t110;
t100 = mrSges(4,1) * t63 + pkin(3) * (-t111 + t116) - t107;
t99 = mrSges(4,2) * t65 - mrSges(4,3) * t63 + Ifges(4,1) * qJDD(2) - qJ(4) * t56 + t106;
t98 = -mrSges(3,2) * t68 + qJ(3) * (t105 - t116) + pkin(2) * (-qJDD(2) * mrSges(4,2) + t103) + mrSges(3,1) * t67 + Ifges(3,3) * qJDD(2) + t99;
t97 = mrSges(2,1) * t74 - mrSges(2,2) * t75 + pkin(1) * t46 + t98;
t48 = -mrSges(3,3) * t67 - qJ(3) * t73 + t113 * t96 + (-Ifges(4,4) + Ifges(3,5)) * qJDD(2) + (mrSges(3,2) - t115) * t87 - t102;
t47 = mrSges(3,3) * t68 - pkin(2) * t73 - t113 * qJDD(2) + (Ifges(3,5) - t114) * t96 + (-mrSges(3,1) - t109) * t87 - t100;
t44 = m(2) * t75 + t108;
t43 = m(2) * t74 + t46;
t42 = mrSges(2,2) * t87 - mrSges(2,3) * t74 - pkin(4) * t46 - t94 * t47 + t95 * t48;
t41 = -mrSges(2,1) * t87 + mrSges(2,3) * t75 + t94 * t48 + t95 * t47 - pkin(1) * (m(3) * t87 + t73) + pkin(4) * t108;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t93 * t42 - t92 * t41 - qJ(1) * (t93 * t43 + t92 * t44), t42, t48, t99, t106; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t92 * t42 + t93 * t41 + qJ(1) * (-t92 * t43 + t93 * t44), t41, t47, Ifges(4,4) * qJDD(2) - t96 * Ifges(4,5) + t115 * t87 + t102, -t96 * Ifges(5,5) - t107; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t97, t97, t98, Ifges(4,5) * qJDD(2) + t109 * t87 + t114 * t96 + t100, -mrSges(5,2) * t87 + t110;];
m_new  = t1;
