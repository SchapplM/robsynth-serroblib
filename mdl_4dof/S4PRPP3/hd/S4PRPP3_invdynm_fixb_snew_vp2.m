% Calculate vector of cutting torques with Newton-Euler for
% S4PRPP3
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
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
% Datum: 2019-05-04 18:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRPP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_invdynm_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:55:42
% EndTime: 2019-05-04 18:55:42
% DurationCPUTime: 0.20s
% Computational Cost: add. (1087->106), mult. (1365->104), div. (0->0), fcn. (282->2), ass. (0->42)
t94 = qJD(2) ^ 2;
t121 = -Ifges(5,5) * qJDD(2) - t94 * Ifges(5,6);
t84 = -g(2) + qJDD(1);
t89 = sin(qJ(2));
t90 = cos(qJ(2));
t68 = -t90 * g(1) + t89 * t84;
t104 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t68;
t118 = -pkin(2) - pkin(3);
t57 = t118 * t94 + t104;
t111 = m(5) * t57 + qJDD(2) * mrSges(5,2);
t112 = mrSges(5,3) * t57 + t94 * Ifges(5,5);
t62 = -t94 * pkin(2) + t104;
t83 = g(3) + qJDD(4);
t120 = (m(5) * pkin(3) + mrSges(5,1)) * t83 + mrSges(4,2) * t62 - qJ(4) * (-t94 * mrSges(5,1) + t111) - t112;
t119 = m(3) + m(4);
t117 = m(5) * t83;
t115 = -mrSges(4,1) - mrSges(5,1);
t114 = -mrSges(1,2) - mrSges(2,3);
t113 = Ifges(4,6) - Ifges(5,6);
t67 = t89 * g(1) + t90 * t84;
t99 = -t94 * qJ(3) + qJDD(3) - t67;
t61 = t118 * qJDD(2) + t99;
t109 = mrSges(5,1) * t61 - mrSges(5,2) * t57 - Ifges(5,3) * qJDD(2);
t107 = mrSges(5,2) * t83 - mrSges(5,3) * t61;
t105 = m(4) * t62 + qJDD(2) * mrSges(4,3) + t111;
t46 = m(3) * t68 - qJDD(2) * mrSges(3,2) + (-mrSges(3,1) + t115) * t94 + t105;
t53 = m(5) * t61 - qJDD(2) * mrSges(5,1) - t94 * mrSges(5,2);
t64 = -qJDD(2) * pkin(2) + t99;
t100 = -m(4) * t64 + qJDD(2) * mrSges(4,1) + t94 * mrSges(4,3) - t53;
t47 = m(3) * t67 + qJDD(2) * mrSges(3,1) - t94 * mrSges(3,2) + t100;
t106 = t90 * t46 - t89 * t47;
t69 = -m(4) * g(3) - t117;
t41 = mrSges(3,3) * t68 - pkin(2) * t69 + (Ifges(4,4) + Ifges(3,5)) * t94 + (mrSges(3,1) + mrSges(4,1)) * g(3) + (Ifges(3,6) - t113) * qJDD(2) + t120;
t97 = mrSges(4,2) * t64 + mrSges(4,3) * g(3) + Ifges(4,4) * qJDD(2) + t94 * Ifges(4,6) - qJ(4) * t53 + t107;
t44 = -mrSges(3,2) * g(3) - mrSges(3,3) * t67 - qJ(3) * t69 + (-Ifges(3,6) - Ifges(5,6)) * t94 + (Ifges(3,5) - Ifges(5,5)) * qJDD(2) + t97;
t103 = pkin(4) * t106 + t90 * t41 + t89 * t44 + pkin(1) * (t119 * g(3) + t117) + mrSges(2,2) * g(1) + mrSges(2,1) * g(3);
t39 = t89 * t46 + t90 * t47;
t102 = mrSges(2,2) * t84 - pkin(4) * t39 - t89 * t41 + t90 * t44;
t98 = -mrSges(4,1) * t64 + mrSges(4,3) * t62 + Ifges(4,2) * qJDD(2) - pkin(3) * t53 - t109;
t96 = -mrSges(3,2) * t68 + qJ(3) * (t115 * t94 + t105) + pkin(2) * t100 + mrSges(3,1) * t67 + Ifges(3,3) * qJDD(2) + t98;
t95 = mrSges(2,1) * t84 + pkin(1) * t39 + t96;
t1 = [-qJ(1) * t117 + mrSges(1,3) * g(2) + (qJ(1) * (-m(2) - t119) + t114) * g(3) + t102, -mrSges(2,3) * g(3) + t102, t44, t97 + t121, t107 + t121; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t103, -mrSges(2,3) * g(1) - t95, t41, t98, -mrSges(5,1) * t83 - Ifges(5,6) * qJDD(2) + t112; t95 + (qJ(1) * m(2) - t114) * g(1) - qJ(1) * t106 - mrSges(1,1) * g(2), t103, t96, -mrSges(4,1) * g(3) - t94 * Ifges(4,4) + t113 * qJDD(2) - t120, t109;];
m_new  = t1;
