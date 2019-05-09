% Calculate vector of cutting torques with Newton-Euler for
% S4PRPP2
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
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2019-05-04 18:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRPP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:54:21
% EndTime: 2019-05-04 18:54:21
% DurationCPUTime: 0.25s
% Computational Cost: add. (2392->104), mult. (3115->114), div. (0->0), fcn. (1236->4), ass. (0->41)
t100 = qJD(2) ^ 2;
t91 = -g(2) + qJDD(1);
t96 = sin(qJ(2));
t97 = cos(qJ(2));
t76 = t96 * g(1) + t97 * t91;
t72 = qJDD(2) * pkin(2) + t76;
t77 = -t97 * g(1) + t96 * t91;
t73 = -t100 * pkin(2) + t77;
t94 = sin(pkin(5));
t95 = cos(pkin(5));
t69 = t94 * t72 + t95 * t73;
t64 = -t100 * pkin(3) + qJDD(2) * qJ(4) + (2 * qJD(4) * qJD(2)) + t69;
t115 = mrSges(5,2) * t64;
t114 = -mrSges(4,1) - mrSges(5,1);
t113 = -mrSges(1,2) - mrSges(2,3);
t112 = m(5) * t64 + qJDD(2) * mrSges(5,3);
t56 = m(4) * t69 - qJDD(2) * mrSges(4,2) + t114 * t100 + t112;
t68 = t95 * t72 - t94 * t73;
t66 = -qJDD(2) * pkin(3) - t100 * qJ(4) + qJDD(4) - t68;
t107 = -m(5) * t66 + qJDD(2) * mrSges(5,1) + t100 * mrSges(5,3);
t57 = m(4) * t68 + qJDD(2) * mrSges(4,1) - t100 * mrSges(4,2) + t107;
t50 = t94 * t56 + t95 * t57;
t90 = -g(3) + qJDD(3);
t111 = (m(4) + m(5)) * t90;
t110 = mrSges(5,2) * t66 + Ifges(5,4) * qJDD(2) + t100 * Ifges(5,6);
t47 = m(3) * t76 + qJDD(2) * mrSges(3,1) - t100 * mrSges(3,2) + t50;
t108 = t95 * t56 - t94 * t57;
t48 = m(3) * t77 - t100 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t108;
t109 = -t96 * t47 + t97 * t48;
t106 = -mrSges(5,1) * t66 + mrSges(5,3) * t64 + Ifges(5,2) * qJDD(2);
t51 = t115 + mrSges(4,3) * t69 + (Ifges(5,4) + Ifges(4,5)) * t100 + (Ifges(4,6) - Ifges(5,6)) * qJDD(2) + (-m(5) * pkin(3) + t114) * t90;
t52 = -mrSges(4,3) * t68 + Ifges(4,5) * qJDD(2) - t100 * Ifges(4,6) + (-m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3)) * t90 + t110;
t40 = mrSges(3,1) * g(3) + mrSges(3,3) * t77 + t100 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t111 + qJ(3) * t108 + t95 * t51 + t94 * t52;
t43 = -mrSges(3,2) * g(3) - mrSges(3,3) * t76 + Ifges(3,5) * qJDD(2) - t100 * Ifges(3,6) - qJ(3) * t50 - t94 * t51 + t95 * t52;
t105 = pkin(4) * t109 + t97 * t40 + t96 * t43 + pkin(1) * (m(3) * g(3) - t111) + mrSges(2,2) * g(1) + mrSges(2,1) * g(3);
t45 = t97 * t47 + t96 * t48;
t104 = mrSges(2,2) * t91 - pkin(4) * t45 - t96 * t40 + t97 * t43;
t103 = -mrSges(4,2) * t69 + t106 + qJ(4) * (-t100 * mrSges(5,1) + t112) + pkin(3) * t107 + mrSges(4,1) * t68 + Ifges(4,3) * qJDD(2);
t102 = mrSges(3,1) * t76 - mrSges(3,2) * t77 + Ifges(3,3) * qJDD(2) + pkin(2) * t50 + t103;
t101 = mrSges(2,1) * t91 + pkin(1) * t45 + t102;
t1 = [mrSges(1,3) * g(2) + qJ(1) * t111 + (qJ(1) * (-m(2) - m(3)) + t113) * g(3) + t104, -mrSges(2,3) * g(3) + t104, t43, t52, -mrSges(5,3) * t90 + t110; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t105, -mrSges(2,3) * g(1) - t101, t40, t51, t106; t101 - qJ(1) * t109 + (qJ(1) * m(2) - t113) * g(1) - mrSges(1,1) * g(2), t105, t102, t103, mrSges(5,1) * t90 - t100 * Ifges(5,4) + Ifges(5,6) * qJDD(2) - t115;];
m_new  = t1;
