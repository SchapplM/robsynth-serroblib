% Calculate vector of cutting torques with Newton-Euler for
% S4PRRP1
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-05-04 19:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:02:13
% EndTime: 2019-05-04 19:02:13
% DurationCPUTime: 0.39s
% Computational Cost: add. (4744->109), mult. (6612->127), div. (0->0), fcn. (3974->6), ass. (0->48)
t100 = (qJD(2) + qJD(3));
t105 = sin(qJ(3));
t107 = cos(qJ(3));
t106 = sin(qJ(2));
t108 = cos(qJ(2));
t103 = sin(pkin(6));
t104 = cos(pkin(6));
t87 = t103 * g(1) - t104 * g(2);
t88 = -t104 * g(1) - t103 * g(2);
t81 = -t106 * t88 + t108 * t87;
t78 = qJDD(2) * pkin(2) + t81;
t109 = qJD(2) ^ 2;
t82 = t106 * t87 + t108 * t88;
t79 = -t109 * pkin(2) + t82;
t75 = t105 * t78 + t107 * t79;
t98 = t100 ^ 2;
t99 = qJDD(2) + qJDD(3);
t70 = -t98 * pkin(3) + t99 * qJ(4) + (2 * qJD(4) * t100) + t75;
t121 = mrSges(5,2) * t70;
t120 = -mrSges(4,1) - mrSges(5,1);
t119 = m(5) * t70 + t99 * mrSges(5,3);
t63 = m(4) * t75 - t99 * mrSges(4,2) + t120 * t98 + t119;
t74 = -t105 * t79 + t107 * t78;
t72 = -t99 * pkin(3) - t98 * qJ(4) + qJDD(4) - t74;
t114 = -m(5) * t72 + t99 * mrSges(5,1) + t98 * mrSges(5,3);
t65 = m(4) * t74 + t99 * mrSges(4,1) - t98 * mrSges(4,2) + t114;
t58 = t105 * t63 + t107 * t65;
t55 = m(3) * t81 + qJDD(2) * mrSges(3,1) - t109 * mrSges(3,2) + t58;
t116 = -t105 * t65 + t107 * t63;
t56 = m(3) * t82 - t109 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t116;
t50 = t106 * t56 + t108 * t55;
t102 = -g(3) + qJDD(1);
t118 = (m(4) + m(5)) * t102;
t117 = mrSges(5,2) * t72 + Ifges(5,4) * t99 + t98 * Ifges(5,6);
t115 = -t106 * t55 + t108 * t56;
t113 = -mrSges(5,1) * t72 + mrSges(5,3) * t70 + Ifges(5,2) * t99;
t112 = -mrSges(4,2) * t75 + t113 + qJ(4) * (-t98 * mrSges(5,1) + t119) + pkin(3) * t114 + mrSges(4,1) * t74 + Ifges(4,3) * t99;
t111 = mrSges(3,1) * t81 - mrSges(3,2) * t82 + Ifges(3,3) * qJDD(2) + pkin(2) * t58 + t112;
t110 = mrSges(2,1) * t87 - mrSges(2,2) * t88 + pkin(1) * t50 + t111;
t60 = -mrSges(4,3) * t74 + Ifges(4,5) * t99 - t98 * Ifges(4,6) + (-m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3)) * t102 + t117;
t59 = t121 + mrSges(4,3) * t75 + (Ifges(4,6) - Ifges(5,6)) * t99 + (Ifges(5,4) + Ifges(4,5)) * t98 + (-m(5) * pkin(3) + t120) * t102;
t51 = mrSges(3,2) * t102 - mrSges(3,3) * t81 + Ifges(3,5) * qJDD(2) - t109 * Ifges(3,6) - pkin(5) * t58 - t105 * t59 + t107 * t60;
t48 = -mrSges(3,1) * t102 + mrSges(3,3) * t82 + t109 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t118 + pkin(5) * t116 + t105 * t60 + t107 * t59;
t47 = m(2) * t88 + t115;
t46 = m(2) * t87 + t50;
t45 = mrSges(2,2) * t102 - mrSges(2,3) * t87 - pkin(4) * t50 - t106 * t48 + t108 * t51;
t44 = -mrSges(2,1) * t102 + mrSges(2,3) * t88 + t106 * t51 + t108 * t48 - pkin(1) * (m(3) * t102 + t118) + pkin(4) * t115;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t104 * t45 - t103 * t44 - qJ(1) * (t103 * t47 + t104 * t46), t45, t51, t60, -mrSges(5,3) * t102 + t117; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t103 * t45 + t104 * t44 + qJ(1) * (-t103 * t46 + t104 * t47), t44, t48, t59, t113; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t110, t110, t111, t112, mrSges(5,1) * t102 - t98 * Ifges(5,4) + Ifges(5,6) * t99 - t121;];
m_new  = t1;
