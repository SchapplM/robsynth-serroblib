% Calculate vector of cutting forces with Newton-Euler
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:26:29
% EndTime: 2019-05-05 14:26:31
% DurationCPUTime: 0.61s
% Computational Cost: add. (3834->155), mult. (7295->177), div. (0->0), fcn. (3177->6), ass. (0->75)
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t110 = t76 * g(1) - t79 * g(2);
t81 = qJD(1) ^ 2;
t36 = -qJDD(1) * pkin(1) - t81 * qJ(2) + qJDD(2) - t110;
t122 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t36;
t97 = -t79 * g(1) - t76 * g(2);
t121 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t97;
t120 = -2 * qJD(5);
t119 = -m(3) - m(4);
t118 = pkin(4) + pkin(8);
t117 = pkin(8) * t81;
t75 = sin(qJ(4));
t116 = t75 * g(3);
t115 = t81 * pkin(7);
t85 = qJDD(3) + (-pkin(1) - qJ(3)) * t81 + t121;
t28 = -qJDD(1) * pkin(7) + t85;
t78 = cos(qJ(4));
t114 = t78 * t28;
t113 = mrSges(3,2) - mrSges(4,3);
t112 = -mrSges(4,2) - mrSges(3,3);
t111 = -mrSges(5,3) - mrSges(6,1);
t109 = qJD(1) * t75;
t108 = t78 * qJD(1);
t107 = -m(2) + t119;
t106 = qJD(1) * qJD(4);
t103 = mrSges(2,1) - t113;
t63 = t75 * t106;
t53 = t78 * qJDD(1) - t63;
t54 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t109;
t56 = mrSges(6,1) * t109 - (qJD(4) * mrSges(6,3));
t100 = t78 * t106;
t52 = t75 * qJDD(1) + t100;
t58 = pkin(5) * t108 - qJD(4) * pkin(8);
t72 = t75 ^ 2;
t83 = pkin(4) * t100 + t108 * t120 - t122 + (-t53 + t63) * qJ(5);
t12 = t83 - t58 * t108 + (-pkin(5) * t72 - pkin(7)) * t81 + t118 * t52;
t49 = (pkin(4) * t75 - qJ(5) * t78) * qJD(1);
t80 = qJD(4) ^ 2;
t90 = -t80 * qJ(5) + t49 * t108 + qJDD(5) - t114;
t15 = t53 * pkin(5) - t118 * qJDD(4) + (pkin(5) * t106 + t78 * t117 - g(3)) * t75 + t90;
t74 = sin(qJ(6));
t77 = cos(qJ(6));
t47 = -t74 * qJD(4) + t77 * t109;
t23 = t47 * qJD(6) + t77 * qJDD(4) + t74 * t52;
t48 = t77 * qJD(4) + t74 * t109;
t25 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t61 = qJD(6) + t108;
t33 = -t61 * mrSges(7,2) + t47 * mrSges(7,3);
t46 = qJDD(6) + t53;
t10 = m(7) * (-t74 * t12 + t77 * t15) - t23 * mrSges(7,3) + t46 * mrSges(7,1) - t48 * t25 + t61 * t33;
t22 = -t48 * qJD(6) - t74 * qJDD(4) + t77 * t52;
t34 = t61 * mrSges(7,1) - t48 * mrSges(7,3);
t11 = m(7) * (t77 * t12 + t74 * t15) + t22 * mrSges(7,3) - t46 * mrSges(7,2) + t47 * t25 - t61 * t34;
t91 = -m(6) * (-qJDD(4) * pkin(4) - t116 + t90) - t77 * t10 - t74 * t11;
t50 = (-mrSges(6,2) * t75 - mrSges(6,3) * t78) * qJD(1);
t98 = qJD(1) * (-t50 - (mrSges(5,1) * t75 + mrSges(5,2) * t78) * qJD(1));
t4 = m(5) * (t114 + t116) + t111 * t53 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t54 - t56) * qJD(4) + t78 * t98 + t91;
t101 = -t78 * g(3) + t75 * t28;
t55 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t108;
t57 = mrSges(6,1) * t108 + qJD(4) * mrSges(6,2);
t84 = -t80 * pkin(4) + qJDD(4) * qJ(5) - t49 * t109 + t101;
t89 = -t22 * mrSges(7,1) - t47 * t33 + m(7) * (-t72 * t117 - t52 * pkin(5) + ((2 * qJD(5)) + t58) * qJD(4) + t84) + t23 * mrSges(7,2) + t48 * t34;
t88 = -m(6) * ((qJD(4) * t120) - t84) + t89;
t7 = m(5) * t101 + t111 * t52 + (-mrSges(5,2) + mrSges(6,3)) * qJDD(4) + (-t55 + t57) * qJD(4) + t75 * t98 + t88;
t102 = -t75 * t4 + t78 * t7;
t99 = m(4) * t85 + qJDD(1) * mrSges(4,2) + t78 * t4 + t75 * t7;
t95 = m(3) * (t81 * pkin(1) - t121) - t99;
t93 = -t74 * t10 + t77 * t11 + m(6) * (t52 * pkin(4) - t115 + t83) - t56 * t109 - t57 * t108 - t52 * mrSges(6,2) - t53 * mrSges(6,3);
t87 = m(5) * (-t122 - t115) + t53 * mrSges(5,2) + t52 * mrSges(5,1) + t55 * t108 + t54 * t109 + t93;
t86 = m(4) * t122 - t87;
t82 = m(3) * t36 + t86;
t2 = -t82 + m(2) * t110 + (-mrSges(2,2) - t112) * t81 + t103 * qJDD(1);
t1 = m(2) * t97 + (-mrSges(2,2) + mrSges(3,3)) * qJDD(1) - t103 * t81 - t95;
t3 = [-m(1) * g(1) + t79 * t1 - t76 * t2, t1, t119 * g(3) + t102, -m(4) * g(3) + t102, t7, t93, t11; -m(1) * g(2) + t76 * t1 + t79 * t2, t2, -qJDD(1) * mrSges(3,3) - t113 * t81 + t95, -t81 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t86, t4, t52 * mrSges(6,1) - qJDD(4) * mrSges(6,3) - qJD(4) * t57 + t50 * t109 - t88, t10; (-m(1) + t107) * g(3) + t102, t107 * g(3) + t102, t113 * qJDD(1) + t112 * t81 + t82, -t81 * mrSges(4,3) + t99, t87, t53 * mrSges(6,1) + qJDD(4) * mrSges(6,2) + qJD(4) * t56 + t50 * t108 - t91, t89;];
f_new  = t3;
