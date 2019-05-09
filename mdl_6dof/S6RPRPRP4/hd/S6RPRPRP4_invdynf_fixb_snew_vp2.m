% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-05-05 17:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:42:34
% EndTime: 2019-05-05 17:42:37
% DurationCPUTime: 1.00s
% Computational Cost: add. (8444->174), mult. (16431->209), div. (0->0), fcn. (8587->8), ass. (0->85)
t133 = -2 * qJD(4);
t84 = sin(qJ(3));
t112 = t84 * qJD(1);
t67 = mrSges(5,1) * t112 + qJD(3) * mrSges(5,2);
t114 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t112 - t67;
t122 = mrSges(4,1) - mrSges(5,2);
t90 = qJD(1) ^ 2;
t125 = t90 * pkin(7);
t111 = qJD(1) * qJD(3);
t87 = cos(qJ(3));
t105 = t87 * t111;
t61 = t84 * qJDD(1) + t105;
t106 = t84 * t111;
t62 = t87 * qJDD(1) - t106;
t113 = qJD(1) * t87;
t65 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t113;
t85 = sin(qJ(1));
t88 = cos(qJ(1));
t107 = t85 * g(1) - t88 * g(2);
t56 = qJDD(1) * pkin(1) + t107;
t103 = -t88 * g(1) - t85 * g(2);
t60 = -t90 * pkin(1) + t103;
t81 = sin(pkin(9));
t82 = cos(pkin(9));
t104 = t82 * t56 - t81 * t60;
t98 = -qJDD(1) * pkin(2) - t104;
t128 = -pkin(3) - pkin(8);
t68 = pkin(4) * t112 - qJD(3) * pkin(8);
t79 = t87 ^ 2;
t93 = pkin(3) * t106 + t112 * t133 + (-t61 - t105) * qJ(4) + t98;
t17 = -t68 * t112 + (-pkin(4) * t79 - pkin(7)) * t90 + t128 * t62 + t93;
t116 = t81 * t56 + t82 * t60;
t29 = -t90 * pkin(2) + qJDD(1) * pkin(7) + t116;
t26 = t84 * t29;
t57 = (-pkin(3) * t87 - qJ(4) * t84) * qJD(1);
t89 = qJD(3) ^ 2;
t101 = -t89 * qJ(4) + t57 * t112 + qJDD(4) + t26;
t126 = pkin(8) * t90;
t80 = -g(3) + qJDD(2);
t21 = t61 * pkin(4) + t128 * qJDD(3) + (-pkin(4) * t111 - t84 * t126 - t80) * t87 + t101;
t83 = sin(qJ(5));
t86 = cos(qJ(5));
t119 = t86 * t17 + t83 * t21;
t54 = t83 * qJD(3) + t86 * t113;
t55 = t86 * qJD(3) - t83 * t113;
t37 = t54 * pkin(5) - t55 * qJ(6);
t71 = qJD(5) + t112;
t43 = -t71 * mrSges(7,1) + t55 * mrSges(7,2);
t53 = qJDD(5) + t61;
t69 = t71 ^ 2;
t108 = m(7) * (-t69 * pkin(5) + t53 * qJ(6) + 0.2e1 * qJD(6) * t71 - t54 * t37 + t119) + t71 * t43 + t53 * mrSges(7,3);
t38 = t54 * mrSges(7,1) - t55 * mrSges(7,3);
t117 = -t54 * mrSges(6,1) - t55 * mrSges(6,2) - t38;
t120 = -mrSges(6,3) - mrSges(7,2);
t34 = t55 * qJD(5) + t83 * qJDD(3) + t86 * t62;
t42 = t71 * mrSges(6,1) - t55 * mrSges(6,3);
t10 = m(6) * t119 - t53 * mrSges(6,2) + t117 * t54 + t120 * t34 - t71 * t42 + t108;
t102 = -t83 * t17 + t86 * t21;
t127 = m(7) * (-t53 * pkin(5) - t69 * qJ(6) + t55 * t37 + qJDD(6) - t102);
t35 = -t54 * qJD(5) + t86 * qJDD(3) - t83 * t62;
t40 = -t71 * mrSges(6,2) - t54 * mrSges(6,3);
t41 = -t54 * mrSges(7,2) + t71 * mrSges(7,3);
t11 = m(6) * t102 - t127 + (t40 + t41) * t71 + t117 * t55 + (mrSges(6,1) + mrSges(7,1)) * t53 + t120 * t35;
t66 = -mrSges(5,1) * t113 - qJD(3) * mrSges(5,3);
t99 = -t86 * t10 + t83 * t11 - m(5) * (-t62 * pkin(3) - t125 + t93) - t66 * t113 + t61 * mrSges(5,3);
t132 = (t114 * t84 - t87 * t65) * qJD(1) - t122 * t62 + m(4) * (t98 - t125) + t61 * mrSges(4,2) - t99;
t118 = t87 * t29 + t84 * t80;
t131 = t89 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t133 - t57 * t113 - t118;
t123 = t87 * t80;
t121 = mrSges(4,3) + mrSges(5,1);
t58 = (mrSges(5,2) * t87 - mrSges(5,3) * t84) * qJD(1);
t115 = t58 + (-mrSges(4,1) * t87 + mrSges(4,2) * t84) * qJD(1);
t97 = -m(5) * (-qJDD(3) * pkin(3) + t101 - t123) - t83 * t10 - t86 * t11;
t6 = m(4) * (-t26 + t123) - t121 * t61 + t122 * qJDD(3) + (t65 - t66) * qJD(3) - t115 * t112 + t97;
t94 = t62 * pkin(4) + qJD(3) * t68 - t79 * t126 - t131;
t96 = -t35 * mrSges(7,3) - t55 * t43 + m(7) * (-0.2e1 * qJD(6) * t55 + (t54 * t71 - t35) * qJ(6) + (t55 * t71 + t34) * pkin(5) + t94) + t34 * mrSges(7,1) + t54 * t41;
t92 = m(6) * t94 + t34 * mrSges(6,1) + t35 * mrSges(6,2) + t54 * t40 + t55 * t42 + t96;
t91 = -m(5) * t131 + t92;
t8 = t91 + t121 * t62 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) - t114 * qJD(3) + m(4) * t118 + t115 * t113;
t109 = m(3) * t80 + t87 * t6 + t84 * t8;
t4 = m(3) * t104 + qJDD(1) * mrSges(3,1) - t90 * mrSges(3,2) - t132;
t3 = m(3) * t116 - t90 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t84 * t6 + t87 * t8;
t2 = m(2) * t103 - t90 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t82 * t3 - t81 * t4;
t1 = m(2) * t107 + qJDD(1) * mrSges(2,1) - t90 * mrSges(2,2) + t81 * t3 + t82 * t4;
t5 = [-m(1) * g(1) - t85 * t1 + t88 * t2, t2, t3, t8, t62 * mrSges(5,2) - t67 * t112 - t99, t10, -t34 * mrSges(7,2) - t54 * t38 + t108; -m(1) * g(2) + t88 * t1 + t85 * t2, t1, t4, t6, -t62 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t67 - t58 * t113 - t91, t11, t96; (-m(1) - m(2)) * g(3) + t109, -m(2) * g(3) + t109, t109, t132, t61 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t66 + t58 * t112 - t97, t92, -t53 * mrSges(7,1) + t35 * mrSges(7,2) + t55 * t38 - t71 * t41 + t127;];
f_new  = t5;
