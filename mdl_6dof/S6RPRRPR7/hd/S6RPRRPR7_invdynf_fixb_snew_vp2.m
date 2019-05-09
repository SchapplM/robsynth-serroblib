% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-05-05 23:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:00:35
% EndTime: 2019-05-05 23:00:42
% DurationCPUTime: 2.64s
% Computational Cost: add. (34950->179), mult. (75226->233), div. (0->0), fcn. (52345->10), ass. (0->93)
t92 = sin(qJ(1));
t96 = cos(qJ(1));
t110 = -t96 * g(1) - t92 * g(2);
t107 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t110;
t127 = 2 * qJD(5);
t126 = -m(2) - m(3);
t125 = -pkin(1) - pkin(7);
t124 = (mrSges(2,1) - mrSges(3,2));
t123 = -mrSges(2,2) + mrSges(3,3);
t118 = qJD(1) * qJD(3);
t91 = sin(qJ(3));
t112 = t91 * t118;
t114 = t92 * g(1) - t96 * g(2);
t97 = qJD(1) ^ 2;
t105 = -t97 * qJ(2) + qJDD(2) - t114;
t63 = t125 * qJDD(1) + t105;
t95 = cos(qJ(3));
t121 = t91 * g(3) + t95 * t63;
t75 = qJDD(1) * t95 - t112;
t38 = (-t75 - t112) * pkin(8) + (-t91 * t95 * t97 + qJDD(3)) * pkin(3) + t121;
t113 = -g(3) * t95 + t91 * t63;
t74 = -qJDD(1) * t91 - t95 * t118;
t119 = qJD(1) * t95;
t78 = qJD(3) * pkin(3) - pkin(8) * t119;
t86 = t91 ^ 2;
t39 = -pkin(3) * t86 * t97 + pkin(8) * t74 - qJD(3) * t78 + t113;
t90 = sin(qJ(4));
t94 = cos(qJ(4));
t122 = t90 * t38 + t94 * t39;
t120 = qJD(1) * t91;
t111 = t94 * t38 - t90 * t39;
t67 = (-t90 * t95 - t91 * t94) * qJD(1);
t46 = t67 * qJD(4) + t74 * t90 + t75 * t94;
t68 = (-t90 * t91 + t94 * t95) * qJD(1);
t83 = qJDD(3) + qJDD(4);
t84 = qJD(3) + qJD(4);
t18 = (t67 * t84 - t46) * qJ(5) + (t67 * t68 + t83) * pkin(4) + t111;
t45 = -qJD(4) * t68 + t74 * t94 - t75 * t90;
t60 = pkin(4) * t84 - qJ(5) * t68;
t66 = t67 ^ 2;
t20 = -t66 * pkin(4) + t45 * qJ(5) - t60 * t84 + t122;
t87 = sin(pkin(10));
t88 = cos(pkin(10));
t53 = t67 * t88 - t68 * t87;
t116 = t53 * t127 + t87 * t18 + t88 * t20;
t55 = -t67 * mrSges(5,1) + mrSges(5,2) * t68;
t59 = -mrSges(5,2) * t84 + t67 * mrSges(5,3);
t54 = t67 * t87 + t68 * t88;
t34 = -pkin(5) * t53 - pkin(9) * t54;
t82 = t84 ^ 2;
t15 = -pkin(5) * t82 + pkin(9) * t83 + t53 * t34 + t116;
t102 = -t74 * pkin(3) + t78 * t119 + (-pkin(8) * t86 + t125) * t97 + t107;
t100 = -t45 * pkin(4) - t66 * qJ(5) + t68 * t60 + qJDD(5) + t102;
t27 = t45 * t88 - t46 * t87;
t28 = t45 * t87 + t46 * t88;
t16 = t100 + (-t53 * t84 - t28) * pkin(9) + (t54 * t84 - t27) * pkin(5);
t89 = sin(qJ(6));
t93 = cos(qJ(6));
t43 = -t54 * t89 + t84 * t93;
t22 = t43 * qJD(6) + t28 * t93 + t83 * t89;
t26 = qJDD(6) - t27;
t44 = t54 * t93 + t84 * t89;
t29 = -mrSges(7,1) * t43 + mrSges(7,2) * t44;
t50 = qJD(6) - t53;
t30 = -mrSges(7,2) * t50 + mrSges(7,3) * t43;
t12 = m(7) * (-t15 * t89 + t16 * t93) - t22 * mrSges(7,3) + t26 * mrSges(7,1) - t44 * t29 + t50 * t30;
t21 = -t44 * qJD(6) - t28 * t89 + t83 * t93;
t31 = mrSges(7,1) * t50 - mrSges(7,3) * t44;
t13 = m(7) * (t15 * t93 + t16 * t89) + t21 * mrSges(7,3) - t26 * mrSges(7,2) + t43 * t29 - t50 * t31;
t33 = -mrSges(6,1) * t53 + mrSges(6,2) * t54;
t48 = mrSges(6,1) * t84 - t54 * mrSges(6,3);
t8 = m(6) * t116 - t83 * mrSges(6,2) + t27 * mrSges(6,3) - t89 * t12 + t93 * t13 + t53 * t33 - t84 * t48;
t109 = -t88 * t18 + t87 * t20;
t103 = m(7) * (-t83 * pkin(5) - t82 * pkin(9) + (t127 + t34) * t54 + t109) - t21 * mrSges(7,1) + t22 * mrSges(7,2) - t43 * t30 + t44 * t31;
t47 = -mrSges(6,2) * t84 + t53 * mrSges(6,3);
t9 = m(6) * (-0.2e1 * qJD(5) * t54 - t109) - t28 * mrSges(6,3) + t83 * mrSges(6,1) - t54 * t33 + t84 * t47 - t103;
t5 = m(5) * t111 + t83 * mrSges(5,1) - t46 * mrSges(5,3) - t68 * t55 + t84 * t59 + t87 * t8 + t88 * t9;
t61 = mrSges(5,1) * t84 - mrSges(5,3) * t68;
t6 = m(5) * t122 - t83 * mrSges(5,2) + t45 * mrSges(5,3) + t67 * t55 - t84 * t61 + t88 * t8 - t87 * t9;
t73 = (mrSges(4,1) * t91 + mrSges(4,2) * t95) * qJD(1);
t76 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t120;
t3 = m(4) * t121 + qJDD(3) * mrSges(4,1) - t75 * mrSges(4,3) + qJD(3) * t76 - t73 * t119 + t94 * t5 + t90 * t6;
t77 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t119;
t4 = m(4) * t113 - qJDD(3) * mrSges(4,2) + t74 * mrSges(4,3) - qJD(3) * t77 - t73 * t120 - t90 * t5 + t94 * t6;
t115 = -t91 * t3 + t95 * t4;
t106 = -m(3) * (-qJDD(1) * pkin(1) + t105) - t95 * t3 - t91 * t4;
t104 = m(6) * t100 - t27 * mrSges(6,1) + t28 * mrSges(6,2) + t93 * t12 + t89 * t13 - t53 * t47 + t54 * t48;
t101 = m(5) * t102 - t45 * mrSges(5,1) + t46 * mrSges(5,2) - t67 * t59 + t68 * t61 + t104;
t99 = -t74 * mrSges(4,1) + t101 + m(4) * (t125 * t97 + t107) + t76 * t120 + t77 * t119 + t75 * mrSges(4,2);
t98 = -m(3) * (t97 * pkin(1) - t107) + t99;
t7 = m(2) * t110 + t123 * qJDD(1) - (t124 * t97) + t98;
t1 = m(2) * t114 + t124 * qJDD(1) + t123 * t97 + t106;
t2 = [-m(1) * g(1) - t1 * t92 + t7 * t96, t7, -m(3) * g(3) + t115, t4, t6, t8, t13; -m(1) * g(2) + t1 * t96 + t7 * t92, t1, -(t97 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t98, t3, t5, t9, t12; (-m(1) + t126) * g(3) + t115, t126 * g(3) + t115, qJDD(1) * mrSges(3,2) - t97 * mrSges(3,3) - t106, t99, t101, t104, t103;];
f_new  = t2;
