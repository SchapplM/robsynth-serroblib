% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-05-05 04:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:06:55
% EndTime: 2019-05-05 04:06:59
% DurationCPUTime: 1.19s
% Computational Cost: add. (11669->173), mult. (22721->214), div. (0->0), fcn. (13458->10), ass. (0->90)
t140 = -2 * qJD(4);
t84 = sin(pkin(10));
t86 = cos(pkin(10));
t66 = t84 * g(1) - t86 * g(2);
t87 = cos(pkin(6));
t129 = t66 * t87;
t67 = -t86 * g(1) - t84 * g(2);
t83 = -g(3) + qJDD(1);
t85 = sin(pkin(6));
t90 = sin(qJ(2));
t93 = cos(qJ(2));
t137 = (t83 * t85 + t129) * t93 - t90 * t67;
t101 = -qJDD(2) * pkin(2) - t137;
t92 = cos(qJ(3));
t117 = qJD(2) * t92;
t89 = sin(qJ(3));
t116 = t89 * qJD(2);
t128 = t85 * t90;
t112 = t83 * t128 + t90 * t129 + t93 * t67;
t95 = qJD(2) ^ 2;
t31 = -t95 * pkin(2) + qJDD(2) * pkin(8) + t112;
t28 = t89 * t31;
t60 = (-pkin(3) * t92 - qJ(4) * t89) * qJD(2);
t94 = qJD(3) ^ 2;
t106 = -t94 * qJ(4) + t60 * t116 + qJDD(4) + t28;
t115 = qJD(2) * qJD(3);
t132 = pkin(9) * t95;
t134 = -pkin(3) - pkin(9);
t47 = -t85 * t66 + t87 * t83;
t109 = t92 * t115;
t63 = t89 * qJDD(2) + t109;
t21 = t63 * pkin(4) + t134 * qJDD(3) + (-pkin(4) * t115 - t89 * t132 - t47) * t92 + t106;
t110 = t89 * t115;
t64 = t92 * qJDD(2) - t110;
t73 = pkin(4) * t116 - qJD(3) * pkin(9);
t82 = t92 ^ 2;
t96 = pkin(3) * t110 + t116 * t140 + (-t63 - t109) * qJ(4) + t101;
t23 = -t73 * t116 + (-pkin(4) * t82 - pkin(8)) * t95 + t134 * t64 + t96;
t88 = sin(qJ(5));
t91 = cos(qJ(5));
t122 = t88 * t21 + t91 * t23;
t58 = t88 * qJD(3) + t91 * t117;
t59 = t91 * qJD(3) - t88 * t117;
t39 = t58 * pkin(5) - t59 * qJ(6);
t76 = qJD(5) + t116;
t46 = -t76 * mrSges(7,1) + t59 * mrSges(7,2);
t55 = qJDD(5) + t63;
t74 = t76 ^ 2;
t113 = m(7) * (-t74 * pkin(5) + t55 * qJ(6) + 0.2e1 * qJD(6) * t76 - t58 * t39 + t122) + t76 * t46 + t55 * mrSges(7,3);
t40 = t58 * mrSges(7,1) - t59 * mrSges(7,3);
t120 = -t58 * mrSges(6,1) - t59 * mrSges(6,2) - t40;
t123 = -mrSges(6,3) - mrSges(7,2);
t36 = t59 * qJD(5) + t88 * qJDD(3) + t91 * t64;
t45 = t76 * mrSges(6,1) - t59 * mrSges(6,3);
t12 = m(6) * t122 - t55 * mrSges(6,2) + t120 * t58 + t123 * t36 - t76 * t45 + t113;
t108 = t91 * t21 - t88 * t23;
t133 = m(7) * (-t55 * pkin(5) - t74 * qJ(6) + t59 * t39 + qJDD(6) - t108);
t37 = -t58 * qJD(5) + t91 * qJDD(3) - t88 * t64;
t43 = -t76 * mrSges(6,2) - t58 * mrSges(6,3);
t44 = -t58 * mrSges(7,2) + t76 * mrSges(7,3);
t13 = m(6) * t108 - t133 + (t43 + t44) * t76 + t120 * t59 + (mrSges(6,1) + mrSges(7,1)) * t55 + t123 * t37;
t130 = t95 * pkin(8);
t70 = -mrSges(5,1) * t117 - qJD(3) * mrSges(5,3);
t104 = -t91 * t12 + t88 * t13 - m(5) * (-t64 * pkin(3) - t130 + t96) - t70 * t117 + t63 * mrSges(5,3);
t71 = mrSges(5,1) * t116 + qJD(3) * mrSges(5,2);
t118 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t116 - t71;
t125 = mrSges(4,1) - mrSges(5,2);
t69 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t117;
t139 = (t118 * t89 - t92 * t69) * qJD(2) - t125 * t64 + m(4) * (t101 - t130) + t63 * mrSges(4,2) - t104;
t121 = t92 * t31 + t89 * t47;
t138 = t94 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t140 - t60 * t117 - t121;
t8 = m(3) * t137 + qJDD(2) * mrSges(3,1) - t95 * mrSges(3,2) - t139;
t131 = t8 * t93;
t127 = t92 * t47;
t124 = mrSges(4,3) + mrSges(5,1);
t61 = (mrSges(5,2) * t92 - mrSges(5,3) * t89) * qJD(2);
t119 = t61 + (-mrSges(4,1) * t92 + mrSges(4,2) * t89) * qJD(2);
t99 = t64 * pkin(4) + qJD(3) * t73 - t82 * t132 - t138;
t102 = -t37 * mrSges(7,3) - t59 * t46 + m(7) * (-0.2e1 * qJD(6) * t59 + (t58 * t76 - t37) * qJ(6) + (t59 * t76 + t36) * pkin(5) + t99) + t36 * mrSges(7,1) + t58 * t44;
t98 = m(6) * t99 + t36 * mrSges(6,1) + t37 * mrSges(6,2) + t58 * t43 + t59 * t45 + t102;
t97 = -m(5) * t138 + t98;
t10 = t97 + t124 * t64 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) - t118 * qJD(3) + m(4) * t121 + t119 * t117;
t103 = -m(5) * (-qJDD(3) * pkin(3) + t106 - t127) - t88 * t12 - t91 * t13;
t9 = m(4) * (-t28 + t127) - t124 * t63 + t125 * qJDD(3) + (t69 - t70) * qJD(3) - t119 * t116 + t103;
t4 = m(3) * t112 - t95 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t92 * t10 - t89 * t9;
t6 = m(3) * t47 + t89 * t10 + t92 * t9;
t111 = m(2) * t83 + t4 * t128 + t85 * t131 + t87 * t6;
t2 = m(2) * t67 + t93 * t4 - t90 * t8;
t1 = m(2) * t66 - t85 * t6 + (t4 * t90 + t131) * t87;
t3 = [-m(1) * g(1) - t84 * t1 + t86 * t2, t2, t4, t10, t64 * mrSges(5,2) - t71 * t116 - t104, t12, -t36 * mrSges(7,2) - t58 * t40 + t113; -m(1) * g(2) + t86 * t1 + t84 * t2, t1, t8, t9, -t64 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t71 - t61 * t117 - t97, t13, t102; -m(1) * g(3) + t111, t111, t6, t139, t63 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t70 + t61 * t116 - t103, t98, -t55 * mrSges(7,1) + t37 * mrSges(7,2) + t59 * t40 - t76 * t44 + t133;];
f_new  = t3;
