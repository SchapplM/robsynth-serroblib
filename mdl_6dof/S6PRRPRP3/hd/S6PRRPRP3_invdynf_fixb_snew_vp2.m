% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRP3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:51:49
% EndTime: 2019-05-05 03:51:54
% DurationCPUTime: 2.52s
% Computational Cost: add. (33951->171), mult. (69689->226), div. (0->0), fcn. (48629->12), ass. (0->91)
t88 = sin(pkin(10));
t91 = cos(pkin(10));
t77 = t88 * g(1) - t91 * g(2);
t92 = cos(pkin(6));
t123 = t77 * t92;
t78 = -t91 * g(1) - t88 * g(2);
t86 = -g(3) + qJDD(1);
t89 = sin(pkin(6));
t95 = sin(qJ(2));
t97 = cos(qJ(2));
t128 = (t86 * t89 + t123) * t97 - t95 * t78;
t96 = cos(qJ(3));
t116 = t96 * qJD(2);
t122 = t89 * t95;
t112 = t86 * t122 + t95 * t123 + t97 * t78;
t99 = qJD(2) ^ 2;
t44 = -t99 * pkin(2) + qJDD(2) * pkin(8) + t112;
t60 = -t89 * t77 + t92 * t86;
t94 = sin(qJ(3));
t118 = t96 * t44 + t94 * t60;
t73 = (-pkin(3) * t96 - qJ(4) * t94) * qJD(2);
t98 = qJD(3) ^ 2;
t27 = -t98 * pkin(3) + qJDD(3) * qJ(4) + t73 * t116 + t118;
t115 = qJD(2) * qJD(3);
t110 = t96 * t115;
t43 = -qJDD(2) * pkin(2) - t99 * pkin(8) - t128;
t75 = t94 * qJDD(2) + t110;
t84 = t94 * t115;
t76 = t96 * qJDD(2) - t84;
t32 = (-t75 - t110) * qJ(4) + (-t76 + t84) * pkin(3) + t43;
t117 = qJD(2) * t94;
t87 = sin(pkin(11));
t90 = cos(pkin(11));
t69 = t87 * qJD(3) + t90 * t117;
t108 = -0.2e1 * qJD(4) * t69 - t87 * t27 + t90 * t32;
t124 = cos(qJ(5));
t58 = t87 * qJDD(3) + t90 * t75;
t68 = t90 * qJD(3) - t87 * t117;
t20 = (-t68 * t116 - t58) * pkin(9) + (t68 * t69 - t76) * pkin(4) + t108;
t113 = 0.2e1 * qJD(4) * t68 + t90 * t27 + t87 * t32;
t57 = t90 * qJDD(3) - t87 * t75;
t59 = -pkin(4) * t116 - t69 * pkin(9);
t67 = t68 ^ 2;
t22 = -t67 * pkin(4) + t57 * pkin(9) + t59 * t116 + t113;
t93 = sin(qJ(5));
t120 = t124 * t22 + t93 * t20;
t50 = -t124 * t68 + t93 * t69;
t51 = t124 * t69 + t93 * t68;
t37 = t50 * pkin(5) - t51 * qJ(6);
t83 = qJD(5) - t116;
t47 = -t83 * mrSges(7,1) + t51 * mrSges(7,2);
t70 = qJDD(5) - t76;
t82 = t83 ^ 2;
t114 = m(7) * (-t82 * pkin(5) + t70 * qJ(6) + 0.2e1 * qJD(6) * t83 - t50 * t37 + t120) + t83 * t47 + t70 * mrSges(7,3);
t38 = t50 * mrSges(7,1) - t51 * mrSges(7,3);
t119 = -t50 * mrSges(6,1) - t51 * mrSges(6,2) - t38;
t121 = -mrSges(6,3) - mrSges(7,2);
t33 = t51 * qJD(5) - t124 * t57 + t93 * t58;
t46 = t83 * mrSges(6,1) - t51 * mrSges(6,3);
t12 = m(6) * t120 - t70 * mrSges(6,2) + t119 * t50 + t121 * t33 - t83 * t46 + t114;
t105 = t124 * t20 - t93 * t22;
t126 = m(7) * (-t70 * pkin(5) - t82 * qJ(6) + t51 * t37 + qJDD(6) - t105);
t34 = -t50 * qJD(5) + t124 * t58 + t93 * t57;
t45 = -t83 * mrSges(6,2) - t50 * mrSges(6,3);
t48 = -t50 * mrSges(7,2) + t83 * mrSges(7,3);
t13 = m(6) * t105 - t126 + (t45 + t48) * t83 + (mrSges(6,1) + mrSges(7,1)) * t70 + t119 * t51 + t121 * t34;
t52 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t55 = mrSges(5,2) * t116 + t68 * mrSges(5,3);
t10 = m(5) * t108 - t76 * mrSges(5,1) - t58 * mrSges(5,3) - t55 * t116 + t93 * t12 + t124 * t13 - t69 * t52;
t56 = -mrSges(5,1) * t116 - t69 * mrSges(5,3);
t11 = m(5) * t113 + t76 * mrSges(5,2) + t57 * mrSges(5,3) + t56 * t116 + t124 * t12 - t93 * t13 + t68 * t52;
t79 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t117;
t80 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t116;
t127 = m(4) * t43 - t76 * mrSges(4,1) + t75 * mrSges(4,2) + t90 * t10 + t87 * t11 + (t94 * t79 - t96 * t80) * qJD(2);
t8 = m(3) * t128 + qJDD(2) * mrSges(3,1) - t99 * mrSges(3,2) - t127;
t125 = t8 * t97;
t109 = -t94 * t44 + t96 * t60;
t26 = -qJDD(3) * pkin(3) - t98 * qJ(4) + t73 * t117 + qJDD(4) - t109;
t101 = -t57 * pkin(4) - t67 * pkin(9) + t69 * t59 + t26;
t104 = t34 * mrSges(7,3) + t51 * t47 - m(7) * (-0.2e1 * qJD(6) * t51 + (t50 * t83 - t34) * qJ(6) + (t51 * t83 + t33) * pkin(5) + t101) - t33 * mrSges(7,1) - t50 * t48;
t102 = m(6) * t101 + t33 * mrSges(6,1) + t34 * mrSges(6,2) + t50 * t45 + t51 * t46 - t104;
t100 = m(5) * t26 - t57 * mrSges(5,1) + t58 * mrSges(5,2) - t68 * t55 + t69 * t56 + t102;
t74 = (-mrSges(4,1) * t96 + mrSges(4,2) * t94) * qJD(2);
t14 = m(4) * t109 + qJDD(3) * mrSges(4,1) - t75 * mrSges(4,3) + qJD(3) * t80 - t74 * t117 - t100;
t9 = m(4) * t118 - qJDD(3) * mrSges(4,2) + t76 * mrSges(4,3) - qJD(3) * t79 - t87 * t10 + t90 * t11 + t74 * t116;
t4 = m(3) * t112 - t99 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t94 * t14 + t96 * t9;
t6 = m(3) * t60 + t96 * t14 + t94 * t9;
t111 = m(2) * t86 + t4 * t122 + t89 * t125 + t92 * t6;
t2 = m(2) * t78 + t97 * t4 - t95 * t8;
t1 = m(2) * t77 - t89 * t6 + (t4 * t95 + t125) * t92;
t3 = [-m(1) * g(1) - t88 * t1 + t91 * t2, t2, t4, t9, t11, t12, -t33 * mrSges(7,2) - t50 * t38 + t114; -m(1) * g(2) + t91 * t1 + t88 * t2, t1, t8, t14, t10, t13, -t104; -m(1) * g(3) + t111, t111, t6, t127, t100, t102, -t70 * mrSges(7,1) + t34 * mrSges(7,2) + t51 * t38 - t83 * t48 + t126;];
f_new  = t3;
