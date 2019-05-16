% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:50:51
% EndTime: 2019-05-05 09:50:57
% DurationCPUTime: 2.66s
% Computational Cost: add. (37554->172), mult. (72147->224), div. (0->0), fcn. (50846->12), ass. (0->93)
t87 = sin(pkin(11));
t89 = cos(pkin(11));
t76 = t87 * g(1) - t89 * g(2);
t90 = cos(pkin(6));
t123 = t76 * t90;
t77 = -t89 * g(1) - t87 * g(2);
t86 = -g(3) + qJDD(1);
t88 = sin(pkin(6));
t94 = sin(qJ(2));
t97 = cos(qJ(2));
t128 = (t86 * t88 + t123) * t97 - t94 * t77;
t96 = cos(qJ(3));
t115 = t96 * qJD(2);
t122 = t88 * t94;
t112 = t86 * t122 + t94 * t123 + t97 * t77;
t99 = qJD(2) ^ 2;
t44 = -t99 * pkin(2) + qJDD(2) * pkin(8) + t112;
t60 = -t88 * t76 + t90 * t86;
t93 = sin(qJ(3));
t117 = t96 * t44 + t93 * t60;
t73 = (-pkin(3) * t96 - pkin(9) * t93) * qJD(2);
t98 = qJD(3) ^ 2;
t31 = -t98 * pkin(3) + qJDD(3) * pkin(9) + t73 * t115 + t117;
t114 = qJD(2) * qJD(3);
t110 = t96 * t114;
t43 = -qJDD(2) * pkin(2) - t99 * pkin(8) - t128;
t74 = t93 * qJDD(2) + t110;
t84 = t93 * t114;
t75 = t96 * qJDD(2) - t84;
t34 = (-t74 - t110) * pkin(9) + (-t75 + t84) * pkin(3) + t43;
t92 = sin(qJ(4));
t95 = cos(qJ(4));
t109 = -t92 * t31 + t95 * t34;
t124 = cos(qJ(5));
t116 = qJD(2) * t93;
t70 = t95 * qJD(3) - t92 * t116;
t52 = t70 * qJD(4) + t92 * qJDD(3) + t95 * t74;
t67 = qJDD(4) - t75;
t71 = t92 * qJD(3) + t95 * t116;
t83 = qJD(4) - t115;
t20 = (t70 * t83 - t52) * pkin(10) + (t70 * t71 + t67) * pkin(4) + t109;
t119 = t95 * t31 + t92 * t34;
t51 = -t71 * qJD(4) + t95 * qJDD(3) - t92 * t74;
t59 = t83 * pkin(4) - t71 * pkin(10);
t66 = t70 ^ 2;
t22 = -t66 * pkin(4) + t51 * pkin(10) - t83 * t59 + t119;
t91 = sin(qJ(5));
t120 = t124 * t22 + t91 * t20;
t53 = -t124 * t70 + t91 * t71;
t54 = t124 * t71 + t91 * t70;
t37 = t53 * pkin(5) - t54 * qJ(6);
t82 = qJD(5) + t83;
t48 = -t82 * mrSges(7,1) + t54 * mrSges(7,2);
t65 = qJDD(5) + t67;
t81 = t82 ^ 2;
t113 = m(7) * (-t81 * pkin(5) + t65 * qJ(6) + 0.2e1 * qJD(6) * t82 - t53 * t37 + t120) + t82 * t48 + t65 * mrSges(7,3);
t38 = t53 * mrSges(7,1) - t54 * mrSges(7,3);
t118 = -t53 * mrSges(6,1) - t54 * mrSges(6,2) - t38;
t121 = -mrSges(6,3) - mrSges(7,2);
t27 = t54 * qJD(5) - t124 * t51 + t91 * t52;
t47 = t82 * mrSges(6,1) - t54 * mrSges(6,3);
t13 = m(6) * t120 - t65 * mrSges(6,2) + t118 * t53 + t121 * t27 - t82 * t47 + t113;
t105 = t124 * t20 - t91 * t22;
t126 = m(7) * (-t65 * pkin(5) - t81 * qJ(6) + t54 * t37 + qJDD(6) - t105);
t28 = -t53 * qJD(5) + t124 * t52 + t91 * t51;
t45 = -t53 * mrSges(7,2) + t82 * mrSges(7,3);
t46 = -t82 * mrSges(6,2) - t53 * mrSges(6,3);
t14 = m(6) * t105 - t126 + (t46 + t45) * t82 + (mrSges(6,1) + mrSges(7,1)) * t65 + t118 * t54 + t121 * t28;
t55 = -t70 * mrSges(5,1) + t71 * mrSges(5,2);
t57 = -t83 * mrSges(5,2) + t70 * mrSges(5,3);
t10 = m(5) * t109 + t67 * mrSges(5,1) - t52 * mrSges(5,3) + t124 * t14 + t91 * t13 - t71 * t55 + t83 * t57;
t58 = t83 * mrSges(5,1) - t71 * mrSges(5,3);
t11 = m(5) * t119 - t67 * mrSges(5,2) + t51 * mrSges(5,3) + t124 * t13 - t91 * t14 + t70 * t55 - t83 * t58;
t78 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t116;
t79 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t115;
t127 = m(4) * t43 - t75 * mrSges(4,1) + t74 * mrSges(4,2) + t95 * t10 + t92 * t11 + (t93 * t78 - t96 * t79) * qJD(2);
t8 = m(3) * t128 + qJDD(2) * mrSges(3,1) - t99 * mrSges(3,2) - t127;
t125 = t8 * t97;
t108 = -t93 * t44 + t96 * t60;
t30 = -qJDD(3) * pkin(3) - t98 * pkin(9) + t73 * t116 - t108;
t102 = -t51 * pkin(4) - t66 * pkin(10) + t71 * t59 + t30;
t104 = t28 * mrSges(7,3) + t54 * t48 - m(7) * (-0.2e1 * qJD(6) * t54 + (t53 * t82 - t28) * qJ(6) + (t54 * t82 + t27) * pkin(5) + t102) - t27 * mrSges(7,1) - t53 * t45;
t101 = m(6) * t102 + t27 * mrSges(6,1) + t28 * mrSges(6,2) + t53 * t46 + t54 * t47 - t104;
t100 = m(5) * t30 - t51 * mrSges(5,1) + t52 * mrSges(5,2) - t70 * t57 + t71 * t58 + t101;
t72 = (-mrSges(4,1) * t96 + mrSges(4,2) * t93) * qJD(2);
t12 = m(4) * t108 + qJDD(3) * mrSges(4,1) - t74 * mrSges(4,3) + qJD(3) * t79 - t72 * t116 - t100;
t9 = m(4) * t117 - qJDD(3) * mrSges(4,2) + t75 * mrSges(4,3) - qJD(3) * t78 - t92 * t10 + t95 * t11 + t72 * t115;
t4 = m(3) * t112 - t99 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t93 * t12 + t96 * t9;
t6 = m(3) * t60 + t96 * t12 + t93 * t9;
t111 = m(2) * t86 + t4 * t122 + t88 * t125 + t90 * t6;
t2 = m(2) * t77 + t97 * t4 - t94 * t8;
t1 = m(2) * t76 - t88 * t6 + (t4 * t94 + t125) * t90;
t3 = [-m(1) * g(1) - t87 * t1 + t89 * t2, t2, t4, t9, t11, t13, -t27 * mrSges(7,2) - t53 * t38 + t113; -m(1) * g(2) + t89 * t1 + t87 * t2, t1, t8, t12, t10, t14, -t104; -m(1) * g(3) + t111, t111, t6, t127, t100, t101, -t65 * mrSges(7,1) + t28 * mrSges(7,2) + t54 * t38 - t82 * t45 + t126;];
f_new  = t3;
