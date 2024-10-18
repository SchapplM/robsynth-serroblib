% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR15_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR15_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:24:24
% EndTime: 2024-09-27 22:24:31
% DurationCPUTime: 3.70s
% Computational Cost: add. (116632->182), mult. (326996->263), div. (0->0), fcn. (264472->14), ass. (0->107)
t100 = cos(qJ(4));
t101 = cos(qJ(3));
t102 = cos(qJ(2));
t104 = qJD(1) ^ 2;
t103 = cos(qJ(1));
t98 = sin(qJ(1));
t118 = t98 * g(1) - g(2) * t103;
t91 = sin(pkin(5));
t135 = pkin(8) * t91;
t73 = qJDD(1) * pkin(1) + t104 * t135 + t118;
t93 = cos(pkin(5));
t132 = t73 * t93;
t114 = -g(1) * t103 - g(2) * t98;
t74 = -pkin(1) * t104 + qJDD(1) * t135 + t114;
t97 = sin(qJ(2));
t115 = t102 * t132 - t97 * t74;
t125 = t104 * t91 ^ 2;
t122 = qJD(1) * qJD(2);
t77 = (qJDD(1) * t97 + t102 * t122) * t91;
t83 = t93 * qJDD(1) + qJDD(2);
t84 = t93 * qJD(1) + qJD(2);
t43 = t83 * pkin(2) - t77 * pkin(9) + (pkin(2) * t97 * t125 + (pkin(9) * qJD(1) * t84 - g(3)) * t91) * t102 + t115;
t131 = t91 * t97;
t109 = -g(3) * t131 + t102 * t74 + t97 * t132;
t121 = t102 ^ 2 * t125;
t124 = qJD(1) * t91;
t120 = t97 * t124;
t76 = pkin(2) * t84 - pkin(9) * t120;
t78 = (qJDD(1) * t102 - t122 * t97) * t91;
t44 = -pkin(2) * t121 + pkin(9) * t78 - t76 * t84 + t109;
t96 = sin(qJ(3));
t116 = t101 * t43 - t96 * t44;
t69 = (t101 * t102 - t96 * t97) * t124;
t53 = qJD(3) * t69 + t101 * t77 + t78 * t96;
t70 = (t101 * t97 + t102 * t96) * t124;
t81 = qJDD(3) + t83;
t82 = qJD(3) + t84;
t26 = (t69 * t82 - t53) * pkin(10) + (t69 * t70 + t81) * pkin(3) + t116;
t127 = t101 * t44 + t96 * t43;
t52 = -qJD(3) * t70 + t101 * t78 - t77 * t96;
t63 = pkin(3) * t82 - pkin(10) * t70;
t67 = t69 ^ 2;
t28 = -pkin(3) * t67 + pkin(10) * t52 - t63 * t82 + t127;
t95 = sin(qJ(4));
t117 = t100 * t26 - t28 * t95;
t92 = cos(pkin(6));
t133 = pkin(11) * t92;
t58 = t100 * t69 - t70 * t95;
t36 = qJD(4) * t58 + t100 * t53 + t52 * t95;
t90 = sin(pkin(6));
t134 = pkin(11) * t90;
t59 = t100 * t70 + t69 * t95;
t45 = -pkin(4) * t58 - t134 * t59;
t80 = qJD(4) + t82;
t110 = t58 * t92 + t80 * t90;
t48 = t110 * pkin(11);
t79 = qJDD(4) + t81;
t19 = pkin(4) * t79 - t133 * t36 - t45 * t59 + t48 * t80 + t117;
t113 = -g(3) * t93 - t73 * t91;
t107 = -pkin(2) * t78 - pkin(9) * t121 + t76 * t120 + t113;
t106 = -pkin(3) * t52 - pkin(10) * t67 + t70 * t63 + t107;
t35 = -qJD(4) * t59 + t100 * t52 - t53 * t95;
t51 = pkin(4) * t80 - t133 * t59;
t21 = -pkin(4) * t35 - t134 * t36 - t48 * t58 + t51 * t59 + t106;
t112 = t19 * t92 + t21 * t90;
t111 = t35 * t92 + t79 * t90;
t128 = t100 * t28 + t95 * t26;
t20 = pkin(11) * t111 + t58 * t45 - t80 * t51 + t128;
t94 = sin(qJ(5));
t99 = cos(qJ(5));
t37 = t110 * t99 - t94 * t59;
t23 = t37 * qJD(5) + t111 * t94 + t99 * t36;
t38 = t110 * t94 + t99 * t59;
t29 = -mrSges(6,1) * t37 + mrSges(6,2) * t38;
t31 = -t35 * t90 + t79 * t92 + qJDD(5);
t49 = -t58 * t90 + t80 * t92 + qJD(5);
t32 = -mrSges(6,2) * t49 + mrSges(6,3) * t37;
t15 = m(6) * (t112 * t99 - t94 * t20) - t23 * mrSges(6,3) + t31 * mrSges(6,1) - t38 * t29 + t49 * t32;
t22 = -t38 * qJD(5) + t111 * t99 - t94 * t36;
t33 = mrSges(6,1) * t49 - mrSges(6,3) * t38;
t16 = m(6) * (t112 * t94 + t99 * t20) + t22 * mrSges(6,3) - t31 * mrSges(6,2) + t37 * t29 - t49 * t33;
t136 = t99 * t15 + t94 * t16;
t126 = t102 * t91;
t18 = m(6) * (-t19 * t90 + t21 * t92) + t23 * mrSges(6,2) - t22 * mrSges(6,1) + t38 * t33 - t37 * t32;
t54 = -mrSges(5,2) * t80 + mrSges(5,3) * t58;
t55 = mrSges(5,1) * t80 - mrSges(5,3) * t59;
t108 = m(5) * t106 - t35 * mrSges(5,1) + t36 * mrSges(5,2) + t136 * t90 + t92 * t18 - t58 * t54 + t59 * t55;
t61 = -mrSges(4,2) * t82 + mrSges(4,3) * t69;
t62 = mrSges(4,1) * t82 - mrSges(4,3) * t70;
t105 = m(4) * t107 - t52 * mrSges(4,1) + t53 * mrSges(4,2) - t69 * t61 + t70 * t62 + t108;
t71 = mrSges(3,1) * t84 - mrSges(3,3) * t120;
t119 = t102 * t124;
t72 = -mrSges(3,2) * t84 + mrSges(3,3) * t119;
t10 = (-t102 * t72 + t71 * t97) * t124 + m(3) * t113 + t77 * mrSges(3,2) - t78 * mrSges(3,1) + t105;
t46 = -mrSges(5,1) * t58 + mrSges(5,2) * t59;
t11 = m(5) * t117 + t79 * mrSges(5,1) - t36 * mrSges(5,3) + t136 * t92 - t90 * t18 - t59 * t46 + t80 * t54;
t12 = m(5) * t128 - mrSges(5,2) * t79 + mrSges(5,3) * t35 - t15 * t94 + t16 * t99 + t46 * t58 - t55 * t80;
t60 = -mrSges(4,1) * t69 + mrSges(4,2) * t70;
t7 = m(4) * t116 + mrSges(4,1) * t81 - mrSges(4,3) * t53 + t100 * t11 + t12 * t95 - t60 * t70 + t61 * t82;
t75 = (-mrSges(3,1) * t102 + mrSges(3,2) * t97) * t124;
t8 = m(4) * t127 - mrSges(4,2) * t81 + mrSges(4,3) * t52 + t100 * t12 - t11 * t95 + t60 * t69 - t62 * t82;
t5 = m(3) * (-g(3) * t126 + t115) - t77 * mrSges(3,3) + t83 * mrSges(3,1) - t75 * t120 + t84 * t72 + t96 * t8 + t101 * t7;
t6 = m(3) * t109 - t83 * mrSges(3,2) + t78 * mrSges(3,3) + t101 * t8 + t119 * t75 - t96 * t7 - t84 * t71;
t123 = t10 * t93 + t126 * t5 + t131 * t6;
t2 = m(2) * t114 - mrSges(2,1) * t104 - qJDD(1) * mrSges(2,2) + t102 * t6 - t5 * t97;
t1 = m(2) * t118 + qJDD(1) * mrSges(2,1) - t104 * mrSges(2,2) - t91 * t10 + (t102 * t5 + t6 * t97) * t93;
t3 = [-m(1) * g(1) - t1 * t98 + t103 * t2, t2, t6, t8, t12, t16; -m(1) * g(2) + t1 * t103 + t2 * t98, t1, t5, t7, t11, t15; (-m(1) - m(2)) * g(3) + t123, -m(2) * g(3) + t123, t10, t105, t108, t18;];
f_new = t3;
