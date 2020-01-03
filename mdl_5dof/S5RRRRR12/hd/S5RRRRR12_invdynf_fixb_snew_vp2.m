% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR12
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR12_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:47:05
% EndTime: 2019-12-31 22:47:18
% DurationCPUTime: 6.28s
% Computational Cost: add. (106815->187), mult. (265459->273), div. (0->0), fcn. (219426->14), ass. (0->110)
t81 = sin(pkin(5));
t116 = qJD(1) * t81;
t82 = cos(pkin(6));
t92 = cos(qJ(2));
t119 = t82 * t92;
t80 = sin(pkin(6));
t86 = sin(qJ(3));
t122 = t80 * t86;
t83 = cos(pkin(5));
t78 = t83 * qJD(1) + qJD(2);
t87 = sin(qJ(2));
t91 = cos(qJ(3));
t60 = t78 * t122 + (t86 * t119 + t87 * t91) * t116;
t112 = qJD(1) * qJD(2);
t75 = (qJDD(1) * t87 + t92 * t112) * t81;
t76 = (qJDD(1) * t92 - t87 * t112) * t81;
t77 = t83 * qJDD(1) + qJDD(2);
t99 = t76 * t82 + t77 * t80;
t48 = -t60 * qJD(3) - t86 * t75 + t99 * t91;
t88 = sin(qJ(1));
t93 = cos(qJ(1));
t106 = t88 * g(1) - t93 * g(2);
t129 = pkin(8) * t81;
t94 = qJD(1) ^ 2;
t72 = qJDD(1) * pkin(1) + t94 * t129 + t106;
t124 = t72 * t83;
t103 = -t93 * g(1) - t88 * g(2);
t73 = -t94 * pkin(1) + qJDD(1) * t129 + t103;
t104 = t92 * t124 - t87 * t73;
t115 = qJD(1) * t87;
t127 = pkin(9) * t82;
t114 = qJD(1) * t92;
t108 = t81 * t114;
t123 = t78 * t80;
t65 = (t82 * t108 + t123) * pkin(9);
t128 = pkin(9) * t80;
t69 = (-pkin(2) * t92 - t87 * t128) * t116;
t35 = -t75 * t127 + t77 * pkin(2) + t78 * t65 + (-g(3) * t92 - t69 * t115) * t81 + t104;
t125 = t35 * t82;
t117 = t87 * t124 + t92 * t73;
t109 = t81 * t115;
t68 = t78 * pkin(2) - t109 * t127;
t36 = -t78 * t68 + (-g(3) * t87 + t69 * t114) * t81 + t99 * pkin(9) + t117;
t126 = t83 * g(3);
t41 = -t75 * t128 - t76 * pkin(2) - t126 + (-t72 + (-t65 * t92 + t68 * t87) * qJD(1)) * t81;
t130 = (t41 * t80 + t125) * t91 - t86 * t36;
t59 = (t91 * t119 - t86 * t87) * t116 + t91 * t123;
t121 = t81 * t87;
t120 = t81 * t92;
t110 = t41 * t122 + t86 * t125 + t91 * t36;
t51 = -t59 * pkin(3) - t60 * pkin(10);
t61 = -t80 * t76 + t82 * t77 + qJDD(3);
t66 = -t80 * t108 + t82 * t78 + qJD(3);
t64 = t66 ^ 2;
t21 = -t64 * pkin(3) + t61 * pkin(10) + t59 * t51 + t110;
t105 = -t80 * t35 + t82 * t41;
t49 = t59 * qJD(3) + t91 * t75 + t99 * t86;
t23 = (-t59 * t66 - t49) * pkin(10) + (t60 * t66 - t48) * pkin(3) + t105;
t85 = sin(qJ(4));
t90 = cos(qJ(4));
t118 = t90 * t21 + t85 * t23;
t53 = -t85 * t60 + t90 * t66;
t54 = t90 * t60 + t85 * t66;
t38 = -t53 * pkin(4) - t54 * pkin(11);
t47 = qJDD(4) - t48;
t58 = qJD(4) - t59;
t57 = t58 ^ 2;
t17 = -t57 * pkin(4) + t47 * pkin(11) + t53 * t38 + t118;
t20 = -t61 * pkin(3) - t64 * pkin(10) + t60 * t51 - t130;
t28 = -t54 * qJD(4) - t85 * t49 + t90 * t61;
t29 = t53 * qJD(4) + t90 * t49 + t85 * t61;
t18 = (-t53 * t58 - t29) * pkin(11) + (t54 * t58 - t28) * pkin(4) + t20;
t84 = sin(qJ(5));
t89 = cos(qJ(5));
t43 = -t84 * t54 + t89 * t58;
t25 = t43 * qJD(5) + t89 * t29 + t84 * t47;
t44 = t89 * t54 + t84 * t58;
t26 = -t43 * mrSges(6,1) + t44 * mrSges(6,2);
t27 = qJDD(5) - t28;
t52 = qJD(5) - t53;
t30 = -t52 * mrSges(6,2) + t43 * mrSges(6,3);
t14 = m(6) * (-t84 * t17 + t89 * t18) - t25 * mrSges(6,3) + t27 * mrSges(6,1) - t44 * t26 + t52 * t30;
t24 = -t44 * qJD(5) - t84 * t29 + t89 * t47;
t31 = t52 * mrSges(6,1) - t44 * mrSges(6,3);
t15 = m(6) * (t89 * t17 + t84 * t18) + t24 * mrSges(6,3) - t27 * mrSges(6,2) + t43 * t26 - t52 * t31;
t37 = -t53 * mrSges(5,1) + t54 * mrSges(5,2);
t46 = t58 * mrSges(5,1) - t54 * mrSges(5,3);
t12 = m(5) * t118 - t47 * mrSges(5,2) + t28 * mrSges(5,3) - t84 * t14 + t89 * t15 + t53 * t37 - t58 * t46;
t101 = -t85 * t21 + t90 * t23;
t45 = -t58 * mrSges(5,2) + t53 * mrSges(5,3);
t96 = m(6) * (-t47 * pkin(4) - t57 * pkin(11) + t54 * t38 - t101) - t24 * mrSges(6,1) + t25 * mrSges(6,2) - t43 * t30 + t44 * t31;
t13 = m(5) * t101 + t47 * mrSges(5,1) - t29 * mrSges(5,3) - t54 * t37 + t58 * t45 - t96;
t55 = -t66 * mrSges(4,2) + t59 * mrSges(4,3);
t56 = t66 * mrSges(4,1) - t60 * mrSges(4,3);
t10 = m(4) * t105 - t48 * mrSges(4,1) + t49 * mrSges(4,2) + t85 * t12 + t90 * t13 - t59 * t55 + t60 * t56;
t50 = -t59 * mrSges(4,1) + t60 * mrSges(4,2);
t95 = m(5) * t20 - t28 * mrSges(5,1) + t29 * mrSges(5,2) + t89 * t14 + t84 * t15 - t53 * t45 + t54 * t46;
t11 = m(4) * t130 + t61 * mrSges(4,1) - t49 * mrSges(4,3) - t60 * t50 + t66 * t55 - t95;
t9 = m(4) * t110 - t61 * mrSges(4,2) + t48 * mrSges(4,3) + t90 * t12 - t85 * t13 + t59 * t50 - t66 * t56;
t102 = t91 * t11 + t86 * t9;
t107 = (-mrSges(3,1) * t92 + mrSges(3,2) * t87) * t116 ^ 2;
t71 = -t78 * mrSges(3,2) + mrSges(3,3) * t108;
t4 = m(3) * (-g(3) * t120 + t104) - t75 * mrSges(3,3) + t77 * mrSges(3,1) - t87 * t107 + t78 * t71 - t80 * t10 + t102 * t82;
t70 = t78 * mrSges(3,1) - mrSges(3,3) * t109;
t6 = m(3) * (-t81 * t72 - t126) + t75 * mrSges(3,2) - t76 * mrSges(3,1) + t82 * t10 + t102 * t80 + (t70 * t87 - t71 * t92) * t116;
t8 = m(3) * (-g(3) * t121 + t117) + t76 * mrSges(3,3) - t77 * mrSges(3,2) + t92 * t107 - t78 * t70 + t91 * t9 - t86 * t11;
t113 = t4 * t120 + t8 * t121 + t83 * t6;
t2 = m(2) * t103 - t94 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t87 * t4 + t92 * t8;
t1 = m(2) * t106 + qJDD(1) * mrSges(2,1) - t94 * mrSges(2,2) - t81 * t6 + (t92 * t4 + t87 * t8) * t83;
t3 = [-m(1) * g(1) - t88 * t1 + t93 * t2, t2, t8, t9, t12, t15; -m(1) * g(2) + t93 * t1 + t88 * t2, t1, t4, t11, t13, t14; (-m(1) - m(2)) * g(3) + t113, -m(2) * g(3) + t113, t6, t10, t95, t96;];
f_new = t3;
