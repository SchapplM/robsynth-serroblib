% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR14
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR14_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR14_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:17:16
% EndTime: 2019-12-31 19:17:26
% DurationCPUTime: 5.50s
% Computational Cost: add. (79257->178), mult. (247111->268), div. (0->0), fcn. (206356->14), ass. (0->107)
t81 = cos(pkin(11));
t82 = cos(pkin(6));
t124 = t81 * t82;
t79 = sin(pkin(6));
t83 = cos(pkin(5));
t127 = t79 * t83;
t78 = sin(pkin(11));
t80 = sin(pkin(5));
t86 = sin(qJ(3));
t90 = cos(qJ(3));
t95 = (t124 * t90 - t78 * t86) * t80 + t90 * t127;
t59 = t95 * qJD(1);
t123 = t82 * t86;
t126 = t79 * t86;
t97 = t83 * t126 + (t123 * t81 + t78 * t90) * t80;
t60 = t97 * qJD(1);
t50 = -t60 * qJD(3) + qJDD(1) * t95;
t119 = qJD(1) * t80;
t112 = qJD(2) * t119;
t125 = t80 * t81;
t87 = sin(qJ(1));
t91 = cos(qJ(1));
t113 = t87 * g(1) - t91 * g(2);
t121 = qJ(2) * t80;
t92 = qJD(1) ^ 2;
t70 = qJDD(1) * pkin(1) + t121 * t92 + t113;
t129 = t70 * t83;
t104 = -g(3) * t125 - 0.2e1 * t78 * t112 + t81 * t129;
t67 = (t124 * t80 + t127) * qJD(1) * pkin(8);
t109 = -t91 * g(1) - t87 * g(2);
t71 = -t92 * pkin(1) + qJDD(1) * t121 + t109;
t103 = -pkin(8) * t78 * t79 - pkin(2) * t81;
t120 = pkin(8) * qJDD(1);
t98 = qJD(1) * t103 * t119 + t120 * t82;
t35 = (pkin(2) * qJDD(1) + qJD(1) * t67) * t83 + (-t80 * t98 - t71) * t78 + t104;
t114 = t78 * t129 + (0.2e1 * t112 + t71) * t81;
t128 = t78 * t80;
t72 = (-pkin(8) * t128 * t82 + pkin(2) * t83) * qJD(1);
t36 = (-qJD(1) * t72 + t120 * t79) * t83 + (-g(3) * t78 + t81 * t98) * t80 + t114;
t110 = -t83 * g(3) + qJDD(2);
t44 = (-t70 + t103 * qJDD(1) + (-t67 * t81 + t72 * t78) * qJD(1)) * t80 + t110;
t130 = (t35 * t82 + t44 * t79) * t90 - t86 * t36;
t115 = t35 * t123 + t44 * t126 + t90 * t36;
t49 = -t59 * pkin(3) - t60 * pkin(9);
t99 = -t125 * t79 + t82 * t83;
t68 = qJD(1) * t99 + qJD(3);
t64 = t68 ^ 2;
t65 = qJDD(1) * t99 + qJDD(3);
t21 = -t64 * pkin(3) + t65 * pkin(9) + t59 * t49 + t115;
t111 = -t79 * t35 + t82 * t44;
t51 = t59 * qJD(3) + qJDD(1) * t97;
t23 = (-t59 * t68 - t51) * pkin(9) + (t60 * t68 - t50) * pkin(3) + t111;
t85 = sin(qJ(4));
t89 = cos(qJ(4));
t122 = t89 * t21 + t85 * t23;
t53 = -t85 * t60 + t89 * t68;
t54 = t89 * t60 + t85 * t68;
t38 = -t53 * pkin(4) - t54 * pkin(10);
t47 = qJDD(4) - t50;
t58 = qJD(4) - t59;
t57 = t58 ^ 2;
t17 = -t57 * pkin(4) + t47 * pkin(10) + t53 * t38 + t122;
t20 = -t65 * pkin(3) - t64 * pkin(9) + t60 * t49 - t130;
t30 = -t54 * qJD(4) - t85 * t51 + t89 * t65;
t31 = t53 * qJD(4) + t89 * t51 + t85 * t65;
t18 = (-t53 * t58 - t31) * pkin(10) + (t54 * t58 - t30) * pkin(4) + t20;
t84 = sin(qJ(5));
t88 = cos(qJ(5));
t42 = -t84 * t54 + t88 * t58;
t25 = t42 * qJD(5) + t88 * t31 + t84 * t47;
t43 = t88 * t54 + t84 * t58;
t26 = -t42 * mrSges(6,1) + t43 * mrSges(6,2);
t52 = qJD(5) - t53;
t27 = -t52 * mrSges(6,2) + t42 * mrSges(6,3);
t29 = qJDD(5) - t30;
t14 = m(6) * (-t84 * t17 + t88 * t18) - t25 * mrSges(6,3) + t29 * mrSges(6,1) - t43 * t26 + t52 * t27;
t24 = -t43 * qJD(5) - t84 * t31 + t88 * t47;
t28 = t52 * mrSges(6,1) - t43 * mrSges(6,3);
t15 = m(6) * (t88 * t17 + t84 * t18) + t24 * mrSges(6,3) - t29 * mrSges(6,2) + t42 * t26 - t52 * t28;
t37 = -t53 * mrSges(5,1) + t54 * mrSges(5,2);
t46 = t58 * mrSges(5,1) - t54 * mrSges(5,3);
t12 = m(5) * t122 - t47 * mrSges(5,2) + t30 * mrSges(5,3) - t84 * t14 + t88 * t15 + t53 * t37 - t58 * t46;
t106 = -t85 * t21 + t89 * t23;
t45 = -t58 * mrSges(5,2) + t53 * mrSges(5,3);
t94 = m(6) * (-t47 * pkin(4) - t57 * pkin(10) + t54 * t38 - t106) - t24 * mrSges(6,1) + t25 * mrSges(6,2) - t42 * t27 + t43 * t28;
t13 = m(5) * t106 + t47 * mrSges(5,1) - t31 * mrSges(5,3) - t54 * t37 + t58 * t45 - t94;
t55 = -t68 * mrSges(4,2) + t59 * mrSges(4,3);
t56 = t68 * mrSges(4,1) - t60 * mrSges(4,3);
t10 = m(4) * t111 - t50 * mrSges(4,1) + t51 * mrSges(4,2) + t85 * t12 + t89 * t13 - t59 * t55 + t60 * t56;
t102 = t83 * mrSges(3,1) - mrSges(3,3) * t128;
t48 = -t59 * mrSges(4,1) + t60 * mrSges(4,2);
t93 = m(5) * t20 - t30 * mrSges(5,1) + t31 * mrSges(5,2) + t88 * t14 + t84 * t15 - t53 * t45 + t54 * t46;
t11 = m(4) * t130 + t65 * mrSges(4,1) - t51 * mrSges(4,3) - t60 * t48 + t68 * t55 - t93;
t9 = m(4) * t115 - t65 * mrSges(4,2) + t50 * mrSges(4,3) + t89 * t12 - t85 * t13 + t59 * t48 - t68 * t56;
t108 = t90 * t11 + t86 * t9;
t107 = -mrSges(3,1) * t81 + mrSges(3,2) * t78;
t69 = t107 * t119;
t101 = -t83 * mrSges(3,2) + mrSges(3,3) * t125;
t74 = t101 * qJD(1);
t4 = m(3) * (-t78 * t71 + t104) - t79 * t10 + t108 * t82 + t102 * qJDD(1) + (-t128 * t69 + t83 * t74) * qJD(1);
t73 = t102 * qJD(1);
t6 = m(3) * t110 + t82 * t10 + t108 * t79 + (-m(3) * t70 + t107 * qJDD(1) + (t73 * t78 - t74 * t81) * qJD(1)) * t80;
t8 = m(3) * (-g(3) * t128 + t114) + t90 * t9 - t86 * t11 + t101 * qJDD(1) + (t125 * t69 - t83 * t73) * qJD(1);
t117 = t4 * t125 + t8 * t128 + t83 * t6;
t2 = m(2) * t109 - t92 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t78 * t4 + t81 * t8;
t1 = m(2) * t113 + qJDD(1) * mrSges(2,1) - t92 * mrSges(2,2) - t80 * t6 + (t81 * t4 + t78 * t8) * t83;
t3 = [-m(1) * g(1) - t87 * t1 + t91 * t2, t2, t8, t9, t12, t15; -m(1) * g(2) + t91 * t1 + t87 * t2, t1, t4, t11, t13, t14; (-m(1) - m(2)) * g(3) + t117, -m(2) * g(3) + t117, t6, t10, t93, t94;];
f_new = t3;
