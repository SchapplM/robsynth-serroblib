% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 06:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:57:01
% EndTime: 2019-05-08 05:57:20
% DurationCPUTime: 5.73s
% Computational Cost: add. (101666->211), mult. (216409->280), div. (0->0), fcn. (173029->12), ass. (0->108)
t108 = sin(qJ(3));
t113 = cos(qJ(3));
t105 = cos(pkin(6));
t100 = t105 * qJDD(1) + qJDD(2);
t104 = sin(pkin(6));
t109 = sin(qJ(2));
t114 = cos(qJ(2));
t137 = qJD(1) * t114;
t116 = qJD(1) ^ 2;
t110 = sin(qJ(1));
t115 = cos(qJ(1));
t127 = t110 * g(1) - t115 * g(2);
t149 = pkin(8) * t104;
t90 = qJDD(1) * pkin(1) + t116 * t149 + t127;
t141 = t105 * t90;
t123 = -t115 * g(1) - t110 * g(2);
t91 = -t116 * pkin(1) + qJDD(1) * t149 + t123;
t142 = t109 * t141 + t114 * t91;
t138 = qJD(1) * t104;
t93 = (-pkin(2) * t114 - pkin(9) * t109) * t138;
t101 = t105 * qJD(1) + qJD(2);
t99 = t101 ^ 2;
t55 = -t99 * pkin(2) + t100 * pkin(9) + (-g(3) * t109 + t93 * t137) * t104 + t142;
t148 = t105 * g(3);
t135 = qJD(1) * qJD(2);
t94 = (qJDD(1) * t109 + t114 * t135) * t104;
t95 = (-qJDD(1) * t114 + t109 * t135) * t104;
t56 = t95 * pkin(2) - t94 * pkin(9) - t148 + (-t90 + (pkin(2) * t109 - pkin(9) * t114) * t101 * qJD(1)) * t104;
t124 = -t108 * t55 + t113 * t56;
t129 = t109 * t138;
t82 = t113 * t101 - t108 * t129;
t83 = t108 * t101 + t113 * t129;
t71 = -t82 * pkin(3) - t83 * pkin(10);
t87 = qJDD(3) + t95;
t128 = t104 * t137;
t98 = qJD(3) - t128;
t96 = t98 ^ 2;
t32 = -t87 * pkin(3) - t96 * pkin(10) + t83 * t71 - t124;
t107 = sin(qJ(4));
t112 = cos(qJ(4));
t69 = t82 * qJD(3) + t108 * t100 + t113 * t94;
t74 = t107 * t98 + t112 * t83;
t43 = -t74 * qJD(4) - t107 * t69 + t112 * t87;
t81 = qJD(4) - t82;
t64 = t81 * pkin(4) - t74 * pkin(11);
t73 = -t107 * t83 + t112 * t98;
t72 = t73 ^ 2;
t119 = -t43 * pkin(4) - t72 * pkin(11) + t74 * t64 + t32;
t106 = sin(qJ(5));
t111 = cos(qJ(5));
t44 = t73 * qJD(4) + t107 * t87 + t112 * t69;
t59 = t106 * t73 + t111 * t74;
t29 = -t59 * qJD(5) - t106 * t44 + t111 * t43;
t58 = -t106 * t74 + t111 * t73;
t30 = t58 * qJD(5) + t106 * t43 + t111 * t44;
t79 = qJD(5) + t81;
t47 = t79 * pkin(5) - t59 * qJ(6);
t48 = t79 * mrSges(7,1) - t59 * mrSges(7,3);
t57 = t58 ^ 2;
t132 = m(7) * (-t29 * pkin(5) - t57 * qJ(6) + t59 * t47 + qJDD(6) + t119) + t30 * mrSges(7,2) + t59 * t48;
t45 = -t79 * mrSges(7,2) + t58 * mrSges(7,3);
t46 = -t79 * mrSges(6,2) + t58 * mrSges(6,3);
t49 = t79 * mrSges(6,1) - t59 * mrSges(6,3);
t152 = m(6) * t119 + t30 * mrSges(6,2) + t59 * t49 + t132 - (t46 + t45) * t58 - (mrSges(6,1) + mrSges(7,1)) * t29;
t62 = -t81 * mrSges(5,2) + t73 * mrSges(5,3);
t63 = t81 * mrSges(5,1) - t74 * mrSges(5,3);
t151 = m(5) * t32 - t43 * mrSges(5,1) + t44 * mrSges(5,2) - t73 * t62 + t74 * t63 + t152;
t143 = t108 * t56 + t113 * t55;
t33 = -t96 * pkin(3) + t87 * pkin(10) + t82 * t71 + t143;
t139 = t104 * t114;
t121 = -g(3) * t139 - t109 * t91 + t114 * t141;
t54 = -t100 * pkin(2) - t99 * pkin(9) + t93 * t129 - t121;
t68 = -t83 * qJD(3) + t113 * t100 - t108 * t94;
t36 = (-t82 * t98 - t69) * pkin(10) + (t83 * t98 - t68) * pkin(3) + t54;
t125 = -t107 * t33 + t112 * t36;
t67 = qJDD(4) - t68;
t21 = (t73 * t81 - t44) * pkin(11) + (t73 * t74 + t67) * pkin(4) + t125;
t145 = t107 * t36 + t112 * t33;
t23 = -t72 * pkin(4) + t43 * pkin(11) - t81 * t64 + t145;
t146 = t106 * t21 + t111 * t23;
t140 = t104 * t109;
t70 = -t82 * mrSges(4,1) + t83 * mrSges(4,2);
t75 = -t98 * mrSges(4,2) + t82 * mrSges(4,3);
t14 = m(4) * t124 + t87 * mrSges(4,1) - t69 * mrSges(4,3) - t83 * t70 + t98 * t75 - t151;
t88 = t101 * mrSges(3,1) - mrSges(3,3) * t129;
t126 = -t106 * t23 + t111 * t21;
t66 = qJDD(5) + t67;
t134 = m(7) * (-0.2e1 * qJD(6) * t59 + (t58 * t79 - t30) * qJ(6) + (t58 * t59 + t66) * pkin(5) + t126) + t79 * t45 + t66 * mrSges(7,1);
t40 = -t58 * mrSges(7,1) + t59 * mrSges(7,2);
t41 = -t58 * mrSges(6,1) + t59 * mrSges(6,2);
t12 = m(6) * t126 + t66 * mrSges(6,1) + t79 * t46 + (-t41 - t40) * t59 + (-mrSges(6,3) - mrSges(7,3)) * t30 + t134;
t133 = m(7) * (-t57 * pkin(5) + t29 * qJ(6) + 0.2e1 * qJD(6) * t58 - t79 * t47 + t146) + t29 * mrSges(7,3) + t58 * t40;
t13 = m(6) * t146 + t29 * mrSges(6,3) + t58 * t41 + (-t49 - t48) * t79 + (-mrSges(6,2) - mrSges(7,2)) * t66 + t133;
t60 = -t73 * mrSges(5,1) + t74 * mrSges(5,2);
t10 = m(5) * t125 + t67 * mrSges(5,1) - t44 * mrSges(5,3) + t106 * t13 + t111 * t12 - t74 * t60 + t81 * t62;
t11 = m(5) * t145 - t67 * mrSges(5,2) + t43 * mrSges(5,3) - t106 * t12 + t111 * t13 + t73 * t60 - t81 * t63;
t76 = t98 * mrSges(4,1) - t83 * mrSges(4,3);
t9 = m(4) * t143 - t87 * mrSges(4,2) + t68 * mrSges(4,3) - t107 * t10 + t112 * t11 + t82 * t70 - t98 * t76;
t92 = (-mrSges(3,1) * t114 + mrSges(3,2) * t109) * t138;
t4 = m(3) * (-g(3) * t140 + t142) - t95 * mrSges(3,3) - t100 * mrSges(3,2) + t92 * t128 - t101 * t88 + t113 * t9 - t108 * t14;
t89 = -t101 * mrSges(3,2) + mrSges(3,3) * t128;
t6 = m(3) * (-t104 * t90 - t148) + t94 * mrSges(3,2) + t95 * mrSges(3,1) + t108 * t9 + t113 * t14 + (t109 * t88 - t114 * t89) * t138;
t118 = m(4) * t54 - t68 * mrSges(4,1) + t69 * mrSges(4,2) + t112 * t10 + t107 * t11 - t82 * t75 + t83 * t76;
t8 = m(3) * t121 + t100 * mrSges(3,1) - t94 * mrSges(3,3) + t101 * t89 - t92 * t129 - t118;
t136 = t105 * t6 + t8 * t139 + t4 * t140;
t2 = m(2) * t123 - t116 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t109 * t8 + t114 * t4;
t1 = m(2) * t127 + qJDD(1) * mrSges(2,1) - t116 * mrSges(2,2) - t104 * t6 + (t109 * t4 + t114 * t8) * t105;
t3 = [-m(1) * g(1) - t110 * t1 + t115 * t2, t2, t4, t9, t11, t13, -t66 * mrSges(7,2) - t79 * t48 + t133; -m(1) * g(2) + t115 * t1 + t110 * t2, t1, t8, t14, t10, t12, -t30 * mrSges(7,3) - t59 * t40 + t134; (-m(1) - m(2)) * g(3) + t136, -m(2) * g(3) + t136, t6, t118, t151, t152, -t29 * mrSges(7,1) - t58 * t45 + t132;];
f_new  = t3;
