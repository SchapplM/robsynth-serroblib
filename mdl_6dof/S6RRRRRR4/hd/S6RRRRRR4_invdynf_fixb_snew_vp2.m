% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 09:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 09:21:50
% EndTime: 2019-05-08 09:22:16
% DurationCPUTime: 8.52s
% Computational Cost: add. (162562->207), mult. (336388->268), div. (0->0), fcn. (250326->12), ass. (0->109)
t106 = sin(qJ(3));
t107 = sin(qJ(2));
t112 = cos(qJ(3));
t113 = cos(qJ(2));
t105 = sin(qJ(4));
t111 = cos(qJ(4));
t132 = qJD(1) * qJD(2);
t100 = t107 * t132;
t130 = t113 * t132;
t116 = qJD(1) ^ 2;
t108 = sin(qJ(1));
t114 = cos(qJ(1));
t129 = g(1) * t108 - t114 * g(2);
t83 = -qJDD(1) * pkin(1) - pkin(7) * t116 - t129;
t92 = qJDD(1) * t107 + t130;
t93 = qJDD(1) * t113 - t100;
t61 = (-t92 - t130) * pkin(8) + (-t93 + t100) * pkin(2) + t83;
t115 = qJD(2) ^ 2;
t125 = -g(1) * t114 - g(2) * t108;
t84 = -pkin(1) * t116 + qJDD(1) * pkin(7) + t125;
t131 = -g(3) * t107 + t113 * t84;
t133 = qJD(1) * t113;
t91 = (-pkin(2) * t113 - pkin(8) * t107) * qJD(1);
t64 = -pkin(2) * t115 + qJDD(2) * pkin(8) + t91 * t133 + t131;
t126 = -t106 * t64 + t112 * t61;
t103 = sin(qJ(6));
t109 = cos(qJ(6));
t104 = sin(qJ(5));
t110 = cos(qJ(5));
t134 = qJD(1) * t107;
t88 = qJD(2) * t112 - t106 * t134;
t70 = qJD(3) * t88 + qJDD(2) * t106 + t112 * t92;
t87 = qJDD(3) - t93;
t89 = qJD(2) * t106 + t112 * t134;
t99 = qJD(3) - t133;
t36 = (t88 * t99 - t70) * pkin(9) + (t88 * t89 + t87) * pkin(3) + t126;
t136 = t106 * t61 + t112 * t64;
t69 = -qJD(3) * t89 + qJDD(2) * t112 - t106 * t92;
t77 = pkin(3) * t99 - pkin(9) * t89;
t86 = t88 ^ 2;
t38 = -pkin(3) * t86 + pkin(9) * t69 - t77 * t99 + t136;
t127 = -t105 * t38 + t111 * t36;
t72 = -t105 * t89 + t111 * t88;
t49 = qJD(4) * t72 + t105 * t69 + t111 * t70;
t73 = t105 * t88 + t111 * t89;
t85 = qJDD(4) + t87;
t98 = qJD(4) + t99;
t24 = (t72 * t98 - t49) * pkin(10) + (t72 * t73 + t85) * pkin(4) + t127;
t137 = t105 * t36 + t111 * t38;
t48 = -qJD(4) * t73 - t105 * t70 + t111 * t69;
t67 = pkin(4) * t98 - pkin(10) * t73;
t71 = t72 ^ 2;
t26 = -pkin(4) * t71 + pkin(10) * t48 - t67 * t98 + t137;
t128 = -t104 * t26 + t110 * t24;
t56 = -t104 * t73 + t110 * t72;
t33 = qJD(5) * t56 + t104 * t48 + t110 * t49;
t57 = t104 * t72 + t110 * t73;
t82 = qJDD(5) + t85;
t95 = qJD(5) + t98;
t15 = (t56 * t95 - t33) * pkin(11) + (t56 * t57 + t82) * pkin(5) + t128;
t138 = t104 * t24 + t110 * t26;
t32 = -qJD(5) * t57 - t104 * t49 + t110 * t48;
t52 = pkin(5) * t95 - pkin(11) * t57;
t55 = t56 ^ 2;
t16 = -pkin(5) * t55 + pkin(11) * t32 - t52 * t95 + t138;
t42 = -t103 * t57 + t109 * t56;
t21 = qJD(6) * t42 + t103 * t32 + t109 * t33;
t43 = t103 * t56 + t109 * t57;
t30 = -mrSges(7,1) * t42 + mrSges(7,2) * t43;
t94 = qJD(6) + t95;
t39 = -mrSges(7,2) * t94 + mrSges(7,3) * t42;
t80 = qJDD(6) + t82;
t13 = m(7) * (-t103 * t16 + t109 * t15) - t21 * mrSges(7,3) + t80 * mrSges(7,1) - t43 * t30 + t94 * t39;
t20 = -qJD(6) * t43 - t103 * t33 + t109 * t32;
t40 = mrSges(7,1) * t94 - mrSges(7,3) * t43;
t14 = m(7) * (t103 * t15 + t109 * t16) + t20 * mrSges(7,3) - t80 * mrSges(7,2) + t42 * t30 - t94 * t40;
t44 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t51 = mrSges(6,1) * t95 - mrSges(6,3) * t57;
t10 = m(6) * t138 - t82 * mrSges(6,2) + t32 * mrSges(6,3) - t103 * t13 + t109 * t14 + t56 * t44 - t95 * t51;
t58 = -mrSges(5,1) * t72 + mrSges(5,2) * t73;
t65 = -mrSges(5,2) * t98 + mrSges(5,3) * t72;
t50 = -mrSges(6,2) * t95 + mrSges(6,3) * t56;
t9 = m(6) * t128 + t82 * mrSges(6,1) - t33 * mrSges(6,3) + t103 * t14 + t109 * t13 - t57 * t44 + t95 * t50;
t7 = m(5) * t127 + t85 * mrSges(5,1) - t49 * mrSges(5,3) + t104 * t10 + t110 * t9 - t73 * t58 + t98 * t65;
t74 = -mrSges(4,1) * t88 + mrSges(4,2) * t89;
t75 = -mrSges(4,2) * t99 + mrSges(4,3) * t88;
t66 = mrSges(5,1) * t98 - mrSges(5,3) * t73;
t8 = m(5) * t137 - t85 * mrSges(5,2) + t48 * mrSges(5,3) + t110 * t10 - t104 * t9 + t72 * t58 - t98 * t66;
t5 = m(4) * t126 + t87 * mrSges(4,1) - t70 * mrSges(4,3) + t105 * t8 + t111 * t7 - t89 * t74 + t99 * t75;
t76 = mrSges(4,1) * t99 - mrSges(4,3) * t89;
t6 = m(4) * t136 - t87 * mrSges(4,2) + t69 * mrSges(4,3) - t105 * t7 + t111 * t8 + t88 * t74 - t99 * t76;
t96 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t134;
t97 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133;
t140 = m(3) * t83 - t93 * mrSges(3,1) + t92 * mrSges(3,2) + t106 * t6 + t112 * t5 + (t107 * t96 - t113 * t97) * qJD(1);
t135 = -t113 * g(3) - t107 * t84;
t63 = -qJDD(2) * pkin(2) - pkin(8) * t115 + t91 * t134 - t135;
t121 = -pkin(3) * t69 - pkin(9) * t86 + t89 * t77 + t63;
t119 = -pkin(4) * t48 - pkin(10) * t71 + t73 * t67 + t121;
t123 = t20 * mrSges(7,1) + t42 * t39 - m(7) * (-pkin(5) * t32 - pkin(11) * t55 + t52 * t57 + t119) - t21 * mrSges(7,2) - t43 * t40;
t120 = -m(6) * t119 + t32 * mrSges(6,1) - t33 * mrSges(6,2) + t56 * t50 - t57 * t51 + t123;
t118 = -m(5) * t121 + t48 * mrSges(5,1) - t49 * mrSges(5,2) + t72 * t65 - t73 * t66 + t120;
t117 = m(4) * t63 - t69 * mrSges(4,1) + t70 * mrSges(4,2) - t88 * t75 + t89 * t76 - t118;
t90 = (-mrSges(3,1) * t113 + mrSges(3,2) * t107) * qJD(1);
t12 = m(3) * t135 + qJDD(2) * mrSges(3,1) - t92 * mrSges(3,3) + qJD(2) * t97 - t90 * t134 - t117;
t4 = m(3) * t131 - qJDD(2) * mrSges(3,2) + t93 * mrSges(3,3) - qJD(2) * t96 - t106 * t5 + t112 * t6 + t90 * t133;
t139 = t107 * t4 + t113 * t12;
t2 = m(2) * t129 + qJDD(1) * mrSges(2,1) - t116 * mrSges(2,2) - t140;
t1 = m(2) * t125 - t116 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t107 * t12 + t113 * t4;
t3 = [-m(1) * g(1) + t1 * t114 - t108 * t2, t1, t4, t6, t8, t10, t14; -m(1) * g(2) + t1 * t108 + t114 * t2, t2, t12, t5, t7, t9, t13; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t140, t117, -t118, -t120, -t123;];
f_new  = t3;
