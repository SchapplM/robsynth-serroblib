% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 11:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:56:10
% EndTime: 2019-05-07 10:56:46
% DurationCPUTime: 8.15s
% Computational Cost: add. (151344->206), mult. (320471->270), div. (0->0), fcn. (235294->12), ass. (0->107)
t107 = sin(qJ(3));
t108 = sin(qJ(2));
t112 = cos(qJ(3));
t113 = cos(qJ(2));
t103 = sin(pkin(11));
t104 = cos(pkin(11));
t133 = qJD(1) * qJD(2);
t100 = t108 * t133;
t130 = t113 * t133;
t116 = qJD(1) ^ 2;
t109 = sin(qJ(1));
t114 = cos(qJ(1));
t129 = t109 * g(1) - t114 * g(2);
t84 = -qJDD(1) * pkin(1) - t116 * pkin(7) - t129;
t93 = qJDD(1) * t108 + t130;
t94 = qJDD(1) * t113 - t100;
t61 = (-t93 - t130) * pkin(8) + (-t94 + t100) * pkin(2) + t84;
t115 = qJD(2) ^ 2;
t125 = -g(1) * t114 - g(2) * t109;
t85 = -pkin(1) * t116 + qJDD(1) * pkin(7) + t125;
t131 = -g(3) * t108 + t113 * t85;
t134 = qJD(1) * t113;
t92 = (-pkin(2) * t113 - pkin(8) * t108) * qJD(1);
t64 = -pkin(2) * t115 + qJDD(2) * pkin(8) + t92 * t134 + t131;
t127 = -t107 * t64 + t112 * t61;
t105 = sin(qJ(6));
t110 = cos(qJ(6));
t106 = sin(qJ(5));
t111 = cos(qJ(5));
t135 = qJD(1) * t108;
t89 = qJD(2) * t112 - t107 * t135;
t72 = qJD(3) * t89 + qJDD(2) * t107 + t112 * t93;
t88 = qJDD(3) - t94;
t90 = qJD(2) * t107 + t112 * t135;
t99 = qJD(3) - t134;
t36 = (t89 * t99 - t72) * qJ(4) + (t89 * t90 + t88) * pkin(3) + t127;
t137 = t107 * t61 + t112 * t64;
t71 = -qJD(3) * t90 + qJDD(2) * t112 - t107 * t93;
t78 = pkin(3) * t99 - qJ(4) * t90;
t87 = t89 ^ 2;
t38 = -pkin(3) * t87 + qJ(4) * t71 - t78 * t99 + t137;
t75 = t103 * t89 + t104 * t90;
t126 = -0.2e1 * qJD(4) * t75 - t103 * t38 + t104 * t36;
t52 = t103 * t71 + t104 * t72;
t74 = -t103 * t90 + t104 * t89;
t21 = (t74 * t99 - t52) * pkin(9) + (t74 * t75 + t88) * pkin(4) + t126;
t132 = 0.2e1 * qJD(4) * t74 + t103 * t36 + t104 * t38;
t51 = -t103 * t72 + t104 * t71;
t67 = pkin(4) * t99 - pkin(9) * t75;
t73 = t74 ^ 2;
t26 = -pkin(4) * t73 + pkin(9) * t51 - t67 * t99 + t132;
t128 = -t106 * t26 + t111 * t21;
t56 = -t106 * t75 + t111 * t74;
t33 = qJD(5) * t56 + t106 * t51 + t111 * t52;
t57 = t106 * t74 + t111 * t75;
t86 = qJDD(5) + t88;
t98 = qJD(5) + t99;
t15 = (t56 * t98 - t33) * pkin(10) + (t56 * t57 + t86) * pkin(5) + t128;
t138 = t106 * t21 + t111 * t26;
t32 = -qJD(5) * t57 - t106 * t52 + t111 * t51;
t49 = pkin(5) * t98 - pkin(10) * t57;
t55 = t56 ^ 2;
t16 = -pkin(5) * t55 + pkin(10) * t32 - t49 * t98 + t138;
t42 = -t105 * t57 + t110 * t56;
t24 = qJD(6) * t42 + t105 * t32 + t110 * t33;
t43 = t105 * t56 + t110 * t57;
t30 = -mrSges(7,1) * t42 + mrSges(7,2) * t43;
t95 = qJD(6) + t98;
t39 = -mrSges(7,2) * t95 + mrSges(7,3) * t42;
t83 = qJDD(6) + t86;
t11 = m(7) * (-t105 * t16 + t110 * t15) - t24 * mrSges(7,3) + t83 * mrSges(7,1) - t43 * t30 + t95 * t39;
t23 = -qJD(6) * t43 - t105 * t33 + t110 * t32;
t40 = mrSges(7,1) * t95 - mrSges(7,3) * t43;
t12 = m(7) * (t105 * t15 + t110 * t16) + t23 * mrSges(7,3) - t83 * mrSges(7,2) + t42 * t30 - t95 * t40;
t44 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t48 = mrSges(6,1) * t98 - mrSges(6,3) * t57;
t10 = m(6) * t138 - t86 * mrSges(6,2) + t32 * mrSges(6,3) - t105 * t11 + t110 * t12 + t56 * t44 - t98 * t48;
t58 = -mrSges(5,1) * t74 + mrSges(5,2) * t75;
t65 = -mrSges(5,2) * t99 + mrSges(5,3) * t74;
t47 = -mrSges(6,2) * t98 + mrSges(6,3) * t56;
t9 = m(6) * t128 + t86 * mrSges(6,1) - t33 * mrSges(6,3) + t105 * t12 + t110 * t11 - t57 * t44 + t98 * t47;
t7 = m(5) * t126 + t88 * mrSges(5,1) - t52 * mrSges(5,3) + t106 * t10 + t111 * t9 - t75 * t58 + t99 * t65;
t76 = -mrSges(4,1) * t89 + mrSges(4,2) * t90;
t77 = -mrSges(4,2) * t99 + mrSges(4,3) * t89;
t66 = mrSges(5,1) * t99 - mrSges(5,3) * t75;
t8 = m(5) * t132 - t88 * mrSges(5,2) + t51 * mrSges(5,3) + t111 * t10 - t106 * t9 + t74 * t58 - t99 * t66;
t5 = m(4) * t127 + t88 * mrSges(4,1) - t72 * mrSges(4,3) + t103 * t8 + t104 * t7 - t90 * t76 + t99 * t77;
t79 = mrSges(4,1) * t99 - mrSges(4,3) * t90;
t6 = m(4) * t137 - t88 * mrSges(4,2) + t71 * mrSges(4,3) - t103 * t7 + t104 * t8 + t89 * t76 - t99 * t79;
t96 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t135;
t97 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t134;
t140 = m(3) * t84 - t94 * mrSges(3,1) + t93 * mrSges(3,2) + t107 * t6 + t112 * t5 + (t108 * t96 - t113 * t97) * qJD(1);
t136 = -t113 * g(3) - t108 * t85;
t63 = -qJDD(2) * pkin(2) - pkin(8) * t115 + t92 * t135 - t136;
t121 = -pkin(3) * t71 - qJ(4) * t87 + t90 * t78 + qJDD(4) + t63;
t119 = -pkin(4) * t51 - pkin(9) * t73 + t75 * t67 + t121;
t123 = t23 * mrSges(7,1) + t42 * t39 - m(7) * (-pkin(5) * t32 - pkin(10) * t55 + t49 * t57 + t119) - t24 * mrSges(7,2) - t43 * t40;
t120 = -m(6) * t119 + t32 * mrSges(6,1) - t33 * mrSges(6,2) + t56 * t47 - t57 * t48 + t123;
t118 = -m(5) * t121 + t51 * mrSges(5,1) - t52 * mrSges(5,2) + t74 * t65 - t75 * t66 + t120;
t117 = m(4) * t63 - t71 * mrSges(4,1) + t72 * mrSges(4,2) - t89 * t77 + t90 * t79 - t118;
t91 = (-mrSges(3,1) * t113 + mrSges(3,2) * t108) * qJD(1);
t14 = m(3) * t136 + qJDD(2) * mrSges(3,1) - t93 * mrSges(3,3) + qJD(2) * t97 - t91 * t135 - t117;
t4 = m(3) * t131 - qJDD(2) * mrSges(3,2) + t94 * mrSges(3,3) - qJD(2) * t96 - t107 * t5 + t112 * t6 + t91 * t134;
t139 = t108 * t4 + t113 * t14;
t2 = m(2) * t129 + qJDD(1) * mrSges(2,1) - t116 * mrSges(2,2) - t140;
t1 = m(2) * t125 - t116 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t108 * t14 + t113 * t4;
t3 = [-m(1) * g(1) + t1 * t114 - t109 * t2, t1, t4, t6, t8, t10, t12; -m(1) * g(2) + t1 * t109 + t114 * t2, t2, t14, t5, t7, t9, t11; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t140, t117, -t118, -t120, -t123;];
f_new  = t3;
