% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 20:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:45:02
% EndTime: 2019-05-07 20:45:25
% DurationCPUTime: 8.21s
% Computational Cost: add. (157337->206), mult. (328350->270), div. (0->0), fcn. (242262->12), ass. (0->107)
t107 = sin(qJ(3));
t108 = sin(qJ(2));
t112 = cos(qJ(3));
t113 = cos(qJ(2));
t106 = sin(qJ(4));
t111 = cos(qJ(4));
t133 = qJD(1) * qJD(2);
t100 = t108 * t133;
t130 = t113 * t133;
t116 = qJD(1) ^ 2;
t109 = sin(qJ(1));
t114 = cos(qJ(1));
t129 = t109 * g(1) - t114 * g(2);
t84 = -qJDD(1) * pkin(1) - t116 * pkin(7) - t129;
t93 = t108 * qJDD(1) + t130;
t94 = t113 * qJDD(1) - t100;
t63 = (-t93 - t130) * pkin(8) + (-t94 + t100) * pkin(2) + t84;
t115 = qJD(2) ^ 2;
t125 = -t114 * g(1) - t109 * g(2);
t85 = -t116 * pkin(1) + qJDD(1) * pkin(7) + t125;
t131 = -t108 * g(3) + t113 * t85;
t134 = t113 * qJD(1);
t92 = (-pkin(2) * t113 - pkin(8) * t108) * qJD(1);
t66 = -t115 * pkin(2) + qJDD(2) * pkin(8) + t92 * t134 + t131;
t127 = -t107 * t66 + t112 * t63;
t105 = sin(qJ(6));
t110 = cos(qJ(6));
t103 = sin(pkin(11));
t104 = cos(pkin(11));
t135 = qJD(1) * t108;
t89 = t112 * qJD(2) - t107 * t135;
t72 = t89 * qJD(3) + t107 * qJDD(2) + t112 * t93;
t88 = qJDD(3) - t94;
t90 = t107 * qJD(2) + t112 * t135;
t99 = qJD(3) - t134;
t36 = (t89 * t99 - t72) * pkin(9) + (t89 * t90 + t88) * pkin(3) + t127;
t137 = t107 * t63 + t112 * t66;
t71 = -t90 * qJD(3) + t112 * qJDD(2) - t107 * t93;
t79 = t99 * pkin(3) - t90 * pkin(9);
t87 = t89 ^ 2;
t38 = -t87 * pkin(3) + t71 * pkin(9) - t99 * t79 + t137;
t128 = -t106 * t38 + t111 * t36;
t74 = -t106 * t90 + t111 * t89;
t49 = t74 * qJD(4) + t106 * t71 + t111 * t72;
t75 = t106 * t89 + t111 * t90;
t86 = qJDD(4) + t88;
t98 = qJD(4) + t99;
t21 = (t74 * t98 - t49) * qJ(5) + (t74 * t75 + t86) * pkin(4) + t128;
t138 = t106 * t36 + t111 * t38;
t48 = -t75 * qJD(4) - t106 * t72 + t111 * t71;
t68 = t98 * pkin(4) - t75 * qJ(5);
t73 = t74 ^ 2;
t26 = -t73 * pkin(4) + t48 * qJ(5) - t98 * t68 + t138;
t59 = t103 * t74 + t104 * t75;
t126 = -0.2e1 * qJD(5) * t59 - t103 * t26 + t104 * t21;
t33 = t103 * t48 + t104 * t49;
t58 = -t103 * t75 + t104 * t74;
t15 = (t58 * t98 - t33) * pkin(10) + (t58 * t59 + t86) * pkin(5) + t126;
t132 = 0.2e1 * qJD(5) * t58 + t103 * t21 + t104 * t26;
t32 = -t103 * t49 + t104 * t48;
t52 = t98 * pkin(5) - t59 * pkin(10);
t57 = t58 ^ 2;
t16 = -t57 * pkin(5) + t32 * pkin(10) - t98 * t52 + t132;
t42 = -t105 * t59 + t110 * t58;
t24 = t42 * qJD(6) + t105 * t32 + t110 * t33;
t43 = t105 * t58 + t110 * t59;
t30 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t95 = qJD(6) + t98;
t39 = -t95 * mrSges(7,2) + t42 * mrSges(7,3);
t83 = qJDD(6) + t86;
t11 = m(7) * (-t105 * t16 + t110 * t15) - t24 * mrSges(7,3) + t83 * mrSges(7,1) - t43 * t30 + t95 * t39;
t23 = -t43 * qJD(6) - t105 * t33 + t110 * t32;
t40 = t95 * mrSges(7,1) - t43 * mrSges(7,3);
t12 = m(7) * (t105 * t15 + t110 * t16) + t23 * mrSges(7,3) - t83 * mrSges(7,2) + t42 * t30 - t95 * t40;
t44 = -t58 * mrSges(6,1) + t59 * mrSges(6,2);
t51 = t98 * mrSges(6,1) - t59 * mrSges(6,3);
t10 = m(6) * t132 - t86 * mrSges(6,2) + t32 * mrSges(6,3) - t105 * t11 + t110 * t12 + t58 * t44 - t98 * t51;
t60 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t67 = -t98 * mrSges(5,2) + t74 * mrSges(5,3);
t50 = -t98 * mrSges(6,2) + t58 * mrSges(6,3);
t9 = m(6) * t126 + t86 * mrSges(6,1) - t33 * mrSges(6,3) + t105 * t12 + t110 * t11 - t59 * t44 + t98 * t50;
t7 = m(5) * t128 + t86 * mrSges(5,1) - t49 * mrSges(5,3) + t103 * t10 + t104 * t9 - t75 * t60 + t98 * t67;
t76 = -t89 * mrSges(4,1) + t90 * mrSges(4,2);
t77 = -t99 * mrSges(4,2) + t89 * mrSges(4,3);
t69 = t98 * mrSges(5,1) - t75 * mrSges(5,3);
t8 = m(5) * t138 - t86 * mrSges(5,2) + t48 * mrSges(5,3) + t104 * t10 - t103 * t9 + t74 * t60 - t98 * t69;
t5 = m(4) * t127 + t88 * mrSges(4,1) - t72 * mrSges(4,3) + t106 * t8 + t111 * t7 - t90 * t76 + t99 * t77;
t78 = t99 * mrSges(4,1) - t90 * mrSges(4,3);
t6 = m(4) * t137 - t88 * mrSges(4,2) + t71 * mrSges(4,3) - t106 * t7 + t111 * t8 + t89 * t76 - t99 * t78;
t96 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t135;
t97 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t134;
t140 = m(3) * t84 - t94 * mrSges(3,1) + t93 * mrSges(3,2) + t107 * t6 + t112 * t5 + (t108 * t96 - t113 * t97) * qJD(1);
t136 = -t113 * g(3) - t108 * t85;
t65 = -qJDD(2) * pkin(2) - t115 * pkin(8) + t92 * t135 - t136;
t121 = -t71 * pkin(3) - t87 * pkin(9) + t90 * t79 + t65;
t119 = -t48 * pkin(4) - t73 * qJ(5) + t75 * t68 + qJDD(5) + t121;
t123 = t23 * mrSges(7,1) + t42 * t39 - m(7) * (-t32 * pkin(5) - t57 * pkin(10) + t59 * t52 + t119) - t24 * mrSges(7,2) - t43 * t40;
t120 = -m(6) * t119 + t32 * mrSges(6,1) - t33 * mrSges(6,2) + t58 * t50 - t59 * t51 + t123;
t118 = -m(5) * t121 + t48 * mrSges(5,1) - t49 * mrSges(5,2) + t74 * t67 - t75 * t69 + t120;
t117 = m(4) * t65 - t71 * mrSges(4,1) + t72 * mrSges(4,2) - t89 * t77 + t90 * t78 - t118;
t91 = (-mrSges(3,1) * t113 + mrSges(3,2) * t108) * qJD(1);
t14 = m(3) * t136 + qJDD(2) * mrSges(3,1) - t93 * mrSges(3,3) + qJD(2) * t97 - t91 * t135 - t117;
t4 = m(3) * t131 - qJDD(2) * mrSges(3,2) + t94 * mrSges(3,3) - qJD(2) * t96 - t107 * t5 + t112 * t6 + t91 * t134;
t139 = t108 * t4 + t113 * t14;
t2 = m(2) * t129 + qJDD(1) * mrSges(2,1) - t116 * mrSges(2,2) - t140;
t1 = m(2) * t125 - t116 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t108 * t14 + t113 * t4;
t3 = [-m(1) * g(1) + t114 * t1 - t109 * t2, t1, t4, t6, t8, t10, t12; -m(1) * g(2) + t109 * t1 + t114 * t2, t2, t14, t5, t7, t9, t11; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t140, t117, -t118, -t120, -t123;];
f_new  = t3;
