% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:37:03
% EndTime: 2019-05-07 07:37:13
% DurationCPUTime: 3.70s
% Computational Cost: add. (50592->201), mult. (104804->261), div. (0->0), fcn. (74708->10), ass. (0->97)
t104 = sin(qJ(2));
t106 = cos(qJ(2));
t108 = qJD(1) ^ 2;
t102 = sin(qJ(5));
t135 = cos(qJ(5));
t100 = sin(pkin(10));
t101 = cos(pkin(10));
t105 = sin(qJ(1));
t107 = cos(qJ(1));
t123 = t105 * g(1) - t107 * g(2);
t117 = -qJDD(1) * pkin(1) - t123;
t128 = qJD(1) * t104;
t126 = qJD(1) * qJD(2);
t90 = t106 * qJDD(1) - t104 * t126;
t93 = qJD(2) * pkin(2) - pkin(8) * t128;
t99 = t106 ^ 2;
t113 = -t90 * pkin(2) + t93 * t128 + (-pkin(8) * t99 - pkin(7)) * t108 + t117;
t103 = sin(qJ(3));
t136 = cos(qJ(3));
t83 = (t103 * t106 + t136 * t104) * qJD(1);
t89 = t104 * qJDD(1) + t106 * t126;
t62 = t83 * qJD(3) + t103 * t89 - t136 * t90;
t127 = qJD(1) * t106;
t82 = t103 * t128 - t136 * t127;
t63 = -t82 * qJD(3) + t103 * t90 + t136 * t89;
t98 = qJD(2) + qJD(3);
t31 = (t82 * t98 - t63) * qJ(4) + (t83 * t98 + t62) * pkin(3) + t113;
t119 = -t107 * g(1) - t105 * g(2);
t85 = -t108 * pkin(1) + qJDD(1) * pkin(7) + t119;
t129 = t104 * t85;
t134 = pkin(2) * t108;
t55 = qJDD(2) * pkin(2) - t89 * pkin(8) - t129 + (pkin(8) * t126 + t104 * t134 - g(3)) * t106;
t122 = -t104 * g(3) + t106 * t85;
t56 = t90 * pkin(8) - qJD(2) * t93 - t99 * t134 + t122;
t130 = t103 * t55 + t136 * t56;
t69 = t82 * pkin(3) - t83 * qJ(4);
t96 = t98 ^ 2;
t97 = qJDD(2) + qJDD(3);
t34 = -t96 * pkin(3) + t97 * qJ(4) - t82 * t69 + t130;
t76 = t100 * t98 + t101 * t83;
t121 = -0.2e1 * qJD(4) * t76 - t100 * t34 + t101 * t31;
t50 = t100 * t97 + t101 * t63;
t75 = -t100 * t83 + t101 * t98;
t20 = (t75 * t82 - t50) * pkin(9) + (t75 * t76 + t62) * pkin(4) + t121;
t124 = 0.2e1 * qJD(4) * t75 + t100 * t31 + t101 * t34;
t49 = -t100 * t63 + t101 * t97;
t67 = t82 * pkin(4) - t76 * pkin(9);
t74 = t75 ^ 2;
t22 = -t74 * pkin(4) + t49 * pkin(9) - t82 * t67 + t124;
t132 = t102 * t20 + t135 * t22;
t47 = t102 * t76 - t135 * t75;
t48 = t102 * t75 + t135 * t76;
t37 = t47 * pkin(5) - t48 * qJ(6);
t81 = qJD(5) + t82;
t44 = -t81 * mrSges(7,1) + t48 * mrSges(7,2);
t60 = qJDD(5) + t62;
t80 = t81 ^ 2;
t125 = m(7) * (-t80 * pkin(5) + t60 * qJ(6) + 0.2e1 * qJD(6) * t81 - t47 * t37 + t132) + t81 * t44 + t60 * mrSges(7,3);
t38 = t47 * mrSges(7,1) - t48 * mrSges(7,3);
t131 = -t47 * mrSges(6,1) - t48 * mrSges(6,2) - t38;
t133 = -mrSges(6,3) - mrSges(7,2);
t29 = t48 * qJD(5) + t102 * t50 - t135 * t49;
t43 = t81 * mrSges(6,1) - t48 * mrSges(6,3);
t12 = m(6) * t132 - t60 * mrSges(6,2) + t131 * t47 + t133 * t29 - t81 * t43 + t125;
t116 = -t102 * t22 + t135 * t20;
t137 = m(7) * (-t60 * pkin(5) - t80 * qJ(6) + t48 * t37 + qJDD(6) - t116);
t30 = -t47 * qJD(5) + t102 * t49 + t135 * t50;
t41 = -t47 * mrSges(7,2) + t81 * mrSges(7,3);
t42 = -t81 * mrSges(6,2) - t47 * mrSges(6,3);
t13 = m(6) * t116 - t137 + (t42 + t41) * t81 + (mrSges(6,1) + mrSges(7,1)) * t60 + t131 * t48 + t133 * t30;
t52 = -t75 * mrSges(5,1) + t76 * mrSges(5,2);
t65 = -t82 * mrSges(5,2) + t75 * mrSges(5,3);
t10 = m(5) * t121 + t62 * mrSges(5,1) - t50 * mrSges(5,3) + t102 * t12 + t135 * t13 - t76 * t52 + t82 * t65;
t66 = t82 * mrSges(5,1) - t76 * mrSges(5,3);
t11 = m(5) * t124 - t62 * mrSges(5,2) + t49 * mrSges(5,3) - t102 * t13 + t135 * t12 + t75 * t52 - t82 * t66;
t77 = -t98 * mrSges(4,2) - t82 * mrSges(4,3);
t78 = t98 * mrSges(4,1) - t83 * mrSges(4,3);
t114 = m(4) * t113 + t62 * mrSges(4,1) + t63 * mrSges(4,2) + t101 * t10 + t100 * t11 + t82 * t77 + t83 * t78;
t91 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t128;
t92 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t127;
t139 = (t104 * t91 - t106 * t92) * qJD(1) + m(3) * (-t108 * pkin(7) + t117) - t90 * mrSges(3,1) + t89 * mrSges(3,2) + t114;
t120 = -t103 * t56 + t136 * t55;
t33 = -t97 * pkin(3) - t96 * qJ(4) + t83 * t69 + qJDD(4) - t120;
t110 = -t49 * pkin(4) - t74 * pkin(9) + t76 * t67 + t33;
t115 = t30 * mrSges(7,3) + t48 * t44 - m(7) * (-0.2e1 * qJD(6) * t48 + (t47 * t81 - t30) * qJ(6) + (t48 * t81 + t29) * pkin(5) + t110) - t29 * mrSges(7,1) - t47 * t41;
t112 = m(6) * t110 + t29 * mrSges(6,1) + t30 * mrSges(6,2) + t47 * t42 + t48 * t43 - t115;
t109 = m(5) * t33 - t49 * mrSges(5,1) + t50 * mrSges(5,2) - t75 * t65 + t76 * t66 + t112;
t70 = t82 * mrSges(4,1) + t83 * mrSges(4,2);
t14 = m(4) * t120 + t97 * mrSges(4,1) - t63 * mrSges(4,3) - t83 * t70 + t98 * t77 - t109;
t7 = m(4) * t130 - t97 * mrSges(4,2) - t62 * mrSges(4,3) - t100 * t10 + t101 * t11 - t82 * t70 - t98 * t78;
t88 = (-mrSges(3,1) * t106 + mrSges(3,2) * t104) * qJD(1);
t4 = m(3) * (-t106 * g(3) - t129) - t89 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t88 * t128 + qJD(2) * t92 + t103 * t7 + t136 * t14;
t5 = m(3) * t122 - qJDD(2) * mrSges(3,2) + t90 * mrSges(3,3) - qJD(2) * t91 - t103 * t14 + t88 * t127 + t136 * t7;
t138 = t104 * t5 + t106 * t4;
t6 = m(2) * t123 + qJDD(1) * mrSges(2,1) - t108 * mrSges(2,2) - t139;
t1 = m(2) * t119 - t108 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t104 * t4 + t106 * t5;
t2 = [-m(1) * g(1) + t107 * t1 - t105 * t6, t1, t5, t7, t11, t12, -t29 * mrSges(7,2) - t47 * t38 + t125; -m(1) * g(2) + t105 * t1 + t107 * t6, t6, t4, t14, t10, t13, -t115; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t139, t114, t109, t112, -t60 * mrSges(7,1) + t30 * mrSges(7,2) + t48 * t38 - t81 * t41 + t137;];
f_new  = t2;
