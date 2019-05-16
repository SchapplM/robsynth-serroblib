% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 10:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:39:22
% EndTime: 2019-05-07 10:39:31
% DurationCPUTime: 2.99s
% Computational Cost: add. (38731->204), mult. (79602->257), div. (0->0), fcn. (55552->10), ass. (0->102)
t150 = -2 * qJD(4);
t104 = sin(qJ(2));
t108 = cos(qJ(2));
t110 = qJD(1) ^ 2;
t100 = t108 ^ 2;
t105 = sin(qJ(1));
t109 = cos(qJ(1));
t128 = t105 * g(1) - t109 * g(2);
t120 = -qJDD(1) * pkin(1) - t128;
t131 = qJD(1) * t104;
t129 = qJD(1) * qJD(2);
t89 = t108 * qJDD(1) - t104 * t129;
t92 = qJD(2) * pkin(2) - pkin(8) * t131;
t115 = -t89 * pkin(2) + t92 * t131 + (-pkin(8) * t100 - pkin(7)) * t110 + t120;
t101 = sin(qJ(6));
t106 = cos(qJ(6));
t102 = sin(qJ(5));
t107 = cos(qJ(5));
t103 = sin(qJ(3));
t130 = qJD(1) * t108;
t144 = cos(qJ(3));
t81 = t103 * t131 - t144 * t130;
t99 = qJD(2) + qJD(3);
t142 = t81 * t99;
t88 = t104 * qJDD(1) + t108 * t129;
t56 = -t81 * qJD(3) + t103 * t89 + t144 * t88;
t82 = (t103 * t108 + t144 * t104) * qJD(1);
t111 = (-t56 + t142) * qJ(4) + t115 + (t99 * pkin(3) + t150) * t82;
t55 = t82 * qJD(3) + t103 * t88 - t144 * t89;
t74 = t82 * pkin(4) - t99 * pkin(9);
t80 = t81 ^ 2;
t19 = -t80 * pkin(4) - t82 * t74 + (pkin(3) + pkin(9)) * t55 + t111;
t124 = -t109 * g(1) - t105 * g(2);
t84 = -t110 * pkin(1) + qJDD(1) * pkin(7) + t124;
t133 = t104 * t84;
t143 = pkin(2) * t110;
t46 = qJDD(2) * pkin(2) - t88 * pkin(8) - t133 + (pkin(8) * t129 + t104 * t143 - g(3)) * t108;
t127 = -t104 * g(3) + t108 * t84;
t47 = t89 * pkin(8) - qJD(2) * t92 - t100 * t143 + t127;
t125 = -t103 * t47 + t144 * t46;
t63 = t81 * pkin(3) - t82 * qJ(4);
t97 = t99 ^ 2;
t98 = qJDD(2) + qJDD(3);
t31 = -t98 * pkin(3) - t97 * qJ(4) + t82 * t63 + qJDD(4) - t125;
t22 = (t81 * t82 - t98) * pkin(9) + (t56 + t142) * pkin(4) + t31;
t126 = -t102 * t19 + t107 * t22;
t68 = -t102 * t99 + t107 * t81;
t36 = t68 * qJD(5) + t102 * t55 + t107 * t98;
t53 = qJDD(5) + t56;
t69 = t102 * t81 + t107 * t99;
t79 = qJD(5) + t82;
t14 = (t68 * t79 - t36) * pkin(10) + (t68 * t69 + t53) * pkin(5) + t126;
t137 = t102 * t22 + t107 * t19;
t35 = -t69 * qJD(5) - t102 * t98 + t107 * t55;
t60 = t79 * pkin(5) - t69 * pkin(10);
t67 = t68 ^ 2;
t15 = -t67 * pkin(5) + t35 * pkin(10) - t79 * t60 + t137;
t40 = -t101 * t69 + t106 * t68;
t27 = t40 * qJD(6) + t101 * t35 + t106 * t36;
t41 = t101 * t68 + t106 * t69;
t33 = -t40 * mrSges(7,1) + t41 * mrSges(7,2);
t77 = qJD(6) + t79;
t37 = -t77 * mrSges(7,2) + t40 * mrSges(7,3);
t48 = qJDD(6) + t53;
t12 = m(7) * (-t101 * t15 + t106 * t14) - t27 * mrSges(7,3) + t48 * mrSges(7,1) - t41 * t33 + t77 * t37;
t26 = -t41 * qJD(6) - t101 * t36 + t106 * t35;
t38 = t77 * mrSges(7,1) - t41 * mrSges(7,3);
t13 = m(7) * (t101 * t14 + t106 * t15) + t26 * mrSges(7,3) - t48 * mrSges(7,2) + t40 * t33 - t77 * t38;
t45 = -t68 * mrSges(6,1) + t69 * mrSges(6,2);
t59 = t79 * mrSges(6,1) - t69 * mrSges(6,3);
t10 = m(6) * t137 - t53 * mrSges(6,2) + t35 * mrSges(6,3) - t101 * t12 + t106 * t13 + t68 * t45 - t79 * t59;
t73 = t82 * mrSges(5,1) + t99 * mrSges(5,2);
t58 = -t79 * mrSges(6,2) + t68 * mrSges(6,3);
t9 = m(6) * t126 + t53 * mrSges(6,1) - t36 * mrSges(6,3) + t101 * t13 + t106 * t12 - t69 * t45 + t79 * t58;
t121 = -t107 * t10 + t102 * t9 - m(5) * (t55 * pkin(3) + t111) + t56 * mrSges(5,3) + t82 * t73;
t72 = t81 * mrSges(5,1) - t99 * mrSges(5,3);
t134 = -t99 * mrSges(4,2) - t81 * mrSges(4,3) - t72;
t139 = mrSges(4,1) - mrSges(5,2);
t71 = t99 * mrSges(4,1) - t82 * mrSges(4,3);
t113 = m(4) * t115 + t56 * mrSges(4,2) + t134 * t81 + t139 * t55 + t82 * t71 - t121;
t90 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t131;
t91 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t130;
t149 = t113 + (t104 * t90 - t108 * t91) * qJD(1) - t89 * mrSges(3,1) + t88 * mrSges(3,2) + m(3) * (-t110 * pkin(7) + t120);
t136 = t103 * t46 + t144 * t47;
t148 = t97 * pkin(3) - t98 * qJ(4) + t99 * t150 + t81 * t63 - t136;
t116 = -t55 * pkin(4) - t80 * pkin(9) + t99 * t74 - t148;
t118 = -t26 * mrSges(7,1) - t40 * t37 + m(7) * (-t35 * pkin(5) - t67 * pkin(10) + t69 * t60 + t116) + t27 * mrSges(7,2) + t41 * t38;
t114 = m(6) * t116 - t35 * mrSges(6,1) + t36 * mrSges(6,2) - t68 * t58 + t69 * t59 + t118;
t112 = -m(5) * t148 + t114;
t65 = -t81 * mrSges(5,2) - t82 * mrSges(5,3);
t135 = -t81 * mrSges(4,1) - t82 * mrSges(4,2) - t65;
t138 = -mrSges(4,3) - mrSges(5,1);
t11 = t112 + (-t71 + t73) * t99 + (-mrSges(4,2) + mrSges(5,3)) * t98 + t135 * t81 + t138 * t55 + m(4) * t136;
t119 = -m(5) * t31 - t102 * t10 - t107 * t9;
t7 = m(4) * t125 + t134 * t99 + t135 * t82 + t138 * t56 + t139 * t98 + t119;
t87 = (-mrSges(3,1) * t108 + mrSges(3,2) * t104) * qJD(1);
t4 = m(3) * (-t108 * g(3) - t133) - t88 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t87 * t131 + qJD(2) * t91 + t103 * t11 + t144 * t7;
t5 = m(3) * t127 - qJDD(2) * mrSges(3,2) + t89 * mrSges(3,3) - qJD(2) * t90 - t103 * t7 + t144 * t11 + t87 * t130;
t146 = t104 * t5 + t108 * t4;
t6 = m(2) * t128 + qJDD(1) * mrSges(2,1) - t110 * mrSges(2,2) - t149;
t1 = m(2) * t124 - t110 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t104 * t4 + t108 * t5;
t2 = [-m(1) * g(1) + t109 * t1 - t105 * t6, t1, t5, t11, -t55 * mrSges(5,2) - t81 * t72 - t121, t10, t13; -m(1) * g(2) + t105 * t1 + t109 * t6, t6, t4, t7, t55 * mrSges(5,1) - t98 * mrSges(5,3) + t81 * t65 - t99 * t73 - t112, t9, t12; (-m(1) - m(2)) * g(3) + t146, -m(2) * g(3) + t146, t149, t113, t56 * mrSges(5,1) + t98 * mrSges(5,2) + t82 * t65 + t99 * t72 - t119, t114, t118;];
f_new  = t2;
