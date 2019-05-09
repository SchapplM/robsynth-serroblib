% Calculate vector of cutting forces with Newton-Euler
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 04:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:19:51
% EndTime: 2019-05-07 04:19:59
% DurationCPUTime: 3.23s
% Computational Cost: add. (39023->207), mult. (88519->262), div. (0->0), fcn. (63622->10), ass. (0->103)
t129 = cos(pkin(10));
t100 = sin(qJ(3));
t101 = sin(qJ(2));
t104 = cos(qJ(3));
t105 = cos(qJ(2));
t82 = (-t100 * t101 + t104 * t105) * qJD(1);
t83 = (t100 * t105 + t101 * t104) * qJD(1);
t98 = sin(pkin(10));
t73 = -t129 * t82 + t98 * t83;
t74 = t129 * t83 + t98 * t82;
t44 = t73 * pkin(4) - t74 * qJ(5);
t150 = (2 * qJD(4)) + t44;
t107 = qJD(1) ^ 2;
t102 = sin(qJ(1));
t106 = cos(qJ(1));
t125 = t102 * g(1) - t106 * g(2);
t119 = -qJDD(1) * pkin(1) - t125;
t128 = qJD(1) * t101;
t126 = qJD(1) * qJD(2);
t89 = t105 * qJDD(1) - t101 * t126;
t92 = qJD(2) * pkin(2) - pkin(8) * t128;
t97 = t105 ^ 2;
t112 = -t89 * pkin(2) + t92 * t128 + (-pkin(8) * t97 - pkin(7)) * t107 + t119;
t88 = t101 * qJDD(1) + t105 * t126;
t57 = -t83 * qJD(3) - t100 * t88 + t104 * t89;
t96 = qJD(2) + qJD(3);
t77 = t96 * pkin(3) - t83 * qJ(4);
t81 = t82 ^ 2;
t110 = -t57 * pkin(3) - t81 * qJ(4) + t83 * t77 + qJDD(4) + t112;
t103 = cos(qJ(6));
t140 = t73 * t96;
t145 = -2 * qJD(5);
t58 = t82 * qJD(3) + t100 * t89 + t104 * t88;
t38 = t129 * t58 + t98 * t57;
t108 = (-t38 + t140) * qJ(5) + t110 + (t96 * pkin(4) + t145) * t74;
t122 = -t106 * g(1) - t102 * g(2);
t85 = -t107 * pkin(1) + qJDD(1) * pkin(7) + t122;
t130 = t101 * t85;
t142 = pkin(2) * t107;
t50 = qJDD(2) * pkin(2) - t88 * pkin(8) - t130 + (pkin(8) * t126 + t101 * t142 - g(3)) * t105;
t124 = -t101 * g(3) + t105 * t85;
t51 = t89 * pkin(8) - qJD(2) * t92 - t97 * t142 + t124;
t123 = -t100 * t51 + t104 * t50;
t95 = qJDD(2) + qJDD(3);
t23 = (t82 * t96 - t58) * qJ(4) + (t82 * t83 + t95) * pkin(3) + t123;
t132 = t100 * t50 + t104 * t51;
t26 = -t81 * pkin(3) + t57 * qJ(4) - t96 * t77 + t132;
t121 = t129 * t23 - t98 * t26;
t94 = t96 ^ 2;
t19 = -t95 * pkin(4) - t94 * qJ(5) + t150 * t74 + qJDD(5) - t121;
t14 = (t73 * t74 - t95) * pkin(9) + (t38 + t140) * pkin(5) + t19;
t37 = -t129 * t57 + t98 * t58;
t64 = t74 * pkin(5) - t96 * pkin(9);
t70 = t73 ^ 2;
t17 = -t74 * t64 - t70 * pkin(5) + t108 + (pkin(4) + pkin(9)) * t37;
t99 = sin(qJ(6));
t54 = t103 * t73 - t99 * t96;
t29 = t54 * qJD(6) + t103 * t95 + t99 * t37;
t36 = qJDD(6) + t38;
t55 = t103 * t96 + t99 * t73;
t39 = -t54 * mrSges(7,1) + t55 * mrSges(7,2);
t69 = qJD(6) + t74;
t40 = -t69 * mrSges(7,2) + t54 * mrSges(7,3);
t12 = m(7) * (t103 * t14 - t99 * t17) - t29 * mrSges(7,3) + t36 * mrSges(7,1) - t55 * t39 + t69 * t40;
t28 = -t55 * qJD(6) + t103 * t37 - t99 * t95;
t41 = t69 * mrSges(7,1) - t55 * mrSges(7,3);
t13 = m(7) * (t103 * t17 + t99 * t14) + t28 * mrSges(7,3) - t36 * mrSges(7,2) + t54 * t39 - t69 * t41;
t63 = t74 * mrSges(6,1) + t96 * mrSges(6,2);
t118 = -t103 * t13 + t99 * t12 - m(6) * (t37 * pkin(4) + t108) + t38 * mrSges(6,3) + t74 * t63;
t62 = t73 * mrSges(6,1) - t96 * mrSges(6,3);
t131 = t96 * mrSges(5,2) + t73 * mrSges(5,3) + t62;
t136 = mrSges(6,2) - mrSges(5,1);
t61 = t96 * mrSges(5,1) - t74 * mrSges(5,3);
t111 = m(5) * t110 + t38 * mrSges(5,2) - t131 * t73 - t136 * t37 + t74 * t61 - t118;
t76 = -t96 * mrSges(4,2) + t82 * mrSges(4,3);
t78 = t96 * mrSges(4,1) - t83 * mrSges(4,3);
t109 = m(4) * t112 - t57 * mrSges(4,1) + t58 * mrSges(4,2) - t82 * t76 + t83 * t78 + t111;
t90 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t128;
t127 = qJD(1) * t105;
t91 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t127;
t149 = t109 + (t101 * t90 - t105 * t91) * qJD(1) - t89 * mrSges(3,1) + t88 * mrSges(3,2) + m(3) * (-t107 * pkin(7) + t119);
t147 = -2 * qJD(4);
t134 = t129 * t26 + t98 * t23;
t117 = t94 * pkin(4) - t95 * qJ(5) - t134;
t67 = t73 * t147;
t115 = -t28 * mrSges(7,1) - t54 * t40 + m(7) * (-t37 * pkin(5) - t70 * pkin(9) - t73 * t44 + t67 + ((2 * qJD(5)) + t64) * t96 - t117) + t29 * mrSges(7,2) + t55 * t41;
t113 = -m(6) * (t96 * t145 + t150 * t73 + t117) + t115;
t46 = -t73 * mrSges(6,2) - t74 * mrSges(6,3);
t133 = -t73 * mrSges(5,1) - t74 * mrSges(5,2) - t46;
t135 = -mrSges(5,3) - mrSges(6,1);
t10 = m(5) * (t67 + t134) + (-t61 + t63) * t96 + (-mrSges(5,2) + mrSges(6,3)) * t95 + t133 * t73 + t135 * t37 + t113;
t75 = -t82 * mrSges(4,1) + t83 * mrSges(4,2);
t116 = -m(6) * t19 - t103 * t12 - t99 * t13;
t9 = m(5) * t121 - t131 * t96 - t136 * t95 + (m(5) * t147 + t133) * t74 + t135 * t38 + t116;
t6 = m(4) * t123 + t95 * mrSges(4,1) - t58 * mrSges(4,3) + t98 * t10 + t129 * t9 - t83 * t75 + t96 * t76;
t7 = m(4) * t132 - t95 * mrSges(4,2) + t57 * mrSges(4,3) + t129 * t10 + t82 * t75 - t96 * t78 - t98 * t9;
t87 = (-mrSges(3,1) * t105 + mrSges(3,2) * t101) * qJD(1);
t4 = m(3) * (-t105 * g(3) - t130) - t88 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t87 * t128 + qJD(2) * t91 + t100 * t7 + t104 * t6;
t5 = m(3) * t124 - qJDD(2) * mrSges(3,2) + t89 * mrSges(3,3) - qJD(2) * t90 - t100 * t6 + t104 * t7 + t87 * t127;
t144 = t101 * t5 + t105 * t4;
t8 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t107 * mrSges(2,2) - t149;
t1 = m(2) * t122 - t107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t101 * t4 + t105 * t5;
t2 = [-m(1) * g(1) + t106 * t1 - t102 * t8, t1, t5, t7, t10, -t37 * mrSges(6,2) - t73 * t62 - t118, t13; -m(1) * g(2) + t102 * t1 + t106 * t8, t8, t4, t6, t9, t37 * mrSges(6,1) - t95 * mrSges(6,3) + t73 * t46 - t96 * t63 - t113, t12; (-m(1) - m(2)) * g(3) + t144, -m(2) * g(3) + t144, t149, t109, t111, t38 * mrSges(6,1) + t95 * mrSges(6,2) + t74 * t46 + t96 * t62 - t116, t115;];
f_new  = t2;
