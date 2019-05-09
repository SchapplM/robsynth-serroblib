% Calculate vector of cutting forces with Newton-Euler
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-05-06 08:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:10:26
% EndTime: 2019-05-06 08:10:32
% DurationCPUTime: 2.58s
% Computational Cost: add. (29803->208), mult. (70845->264), div. (0->0), fcn. (48892->10), ass. (0->101)
t105 = sin(pkin(10));
t138 = cos(pkin(10));
t106 = sin(pkin(9));
t108 = sin(qJ(2));
t111 = cos(qJ(2));
t139 = cos(pkin(9));
t89 = (t106 * t111 + t139 * t108) * qJD(1);
t77 = -t138 * qJD(2) + t105 * t89;
t135 = qJD(1) * t111;
t136 = qJD(1) * t108;
t88 = t106 * t136 - t139 * t135;
t146 = t77 * t88;
t134 = qJD(1) * qJD(2);
t95 = t108 * qJDD(1) + t111 * t134;
t96 = t111 * qJDD(1) - t108 * t134;
t73 = t106 * t96 + t139 * t95;
t64 = t105 * qJDD(2) + t138 * t73;
t153 = (-t64 + t146) * qJ(5);
t114 = qJD(1) ^ 2;
t107 = sin(qJ(6));
t110 = cos(qJ(6));
t113 = qJD(2) ^ 2;
t109 = sin(qJ(1));
t112 = cos(qJ(1));
t127 = -t112 * g(1) - t109 * g(2);
t92 = -t114 * pkin(1) + qJDD(1) * pkin(7) + t127;
t140 = t108 * t92;
t147 = pkin(2) * t114;
t49 = qJDD(2) * pkin(2) - t95 * qJ(3) - t140 + (qJ(3) * t134 + t108 * t147 - g(3)) * t111;
t104 = t111 ^ 2;
t131 = -t108 * g(3) + t111 * t92;
t97 = qJD(2) * pkin(2) - qJ(3) * t136;
t53 = t96 * qJ(3) - qJD(2) * t97 - t104 * t147 + t131;
t132 = -0.2e1 * qJD(3) * t88 + t106 * t49 + t139 * t53;
t66 = t88 * pkin(3) - t89 * qJ(4);
t27 = -t113 * pkin(3) + qJDD(2) * qJ(4) - t88 * t66 + t132;
t129 = t109 * g(1) - t112 * g(2);
t124 = -qJDD(1) * pkin(1) - t129;
t115 = -t96 * pkin(2) + qJDD(3) + t97 * t136 + (-qJ(3) * t104 - pkin(7)) * t114 + t124;
t72 = t106 * t95 - t139 * t96;
t29 = (qJD(2) * t88 - t73) * qJ(4) + (qJD(2) * t89 + t72) * pkin(3) + t115;
t126 = -t105 * t27 + t138 * t29;
t78 = t105 * qJD(2) + t138 * t89;
t50 = pkin(4) * t77 - qJ(5) * t78;
t87 = t88 ^ 2;
t21 = -t72 * pkin(4) - t87 * qJ(5) + qJDD(5) - t126 + ((2 * qJD(4)) + t50) * t78;
t16 = (-t64 - t146) * pkin(8) + (t77 * t78 - t72) * pkin(5) + t21;
t150 = -2 * qJD(4);
t133 = t105 * t29 + t138 * t27 + t77 * t150;
t149 = 2 * qJD(5);
t120 = -t87 * pkin(4) + t72 * qJ(5) + t88 * t149 - t77 * t50 + t133;
t61 = -t88 * pkin(5) - t78 * pkin(8);
t63 = -t138 * qJDD(2) + t105 * t73;
t76 = t77 ^ 2;
t17 = -t76 * pkin(5) + t63 * pkin(8) + t88 * t61 + t120;
t45 = -t107 * t78 + t110 * t77;
t33 = t45 * qJD(6) + t107 * t63 + t110 * t64;
t46 = t107 * t77 + t110 * t78;
t36 = -mrSges(7,1) * t45 + mrSges(7,2) * t46;
t86 = qJD(6) - t88;
t39 = -t86 * mrSges(7,2) + t45 * mrSges(7,3);
t71 = qJDD(6) - t72;
t14 = m(7) * (-t107 * t17 + t110 * t16) - t33 * mrSges(7,3) + t71 * mrSges(7,1) - t46 * t36 + t86 * t39;
t32 = -t46 * qJD(6) - t107 * t64 + t110 * t63;
t40 = t86 * mrSges(7,1) - t46 * mrSges(7,3);
t15 = m(7) * (t107 * t16 + t110 * t17) + t32 * mrSges(7,3) - t71 * mrSges(7,2) + t45 * t36 - t86 * t40;
t121 = -m(6) * t21 - t107 * t15 - t110 * t14;
t51 = mrSges(6,1) * t77 - mrSges(6,3) * t78;
t142 = -mrSges(5,1) * t77 - mrSges(5,2) * t78 - t51;
t144 = -mrSges(5,3) - mrSges(6,2);
t57 = -t77 * mrSges(6,2) + t88 * mrSges(6,3);
t58 = -t88 * mrSges(5,2) - t77 * mrSges(5,3);
t11 = m(5) * t126 + (t58 + t57) * t88 + (m(5) * t150 + t142) * t78 + (mrSges(5,1) + mrSges(6,1)) * t72 + t144 * t64 + t121;
t79 = -qJD(2) * mrSges(4,2) - t88 * mrSges(4,3);
t80 = qJD(2) * mrSges(4,1) - t89 * mrSges(4,3);
t60 = -t88 * mrSges(6,1) + t78 * mrSges(6,2);
t123 = m(6) * t120 + t72 * mrSges(6,3) - t107 * t14 + t110 * t15 + t88 * t60;
t59 = t88 * mrSges(5,1) - t78 * mrSges(5,3);
t9 = m(5) * t133 - t72 * mrSges(5,2) + t142 * t77 + t144 * t63 - t88 * t59 + t123;
t119 = m(4) * t115 + t72 * mrSges(4,1) + t73 * mrSges(4,2) + t105 * t9 + t138 * t11 + t88 * t79 + t89 * t80;
t98 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t136;
t99 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t135;
t152 = (t108 * t98 - t111 * t99) * qJD(1) + m(3) * (-t114 * pkin(7) + t124) - t96 * mrSges(3,1) + t95 * mrSges(3,2) + t119;
t143 = -t106 * t53 + t139 * t49;
t118 = qJDD(2) * pkin(3) + t113 * qJ(4) - t89 * t66 - qJDD(4) + t143;
t137 = qJD(3) * t89;
t83 = -0.2e1 * t137;
t128 = m(7) * (-t76 * pkin(8) + t83 + (-pkin(4) - pkin(5)) * t63 - t153 + (-pkin(4) * t88 + t149 + t61) * t78 + t118) + t33 * mrSges(7,2) - t32 * mrSges(7,1) + t46 * t40 - t45 * t39;
t26 = -t118 + 0.2e1 * t137;
t122 = m(6) * (-0.2e1 * qJD(5) * t78 + t153 + (t78 * t88 + t63) * pkin(4) + t26) + t77 * t57 + t63 * mrSges(6,1) - t128;
t151 = m(5) * t26 + t63 * mrSges(5,1) + (t59 - t60) * t78 + (mrSges(5,2) - mrSges(6,3)) * t64 + t77 * t58 + t122;
t67 = t88 * mrSges(4,1) + t89 * mrSges(4,2);
t12 = qJDD(2) * mrSges(4,1) - t73 * mrSges(4,3) - t89 * t67 + qJD(2) * t79 + m(4) * (t83 + t143) - t151;
t7 = m(4) * t132 - qJDD(2) * mrSges(4,2) - t72 * mrSges(4,3) - qJD(2) * t80 - t105 * t11 + t138 * t9 - t88 * t67;
t94 = (-mrSges(3,1) * t111 + mrSges(3,2) * t108) * qJD(1);
t4 = m(3) * (-t111 * g(3) - t140) - t95 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t94 * t136 + qJD(2) * t99 + t106 * t7 + t139 * t12;
t5 = m(3) * t131 - qJDD(2) * mrSges(3,2) + t96 * mrSges(3,3) - qJD(2) * t98 - t106 * t12 + t94 * t135 + t139 * t7;
t148 = t108 * t5 + t111 * t4;
t6 = m(2) * t129 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t152;
t1 = m(2) * t127 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t108 * t4 + t111 * t5;
t2 = [-m(1) * g(1) + t112 * t1 - t109 * t6, t1, t5, t7, t9, -t63 * mrSges(6,2) - t77 * t51 + t123, t15; -m(1) * g(2) + t109 * t1 + t112 * t6, t6, t4, t12, t11, -t64 * mrSges(6,3) - t78 * t60 + t122, t14; (-m(1) - m(2)) * g(3) + t148, -m(2) * g(3) + t148, t152, t119, t151, -t72 * mrSges(6,1) + t64 * mrSges(6,2) + t78 * t51 - t88 * t57 - t121, t128;];
f_new  = t2;
