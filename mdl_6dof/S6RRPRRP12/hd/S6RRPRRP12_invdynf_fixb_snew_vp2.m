% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP12
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 19:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP12_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP12_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP12_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:59:18
% EndTime: 2019-05-06 18:59:23
% DurationCPUTime: 1.84s
% Computational Cost: add. (18932->201), mult. (38590->243), div. (0->0), fcn. (23257->8), ass. (0->94)
t147 = -2 * qJD(3);
t102 = cos(qJ(2));
t100 = sin(qJ(1));
t103 = cos(qJ(1));
t122 = t100 * g(1) - t103 * g(2);
t117 = -qJDD(1) * pkin(1) - t122;
t101 = cos(qJ(4));
t127 = qJD(1) * qJD(2);
t120 = t102 * t127;
t99 = sin(qJ(2));
t123 = t99 * t127;
t77 = t99 * qJDD(1) + t120;
t91 = t99 * qJD(1);
t110 = pkin(2) * t123 + t91 * t147 + (-t77 - t120) * qJ(3) + t117;
t128 = qJD(1) * t102;
t105 = qJD(1) ^ 2;
t139 = t105 * pkin(7);
t140 = cos(qJ(5));
t78 = t102 * qJDD(1) - t123;
t84 = pkin(3) * t91 - qJD(2) * pkin(8);
t96 = t102 ^ 2;
t27 = -t84 * t91 + (-pkin(2) - pkin(8)) * t78 + (-pkin(3) * t96 - pkin(7)) * t105 + t110;
t104 = qJD(2) ^ 2;
t119 = -t103 * g(1) - t100 * g(2);
t66 = -t105 * pkin(1) + qJDD(1) * pkin(7) + t119;
t132 = -t102 * g(3) - t99 * t66;
t74 = (-pkin(2) * t102 - qJ(3) * t99) * qJD(1);
t42 = -qJDD(2) * pkin(2) - t104 * qJ(3) + t74 * t91 + qJDD(3) - t132;
t32 = (-t102 * t105 * t99 - qJDD(2)) * pkin(8) + (t77 - t120) * pkin(3) + t42;
t98 = sin(qJ(4));
t121 = t101 * t32 - t98 * t27;
t72 = -t98 * qJD(2) - t101 * t128;
t52 = t72 * qJD(4) + t101 * qJDD(2) - t98 * t78;
t71 = qJDD(4) + t77;
t73 = t101 * qJD(2) - t98 * t128;
t88 = t91 + qJD(4);
t17 = (t72 * t88 - t52) * pkin(9) + (t72 * t73 + t71) * pkin(4) + t121;
t134 = t101 * t27 + t98 * t32;
t51 = -t73 * qJD(4) - t98 * qJDD(2) - t101 * t78;
t58 = t88 * pkin(4) - t73 * pkin(9);
t70 = t72 ^ 2;
t19 = -t70 * pkin(4) + t51 * pkin(9) - t88 * t58 + t134;
t97 = sin(qJ(5));
t135 = t140 * t19 + t97 * t17;
t53 = -t140 * t72 + t97 * t73;
t54 = t140 * t73 + t97 * t72;
t35 = t53 * pkin(5) - t54 * qJ(6);
t86 = qJD(5) + t88;
t46 = -t86 * mrSges(7,1) + t54 * mrSges(7,2);
t67 = qJDD(5) + t71;
t85 = t86 ^ 2;
t125 = m(7) * (-t85 * pkin(5) + t67 * qJ(6) + 0.2e1 * qJD(6) * t86 - t53 * t35 + t135) + t86 * t46 + t67 * mrSges(7,3);
t36 = t53 * mrSges(7,1) - t54 * mrSges(7,3);
t133 = -t53 * mrSges(6,1) - t54 * mrSges(6,2) - t36;
t136 = -mrSges(6,3) - mrSges(7,2);
t24 = t54 * qJD(5) - t140 * t51 + t97 * t52;
t45 = t86 * mrSges(6,1) - t54 * mrSges(6,3);
t10 = m(6) * t135 - t67 * mrSges(6,2) + t133 * t53 + t136 * t24 - t86 * t45 + t125;
t116 = t140 * t17 - t97 * t19;
t141 = m(7) * (-t67 * pkin(5) - t85 * qJ(6) + t54 * t35 + qJDD(6) - t116);
t25 = -t53 * qJD(5) + t140 * t52 + t97 * t51;
t43 = -t53 * mrSges(7,2) + t86 * mrSges(7,3);
t44 = -t86 * mrSges(6,2) - t53 * mrSges(6,3);
t11 = m(6) * t116 - t141 + (t44 + t43) * t86 + (mrSges(6,1) + mrSges(7,1)) * t67 + t133 * t54 + t136 * t25;
t55 = -t72 * mrSges(5,1) + t73 * mrSges(5,2);
t56 = -t88 * mrSges(5,2) + t72 * mrSges(5,3);
t6 = m(5) * t121 + t71 * mrSges(5,1) - t52 * mrSges(5,3) + t97 * t10 + t140 * t11 - t73 * t55 + t88 * t56;
t57 = t88 * mrSges(5,1) - t73 * mrSges(5,3);
t7 = m(5) * t134 - t71 * mrSges(5,2) + t51 * mrSges(5,3) + t140 * t10 - t97 * t11 + t72 * t55 - t88 * t57;
t82 = -mrSges(4,1) * t128 - qJD(2) * mrSges(4,3);
t118 = -t101 * t7 + t98 * t6 - m(4) * (-t78 * pkin(2) + t110 - t139) - t82 * t128 + t77 * mrSges(4,3);
t83 = mrSges(4,1) * t91 + qJD(2) * mrSges(4,2);
t130 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t91 - t83;
t138 = mrSges(3,1) - mrSges(4,2);
t81 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t128;
t146 = (-t102 * t81 + t130 * t99) * qJD(1) - t138 * t78 + m(3) * (t117 - t139) + t77 * mrSges(3,2) - t118;
t124 = -t99 * g(3) + t102 * t66;
t145 = t104 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t147 - t74 * t128 - t124;
t115 = -m(4) * t42 - t101 * t6 - t98 * t7;
t75 = (mrSges(4,2) * t102 - mrSges(4,3) * t99) * qJD(1);
t131 = t75 + (-mrSges(3,1) * t102 + mrSges(3,2) * t99) * qJD(1);
t137 = -mrSges(3,3) - mrSges(4,1);
t4 = m(3) * t132 + t137 * t77 + t138 * qJDD(2) + (t81 - t82) * qJD(2) - t131 * t91 + t115;
t111 = -t96 * t105 * pkin(8) + t78 * pkin(3) + qJD(2) * t84 - t145;
t108 = -t51 * pkin(4) - t70 * pkin(9) + t73 * t58 + t111;
t113 = t25 * mrSges(7,3) + t54 * t46 - m(7) * (t108 + (t53 * t86 - t25) * qJ(6) + (t54 * t86 + t24) * pkin(5) - 0.2e1 * qJD(6) * t54) - t24 * mrSges(7,1) - t53 * t43;
t109 = m(6) * t108 + t24 * mrSges(6,1) + t25 * mrSges(6,2) + t53 * t44 + t54 * t45 - t113;
t107 = -m(5) * t111 + t51 * mrSges(5,1) - t52 * mrSges(5,2) + t72 * t56 - t73 * t57 - t109;
t106 = m(4) * t145 + t107;
t9 = -t106 - t137 * t78 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t130 * qJD(2) + m(3) * t124 + t131 * t128;
t142 = t102 * t4 + t99 * t9;
t2 = m(2) * t122 + qJDD(1) * mrSges(2,1) - t105 * mrSges(2,2) - t146;
t1 = m(2) * t119 - t105 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t102 * t9 - t99 * t4;
t3 = [-m(1) * g(1) + t103 * t1 - t100 * t2, t1, t9, t78 * mrSges(4,2) - t83 * t91 - t118, t7, t10, -t24 * mrSges(7,2) - t53 * t36 + t125; -m(1) * g(2) + t100 * t1 + t103 * t2, t2, t4, -t78 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t83 - t75 * t128 + t106, t6, t11, -t113; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t146, t77 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t82 + t75 * t91 - t115, -t107, t109, -t67 * mrSges(7,1) + t25 * mrSges(7,2) + t54 * t36 - t86 * t43 + t141;];
f_new  = t3;
