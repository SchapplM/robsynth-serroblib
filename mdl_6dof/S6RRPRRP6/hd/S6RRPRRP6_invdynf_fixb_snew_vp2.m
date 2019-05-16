% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:58:20
% EndTime: 2019-05-06 17:58:34
% DurationCPUTime: 5.77s
% Computational Cost: add. (70349->209), mult. (183743->281), div. (0->0), fcn. (144017->12), ass. (0->108)
t151 = -2 * qJD(3);
t101 = sin(pkin(11));
t103 = cos(pkin(11));
t107 = sin(qJ(2));
t110 = cos(qJ(2));
t102 = sin(pkin(6));
t134 = qJD(1) * t102;
t84 = (t101 * t107 - t103 * t110) * t134;
t104 = cos(pkin(6));
t112 = qJD(1) ^ 2;
t108 = sin(qJ(1));
t111 = cos(qJ(1));
t125 = t108 * g(1) - t111 * g(2);
t146 = pkin(8) * t102;
t90 = qJDD(1) * pkin(1) + t112 * t146 + t125;
t139 = t104 * t90;
t121 = -t111 * g(1) - t108 * g(2);
t91 = -t112 * pkin(1) + qJDD(1) * t146 + t121;
t123 = -t107 * t91 + t110 * t139;
t138 = t112 * t102 ^ 2;
t132 = qJD(1) * qJD(2);
t93 = (qJDD(1) * t107 + t110 * t132) * t102;
t96 = t104 * qJDD(1) + qJDD(2);
t97 = t104 * qJD(1) + qJD(2);
t44 = t96 * pkin(2) - t93 * qJ(3) + (pkin(2) * t107 * t138 + (qJ(3) * qJD(1) * t97 - g(3)) * t102) * t110 + t123;
t136 = t102 * t107;
t119 = -g(3) * t136 + t107 * t139 + t110 * t91;
t128 = t110 ^ 2 * t138;
t127 = t107 * t134;
t87 = t97 * pkin(2) - qJ(3) * t127;
t94 = (qJDD(1) * t110 - t107 * t132) * t102;
t49 = -pkin(2) * t128 + t94 * qJ(3) - t97 * t87 + t119;
t85 = (t101 * t110 + t103 * t107) * t134;
t150 = -t101 * t49 + t103 * t44 + t85 * t151;
t106 = sin(qJ(4));
t109 = cos(qJ(4));
t129 = t101 * t44 + t103 * t49 + t84 * t151;
t65 = t84 * pkin(3) - t85 * pkin(9);
t95 = t97 ^ 2;
t29 = -t95 * pkin(3) + t96 * pkin(9) - t84 * t65 + t129;
t120 = -t104 * g(3) - t102 * t90;
t114 = -t94 * pkin(2) - qJ(3) * t128 + t87 * t127 + qJDD(3) + t120;
t68 = -t101 * t93 + t103 * t94;
t69 = t101 * t94 + t103 * t93;
t31 = (t84 * t97 - t69) * pkin(9) + (t85 * t97 - t68) * pkin(3) + t114;
t122 = -t106 * t29 + t109 * t31;
t72 = -t106 * t85 + t109 * t97;
t73 = t106 * t97 + t109 * t85;
t56 = -t72 * pkin(4) - t73 * pkin(10);
t67 = qJDD(4) - t68;
t83 = qJD(4) + t84;
t82 = t83 ^ 2;
t22 = -t67 * pkin(4) - t82 * pkin(10) + t73 * t56 - t122;
t105 = sin(qJ(5));
t147 = cos(qJ(5));
t53 = t72 * qJD(4) + t106 * t96 + t109 * t69;
t59 = t105 * t83 + t147 * t73;
t33 = t59 * qJD(5) + t105 * t53 - t147 * t67;
t58 = t105 * t73 - t147 * t83;
t34 = -t58 * qJD(5) + t105 * t67 + t147 * t53;
t71 = qJD(5) - t72;
t45 = -t58 * mrSges(7,2) + t71 * mrSges(7,3);
t130 = m(7) * (-0.2e1 * qJD(6) * t59 + (t58 * t71 - t34) * qJ(6) + (t59 * t71 + t33) * pkin(5) + t22) + t33 * mrSges(7,1) + t58 * t45;
t46 = -t71 * mrSges(6,2) - t58 * mrSges(6,3);
t47 = t71 * mrSges(6,1) - t59 * mrSges(6,3);
t48 = -t71 * mrSges(7,1) + t59 * mrSges(7,2);
t149 = m(6) * t22 + t33 * mrSges(6,1) + (t47 - t48) * t59 + (mrSges(6,2) - mrSges(7,3)) * t34 + t58 * t46 + t130;
t142 = t106 * t31 + t109 * t29;
t23 = -t82 * pkin(4) + t67 * pkin(10) + t72 * t56 + t142;
t28 = -t96 * pkin(3) - t95 * pkin(9) + t85 * t65 - t150;
t52 = -t73 * qJD(4) - t106 * t69 + t109 * t96;
t25 = (-t72 * t83 - t53) * pkin(10) + (t73 * t83 - t52) * pkin(4) + t28;
t117 = -t105 * t23 + t147 * t25;
t37 = t58 * pkin(5) - t59 * qJ(6);
t51 = qJDD(5) - t52;
t70 = t71 ^ 2;
t148 = m(7) * (-t51 * pkin(5) - t70 * qJ(6) + t59 * t37 + qJDD(6) - t117);
t144 = -mrSges(6,3) - mrSges(7,2);
t143 = t105 * t25 + t147 * t23;
t38 = t58 * mrSges(7,1) - t59 * mrSges(7,3);
t141 = -t58 * mrSges(6,1) - t59 * mrSges(6,2) - t38;
t135 = t102 * t110;
t131 = m(7) * (-t70 * pkin(5) + t51 * qJ(6) + 0.2e1 * qJD(6) * t71 - t58 * t37 + t143) + t71 * t48 + t51 * mrSges(7,3);
t15 = m(6) * t143 - t51 * mrSges(6,2) + t141 * t58 + t144 * t33 - t71 * t47 + t131;
t16 = m(6) * t117 - t148 + (t46 + t45) * t71 + t141 * t59 + (mrSges(6,1) + mrSges(7,1)) * t51 + t144 * t34;
t60 = -t83 * mrSges(5,2) + t72 * mrSges(5,3);
t61 = t83 * mrSges(5,1) - t73 * mrSges(5,3);
t113 = m(5) * t28 - t52 * mrSges(5,1) + t53 * mrSges(5,2) + t105 * t15 + t147 * t16 - t72 * t60 + t73 * t61;
t63 = t84 * mrSges(4,1) + t85 * mrSges(4,2);
t74 = -t97 * mrSges(4,2) - t84 * mrSges(4,3);
t10 = m(4) * t150 + t96 * mrSges(4,1) - t69 * mrSges(4,3) - t85 * t63 + t97 * t74 - t113;
t55 = -t72 * mrSges(5,1) + t73 * mrSges(5,2);
t12 = m(5) * t142 - t67 * mrSges(5,2) + t52 * mrSges(5,3) - t105 * t16 + t147 * t15 + t72 * t55 - t83 * t61;
t14 = m(5) * t122 + t67 * mrSges(5,1) - t53 * mrSges(5,3) - t73 * t55 + t83 * t60 - t149;
t75 = t97 * mrSges(4,1) - t85 * mrSges(4,3);
t7 = m(4) * t129 - t96 * mrSges(4,2) + t68 * mrSges(4,3) - t106 * t14 + t109 * t12 - t84 * t63 - t97 * t75;
t126 = t110 * t134;
t89 = -t97 * mrSges(3,2) + mrSges(3,3) * t126;
t92 = (-mrSges(3,1) * t110 + mrSges(3,2) * t107) * t134;
t5 = m(3) * (-g(3) * t135 + t123) - t93 * mrSges(3,3) + t96 * mrSges(3,1) - t92 * t127 + t97 * t89 + t101 * t7 + t103 * t10;
t88 = t97 * mrSges(3,1) - mrSges(3,3) * t127;
t6 = m(3) * t119 - t96 * mrSges(3,2) + t94 * mrSges(3,3) - t101 * t10 + t103 * t7 + t92 * t126 - t97 * t88;
t116 = m(4) * t114 - t68 * mrSges(4,1) + t69 * mrSges(4,2) + t106 * t12 + t109 * t14 + t84 * t74 + t85 * t75;
t9 = m(3) * t120 + t93 * mrSges(3,2) - t94 * mrSges(3,1) + (t107 * t88 - t110 * t89) * t134 + t116;
t133 = t104 * t9 + t5 * t135 + t6 * t136;
t2 = m(2) * t121 - t112 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t107 * t5 + t110 * t6;
t1 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t112 * mrSges(2,2) - t102 * t9 + (t107 * t6 + t110 * t5) * t104;
t3 = [-m(1) * g(1) - t108 * t1 + t111 * t2, t2, t6, t7, t12, t15, -t33 * mrSges(7,2) - t58 * t38 + t131; -m(1) * g(2) + t111 * t1 + t108 * t2, t1, t5, t10, t14, t16, -t34 * mrSges(7,3) - t59 * t48 + t130; (-m(1) - m(2)) * g(3) + t133, -m(2) * g(3) + t133, t9, t116, t113, t149, -t51 * mrSges(7,1) + t34 * mrSges(7,2) + t59 * t38 - t71 * t45 + t148;];
f_new  = t3;
