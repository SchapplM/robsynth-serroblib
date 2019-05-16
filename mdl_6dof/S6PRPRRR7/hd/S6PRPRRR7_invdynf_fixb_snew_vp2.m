% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 02:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_invdynf_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:49:32
% EndTime: 2019-05-05 01:49:47
% DurationCPUTime: 14.57s
% Computational Cost: add. (256429->186), mult. (724266->280), div. (0->0), fcn. (618043->18), ass. (0->122)
t105 = sin(qJ(4));
t109 = cos(qJ(4));
t101 = cos(pkin(7));
t100 = cos(pkin(8));
t96 = sin(pkin(7));
t144 = t100 * t96;
t95 = sin(pkin(8));
t98 = cos(pkin(14));
t118 = t101 * t95 + t98 * t144;
t93 = sin(pkin(14));
t149 = t93 * t96;
t116 = -t105 * t149 + t118 * t109;
t67 = t116 * qJD(2);
t114 = t118 * t105 + t109 * t149;
t68 = t114 * qJD(2);
t53 = -t68 * qJD(4) + t116 * qJDD(2);
t139 = qJD(2) * t96;
t131 = qJD(3) * t139;
t111 = qJD(2) ^ 2;
t106 = sin(qJ(2));
t110 = cos(qJ(2));
t97 = sin(pkin(6));
t140 = t110 * t97;
t102 = cos(pkin(6));
t94 = sin(pkin(13));
t99 = cos(pkin(13));
t87 = g(1) * t94 - g(2) * t99;
t142 = t102 * t87;
t88 = -g(1) * t99 - g(2) * t94;
t92 = -g(3) + qJDD(1);
t128 = -t106 * t88 + t110 * t142 + t92 * t140;
t146 = qJ(3) * t96;
t63 = qJDD(2) * pkin(2) + t111 * t146 + t128;
t143 = t101 * t63;
t148 = t96 * t98;
t77 = t102 * t92 - t87 * t97;
t135 = -0.2e1 * t93 * t131 + t98 * t143 + t77 * t148;
t141 = t106 * t97;
t134 = t106 * t142 + t110 * t88 + t92 * t141;
t64 = -pkin(2) * t111 + qJDD(2) * t146 + t134;
t151 = pkin(10) * t93;
t122 = -pkin(3) * t98 - t95 * t151;
t74 = t122 * t139;
t117 = pkin(10) * t118;
t75 = qJD(2) * t117;
t33 = (pkin(3) * qJDD(2) + qJD(2) * t75) * t101 + (-t64 + (-pkin(10) * qJDD(2) * t100 - qJD(2) * t74) * t96) * t93 + t135;
t137 = t101 * t77 + qJDD(3);
t79 = (pkin(3) * t101 - t144 * t151) * qJD(2);
t42 = (-t63 + t122 * qJDD(2) + (-t75 * t98 + t79 * t93) * qJD(2)) * t96 + t137;
t153 = t100 * t33 + t42 * t95;
t129 = t93 * t143 + t77 * t149 + (0.2e1 * t131 + t64) * t98;
t34 = (-t101 * t79 + t74 * t148) * qJD(2) + qJDD(2) * t117 + t129;
t152 = -t105 * t34 + t153 * t109;
t104 = sin(qJ(5));
t108 = cos(qJ(5));
t136 = t153 * t105 + t109 * t34;
t52 = -pkin(4) * t67 - pkin(11) * t68;
t119 = t100 * t101 - t95 * t148;
t76 = t119 * qJD(2) + qJD(4);
t72 = t76 ^ 2;
t73 = t119 * qJDD(2) + qJDD(4);
t24 = -pkin(4) * t72 + pkin(11) * t73 + t52 * t67 + t136;
t130 = t100 * t42 - t95 * t33;
t54 = t67 * qJD(4) + t114 * qJDD(2);
t26 = (-t67 * t76 - t54) * pkin(11) + (t68 * t76 - t53) * pkin(4) + t130;
t147 = t104 * t26 + t108 * t24;
t103 = sin(qJ(6));
t107 = cos(qJ(6));
t56 = -t104 * t68 + t108 * t76;
t57 = t104 * t76 + t108 * t68;
t44 = -pkin(5) * t56 - pkin(12) * t57;
t50 = qJDD(5) - t53;
t66 = qJD(5) - t67;
t65 = t66 ^ 2;
t20 = -pkin(5) * t65 + pkin(12) * t50 + t44 * t56 + t147;
t23 = -t73 * pkin(4) - t72 * pkin(11) + t68 * t52 - t152;
t38 = -qJD(5) * t57 - t104 * t54 + t108 * t73;
t39 = qJD(5) * t56 + t104 * t73 + t108 * t54;
t21 = (-t56 * t66 - t39) * pkin(12) + (t57 * t66 - t38) * pkin(5) + t23;
t46 = -t103 * t57 + t107 * t66;
t28 = qJD(6) * t46 + t103 * t50 + t107 * t39;
t47 = t103 * t66 + t107 * t57;
t32 = -mrSges(7,1) * t46 + mrSges(7,2) * t47;
t55 = qJD(6) - t56;
t35 = -mrSges(7,2) * t55 + mrSges(7,3) * t46;
t37 = qJDD(6) - t38;
t17 = m(7) * (-t103 * t20 + t107 * t21) - t28 * mrSges(7,3) + t37 * mrSges(7,1) - t47 * t32 + t55 * t35;
t27 = -qJD(6) * t47 - t103 * t39 + t107 * t50;
t36 = mrSges(7,1) * t55 - mrSges(7,3) * t47;
t18 = m(7) * (t103 * t21 + t107 * t20) + t27 * mrSges(7,3) - t37 * mrSges(7,2) + t46 * t32 - t55 * t36;
t43 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t49 = mrSges(6,1) * t66 - mrSges(6,3) * t57;
t15 = m(6) * t147 - t50 * mrSges(6,2) + t38 * mrSges(6,3) - t103 * t17 + t107 * t18 + t56 * t43 - t66 * t49;
t124 = -t104 * t24 + t108 * t26;
t113 = m(7) * (-pkin(5) * t50 - pkin(12) * t65 + t44 * t57 - t124) - t27 * mrSges(7,1) + t28 * mrSges(7,2) - t46 * t35 + t47 * t36;
t48 = -mrSges(6,2) * t66 + mrSges(6,3) * t56;
t16 = m(6) * t124 + t50 * mrSges(6,1) - t39 * mrSges(6,3) - t57 * t43 + t66 * t48 - t113;
t51 = -mrSges(5,1) * t67 + mrSges(5,2) * t68;
t59 = mrSges(5,1) * t76 - mrSges(5,3) * t68;
t12 = m(5) * t136 - t73 * mrSges(5,2) + t53 * mrSges(5,3) - t104 * t16 + t108 * t15 + t67 * t51 - t76 * t59;
t112 = m(6) * t23 - t38 * mrSges(6,1) + t39 * mrSges(6,2) + t103 * t18 + t107 * t17 - t56 * t48 + t57 * t49;
t58 = -mrSges(5,2) * t76 + mrSges(5,3) * t67;
t14 = m(5) * t152 + t73 * mrSges(5,1) - t54 * mrSges(5,3) - t68 * t51 + t76 * t58 - t112;
t123 = t105 * t12 + t109 * t14;
t126 = -mrSges(4,1) * t98 + mrSges(4,2) * t93;
t13 = m(5) * t130 - t53 * mrSges(5,1) + t54 * mrSges(5,2) + t104 * t15 + t108 * t16 - t67 * t58 + t68 * t59;
t121 = mrSges(4,1) * t101 - mrSges(4,3) * t149;
t82 = t121 * qJD(2);
t120 = -mrSges(4,2) * t101 + mrSges(4,3) * t148;
t83 = t120 * qJD(2);
t10 = m(4) * t137 + t100 * t13 + t123 * t95 + (-m(4) * t63 + t126 * qJDD(2) + (t82 * t93 - t83 * t98) * qJD(2)) * t96;
t78 = t126 * t139;
t11 = m(4) * t129 + t109 * t12 - t105 * t14 + t120 * qJDD(2) + (-t101 * t82 + t78 * t148) * qJD(2);
t9 = m(4) * (-t93 * t64 + t135) - t95 * t13 + t123 * t100 + t121 * qJDD(2) + (t101 * t83 - t78 * t149) * qJD(2);
t127 = t11 * t93 + t9 * t98;
t4 = m(3) * t128 + qJDD(2) * mrSges(3,1) - t111 * mrSges(3,2) - t96 * t10 + t127 * t101;
t6 = m(3) * t77 + t101 * t10 + t127 * t96;
t8 = m(3) * t134 - t111 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t98 * t11 - t93 * t9;
t132 = m(2) * t92 + t102 * t6 + t4 * t140 + t8 * t141;
t2 = m(2) * t88 - t106 * t4 + t110 * t8;
t1 = m(2) * t87 - t97 * t6 + (t106 * t8 + t110 * t4) * t102;
t3 = [-m(1) * g(1) - t1 * t94 + t2 * t99, t2, t8, t11, t12, t15, t18; -m(1) * g(2) + t1 * t99 + t2 * t94, t1, t4, t9, t14, t16, t17; -m(1) * g(3) + t132, t132, t6, t10, t13, t112, t113;];
f_new  = t3;
