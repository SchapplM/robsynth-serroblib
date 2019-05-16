% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP5
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
% Datum: 2019-05-06 17:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:45:28
% EndTime: 2019-05-06 17:45:44
% DurationCPUTime: 6.16s
% Computational Cost: add. (71486->209), mult. (187063->281), div. (0->0), fcn. (147030->12), ass. (0->107)
t150 = -2 * qJD(3);
t102 = sin(pkin(11));
t104 = cos(pkin(11));
t108 = sin(qJ(2));
t112 = cos(qJ(2));
t103 = sin(pkin(6));
t137 = qJD(1) * t103;
t85 = (t102 * t108 - t104 * t112) * t137;
t105 = cos(pkin(6));
t114 = qJD(1) ^ 2;
t109 = sin(qJ(1));
t113 = cos(qJ(1));
t127 = t109 * g(1) - t113 * g(2);
t147 = pkin(8) * t103;
t91 = qJDD(1) * pkin(1) + t114 * t147 + t127;
t142 = t105 * t91;
t122 = -t113 * g(1) - t109 * g(2);
t92 = -t114 * pkin(1) + qJDD(1) * t147 + t122;
t124 = -t108 * t92 + t112 * t142;
t140 = t103 ^ 2 * t114;
t135 = qJD(1) * qJD(2);
t94 = (qJDD(1) * t108 + t112 * t135) * t103;
t97 = t105 * qJDD(1) + qJDD(2);
t98 = t105 * qJD(1) + qJD(2);
t46 = t97 * pkin(2) - t94 * qJ(3) + (pkin(2) * t108 * t140 + (qJ(3) * qJD(1) * t98 - g(3)) * t103) * t112 + t124;
t139 = t103 * t108;
t120 = -g(3) * t139 + t108 * t142 + t112 * t92;
t130 = t112 ^ 2 * t140;
t129 = t108 * t137;
t88 = t98 * pkin(2) - qJ(3) * t129;
t95 = (qJDD(1) * t112 - t108 * t135) * t103;
t52 = -pkin(2) * t130 + t95 * qJ(3) - t98 * t88 + t120;
t86 = (t102 * t112 + t104 * t108) * t137;
t149 = -t102 * t52 + t104 * t46 + t86 * t150;
t107 = sin(qJ(4));
t111 = cos(qJ(4));
t131 = t102 * t46 + t104 * t52 + t85 * t150;
t68 = t85 * pkin(3) - t86 * pkin(9);
t96 = t98 ^ 2;
t30 = -t96 * pkin(3) + t97 * pkin(9) - t85 * t68 + t131;
t121 = -t105 * g(3) - t103 * t91;
t116 = -t95 * pkin(2) - qJ(3) * t130 + t88 * t129 + qJDD(3) + t121;
t71 = -t102 * t94 + t104 * t95;
t72 = t102 * t95 + t104 * t94;
t32 = (t85 * t98 - t72) * pkin(9) + (t86 * t98 - t71) * pkin(3) + t116;
t123 = -t107 * t30 + t111 * t32;
t74 = -t107 * t86 + t111 * t98;
t75 = t107 * t98 + t111 * t86;
t59 = -t74 * pkin(4) - t75 * pkin(10);
t70 = qJDD(4) - t71;
t84 = qJD(4) + t85;
t83 = t84 ^ 2;
t22 = -t70 * pkin(4) - t83 * pkin(10) + t75 * t59 - t123;
t106 = sin(qJ(5));
t110 = cos(qJ(5));
t56 = t74 * qJD(4) + t107 * t97 + t111 * t72;
t63 = t106 * t84 + t110 * t75;
t35 = -t63 * qJD(5) - t106 * t56 + t110 * t70;
t62 = -t106 * t75 + t110 * t84;
t36 = t62 * qJD(5) + t106 * t70 + t110 * t56;
t73 = qJD(5) - t74;
t49 = t73 * pkin(5) - t63 * qJ(6);
t50 = t73 * mrSges(7,1) - t63 * mrSges(7,3);
t61 = t62 ^ 2;
t132 = m(7) * (-t35 * pkin(5) - t61 * qJ(6) + t63 * t49 + qJDD(6) + t22) + t36 * mrSges(7,2) + t63 * t50;
t47 = -t73 * mrSges(7,2) + t62 * mrSges(7,3);
t48 = -t73 * mrSges(6,2) + t62 * mrSges(6,3);
t51 = t73 * mrSges(6,1) - t63 * mrSges(6,3);
t148 = m(6) * t22 + t36 * mrSges(6,2) - (t48 + t47) * t62 - (mrSges(6,1) + mrSges(7,1)) * t35 + t63 * t51 + t132;
t144 = t107 * t32 + t111 * t30;
t23 = -t83 * pkin(4) + t70 * pkin(10) + t74 * t59 + t144;
t29 = -t97 * pkin(3) - t96 * pkin(9) + t86 * t68 - t149;
t55 = -t75 * qJD(4) - t107 * t72 + t111 * t97;
t26 = (-t74 * t84 - t56) * pkin(10) + (t75 * t84 - t55) * pkin(4) + t29;
t145 = t106 * t26 + t110 * t23;
t138 = t103 * t112;
t125 = -t106 * t23 + t110 * t26;
t54 = qJDD(5) - t55;
t134 = m(7) * (-0.2e1 * qJD(6) * t63 + (t62 * t73 - t36) * qJ(6) + (t62 * t63 + t54) * pkin(5) + t125) + t73 * t47 + t54 * mrSges(7,1);
t40 = -t62 * mrSges(7,1) + t63 * mrSges(7,2);
t41 = -t62 * mrSges(6,1) + t63 * mrSges(6,2);
t13 = m(6) * t125 + t54 * mrSges(6,1) + t73 * t48 + (-t41 - t40) * t63 + (-mrSges(6,3) - mrSges(7,3)) * t36 + t134;
t133 = m(7) * (-t61 * pkin(5) + t35 * qJ(6) + 0.2e1 * qJD(6) * t62 - t73 * t49 + t145) + t35 * mrSges(7,3) + t62 * t40;
t14 = m(6) * t145 + t35 * mrSges(6,3) + t62 * t41 + (-t51 - t50) * t73 + (-mrSges(6,2) - mrSges(7,2)) * t54 + t133;
t64 = -t84 * mrSges(5,2) + t74 * mrSges(5,3);
t65 = t84 * mrSges(5,1) - t75 * mrSges(5,3);
t115 = m(5) * t29 - t55 * mrSges(5,1) + t56 * mrSges(5,2) + t106 * t14 + t110 * t13 - t74 * t64 + t75 * t65;
t67 = t85 * mrSges(4,1) + t86 * mrSges(4,2);
t76 = -t98 * mrSges(4,2) - t85 * mrSges(4,3);
t10 = m(4) * t149 + t97 * mrSges(4,1) - t72 * mrSges(4,3) - t86 * t67 + t98 * t76 - t115;
t58 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t12 = m(5) * t144 - t70 * mrSges(5,2) + t55 * mrSges(5,3) - t106 * t13 + t110 * t14 + t74 * t58 - t84 * t65;
t16 = m(5) * t123 + t70 * mrSges(5,1) - t56 * mrSges(5,3) - t75 * t58 + t84 * t64 - t148;
t77 = t98 * mrSges(4,1) - t86 * mrSges(4,3);
t7 = m(4) * t131 - t97 * mrSges(4,2) + t71 * mrSges(4,3) - t107 * t16 + t111 * t12 - t85 * t67 - t98 * t77;
t128 = t112 * t137;
t90 = -t98 * mrSges(3,2) + mrSges(3,3) * t128;
t93 = (-mrSges(3,1) * t112 + mrSges(3,2) * t108) * t137;
t5 = m(3) * (-g(3) * t138 + t124) - t94 * mrSges(3,3) + t97 * mrSges(3,1) - t93 * t129 + t98 * t90 + t102 * t7 + t104 * t10;
t89 = t98 * mrSges(3,1) - mrSges(3,3) * t129;
t6 = m(3) * t120 - t97 * mrSges(3,2) + t95 * mrSges(3,3) - t102 * t10 + t104 * t7 + t93 * t128 - t98 * t89;
t118 = m(4) * t116 - t71 * mrSges(4,1) + t72 * mrSges(4,2) + t107 * t12 + t111 * t16 + t85 * t76 + t86 * t77;
t9 = m(3) * t121 + t94 * mrSges(3,2) - t95 * mrSges(3,1) + (t108 * t89 - t112 * t90) * t137 + t118;
t136 = t105 * t9 + t5 * t138 + t6 * t139;
t2 = m(2) * t122 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t108 * t5 + t112 * t6;
t1 = m(2) * t127 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t103 * t9 + (t108 * t6 + t112 * t5) * t105;
t3 = [-m(1) * g(1) - t109 * t1 + t113 * t2, t2, t6, t7, t12, t14, -t54 * mrSges(7,2) - t73 * t50 + t133; -m(1) * g(2) + t113 * t1 + t109 * t2, t1, t5, t10, t16, t13, -t36 * mrSges(7,3) - t63 * t40 + t134; (-m(1) - m(2)) * g(3) + t136, -m(2) * g(3) + t136, t9, t118, t115, t148, -t35 * mrSges(7,1) - t62 * t47 + t132;];
f_new  = t3;
