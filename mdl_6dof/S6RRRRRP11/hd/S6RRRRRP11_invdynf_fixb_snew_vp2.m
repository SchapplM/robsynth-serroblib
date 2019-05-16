% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP11
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 07:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:41:25
% EndTime: 2019-05-08 06:41:57
% DurationCPUTime: 11.38s
% Computational Cost: add. (207030->225), mult. (510901->310), div. (0->0), fcn. (429614->14), ass. (0->121)
t100 = sin(pkin(7));
t106 = sin(qJ(3));
t111 = cos(qJ(3));
t102 = cos(pkin(7));
t101 = sin(pkin(6));
t112 = cos(qJ(2));
t107 = sin(qJ(2));
t103 = cos(pkin(6));
t114 = qJD(1) ^ 2;
t108 = sin(qJ(1));
t113 = cos(qJ(1));
t127 = t108 * g(1) - t113 * g(2);
t155 = pkin(9) * t101;
t92 = qJDD(1) * pkin(1) + t114 * t155 + t127;
t144 = t103 * t92;
t122 = -t113 * g(1) - t108 * g(2);
t93 = -t114 * pkin(1) + qJDD(1) * t155 + t122;
t124 = -t107 * t93 + t112 * t144;
t138 = qJD(1) * t107;
t153 = pkin(10) * t102;
t137 = qJD(1) * t112;
t128 = t101 * t137;
t98 = t103 * qJD(1) + qJD(2);
t146 = t100 * t98;
t85 = (t102 * t128 + t146) * pkin(10);
t139 = qJD(1) * t101;
t154 = pkin(10) * t100;
t89 = (-pkin(2) * t112 - t107 * t154) * t139;
t135 = qJD(1) * qJD(2);
t95 = (qJDD(1) * t107 + t112 * t135) * t101;
t97 = t103 * qJDD(1) + qJDD(2);
t54 = -t95 * t153 + t97 * pkin(2) + t98 * t85 + (-g(3) * t112 - t89 * t138) * t101 + t124;
t145 = t102 * t54;
t96 = (qJDD(1) * t112 - t107 * t135) * t101;
t119 = t100 * t97 + t102 * t96;
t147 = t107 * t144 + t112 * t93;
t129 = t101 * t138;
t88 = t98 * pkin(2) - t129 * t153;
t55 = -t98 * t88 + (-g(3) * t107 + t89 * t137) * t101 + t119 * pkin(10) + t147;
t152 = t103 * g(3);
t60 = -t95 * t154 - t96 * pkin(2) - t152 + (-t92 + (t107 * t88 - t112 * t85) * qJD(1)) * t101;
t157 = -t106 * t55 + (t100 * t60 + t145) * t111;
t140 = t102 * t112;
t143 = t100 * t106;
t80 = t98 * t143 + (t106 * t140 + t107 * t111) * t139;
t68 = -t80 * qJD(3) - t106 * t95 + t119 * t111;
t79 = (-t106 * t107 + t111 * t140) * t139 + t111 * t146;
t105 = sin(qJ(4));
t110 = cos(qJ(4));
t131 = t106 * t145 + t111 * t55 + t60 * t143;
t71 = -t79 * pkin(3) - t80 * pkin(11);
t81 = -t100 * t96 + t102 * t97 + qJDD(3);
t86 = -t100 * t128 + t102 * t98 + qJD(3);
t84 = t86 ^ 2;
t29 = -t84 * pkin(3) + t81 * pkin(11) + t79 * t71 + t131;
t126 = -t100 * t54 + t102 * t60;
t69 = t79 * qJD(3) + t119 * t106 + t111 * t95;
t31 = (-t79 * t86 - t69) * pkin(11) + (t80 * t86 - t68) * pkin(3) + t126;
t123 = -t105 * t29 + t110 * t31;
t73 = -t105 * t80 + t110 * t86;
t74 = t105 * t86 + t110 * t80;
t57 = -t73 * pkin(4) - t74 * pkin(12);
t67 = qJDD(4) - t68;
t78 = qJD(4) - t79;
t77 = t78 ^ 2;
t21 = -t67 * pkin(4) - t77 * pkin(12) + t74 * t57 - t123;
t104 = sin(qJ(5));
t109 = cos(qJ(5));
t44 = t73 * qJD(4) + t105 * t81 + t110 * t69;
t64 = t104 * t78 + t109 * t74;
t34 = -t64 * qJD(5) - t104 * t44 + t109 * t67;
t63 = -t104 * t74 + t109 * t78;
t35 = t63 * qJD(5) + t104 * t67 + t109 * t44;
t72 = qJD(5) - t73;
t48 = t72 * pkin(5) - t64 * qJ(6);
t49 = t72 * mrSges(7,1) - t64 * mrSges(7,3);
t62 = t63 ^ 2;
t132 = m(7) * (-t34 * pkin(5) - t62 * qJ(6) + t64 * t48 + qJDD(6) + t21) + t35 * mrSges(7,2) + t64 * t49;
t46 = -t72 * mrSges(7,2) + t63 * mrSges(7,3);
t47 = -t72 * mrSges(6,2) + t63 * mrSges(6,3);
t50 = t72 * mrSges(6,1) - t64 * mrSges(6,3);
t156 = m(6) * t21 + t35 * mrSges(6,2) - (t47 + t46) * t63 - (mrSges(6,1) + mrSges(7,1)) * t34 + t64 * t50 + t132;
t149 = t105 * t31 + t110 * t29;
t22 = -t77 * pkin(4) + t67 * pkin(12) + t73 * t57 + t149;
t28 = -t81 * pkin(3) - t84 * pkin(11) + t80 * t71 - t157;
t43 = -t74 * qJD(4) - t105 * t69 + t110 * t81;
t25 = (-t73 * t78 - t44) * pkin(12) + (t74 * t78 - t43) * pkin(4) + t28;
t150 = t104 * t25 + t109 * t22;
t142 = t101 * t107;
t141 = t101 * t112;
t125 = -t104 * t22 + t109 * t25;
t42 = qJDD(5) - t43;
t134 = m(7) * (-0.2e1 * qJD(6) * t64 + (t63 * t72 - t35) * qJ(6) + (t63 * t64 + t42) * pkin(5) + t125) + t72 * t46 + t42 * mrSges(7,1);
t39 = -t63 * mrSges(7,1) + t64 * mrSges(7,2);
t40 = -t63 * mrSges(6,1) + t64 * mrSges(6,2);
t13 = m(6) * t125 + t42 * mrSges(6,1) + t72 * t47 + (-t40 - t39) * t64 + (-mrSges(6,3) - mrSges(7,3)) * t35 + t134;
t133 = m(7) * (-t62 * pkin(5) + t34 * qJ(6) + 0.2e1 * qJD(6) * t63 - t72 * t48 + t150) + t34 * mrSges(7,3) + t63 * t39;
t14 = m(6) * t150 + t34 * mrSges(6,3) + t63 * t40 + (-t50 - t49) * t72 + (-mrSges(6,2) - mrSges(7,2)) * t42 + t133;
t56 = -t73 * mrSges(5,1) + t74 * mrSges(5,2);
t66 = t78 * mrSges(5,1) - t74 * mrSges(5,3);
t12 = m(5) * t149 - t67 * mrSges(5,2) + t43 * mrSges(5,3) - t104 * t13 + t109 * t14 + t73 * t56 - t78 * t66;
t65 = -t78 * mrSges(5,2) + t73 * mrSges(5,3);
t15 = m(5) * t123 + t67 * mrSges(5,1) - t44 * mrSges(5,3) - t74 * t56 + t78 * t65 - t156;
t75 = -t86 * mrSges(4,2) + t79 * mrSges(4,3);
t76 = t86 * mrSges(4,1) - t80 * mrSges(4,3);
t10 = m(4) * t126 - t68 * mrSges(4,1) + t69 * mrSges(4,2) + t105 * t12 + t110 * t15 - t79 * t75 + t80 * t76;
t115 = m(5) * t28 - t43 * mrSges(5,1) + t44 * mrSges(5,2) + t104 * t14 + t109 * t13 - t73 * t65 + t74 * t66;
t70 = -t79 * mrSges(4,1) + t80 * mrSges(4,2);
t11 = m(4) * t157 + t81 * mrSges(4,1) - t69 * mrSges(4,3) - t80 * t70 + t86 * t75 - t115;
t9 = m(4) * t131 - t81 * mrSges(4,2) + t68 * mrSges(4,3) - t105 * t15 + t110 * t12 + t79 * t70 - t86 * t76;
t121 = t106 * t9 + t111 * t11;
t91 = -t98 * mrSges(3,2) + mrSges(3,3) * t128;
t94 = (-mrSges(3,1) * t112 + mrSges(3,2) * t107) * t139;
t4 = m(3) * (-g(3) * t141 + t124) - t95 * mrSges(3,3) + t97 * mrSges(3,1) - t94 * t129 + t98 * t91 - t100 * t10 + t121 * t102;
t90 = t98 * mrSges(3,1) - mrSges(3,3) * t129;
t6 = m(3) * (-t101 * t92 - t152) + t95 * mrSges(3,2) - t96 * mrSges(3,1) + t102 * t10 + t121 * t100 + (t107 * t90 - t112 * t91) * t139;
t8 = m(3) * (-g(3) * t142 + t147) + t96 * mrSges(3,3) - t97 * mrSges(3,2) + t94 * t128 - t98 * t90 + t111 * t9 - t106 * t11;
t136 = t103 * t6 + t4 * t141 + t8 * t142;
t2 = m(2) * t122 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t107 * t4 + t112 * t8;
t1 = m(2) * t127 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t101 * t6 + (t107 * t8 + t112 * t4) * t103;
t3 = [-m(1) * g(1) - t108 * t1 + t113 * t2, t2, t8, t9, t12, t14, -t42 * mrSges(7,2) - t72 * t49 + t133; -m(1) * g(2) + t113 * t1 + t108 * t2, t1, t4, t11, t15, t13, -t35 * mrSges(7,3) - t64 * t39 + t134; (-m(1) - m(2)) * g(3) + t136, -m(2) * g(3) + t136, t6, t10, t115, t156, -t34 * mrSges(7,1) - t63 * t46 + t132;];
f_new  = t3;
