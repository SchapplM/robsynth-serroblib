% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP12
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
% Datum: 2019-05-08 07:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP12_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP12_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 07:17:52
% EndTime: 2019-05-08 07:18:23
% DurationCPUTime: 11.10s
% Computational Cost: add. (203736->225), mult. (502473->310), div. (0->0), fcn. (421752->14), ass. (0->122)
t105 = sin(qJ(3));
t109 = cos(qJ(3));
t101 = cos(pkin(7));
t100 = sin(pkin(6));
t110 = cos(qJ(2));
t106 = sin(qJ(2));
t102 = cos(pkin(6));
t112 = qJD(1) ^ 2;
t107 = sin(qJ(1));
t111 = cos(qJ(1));
t125 = t107 * g(1) - t111 * g(2);
t153 = pkin(9) * t100;
t91 = qJDD(1) * pkin(1) + t112 * t153 + t125;
t141 = t102 * t91;
t121 = -t111 * g(1) - t107 * g(2);
t92 = -t112 * pkin(1) + qJDD(1) * t153 + t121;
t123 = -t106 * t92 + t110 * t141;
t135 = qJD(1) * t106;
t152 = pkin(10) * t101;
t134 = qJD(1) * t110;
t126 = t100 * t134;
t97 = t102 * qJD(1) + qJD(2);
t99 = sin(pkin(7));
t150 = t97 * t99;
t84 = (t101 * t126 + t150) * pkin(10);
t136 = qJD(1) * t100;
t155 = pkin(10) * t99;
t88 = (-pkin(2) * t110 - t106 * t155) * t136;
t132 = qJD(1) * qJD(2);
t94 = (qJDD(1) * t106 + t110 * t132) * t100;
t96 = t102 * qJDD(1) + qJDD(2);
t51 = -t94 * t152 + t96 * pkin(2) + t97 * t84 + (-g(3) * t110 - t88 * t135) * t100 + t123;
t142 = t101 * t51;
t95 = (qJDD(1) * t110 - t106 * t132) * t100;
t119 = t101 * t95 + t96 * t99;
t143 = t106 * t141 + t110 * t92;
t127 = t100 * t135;
t87 = t97 * pkin(2) - t127 * t152;
t52 = -t97 * t87 + (-g(3) * t106 + t88 * t134) * t100 + t119 * pkin(10) + t143;
t151 = t102 * g(3);
t57 = -t94 * t155 - t95 * pkin(2) - t151 + (-t91 + (t106 * t87 - t110 * t84) * qJD(1)) * t100;
t158 = -t105 * t52 + (t57 * t99 + t142) * t109;
t137 = t101 * t110;
t140 = t105 * t99;
t79 = t97 * t140 + (t105 * t137 + t106 * t109) * t136;
t65 = -t79 * qJD(3) - t105 * t94 + t119 * t109;
t78 = (-t105 * t106 + t109 * t137) * t136 + t109 * t150;
t104 = sin(qJ(4));
t108 = cos(qJ(4));
t129 = t105 * t142 + t109 * t52 + t57 * t140;
t68 = -t78 * pkin(3) - t79 * pkin(11);
t80 = t101 * t96 - t99 * t95 + qJDD(3);
t85 = t101 * t97 - t99 * t126 + qJD(3);
t83 = t85 ^ 2;
t28 = -t83 * pkin(3) + t80 * pkin(11) + t78 * t68 + t129;
t124 = t101 * t57 - t99 * t51;
t66 = t78 * qJD(3) + t119 * t105 + t109 * t94;
t30 = (-t78 * t85 - t66) * pkin(11) + (t79 * t85 - t65) * pkin(3) + t124;
t122 = -t104 * t28 + t108 * t30;
t71 = -t104 * t79 + t108 * t85;
t72 = t104 * t85 + t108 * t79;
t54 = -t71 * pkin(4) - t72 * pkin(12);
t64 = qJDD(4) - t65;
t77 = qJD(4) - t78;
t76 = t77 ^ 2;
t21 = -t64 * pkin(4) - t76 * pkin(12) + t72 * t54 - t122;
t103 = sin(qJ(5));
t154 = cos(qJ(5));
t42 = t71 * qJD(4) + t104 * t80 + t108 * t66;
t60 = t103 * t77 + t154 * t72;
t32 = t60 * qJD(5) + t103 * t42 - t154 * t64;
t59 = t103 * t72 - t154 * t77;
t33 = -t59 * qJD(5) + t103 * t64 + t154 * t42;
t70 = qJD(5) - t71;
t44 = -t59 * mrSges(7,2) + t70 * mrSges(7,3);
t130 = m(7) * (-0.2e1 * qJD(6) * t60 + (t59 * t70 - t33) * qJ(6) + (t60 * t70 + t32) * pkin(5) + t21) + t32 * mrSges(7,1) + t59 * t44;
t45 = -t70 * mrSges(6,2) - t59 * mrSges(6,3);
t46 = t70 * mrSges(6,1) - t60 * mrSges(6,3);
t47 = -t70 * mrSges(7,1) + t60 * mrSges(7,2);
t157 = m(6) * t21 + t32 * mrSges(6,1) + (t46 - t47) * t60 + (mrSges(6,2) - mrSges(7,3)) * t33 + t59 * t45 + t130;
t146 = t104 * t30 + t108 * t28;
t22 = -t76 * pkin(4) + t64 * pkin(12) + t71 * t54 + t146;
t27 = -t80 * pkin(3) - t83 * pkin(11) + t79 * t68 - t158;
t41 = -t72 * qJD(4) - t104 * t66 + t108 * t80;
t24 = (-t71 * t77 - t42) * pkin(12) + (t72 * t77 - t41) * pkin(4) + t27;
t116 = -t103 * t22 + t154 * t24;
t36 = t59 * pkin(5) - t60 * qJ(6);
t40 = qJDD(5) - t41;
t69 = t70 ^ 2;
t156 = m(7) * (-t40 * pkin(5) - t69 * qJ(6) + t60 * t36 + qJDD(6) - t116);
t148 = -mrSges(6,3) - mrSges(7,2);
t147 = t103 * t24 + t154 * t22;
t37 = t59 * mrSges(7,1) - t60 * mrSges(7,3);
t145 = -t59 * mrSges(6,1) - t60 * mrSges(6,2) - t37;
t139 = t100 * t106;
t138 = t100 * t110;
t131 = m(7) * (-t69 * pkin(5) + t40 * qJ(6) + 0.2e1 * qJD(6) * t70 - t59 * t36 + t147) + t70 * t47 + t40 * mrSges(7,3);
t14 = m(6) * t147 - t40 * mrSges(6,2) + t145 * t59 + t148 * t32 - t70 * t46 + t131;
t15 = m(6) * t116 - t156 + (t45 + t44) * t70 + t145 * t60 + (mrSges(6,1) + mrSges(7,1)) * t40 + t148 * t33;
t53 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t62 = t77 * mrSges(5,1) - t72 * mrSges(5,3);
t12 = m(5) * t146 - t64 * mrSges(5,2) + t41 * mrSges(5,3) - t103 * t15 + t154 * t14 + t71 * t53 - t77 * t62;
t61 = -t77 * mrSges(5,2) + t71 * mrSges(5,3);
t13 = m(5) * t122 + t64 * mrSges(5,1) - t42 * mrSges(5,3) - t72 * t53 + t77 * t61 - t157;
t73 = -t85 * mrSges(4,2) + t78 * mrSges(4,3);
t74 = t85 * mrSges(4,1) - t79 * mrSges(4,3);
t10 = m(4) * t124 - t65 * mrSges(4,1) + t66 * mrSges(4,2) + t104 * t12 + t108 * t13 - t78 * t73 + t79 * t74;
t113 = m(5) * t27 - t41 * mrSges(5,1) + t42 * mrSges(5,2) + t103 * t14 + t154 * t15 - t71 * t61 + t72 * t62;
t67 = -t78 * mrSges(4,1) + t79 * mrSges(4,2);
t11 = m(4) * t158 + t80 * mrSges(4,1) - t66 * mrSges(4,3) - t79 * t67 + t85 * t73 - t113;
t9 = m(4) * t129 - t80 * mrSges(4,2) + t65 * mrSges(4,3) - t104 * t13 + t108 * t12 + t78 * t67 - t85 * t74;
t118 = t105 * t9 + t109 * t11;
t90 = -t97 * mrSges(3,2) + mrSges(3,3) * t126;
t93 = (-mrSges(3,1) * t110 + mrSges(3,2) * t106) * t136;
t4 = m(3) * (-g(3) * t138 + t123) - t94 * mrSges(3,3) + t96 * mrSges(3,1) - t93 * t127 + t97 * t90 - t99 * t10 + t118 * t101;
t89 = t97 * mrSges(3,1) - mrSges(3,3) * t127;
t6 = m(3) * (-t100 * t91 - t151) + t94 * mrSges(3,2) - t95 * mrSges(3,1) + t101 * t10 + t118 * t99 + (t106 * t89 - t110 * t90) * t136;
t8 = m(3) * (-g(3) * t139 + t143) + t95 * mrSges(3,3) - t96 * mrSges(3,2) + t93 * t126 - t97 * t89 + t109 * t9 - t105 * t11;
t133 = t102 * t6 + t4 * t138 + t8 * t139;
t2 = m(2) * t121 - t112 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t4 + t110 * t8;
t1 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t112 * mrSges(2,2) - t100 * t6 + (t106 * t8 + t110 * t4) * t102;
t3 = [-m(1) * g(1) - t107 * t1 + t111 * t2, t2, t8, t9, t12, t14, -t32 * mrSges(7,2) - t59 * t37 + t131; -m(1) * g(2) + t111 * t1 + t107 * t2, t1, t4, t11, t13, t15, -t33 * mrSges(7,3) - t60 * t47 + t130; (-m(1) - m(2)) * g(3) + t133, -m(2) * g(3) + t133, t6, t10, t113, t157, -t40 * mrSges(7,1) + t33 * mrSges(7,2) + t60 * t37 - t70 * t44 + t156;];
f_new  = t3;
