% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 23:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:19:01
% EndTime: 2019-05-05 23:19:30
% DurationCPUTime: 18.91s
% Computational Cost: add. (300581->220), mult. (943748->322), div. (0->0), fcn. (811565->16), ass. (0->127)
t103 = sin(pkin(12));
t105 = sin(pkin(6));
t112 = sin(qJ(3));
t116 = cos(qJ(3));
t107 = cos(pkin(12));
t108 = cos(pkin(7));
t151 = t107 * t108;
t104 = sin(pkin(7));
t109 = cos(pkin(6));
t154 = t104 * t109;
t121 = (-t103 * t112 + t116 * t151) * t105 + t116 * t154;
t83 = t121 * qJD(1);
t150 = t108 * t112;
t153 = t104 * t112;
t123 = t109 * t153 + (t103 * t116 + t107 * t150) * t105;
t84 = t123 * qJD(1);
t74 = -t84 * qJD(3) + t121 * qJDD(1);
t131 = -pkin(9) * t103 * t104 - pkin(2) * t107;
t148 = qJD(1) * t105;
t157 = pkin(9) * qJDD(1);
t128 = qJD(1) * t131 * t148 + t108 * t157;
t142 = qJD(2) * t148;
t152 = t105 * t107;
t118 = qJD(1) ^ 2;
t113 = sin(qJ(1));
t117 = cos(qJ(1));
t141 = t113 * g(1) - g(2) * t117;
t156 = qJ(2) * t105;
t94 = qJDD(1) * pkin(1) + t118 * t156 + t141;
t158 = t109 * t94;
t132 = -g(3) * t152 - 0.2e1 * t103 * t142 + t107 * t158;
t91 = (t105 * t151 + t154) * qJD(1) * pkin(9);
t137 = -g(1) * t117 - g(2) * t113;
t95 = -pkin(1) * t118 + qJDD(1) * t156 + t137;
t57 = (pkin(2) * qJDD(1) + qJD(1) * t91) * t109 + (-t128 * t105 - t95) * t103 + t132;
t144 = t103 * t158 + (0.2e1 * t142 + t95) * t107;
t155 = t103 * t105;
t96 = (-pkin(9) * t108 * t155 + pkin(2) * t109) * qJD(1);
t58 = (-qJD(1) * t96 + t104 * t157) * t109 + (-g(3) * t103 + t128 * t107) * t105 + t144;
t138 = -t109 * g(3) + qJDD(2);
t67 = (-t94 + t131 * qJDD(1) + (t103 * t96 - t107 * t91) * qJD(1)) * t105 + t138;
t161 = -t112 * t58 + (t104 * t67 + t108 * t57) * t116;
t160 = 2 * qJD(5);
t111 = sin(qJ(4));
t115 = cos(qJ(4));
t145 = t116 * t58 + t57 * t150 + t67 * t153;
t73 = -pkin(3) * t83 - pkin(10) * t84;
t126 = -t104 * t152 + t108 * t109;
t92 = t126 * qJD(1) + qJD(3);
t88 = t92 ^ 2;
t89 = t126 * qJDD(1) + qJDD(3);
t33 = -pkin(3) * t88 + pkin(10) * t89 + t73 * t83 + t145;
t140 = -t104 * t57 + t108 * t67;
t75 = t83 * qJD(3) + t123 * qJDD(1);
t36 = (-t83 * t92 - t75) * pkin(10) + (t84 * t92 - t74) * pkin(3) + t140;
t159 = t111 * t36 + t115 * t33;
t102 = sin(pkin(13));
t106 = cos(pkin(13));
t110 = sin(qJ(6));
t114 = cos(qJ(6));
t139 = -t111 * t33 + t115 * t36;
t77 = -t111 * t84 + t115 * t92;
t53 = qJD(4) * t77 + t111 * t89 + t115 * t75;
t71 = qJDD(4) - t74;
t78 = t111 * t92 + t115 * t84;
t82 = qJD(4) - t83;
t24 = (t77 * t82 - t53) * qJ(5) + (t77 * t78 + t71) * pkin(4) + t139;
t52 = -qJD(4) * t78 - t111 * t75 + t115 * t89;
t69 = pkin(4) * t82 - qJ(5) * t78;
t76 = t77 ^ 2;
t26 = -pkin(4) * t76 + t52 * qJ(5) - t69 * t82 + t159;
t61 = -t102 * t78 + t106 * t77;
t146 = t102 * t24 + t106 * t26 + t61 * t160;
t62 = t102 * t77 + t106 * t78;
t46 = -pkin(5) * t61 - pkin(11) * t62;
t81 = t82 ^ 2;
t21 = -pkin(5) * t81 + pkin(11) * t71 + t61 * t46 + t146;
t32 = -t89 * pkin(3) - t88 * pkin(10) + t84 * t73 - t161;
t120 = -t52 * pkin(4) - t76 * qJ(5) + t78 * t69 + qJDD(5) + t32;
t40 = -t102 * t53 + t106 * t52;
t41 = t102 * t52 + t106 * t53;
t22 = (-t61 * t82 - t41) * pkin(11) + (t62 * t82 - t40) * pkin(5) + t120;
t47 = -t110 * t62 + t114 * t82;
t30 = t47 * qJD(6) + t110 * t71 + t114 * t41;
t48 = t110 * t82 + t114 * t62;
t37 = -mrSges(7,1) * t47 + mrSges(7,2) * t48;
t39 = qJDD(6) - t40;
t60 = qJD(6) - t61;
t42 = -mrSges(7,2) * t60 + mrSges(7,3) * t47;
t18 = m(7) * (-t110 * t21 + t114 * t22) - t30 * mrSges(7,3) + t39 * mrSges(7,1) - t48 * t37 + t60 * t42;
t29 = -t48 * qJD(6) - t110 * t41 + t114 * t71;
t43 = mrSges(7,1) * t60 - mrSges(7,3) * t48;
t19 = m(7) * (t110 * t22 + t114 * t21) + t29 * mrSges(7,3) - t39 * mrSges(7,2) + t47 * t37 - t60 * t43;
t45 = -mrSges(6,1) * t61 + mrSges(6,2) * t62;
t50 = mrSges(6,1) * t82 - t62 * mrSges(6,3);
t13 = m(6) * t146 - t71 * mrSges(6,2) + t40 * mrSges(6,3) - t110 * t18 + t114 * t19 + t61 * t45 - t82 * t50;
t134 = t102 * t26 - t106 * t24;
t124 = m(7) * (-t71 * pkin(5) - t81 * pkin(11) + (t160 + t46) * t62 + t134) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t47 * t42 + t48 * t43;
t49 = -mrSges(6,2) * t82 + t61 * mrSges(6,3);
t15 = m(6) * (-0.2e1 * qJD(5) * t62 - t134) - t41 * mrSges(6,3) + t71 * mrSges(6,1) - t62 * t45 + t82 * t49 - t124;
t63 = -mrSges(5,1) * t77 + mrSges(5,2) * t78;
t68 = -mrSges(5,2) * t82 + mrSges(5,3) * t77;
t11 = m(5) * t139 + t71 * mrSges(5,1) - t53 * mrSges(5,3) + t102 * t13 + t106 * t15 - t78 * t63 + t82 * t68;
t70 = mrSges(5,1) * t82 - mrSges(5,3) * t78;
t12 = m(5) * t159 - t71 * mrSges(5,2) + t52 * mrSges(5,3) - t102 * t15 + t106 * t13 + t77 * t63 - t82 * t70;
t79 = -mrSges(4,2) * t92 + mrSges(4,3) * t83;
t80 = mrSges(4,1) * t92 - mrSges(4,3) * t84;
t10 = m(4) * t140 - t74 * mrSges(4,1) + t75 * mrSges(4,2) + t115 * t11 + t111 * t12 - t83 * t79 + t84 * t80;
t130 = mrSges(3,1) * t109 - mrSges(3,3) * t155;
t125 = -m(6) * t120 + t40 * mrSges(6,1) - t41 * mrSges(6,2) - t110 * t19 - t114 * t18 + t61 * t49 - t62 * t50;
t119 = m(5) * t32 - t52 * mrSges(5,1) + t53 * mrSges(5,2) - t77 * t68 + t78 * t70 - t125;
t72 = -mrSges(4,1) * t83 + mrSges(4,2) * t84;
t14 = m(4) * t161 + t89 * mrSges(4,1) - t75 * mrSges(4,3) - t84 * t72 + t92 * t79 - t119;
t9 = m(4) * t145 - t89 * mrSges(4,2) + t74 * mrSges(4,3) - t111 * t11 + t115 * t12 + t83 * t72 - t92 * t80;
t136 = t112 * t9 + t116 * t14;
t135 = -mrSges(3,1) * t107 + mrSges(3,2) * t103;
t93 = t135 * t148;
t129 = -mrSges(3,2) * t109 + mrSges(3,3) * t152;
t98 = t129 * qJD(1);
t4 = m(3) * (-t103 * t95 + t132) - t104 * t10 + t136 * t108 + t130 * qJDD(1) + (t109 * t98 - t93 * t155) * qJD(1);
t97 = t130 * qJD(1);
t6 = m(3) * t138 + t108 * t10 + t136 * t104 + (-m(3) * t94 + t135 * qJDD(1) + (t103 * t97 - t107 * t98) * qJD(1)) * t105;
t8 = m(3) * (-g(3) * t155 + t144) + t116 * t9 - t112 * t14 + t129 * qJDD(1) + (-t109 * t97 + t93 * t152) * qJD(1);
t147 = t109 * t6 + t4 * t152 + t8 * t155;
t2 = m(2) * t137 - t118 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t103 * t4 + t107 * t8;
t1 = m(2) * t141 + qJDD(1) * mrSges(2,1) - t118 * mrSges(2,2) - t105 * t6 + (t103 * t8 + t107 * t4) * t109;
t3 = [-m(1) * g(1) - t1 * t113 + t117 * t2, t2, t8, t9, t12, t13, t19; -m(1) * g(2) + t1 * t117 + t113 * t2, t1, t4, t14, t11, t15, t18; (-m(1) - m(2)) * g(3) + t147, -m(2) * g(3) + t147, t6, t10, t119, -t125, t124;];
f_new  = t3;
