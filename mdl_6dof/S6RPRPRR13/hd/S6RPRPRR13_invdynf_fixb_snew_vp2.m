% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 21:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR13_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:53:20
% EndTime: 2019-05-05 20:53:34
% DurationCPUTime: 6.77s
% Computational Cost: add. (102865->218), mult. (325887->305), div. (0->0), fcn. (273573->14), ass. (0->124)
t112 = sin(qJ(3));
t161 = cos(pkin(7));
t144 = t112 * t161;
t107 = sin(pkin(7));
t106 = sin(pkin(12));
t108 = sin(pkin(6));
t109 = cos(pkin(12));
t170 = pkin(9) * t107;
t132 = -pkin(2) * t109 - t106 * t170;
t162 = cos(pkin(6));
t135 = -g(3) * t162 + qJDD(2);
t143 = t162 * t107;
t146 = t108 * t161;
t89 = (t109 * t146 + t143) * qJD(1) * pkin(9);
t117 = qJD(1) ^ 2;
t113 = sin(qJ(1));
t116 = cos(qJ(1));
t149 = t113 * g(1) - g(2) * t116;
t160 = qJ(2) * t108;
t93 = qJDD(1) * pkin(1) + t117 * t160 + t149;
t95 = (-pkin(9) * t106 * t146 + pkin(2) * t162) * qJD(1);
t54 = (-t93 + t132 * qJDD(1) + (t106 * t95 - t109 * t89) * qJD(1)) * t108 + t135;
t163 = t107 * t54;
t171 = cos(qJ(3));
t156 = qJD(1) * t108;
t125 = pkin(9) * qJDD(1) * t161 + qJD(1) * t132 * t156;
t145 = t109 * t162;
t150 = qJD(2) * t156;
t157 = t108 * t109;
t134 = -g(3) * t157 - 0.2e1 * t106 * t150 + t93 * t145;
t141 = qJDD(1) * t162;
t142 = qJD(1) * t162;
t139 = -g(1) * t116 - g(2) * t113;
t94 = -pkin(1) * t117 + qJDD(1) * t160 + t139;
t44 = pkin(2) * t141 + t89 * t142 + (-t108 * t125 - t94) * t106 + t134;
t147 = t106 * t162;
t153 = t93 * t147 + (0.2e1 * t150 + t94) * t109;
t45 = t141 * t170 - t95 * t142 + (-g(3) * t106 + t109 * t125) * t108 + t153;
t154 = t112 * t163 + t44 * t144 + t171 * t45;
t140 = t161 * t171;
t158 = t106 * t108;
t172 = t112 * t158 - t140 * t157 - t171 * t143;
t79 = t172 * qJD(1);
t118 = t112 * t143 + (t171 * t106 + t109 * t144) * t108;
t80 = t118 * qJD(1);
t61 = pkin(3) * t79 - qJ(4) * t80;
t173 = -t107 * t157 + t162 * t161;
t90 = -t173 * qJD(1) - qJD(3);
t86 = t90 ^ 2;
t87 = t173 * qJDD(1) + qJDD(3);
t174 = pkin(3) * t86 - t87 * qJ(4) + 0.2e1 * qJD(4) * t90 + t79 * t61 - t154;
t169 = t79 * t90;
t168 = mrSges(4,1) - mrSges(5,2);
t167 = -mrSges(4,3) - mrSges(5,1);
t111 = sin(qJ(5));
t115 = cos(qJ(5));
t120 = -t112 * t45 + t44 * t140 + t163 * t171;
t30 = -t87 * pkin(3) - t86 * qJ(4) + t80 * t61 + qJDD(4) - t120;
t65 = -t79 * qJD(3) + qJDD(1) * t118;
t23 = (t79 * t80 - t87) * pkin(10) + (t65 - t169) * pkin(4) + t30;
t148 = -t107 * t44 + t161 * t54;
t123 = (-t65 - t169) * qJ(4) + t148 + (-t90 * pkin(3) - 0.2e1 * qJD(4)) * t80;
t64 = qJD(3) * t80 + t172 * qJDD(1);
t73 = pkin(4) * t80 + pkin(10) * t90;
t78 = t79 ^ 2;
t27 = -pkin(4) * t78 - t73 * t80 + (pkin(3) + pkin(10)) * t64 + t123;
t166 = t111 * t23 + t115 * t27;
t63 = -mrSges(5,2) * t79 - mrSges(5,3) * t80;
t165 = -mrSges(4,1) * t79 - mrSges(4,2) * t80 - t63;
t70 = mrSges(5,1) * t79 + mrSges(5,3) * t90;
t164 = mrSges(4,2) * t90 - mrSges(4,3) * t79 - t70;
t110 = sin(qJ(6));
t114 = cos(qJ(6));
t67 = t111 * t90 + t115 * t79;
t68 = t111 * t79 - t115 * t90;
t47 = -pkin(5) * t67 - pkin(11) * t68;
t60 = qJDD(5) + t65;
t77 = qJD(5) + t80;
t76 = t77 ^ 2;
t20 = -pkin(5) * t76 + t60 * pkin(11) + t47 * t67 + t166;
t119 = -t64 * pkin(4) - pkin(10) * t78 - t90 * t73 - t174;
t39 = -qJD(5) * t68 - t111 * t87 + t115 * t64;
t40 = qJD(5) * t67 + t111 * t64 + t115 * t87;
t21 = (-t67 * t77 - t40) * pkin(11) + (t68 * t77 - t39) * pkin(5) + t119;
t52 = -t110 * t68 + t114 * t77;
t32 = t52 * qJD(6) + t110 * t60 + t114 * t40;
t53 = t110 * t77 + t114 * t68;
t33 = -mrSges(7,1) * t52 + mrSges(7,2) * t53;
t66 = qJD(6) - t67;
t34 = -mrSges(7,2) * t66 + mrSges(7,3) * t52;
t37 = qJDD(6) - t39;
t17 = m(7) * (-t110 * t20 + t114 * t21) - t32 * mrSges(7,3) + t37 * mrSges(7,1) - t53 * t33 + t66 * t34;
t31 = -t53 * qJD(6) - t110 * t40 + t114 * t60;
t35 = mrSges(7,1) * t66 - mrSges(7,3) * t53;
t18 = m(7) * (t110 * t21 + t114 * t20) + t31 * mrSges(7,3) - t37 * mrSges(7,2) + t52 * t33 - t66 * t35;
t46 = -mrSges(6,1) * t67 + mrSges(6,2) * t68;
t56 = mrSges(6,1) * t77 - mrSges(6,3) * t68;
t13 = m(6) * t166 - t60 * mrSges(6,2) + t39 * mrSges(6,3) - t110 * t17 + t114 * t18 + t67 * t46 - t77 * t56;
t136 = -t111 * t27 + t115 * t23;
t121 = m(7) * (-t60 * pkin(5) - pkin(11) * t76 + t47 * t68 - t136) - t31 * mrSges(7,1) + t32 * mrSges(7,2) - t52 * t34 + t53 * t35;
t55 = -mrSges(6,2) * t77 + mrSges(6,3) * t67;
t14 = m(6) * t136 + t60 * mrSges(6,1) - t40 * mrSges(6,3) - t68 * t46 + t77 * t55 - t121;
t71 = mrSges(5,1) * t80 - mrSges(5,2) * t90;
t131 = -t111 * t14 + t115 * t13 + m(5) * (t64 * pkin(3) + t123) - t80 * t71 - t65 * mrSges(5,3);
t72 = -mrSges(4,1) * t90 - mrSges(4,3) * t80;
t10 = m(4) * t148 + t65 * mrSges(4,2) + t164 * t79 + t168 * t64 + t80 * t72 + t131;
t124 = m(6) * t119 - t39 * mrSges(6,1) + t40 * mrSges(6,2) + t110 * t18 + t114 * t17 - t67 * t55 + t68 * t56;
t122 = -m(5) * t174 + t124;
t11 = m(4) * t154 + (t72 - t71) * t90 + (-mrSges(4,2) + mrSges(5,3)) * t87 + t165 * t79 + t167 * t64 + t122;
t127 = mrSges(3,1) * t162 - mrSges(3,3) * t158;
t128 = -m(5) * t30 - t111 * t13 - t115 * t14;
t9 = m(4) * t120 - t164 * t90 + t165 * t80 + t167 * t65 + t168 * t87 + t128;
t137 = -mrSges(3,1) * t109 + mrSges(3,2) * t106;
t92 = t137 * t156;
t126 = -mrSges(3,2) * t162 + mrSges(3,3) * t157;
t97 = t126 * qJD(1);
t4 = m(3) * (-t106 * t94 + t134) + t11 * t144 + t9 * t140 - t107 * t10 + t127 * qJDD(1) + (-t158 * t92 + t162 * t97) * qJD(1);
t96 = t127 * qJD(1);
t6 = m(3) * t135 + t161 * t10 + (t112 * t11 + t171 * t9) * t107 + (-m(3) * t93 + t137 * qJDD(1) + (t106 * t96 - t109 * t97) * qJD(1)) * t108;
t8 = m(3) * (-g(3) * t158 + t153) + t171 * t11 - t112 * t9 + t126 * qJDD(1) + (t157 * t92 - t162 * t96) * qJD(1);
t155 = t4 * t157 + t8 * t158 + t162 * t6;
t2 = m(2) * t139 - t117 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t4 + t109 * t8;
t1 = m(2) * t149 + qJDD(1) * mrSges(2,1) - t117 * mrSges(2,2) - t108 * t6 + t145 * t4 + t147 * t8;
t3 = [-m(1) * g(1) - t1 * t113 + t116 * t2, t2, t8, t11, -t64 * mrSges(5,2) - t79 * t70 + t131, t13, t18; -m(1) * g(2) + t1 * t116 + t113 * t2, t1, t4, t9, t64 * mrSges(5,1) - t87 * mrSges(5,3) + t79 * t63 + t90 * t71 - t122, t14, t17; (-m(1) - m(2)) * g(3) + t155, -m(2) * g(3) + t155, t6, t10, t65 * mrSges(5,1) + t87 * mrSges(5,2) + t80 * t63 - t90 * t70 - t128, t124, t121;];
f_new  = t3;
