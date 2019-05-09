% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-06 00:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR12_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:41:16
% EndTime: 2019-05-06 00:41:30
% DurationCPUTime: 7.55s
% Computational Cost: add. (115643->219), mult. (360952->306), div. (0->0), fcn. (305349->14), ass. (0->123)
t101 = sin(pkin(12));
t109 = sin(qJ(3));
t112 = cos(qJ(3));
t104 = cos(pkin(12));
t105 = cos(pkin(7));
t147 = t104 * t105;
t167 = -t101 * t109 + t112 * t147;
t103 = sin(pkin(6));
t102 = sin(pkin(7));
t106 = cos(pkin(6));
t150 = t102 * t106;
t140 = t112 * t150;
t146 = t105 * t109;
t149 = t102 * t109;
t117 = t106 * t149 + (t101 * t112 + t104 * t146) * t103;
t79 = t117 * qJD(1);
t67 = -t79 * qJD(3) + (t167 * t103 + t140) * qJDD(1);
t128 = -pkin(9) * t101 * t102 - pkin(2) * t104;
t144 = qJD(1) * t103;
t154 = pkin(9) * qJDD(1);
t124 = qJD(1) * t128 * t144 + t105 * t154;
t138 = qJD(2) * t144;
t148 = t103 * t104;
t114 = qJD(1) ^ 2;
t110 = sin(qJ(1));
t113 = cos(qJ(1));
t137 = t110 * g(1) - g(2) * t113;
t153 = qJ(2) * t103;
t91 = qJDD(1) * pkin(1) + t114 * t153 + t137;
t155 = t106 * t91;
t129 = -g(3) * t148 - 0.2e1 * t101 * t138 + t104 * t155;
t88 = (t103 * t147 + t150) * qJD(1) * pkin(9);
t133 = -g(1) * t113 - g(2) * t110;
t92 = -pkin(1) * t114 + qJDD(1) * t153 + t133;
t47 = (pkin(2) * qJDD(1) + qJD(1) * t88) * t106 + (-t124 * t103 - t92) * t101 + t129;
t141 = t101 * t155 + (0.2e1 * t138 + t92) * t104;
t152 = t101 * t103;
t93 = (-pkin(9) * t105 * t152 + pkin(2) * t106) * qJD(1);
t48 = (-qJD(1) * t93 + t102 * t154) * t106 + (-g(3) * t101 + t124 * t104) * t103 + t141;
t135 = -t106 * g(3) + qJDD(2);
t57 = (-t91 + t128 * qJDD(1) + (t101 * t93 - t104 * t88) * qJD(1)) * t103 + t135;
t166 = -t109 * t48 + (t102 * t57 + t105 * t47) * t112;
t107 = sin(qJ(6));
t111 = cos(qJ(6));
t108 = sin(qJ(4));
t162 = cos(qJ(4));
t122 = -t102 * t148 + t105 * t106;
t89 = t122 * qJD(1) + qJD(3);
t72 = t108 * t79 - t162 * t89;
t78 = qJD(1) * t140 + t167 * t144;
t77 = qJD(4) - t78;
t161 = t72 * t77;
t163 = -2 * qJD(5);
t66 = -pkin(3) * t78 - pkin(10) * t79;
t85 = t89 ^ 2;
t86 = t122 * qJDD(1) + qJDD(3);
t27 = -t86 * pkin(3) - t85 * pkin(10) + t79 * t66 - t166;
t68 = t78 * qJD(3) + t117 * qJDD(1);
t43 = -t72 * qJD(4) + t108 * t86 + t162 * t68;
t73 = t108 * t89 + t162 * t79;
t115 = (-t43 + t161) * qJ(5) + t27 + (t77 * pkin(4) + t163) * t73;
t142 = t112 * t48 + t47 * t146 + t57 * t149;
t28 = -pkin(3) * t85 + pkin(10) * t86 + t66 * t78 + t142;
t136 = -t102 * t47 + t105 * t57;
t30 = (-t78 * t89 - t68) * pkin(10) + (t79 * t89 - t67) * pkin(3) + t136;
t134 = -t108 * t28 + t162 * t30;
t49 = pkin(4) * t72 - qJ(5) * t73;
t64 = qJDD(4) - t67;
t76 = t77 ^ 2;
t22 = -t64 * pkin(4) - t76 * qJ(5) + t73 * t49 + qJDD(5) - t134;
t17 = (t72 * t73 - t64) * pkin(11) + (t43 + t161) * pkin(5) + t22;
t42 = qJD(4) * t73 + t108 * t68 - t162 * t86;
t62 = pkin(5) * t73 - pkin(11) * t77;
t71 = t72 ^ 2;
t20 = -pkin(5) * t71 - t62 * t73 + (pkin(4) + pkin(11)) * t42 + t115;
t55 = -t107 * t77 + t111 * t72;
t33 = t55 * qJD(6) + t107 * t42 + t111 * t64;
t56 = t107 * t72 + t111 * t77;
t35 = -mrSges(7,1) * t55 + mrSges(7,2) * t56;
t70 = qJD(6) + t73;
t37 = -mrSges(7,2) * t70 + t55 * mrSges(7,3);
t40 = qJDD(6) + t43;
t15 = m(7) * (-t107 * t20 + t111 * t17) - t33 * mrSges(7,3) + t40 * mrSges(7,1) - t56 * t35 + t70 * t37;
t32 = -t56 * qJD(6) - t107 * t64 + t111 * t42;
t38 = mrSges(7,1) * t70 - t56 * mrSges(7,3);
t16 = m(7) * (t107 * t17 + t111 * t20) + t32 * mrSges(7,3) - t40 * mrSges(7,2) + t55 * t35 - t70 * t38;
t59 = mrSges(6,1) * t73 + mrSges(6,2) * t77;
t125 = t107 * t15 - t111 * t16 - m(6) * (t42 * pkin(4) + t115) + t43 * mrSges(6,3) + t73 * t59;
t58 = mrSges(6,1) * t72 - mrSges(6,3) * t77;
t156 = -mrSges(5,2) * t77 - mrSges(5,3) * t72 - t58;
t160 = mrSges(5,1) - mrSges(6,2);
t61 = mrSges(5,1) * t77 - mrSges(5,3) * t73;
t164 = m(5) * t27 + t43 * mrSges(5,2) + t156 * t72 + t160 * t42 + t73 * t61 - t125;
t159 = -mrSges(5,3) - mrSges(6,1);
t158 = t108 * t30 + t162 * t28;
t51 = -mrSges(6,2) * t72 - mrSges(6,3) * t73;
t157 = -mrSges(5,1) * t72 - mrSges(5,2) * t73 - t51;
t121 = -m(6) * t22 - t107 * t16 - t111 * t15;
t12 = m(5) * t134 + t156 * t77 + t157 * t73 + t159 * t43 + t160 * t64 + t121;
t119 = -pkin(4) * t76 + t64 * qJ(5) - t49 * t72 + t158;
t120 = -t32 * mrSges(7,1) - t55 * t37 + m(7) * (-t42 * pkin(5) - pkin(11) * t71 + ((2 * qJD(5)) + t62) * t77 + t119) + t33 * mrSges(7,2) + t56 * t38;
t118 = -m(6) * (t77 * t163 - t119) + t120;
t13 = m(5) * t158 + (-t61 + t59) * t77 + t157 * t72 + (-mrSges(5,2) + mrSges(6,3)) * t64 + t159 * t42 + t118;
t74 = -mrSges(4,2) * t89 + mrSges(4,3) * t78;
t75 = mrSges(4,1) * t89 - mrSges(4,3) * t79;
t10 = m(4) * t136 - t67 * mrSges(4,1) + t68 * mrSges(4,2) + t108 * t13 + t162 * t12 - t78 * t74 + t79 * t75;
t127 = mrSges(3,1) * t106 - mrSges(3,3) * t152;
t65 = -mrSges(4,1) * t78 + mrSges(4,2) * t79;
t11 = m(4) * t166 + t86 * mrSges(4,1) - t68 * mrSges(4,3) - t79 * t65 + t89 * t74 - t164;
t9 = m(4) * t142 - t86 * mrSges(4,2) + t67 * mrSges(4,3) - t108 * t12 + t162 * t13 + t78 * t65 - t89 * t75;
t132 = t109 * t9 + t11 * t112;
t131 = -mrSges(3,1) * t104 + mrSges(3,2) * t101;
t90 = t131 * t144;
t126 = -mrSges(3,2) * t106 + mrSges(3,3) * t148;
t95 = t126 * qJD(1);
t4 = m(3) * (-t101 * t92 + t129) - t102 * t10 + t132 * t105 + t127 * qJDD(1) + (t106 * t95 - t90 * t152) * qJD(1);
t94 = t127 * qJD(1);
t6 = m(3) * t135 + t105 * t10 + t132 * t102 + (-m(3) * t91 + t131 * qJDD(1) + (t101 * t94 - t104 * t95) * qJD(1)) * t103;
t8 = m(3) * (-g(3) * t152 + t141) + t112 * t9 - t109 * t11 + t126 * qJDD(1) + (-t106 * t94 + t90 * t148) * qJD(1);
t143 = t106 * t6 + t4 * t148 + t8 * t152;
t2 = m(2) * t133 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t101 * t4 + t104 * t8;
t1 = m(2) * t137 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t103 * t6 + (t101 * t8 + t104 * t4) * t106;
t3 = [-m(1) * g(1) - t1 * t110 + t113 * t2, t2, t8, t9, t13, -t42 * mrSges(6,2) - t72 * t58 - t125, t16; -m(1) * g(2) + t1 * t113 + t110 * t2, t1, t4, t11, t12, t42 * mrSges(6,1) - t64 * mrSges(6,3) + t72 * t51 - t77 * t59 - t118, t15; (-m(1) - m(2)) * g(3) + t143, -m(2) * g(3) + t143, t6, t10, t164, t43 * mrSges(6,1) + t64 * mrSges(6,2) + t73 * t51 + t77 * t58 - t121, t120;];
f_new  = t3;
