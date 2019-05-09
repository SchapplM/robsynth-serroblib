% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-08 03:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR15_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR15_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR15_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 03:19:16
% EndTime: 2019-05-08 03:19:40
% DurationCPUTime: 8.82s
% Computational Cost: add. (157017->226), mult. (388511->307), div. (0->0), fcn. (324783->14), ass. (0->124)
t106 = cos(pkin(6));
t101 = qJD(1) * t106 + qJD(2);
t103 = sin(pkin(7));
t105 = cos(pkin(7));
t104 = sin(pkin(6));
t114 = cos(qJ(2));
t139 = qJD(1) * t114;
t134 = t104 * t139;
t163 = t101 * t103 + t105 * t134;
t109 = sin(qJ(3));
t113 = cos(qJ(3));
t100 = qJDD(1) * t106 + qJDD(2);
t110 = sin(qJ(2));
t116 = qJD(1) ^ 2;
t111 = sin(qJ(1));
t115 = cos(qJ(1));
t133 = g(1) * t111 - g(2) * t115;
t158 = pkin(9) * t104;
t93 = qJDD(1) * pkin(1) + t116 * t158 + t133;
t147 = t106 * t93;
t128 = -g(1) * t115 - g(2) * t111;
t94 = -pkin(1) * t116 + qJDD(1) * t158 + t128;
t131 = -t110 * t94 + t114 * t147;
t140 = qJD(1) * t110;
t156 = pkin(10) * t105;
t86 = t163 * pkin(10);
t141 = qJD(1) * t104;
t157 = pkin(10) * t103;
t90 = (-pkin(2) * t114 - t110 * t157) * t141;
t137 = qJD(1) * qJD(2);
t96 = (qJDD(1) * t110 + t114 * t137) * t104;
t47 = -t96 * t156 + t100 * pkin(2) + t101 * t86 + (-g(3) * t114 - t140 * t90) * t104 + t131;
t97 = (qJDD(1) * t114 - t110 * t137) * t104;
t124 = t100 * t103 + t105 * t97;
t148 = t110 * t147 + t114 * t94;
t135 = t104 * t140;
t89 = pkin(2) * t101 - t135 * t156;
t48 = -t101 * t89 + (-g(3) * t110 + t139 * t90) * t104 + t124 * pkin(10) + t148;
t155 = t106 * g(3);
t54 = -t96 * t157 - t97 * pkin(2) - t155 + (-t93 + (t110 * t89 - t114 * t86) * qJD(1)) * t104;
t162 = -t109 * t48 + (t103 * t54 + t105 * t47) * t113;
t142 = t105 * t109;
t145 = t103 * t109;
t80 = t101 * t145 + (t110 * t113 + t114 * t142) * t141;
t65 = -t80 * qJD(3) - t109 * t96 + t113 * t124;
t107 = sin(qJ(6));
t112 = cos(qJ(6));
t108 = sin(qJ(4));
t159 = cos(qJ(4));
t87 = t101 * t105 - t103 * t134 + qJD(3);
t72 = t108 * t80 - t159 * t87;
t79 = -t109 * t135 + t113 * t163;
t77 = qJD(4) - t79;
t154 = t72 * t77;
t160 = -2 * qJD(5);
t68 = -pkin(3) * t79 - pkin(11) * t80;
t81 = t100 * t105 - t103 * t97 + qJDD(3);
t85 = t87 ^ 2;
t27 = -t81 * pkin(3) - t85 * pkin(11) + t68 * t80 - t162;
t66 = t79 * qJD(3) + t109 * t124 + t113 * t96;
t39 = -qJD(4) * t72 + t108 * t81 + t159 * t66;
t73 = t108 * t87 + t159 * t80;
t117 = (-t39 + t154) * qJ(5) + t27 + (pkin(4) * t77 + t160) * t73;
t136 = t113 * t48 + t142 * t47 + t145 * t54;
t28 = -pkin(3) * t85 + pkin(11) * t81 + t68 * t79 + t136;
t132 = -t103 * t47 + t105 * t54;
t30 = (-t79 * t87 - t66) * pkin(11) + (t80 * t87 - t65) * pkin(3) + t132;
t129 = -t108 * t28 + t159 * t30;
t49 = pkin(4) * t72 - qJ(5) * t73;
t64 = qJDD(4) - t65;
t76 = t77 ^ 2;
t22 = -t64 * pkin(4) - t76 * qJ(5) + t49 * t73 + qJDD(5) - t129;
t17 = (t72 * t73 - t64) * pkin(12) + (t39 + t154) * pkin(5) + t22;
t38 = qJD(4) * t73 + t108 * t66 - t159 * t81;
t62 = pkin(5) * t73 - pkin(12) * t77;
t71 = t72 ^ 2;
t20 = -t71 * pkin(5) - t73 * t62 + (pkin(4) + pkin(12)) * t38 + t117;
t56 = -t107 * t77 + t112 * t72;
t33 = qJD(6) * t56 + t107 * t38 + t112 * t64;
t57 = t107 * t72 + t112 * t77;
t35 = -mrSges(7,1) * t56 + mrSges(7,2) * t57;
t37 = qJDD(6) + t39;
t70 = qJD(6) + t73;
t41 = -mrSges(7,2) * t70 + mrSges(7,3) * t56;
t15 = m(7) * (-t107 * t20 + t112 * t17) - t33 * mrSges(7,3) + t37 * mrSges(7,1) - t57 * t35 + t70 * t41;
t32 = -qJD(6) * t57 - t107 * t64 + t112 * t38;
t42 = mrSges(7,1) * t70 - mrSges(7,3) * t57;
t16 = m(7) * (t107 * t17 + t112 * t20) + t32 * mrSges(7,3) - t37 * mrSges(7,2) + t56 * t35 - t70 * t42;
t59 = mrSges(6,1) * t73 + mrSges(6,2) * t77;
t123 = t107 * t15 - t112 * t16 - m(6) * (t38 * pkin(4) + t117) + t39 * mrSges(6,3) + t73 * t59;
t58 = mrSges(6,1) * t72 - mrSges(6,3) * t77;
t149 = -mrSges(5,2) * t77 - mrSges(5,3) * t72 - t58;
t153 = mrSges(5,1) - mrSges(6,2);
t61 = mrSges(5,1) * t77 - mrSges(5,3) * t73;
t161 = m(5) * t27 + t39 * mrSges(5,2) + t149 * t72 + t153 * t38 + t73 * t61 - t123;
t152 = -mrSges(5,3) - mrSges(6,1);
t151 = t108 * t30 + t159 * t28;
t51 = -t72 * mrSges(6,2) - mrSges(6,3) * t73;
t150 = -mrSges(5,1) * t72 - t73 * mrSges(5,2) - t51;
t144 = t104 * t110;
t143 = t104 * t114;
t122 = -m(6) * t22 - t107 * t16 - t112 * t15;
t12 = m(5) * t129 + t149 * t77 + t150 * t73 + t152 * t39 + t153 * t64 + t122;
t120 = -t76 * pkin(4) + t64 * qJ(5) - t72 * t49 + t151;
t121 = -t32 * mrSges(7,1) - t56 * t41 + m(7) * (-t38 * pkin(5) - t71 * pkin(12) + ((2 * qJD(5)) + t62) * t77 + t120) + t33 * mrSges(7,2) + t57 * t42;
t119 = -m(6) * (t160 * t77 - t120) + t121;
t13 = m(5) * t151 + (-t61 + t59) * t77 + t150 * t72 + (-mrSges(5,2) + mrSges(6,3)) * t64 + t152 * t38 + t119;
t74 = -mrSges(4,2) * t87 + mrSges(4,3) * t79;
t75 = mrSges(4,1) * t87 - mrSges(4,3) * t80;
t10 = m(4) * t132 - mrSges(4,1) * t65 + mrSges(4,2) * t66 + t108 * t13 + t12 * t159 - t74 * t79 + t75 * t80;
t67 = -mrSges(4,1) * t79 + mrSges(4,2) * t80;
t11 = m(4) * t162 + t81 * mrSges(4,1) - t66 * mrSges(4,3) - t80 * t67 + t87 * t74 - t161;
t9 = m(4) * t136 - mrSges(4,2) * t81 + mrSges(4,3) * t65 - t108 * t12 + t13 * t159 + t67 * t79 - t75 * t87;
t127 = t109 * t9 + t11 * t113;
t92 = -mrSges(3,2) * t101 + mrSges(3,3) * t134;
t95 = (-mrSges(3,1) * t114 + mrSges(3,2) * t110) * t141;
t4 = m(3) * (-g(3) * t143 + t131) - t96 * mrSges(3,3) + t100 * mrSges(3,1) - t95 * t135 + t101 * t92 - t103 * t10 + t127 * t105;
t91 = mrSges(3,1) * t101 - mrSges(3,3) * t135;
t6 = m(3) * (-t104 * t93 - t155) + t96 * mrSges(3,2) - t97 * mrSges(3,1) + t105 * t10 + t127 * t103 + (t110 * t91 - t114 * t92) * t141;
t8 = m(3) * (-g(3) * t144 + t148) + t97 * mrSges(3,3) - t100 * mrSges(3,2) + t95 * t134 - t101 * t91 + t113 * t9 - t109 * t11;
t138 = t106 * t6 + t143 * t4 + t144 * t8;
t2 = m(2) * t128 - mrSges(2,1) * t116 - qJDD(1) * mrSges(2,2) - t110 * t4 + t114 * t8;
t1 = m(2) * t133 + qJDD(1) * mrSges(2,1) - t116 * mrSges(2,2) - t104 * t6 + (t110 * t8 + t114 * t4) * t106;
t3 = [-m(1) * g(1) - t1 * t111 + t115 * t2, t2, t8, t9, t13, -t38 * mrSges(6,2) - t72 * t58 - t123, t16; -m(1) * g(2) + t1 * t115 + t111 * t2, t1, t4, t11, t12, t38 * mrSges(6,1) - t64 * mrSges(6,3) + t72 * t51 - t77 * t59 - t119, t15; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t6, t10, t161, t39 * mrSges(6,1) + t64 * mrSges(6,2) + t73 * t51 + t77 * t58 - t122, t121;];
f_new  = t3;
