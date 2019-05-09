% Calculate vector of cutting forces with Newton-Euler
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-05-06 08:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:51:43
% EndTime: 2019-05-06 08:51:46
% DurationCPUTime: 1.43s
% Computational Cost: add. (12003->209), mult. (26491->245), div. (0->0), fcn. (16194->8), ass. (0->100)
t156 = 2 * qJD(4);
t108 = qJD(2) ^ 2;
t103 = sin(qJ(2));
t135 = qJD(1) * t103;
t106 = cos(qJ(2));
t109 = qJD(1) ^ 2;
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t124 = -t107 * g(1) - t104 * g(2);
t78 = -t109 * pkin(1) + qJDD(1) * pkin(7) + t124;
t140 = -t106 * g(3) - t103 * t78;
t86 = (-pkin(2) * t106 - qJ(3) * t103) * qJD(1);
t116 = qJDD(2) * pkin(2) + t108 * qJ(3) - t86 * t135 - qJDD(3) + t140;
t134 = t106 * qJD(1);
t101 = sin(pkin(9));
t138 = cos(pkin(9));
t82 = -t138 * qJD(2) + t101 * t135;
t132 = t82 * t134;
t155 = -qJ(4) * t132 - t116;
t136 = t106 ^ 2 * t109;
t133 = qJD(1) * qJD(2);
t127 = t106 * t133;
t130 = t104 * g(1) - t107 * g(2);
t77 = -qJDD(1) * pkin(1) - t109 * pkin(7) - t130;
t88 = t103 * qJDD(1) + t127;
t95 = t103 * t133;
t89 = t106 * qJDD(1) - t95;
t36 = (-t88 - t127) * qJ(3) + (-t89 + t95) * pkin(2) + t77;
t129 = -t103 * g(3) + t106 * t78;
t40 = -t108 * pkin(2) + qJDD(2) * qJ(3) + t86 * t134 + t129;
t144 = t101 * t36 + t138 * t40;
t83 = t101 * qJD(2) + t138 * t135;
t52 = t82 * pkin(3) - t83 * qJ(4);
t154 = pkin(3) * t136 + t89 * qJ(4) + t134 * t156 + t82 * t52 - t144;
t102 = sin(qJ(6));
t105 = cos(qJ(6));
t148 = t82 * t83;
t151 = 2 * qJD(5);
t123 = -t101 * t40 + t138 * t36;
t25 = t89 * pkin(3) - qJ(4) * t136 + qJDD(4) - t123 + ((2 * qJD(3)) + t52) * t83;
t68 = t101 * qJDD(2) + t138 * t88;
t110 = t134 * t151 + (t89 + t148) * qJ(5) + (t68 - t132) * pkin(4) + t25;
t69 = -pkin(5) * t134 - t82 * pkin(8);
t80 = t83 ^ 2;
t14 = -t80 * pkin(5) + t68 * pkin(8) + t69 * t134 + t110;
t60 = t83 * pkin(4) + qJ(5) * t134;
t67 = -t138 * qJDD(2) + t101 * t88;
t137 = qJD(3) * t82;
t75 = -0.2e1 * t137;
t79 = t82 ^ 2;
t112 = -t67 * pkin(4) - t79 * qJ(5) - t60 * t134 + qJDD(5) - t154 + t75;
t131 = t83 * t134;
t15 = t112 + (-t67 - t131) * pkin(8) + (-t89 + t148) * pkin(5);
t49 = -t102 * t82 + t105 * t83;
t31 = t49 * qJD(6) + t102 * t68 + t105 * t67;
t50 = t102 * t83 + t105 * t82;
t34 = -mrSges(7,1) * t49 + mrSges(7,2) * t50;
t92 = qJD(6) - t134;
t41 = -t92 * mrSges(7,2) + t49 * mrSges(7,3);
t85 = qJDD(6) - t89;
t12 = m(7) * (-t102 * t14 + t105 * t15) - t31 * mrSges(7,3) + t85 * mrSges(7,1) - t50 * t34 + t92 * t41;
t30 = -t50 * qJD(6) - t102 * t67 + t105 * t68;
t42 = t92 * mrSges(7,1) - t50 * mrSges(7,3);
t13 = m(7) * (t102 * t15 + t105 * t14) + t30 * mrSges(7,3) - t85 * mrSges(7,2) + t49 * t34 - t92 * t42;
t64 = -mrSges(6,1) * t134 + t82 * mrSges(6,2);
t121 = -m(6) * t110 - t89 * mrSges(6,3) + t102 * t12 - t105 * t13 - t64 * t134;
t118 = m(5) * t25 - t121;
t62 = t82 * mrSges(5,1) + mrSges(5,3) * t134;
t141 = -mrSges(4,2) * t134 + t82 * mrSges(4,3) + t62;
t54 = -t82 * mrSges(5,2) - t83 * mrSges(5,3);
t143 = -t82 * mrSges(4,1) - t83 * mrSges(4,2) - t54;
t145 = -mrSges(4,3) - mrSges(5,1);
t147 = mrSges(4,1) - mrSges(5,2);
t51 = t83 * mrSges(6,1) - t82 * mrSges(6,3);
t5 = m(4) * t123 - t147 * t89 + t141 * t134 + (-0.2e1 * m(4) * qJD(3) + t143 + t51) * t83 + (mrSges(6,2) + t145) * t68 - t118;
t126 = m(6) * t112 + t67 * mrSges(6,2) + t102 * t13 + t105 * t12 + t82 * t51;
t120 = m(5) * (0.2e1 * t137 + t154) - t126;
t61 = -t83 * mrSges(6,2) + mrSges(6,3) * t134;
t63 = t83 * mrSges(5,1) - mrSges(5,2) * t134;
t142 = t61 + t63;
t146 = mrSges(6,1) + mrSges(5,3);
t66 = -mrSges(4,1) * t134 - t83 * mrSges(4,3);
t6 = m(4) * (t75 + t144) + t143 * t82 + t145 * t67 + (mrSges(4,2) - t146) * t89 + (t66 - t142) * t134 - t120;
t90 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t135;
t91 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t134;
t153 = m(3) * t77 - t89 * mrSges(3,1) + t88 * mrSges(3,2) + (t103 * t90 - t106 * t91) * qJD(1) + t101 * t6 + t138 * t5;
t113 = -0.2e1 * qJD(4) * t83 + (t67 - t131) * pkin(3) + t155;
t149 = t79 * pkin(4);
t125 = m(7) * (t113 + (t151 + t69) * t82 + (-pkin(5) - qJ(4)) * t68 - t83 * t60 - t80 * pkin(8) - t149 + t67 * qJ(5)) + t31 * mrSges(7,2) - t30 * mrSges(7,1) + t50 * t42 - t49 * t41;
t139 = t68 * qJ(4);
t117 = m(6) * (t149 + t139 - 0.2e1 * qJD(5) * t82 + (-pkin(3) - qJ(5)) * t67 + (pkin(3) * t134 + t156 + t60) * t83 - t155) - t67 * mrSges(6,3) + t68 * mrSges(6,1) - t82 * t64 + t83 * t61 - t125;
t115 = t83 * t63 + t68 * mrSges(5,3) - m(5) * (t113 - t139) + t117;
t152 = -m(4) * t116 + t68 * mrSges(4,2) - t141 * t82 + t147 * t67 + t83 * t66 - t115;
t87 = (-mrSges(3,1) * t106 + mrSges(3,2) * t103) * qJD(1);
t4 = m(3) * t129 - qJDD(2) * mrSges(3,2) + t89 * mrSges(3,3) - qJD(2) * t90 - t101 * t5 + t87 * t134 + t138 * t6;
t8 = m(3) * t140 + qJDD(2) * mrSges(3,1) - t88 * mrSges(3,3) + qJD(2) * t91 - t87 * t135 - t152;
t150 = t103 * t4 + t106 * t8;
t2 = m(2) * t130 + qJDD(1) * mrSges(2,1) - t109 * mrSges(2,2) - t153;
t1 = m(2) * t124 - t109 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t103 * t8 + t106 * t4;
t3 = [-m(1) * g(1) + t107 * t1 - t104 * t2, t1, t4, t6, -t67 * mrSges(5,2) - t82 * t62 - t115, -t68 * mrSges(6,2) - t83 * t51 - t121, t13; -m(1) * g(2) + t104 * t1 + t107 * t2, t2, t8, t5, t67 * mrSges(5,1) + t142 * t134 + t146 * t89 + t82 * t54 + t120, t117, t12; (-m(1) - m(2)) * g(3) + t150, -m(2) * g(3) + t150, t153, t152, -t62 * t134 - t89 * mrSges(5,2) + (-t51 + t54) * t83 + (mrSges(5,1) - mrSges(6,2)) * t68 + t118, -t89 * mrSges(6,1) - t61 * t134 + t126, t125;];
f_new  = t3;
