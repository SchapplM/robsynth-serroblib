% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 10:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:52:27
% EndTime: 2019-05-06 10:52:33
% DurationCPUTime: 2.86s
% Computational Cost: add. (32606->202), mult. (74406->262), div. (0->0), fcn. (49035->10), ass. (0->99)
t114 = sin(qJ(2));
t118 = cos(qJ(2));
t140 = qJD(1) * qJD(2);
t136 = t118 * t140;
t121 = qJD(1) ^ 2;
t115 = sin(qJ(1));
t119 = cos(qJ(1));
t144 = g(1) * t115 - g(2) * t119;
t78 = -qJDD(1) * pkin(1) - t121 * pkin(7) - t144;
t87 = qJDD(1) * t114 + t136;
t137 = t114 * t140;
t88 = qJDD(1) * t118 - t137;
t131 = -t88 * pkin(2) + t78 + (-t136 - t87) * qJ(3);
t142 = qJD(1) * t114;
t143 = t118 ^ 2 * t121;
t152 = 2 * qJD(3);
t90 = -qJD(2) * pkin(3) - qJ(4) * t142;
t123 = -pkin(2) * t137 + t88 * pkin(3) - qJ(4) * t143 + qJDD(4) - t131 + (t152 + t90) * t142;
t112 = sin(qJ(6));
t116 = cos(qJ(6));
t109 = sin(pkin(10));
t110 = cos(pkin(10));
t63 = -t109 * t87 - t110 * t88;
t77 = (-t109 * t118 + t110 * t114) * qJD(1);
t67 = -qJD(2) * pkin(4) - pkin(8) * t77;
t76 = (-t109 * t114 - t110 * t118) * qJD(1);
t75 = t76 ^ 2;
t122 = -t63 * pkin(4) - t75 * pkin(8) + t67 * t77 + t123;
t104 = -qJD(2) + qJD(5);
t102 = t104 ^ 2;
t103 = -qJDD(2) + qJDD(5);
t113 = sin(qJ(5));
t117 = cos(qJ(5));
t120 = qJD(2) ^ 2;
t134 = -g(1) * t119 - g(2) * t115;
t79 = -pkin(1) * t121 + qJDD(1) * pkin(7) + t134;
t138 = -g(3) * t114 + t118 * t79;
t141 = qJD(1) * t118;
t84 = (-pkin(2) * t118 - qJ(3) * t114) * qJD(1);
t128 = -pkin(2) * t120 + qJDD(2) * qJ(3) + qJD(2) * t152 + t141 * t84 + t138;
t44 = -pkin(3) * t143 - qJ(4) * t88 + qJD(2) * t90 + t128;
t145 = -g(3) * t118 - t114 * t79;
t51 = -qJDD(2) * pkin(2) - t120 * qJ(3) + t142 * t84 + qJDD(3) - t145;
t45 = (-t87 + t136) * qJ(4) + (-t114 * t118 * t121 - qJDD(2)) * pkin(3) + t51;
t135 = -0.2e1 * qJD(4) * t77 - t109 * t44 + t110 * t45;
t64 = -t109 * t88 + t110 * t87;
t20 = (-qJD(2) * t76 - t64) * pkin(8) + (t76 * t77 - qJDD(2)) * pkin(4) + t135;
t139 = 0.2e1 * qJD(4) * t76 + t109 * t45 + t110 * t44;
t22 = -pkin(4) * t75 + pkin(8) * t63 + qJD(2) * t67 + t139;
t148 = t113 * t20 + t117 * t22;
t55 = -t113 * t77 + t117 * t76;
t56 = t113 * t76 + t117 * t77;
t43 = -pkin(5) * t55 - pkin(9) * t56;
t17 = -pkin(5) * t102 + pkin(9) * t103 + t43 * t55 + t148;
t30 = -qJD(5) * t56 - t113 * t64 + t117 * t63;
t31 = qJD(5) * t55 + t113 * t63 + t117 * t64;
t18 = (-t104 * t55 - t31) * pkin(9) + (t104 * t56 - t30) * pkin(5) + t122;
t49 = t104 * t116 - t112 * t56;
t26 = qJD(6) * t49 + t103 * t112 + t116 * t31;
t29 = qJDD(6) - t30;
t50 = t104 * t112 + t116 * t56;
t32 = -mrSges(7,1) * t49 + mrSges(7,2) * t50;
t54 = qJD(6) - t55;
t35 = -mrSges(7,2) * t54 + mrSges(7,3) * t49;
t14 = m(7) * (-t112 * t17 + t116 * t18) - t26 * mrSges(7,3) + t29 * mrSges(7,1) - t50 * t32 + t54 * t35;
t25 = -qJD(6) * t50 + t103 * t116 - t112 * t31;
t36 = mrSges(7,1) * t54 - mrSges(7,3) * t50;
t15 = m(7) * (t112 * t18 + t116 * t17) + t25 * mrSges(7,3) - t29 * mrSges(7,2) + t49 * t32 - t54 * t36;
t52 = -mrSges(6,2) * t104 + mrSges(6,3) * t55;
t53 = mrSges(6,1) * t104 - mrSges(6,3) * t56;
t132 = m(6) * t122 - mrSges(6,1) * t30 + mrSges(6,2) * t31 + t112 * t15 + t116 * t14 - t52 * t55 + t53 * t56;
t65 = qJD(2) * mrSges(5,2) + mrSges(5,3) * t76;
t66 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t77;
t127 = m(5) * t123 - mrSges(5,1) * t63 + mrSges(5,2) * t64 - t65 * t76 + t66 * t77 + t132;
t125 = -t127 + m(4) * ((pkin(2) * qJD(2) - (2 * qJD(3))) * t142 + t131) - t88 * mrSges(4,1);
t94 = mrSges(4,2) * t141 + qJD(2) * mrSges(4,3);
t146 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t141 + t94;
t91 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t142;
t92 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t142;
t155 = (-t146 * t118 + (t91 - t92) * t114) * qJD(1) + m(3) * t78 - t88 * mrSges(3,1) + (mrSges(3,2) - mrSges(4,3)) * t87 + t125;
t42 = -mrSges(6,1) * t55 + mrSges(6,2) * t56;
t10 = m(6) * t148 - mrSges(6,2) * t103 + mrSges(6,3) * t30 - t104 * t53 - t112 * t14 + t116 * t15 + t42 * t55;
t133 = -t113 * t22 + t117 * t20;
t126 = m(7) * (-pkin(5) * t103 - pkin(9) * t102 + t43 * t56 - t133) - t25 * mrSges(7,1) + t26 * mrSges(7,2) - t49 * t35 + t50 * t36;
t11 = m(6) * t133 + mrSges(6,1) * t103 - mrSges(6,3) * t31 + t104 * t52 - t42 * t56 - t126;
t60 = -mrSges(5,1) * t76 + mrSges(5,2) * t77;
t7 = m(5) * t135 - qJDD(2) * mrSges(5,1) - mrSges(5,3) * t64 - qJD(2) * t65 + t10 * t113 + t11 * t117 - t60 * t77;
t8 = m(5) * t139 + qJDD(2) * mrSges(5,2) + mrSges(5,3) * t63 + qJD(2) * t66 + t10 * t117 - t11 * t113 + t60 * t76;
t85 = (-mrSges(4,1) * t118 - mrSges(4,3) * t114) * qJD(1);
t130 = m(4) * t128 + qJDD(2) * mrSges(4,3) + qJD(2) * t92 - t109 * t7 + t110 * t8 + t141 * t85;
t149 = mrSges(3,3) + mrSges(4,2);
t86 = (-mrSges(3,1) * t118 + mrSges(3,2) * t114) * qJD(1);
t4 = m(3) * t138 - qJDD(2) * mrSges(3,2) - qJD(2) * t91 + t141 * t86 + t149 * t88 + t130;
t129 = -m(4) * t51 - t109 * t8 - t110 * t7;
t5 = m(3) * t145 - t149 * t87 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t146 * qJD(2) + (-t85 - t86) * t142 + t129;
t151 = t114 * t4 + t118 * t5;
t9 = m(2) * t144 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t155;
t1 = m(2) * t134 - mrSges(2,1) * t121 - qJDD(1) * mrSges(2,2) - t114 * t5 + t118 * t4;
t2 = [-m(1) * g(1) + t1 * t119 - t115 * t9, t1, t4, t88 * mrSges(4,2) + t130, t8, t10, t15; -m(1) * g(2) + t1 * t115 + t119 * t9, t9, t5, t125 + (-t114 * t92 - t118 * t94) * qJD(1) - t87 * mrSges(4,3), t7, t11, t14; (-m(1) - m(2)) * g(3) + t151, -m(2) * g(3) + t151, t155, -qJDD(2) * mrSges(4,1) + t87 * mrSges(4,2) - qJD(2) * t94 + t142 * t85 - t129, t127, t132, t126;];
f_new  = t2;
