% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 14:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:43:04
% EndTime: 2019-05-06 14:43:11
% DurationCPUTime: 2.96s
% Computational Cost: add. (34357->203), mult. (76129->262), div. (0->0), fcn. (50239->10), ass. (0->100)
t113 = sin(qJ(2));
t117 = cos(qJ(2));
t139 = qJD(1) * qJD(2);
t136 = t117 * t139;
t120 = qJD(1) ^ 2;
t114 = sin(qJ(1));
t118 = cos(qJ(1));
t143 = t114 * g(1) - t118 * g(2);
t77 = -qJDD(1) * pkin(1) - t120 * pkin(7) - t143;
t86 = qJDD(1) * t113 + t136;
t135 = t113 * t139;
t87 = qJDD(1) * t117 - t135;
t130 = -t87 * pkin(2) + t77 + (-t136 - t86) * qJ(3);
t141 = qJD(1) * t113;
t142 = t117 ^ 2 * t120;
t152 = 2 * qJD(3);
t93 = -qJD(2) * pkin(3) - pkin(8) * t141;
t122 = -pkin(2) * t135 + t87 * pkin(3) - pkin(8) * t142 - t130 + (t152 + t93) * t141;
t111 = sin(qJ(6));
t115 = cos(qJ(6));
t112 = sin(qJ(4));
t116 = cos(qJ(4));
t76 = (-t112 * t117 + t113 * t116) * qJD(1);
t54 = -qJD(4) * t76 - t112 * t86 - t116 * t87;
t103 = -qJD(2) + qJD(4);
t67 = pkin(4) * t103 - qJ(5) * t76;
t75 = (-t112 * t113 - t116 * t117) * qJD(1);
t74 = t75 ^ 2;
t121 = -t54 * pkin(4) - t74 * qJ(5) + t76 * t67 + qJDD(5) + t122;
t101 = t103 ^ 2;
t102 = -qJDD(2) + qJDD(4);
t108 = sin(pkin(10));
t109 = cos(pkin(10));
t151 = 2 * qJD(5);
t119 = qJD(2) ^ 2;
t133 = -g(1) * t118 - g(2) * t114;
t78 = -pkin(1) * t120 + qJDD(1) * pkin(7) + t133;
t137 = -g(3) * t113 + t117 * t78;
t140 = qJD(1) * t117;
t83 = (-pkin(2) * t117 - qJ(3) * t113) * qJD(1);
t127 = -pkin(2) * t119 + qJDD(2) * qJ(3) + qJD(2) * t152 + t83 * t140 + t137;
t44 = -pkin(3) * t142 - pkin(8) * t87 + qJD(2) * t93 + t127;
t146 = -t117 * g(3) - t113 * t78;
t53 = -qJDD(2) * pkin(2) - t119 * qJ(3) + t83 * t141 + qJDD(3) - t146;
t45 = (-t86 + t136) * pkin(8) + (-t113 * t117 * t120 - qJDD(2)) * pkin(3) + t53;
t134 = -t112 * t44 + t116 * t45;
t55 = qJD(4) * t75 - t112 * t87 + t116 * t86;
t20 = (t103 * t75 - t55) * qJ(5) + (t75 * t76 + t102) * pkin(4) + t134;
t147 = t112 * t45 + t116 * t44;
t22 = -pkin(4) * t74 + t54 * qJ(5) - t103 * t67 + t147;
t63 = -t108 * t76 + t109 * t75;
t138 = t108 * t20 + t109 * t22 + t63 * t151;
t64 = t108 * t75 + t109 * t76;
t43 = -pkin(5) * t63 - pkin(9) * t64;
t17 = -pkin(5) * t101 + pkin(9) * t102 + t63 * t43 + t138;
t30 = -t108 * t55 + t109 * t54;
t31 = t108 * t54 + t109 * t55;
t18 = t121 + (-t103 * t63 - t31) * pkin(9) + (t103 * t64 - t30) * pkin(5);
t51 = t103 * t115 - t111 * t64;
t26 = t51 * qJD(6) + t102 * t111 + t115 * t31;
t29 = qJDD(6) - t30;
t52 = t103 * t111 + t115 * t64;
t32 = -mrSges(7,1) * t51 + mrSges(7,2) * t52;
t59 = qJD(6) - t63;
t35 = -mrSges(7,2) * t59 + mrSges(7,3) * t51;
t14 = m(7) * (-t111 * t17 + t115 * t18) - t26 * mrSges(7,3) + t29 * mrSges(7,1) - t52 * t32 + t59 * t35;
t25 = -t52 * qJD(6) + t102 * t115 - t111 * t31;
t36 = mrSges(7,1) * t59 - mrSges(7,3) * t52;
t15 = m(7) * (t111 * t18 + t115 * t17) + t25 * mrSges(7,3) - t29 * mrSges(7,2) + t51 * t32 - t59 * t36;
t56 = -mrSges(6,2) * t103 + t63 * mrSges(6,3);
t57 = mrSges(6,1) * t103 - t64 * mrSges(6,3);
t131 = m(6) * t121 - t30 * mrSges(6,1) + t31 * mrSges(6,2) + t111 * t15 + t115 * t14 - t63 * t56 + t64 * t57;
t66 = -mrSges(5,2) * t103 + mrSges(5,3) * t75;
t68 = mrSges(5,1) * t103 - mrSges(5,3) * t76;
t126 = m(5) * t122 - t54 * mrSges(5,1) + t55 * mrSges(5,2) - t75 * t66 + t76 * t68 + t131;
t124 = -t126 - t87 * mrSges(4,1) + m(4) * ((pkin(2) * qJD(2) - (2 * qJD(3))) * t141 + t130);
t92 = mrSges(4,2) * t140 + qJD(2) * mrSges(4,3);
t144 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t140 + t92;
t89 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t141;
t90 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t141;
t155 = ((t89 - t90) * t113 - t144 * t117) * qJD(1) + m(3) * t77 - t87 * mrSges(3,1) + (mrSges(3,2) - mrSges(4,3)) * t86 + t124;
t42 = -mrSges(6,1) * t63 + mrSges(6,2) * t64;
t10 = m(6) * t138 - t102 * mrSges(6,2) + t30 * mrSges(6,3) - t103 * t57 - t111 * t14 + t115 * t15 + t63 * t42;
t132 = t108 * t22 - t109 * t20;
t125 = m(7) * (-t102 * pkin(5) - t101 * pkin(9) + (t151 + t43) * t64 + t132) - t25 * mrSges(7,1) + t26 * mrSges(7,2) - t51 * t35 + t52 * t36;
t11 = m(6) * (-0.2e1 * qJD(5) * t64 - t132) - t31 * mrSges(6,3) + t102 * mrSges(6,1) - t64 * t42 + t103 * t56 - t125;
t65 = -mrSges(5,1) * t75 + mrSges(5,2) * t76;
t7 = m(5) * t134 + t102 * mrSges(5,1) - t55 * mrSges(5,3) + t108 * t10 + t103 * t66 + t109 * t11 - t76 * t65;
t8 = m(5) * t147 - t102 * mrSges(5,2) + t54 * mrSges(5,3) + t109 * t10 - t103 * t68 - t108 * t11 + t75 * t65;
t84 = (-mrSges(4,1) * t117 - mrSges(4,3) * t113) * qJD(1);
t129 = m(4) * t127 + qJDD(2) * mrSges(4,3) + qJD(2) * t90 - t112 * t7 + t116 * t8 + t84 * t140;
t148 = mrSges(3,3) + mrSges(4,2);
t85 = (-mrSges(3,1) * t117 + mrSges(3,2) * t113) * qJD(1);
t4 = m(3) * t137 - qJDD(2) * mrSges(3,2) - qJD(2) * t89 + t85 * t140 + t148 * t87 + t129;
t128 = -m(4) * t53 - t112 * t8 - t116 * t7;
t5 = m(3) * t146 - t148 * t86 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t144 * qJD(2) + (-t84 - t85) * t141 + t128;
t150 = t113 * t4 + t117 * t5;
t9 = m(2) * t143 + qJDD(1) * mrSges(2,1) - t120 * mrSges(2,2) - t155;
t1 = m(2) * t133 - t120 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t113 * t5 + t117 * t4;
t2 = [-m(1) * g(1) + t1 * t118 - t114 * t9, t1, t4, t87 * mrSges(4,2) + t129, t8, t10, t15; -m(1) * g(2) + t1 * t114 + t118 * t9, t9, t5, t124 + (-t113 * t90 - t117 * t92) * qJD(1) - t86 * mrSges(4,3), t7, t11, t14; (-m(1) - m(2)) * g(3) + t150, -m(2) * g(3) + t150, t155, -qJDD(2) * mrSges(4,1) + t86 * mrSges(4,2) - qJD(2) * t92 + t84 * t141 - t128, t126, t131, t125;];
f_new  = t2;
