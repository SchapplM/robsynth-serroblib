% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-06 22:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 21:58:16
% EndTime: 2019-05-06 21:58:23
% DurationCPUTime: 2.99s
% Computational Cost: add. (36692->203), mult. (79037->260), div. (0->0), fcn. (53039->10), ass. (0->101)
t113 = sin(qJ(2));
t118 = cos(qJ(2));
t139 = qJD(1) * qJD(2);
t136 = t118 * t139;
t121 = qJD(1) ^ 2;
t114 = sin(qJ(1));
t119 = cos(qJ(1));
t142 = t114 * g(1) - t119 * g(2);
t76 = -qJDD(1) * pkin(1) - t121 * pkin(7) - t142;
t85 = t113 * qJDD(1) + t136;
t137 = t113 * t139;
t86 = t118 * qJDD(1) - t137;
t131 = -t86 * pkin(2) + t76 + (-t136 - t85) * qJ(3);
t141 = qJD(1) * t113;
t143 = t118 ^ 2 * t121;
t152 = 2 * qJD(3);
t92 = -qJD(2) * pkin(3) - pkin(8) * t141;
t123 = -pkin(2) * t137 + t86 * pkin(3) - pkin(8) * t143 - t131 + (t152 + t92) * t141;
t110 = sin(qJ(6));
t115 = cos(qJ(6));
t112 = sin(qJ(4));
t117 = cos(qJ(4));
t75 = (-t112 * t118 + t113 * t117) * qJD(1);
t54 = -t75 * qJD(4) - t112 * t85 - t117 * t86;
t104 = -qJD(2) + qJD(4);
t67 = t104 * pkin(4) - t75 * pkin(9);
t74 = (-t112 * t113 - t117 * t118) * qJD(1);
t73 = t74 ^ 2;
t122 = -t54 * pkin(4) - t73 * pkin(9) + t75 * t67 + t123;
t111 = sin(qJ(5));
t116 = cos(qJ(5));
t103 = -qJDD(2) + qJDD(4);
t120 = qJD(2) ^ 2;
t134 = -t119 * g(1) - t114 * g(2);
t77 = -t121 * pkin(1) + qJDD(1) * pkin(7) + t134;
t138 = -t113 * g(3) + t118 * t77;
t140 = qJD(1) * t118;
t82 = (-pkin(2) * t118 - qJ(3) * t113) * qJD(1);
t128 = -t120 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t152 + t82 * t140 + t138;
t42 = -pkin(3) * t143 - t86 * pkin(8) + qJD(2) * t92 + t128;
t144 = -t118 * g(3) - t113 * t77;
t53 = -qJDD(2) * pkin(2) - t120 * qJ(3) + t82 * t141 + qJDD(3) - t144;
t44 = (-t85 + t136) * pkin(8) + (-t113 * t118 * t121 - qJDD(2)) * pkin(3) + t53;
t135 = -t112 * t42 + t117 * t44;
t55 = t74 * qJD(4) - t112 * t86 + t117 * t85;
t20 = (t104 * t74 - t55) * pkin(9) + (t74 * t75 + t103) * pkin(4) + t135;
t147 = t112 * t44 + t117 * t42;
t22 = -t73 * pkin(4) + t54 * pkin(9) - t104 * t67 + t147;
t148 = t111 * t20 + t116 * t22;
t62 = -t111 * t75 + t116 * t74;
t63 = t111 * t74 + t116 * t75;
t45 = -t62 * pkin(5) - t63 * pkin(10);
t98 = qJD(5) + t104;
t96 = t98 ^ 2;
t97 = qJDD(5) + t103;
t17 = -t96 * pkin(5) + t97 * pkin(10) + t62 * t45 + t148;
t30 = -t63 * qJD(5) - t111 * t55 + t116 * t54;
t31 = t62 * qJD(5) + t111 * t54 + t116 * t55;
t18 = t122 + (t63 * t98 - t30) * pkin(5) + (-t62 * t98 - t31) * pkin(10);
t49 = -t110 * t63 + t115 * t98;
t24 = t49 * qJD(6) + t110 * t97 + t115 * t31;
t29 = qJDD(6) - t30;
t50 = t110 * t98 + t115 * t63;
t32 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t58 = qJD(6) - t62;
t35 = -t58 * mrSges(7,2) + t49 * mrSges(7,3);
t14 = m(7) * (-t110 * t17 + t115 * t18) - t24 * mrSges(7,3) + t29 * mrSges(7,1) - t50 * t32 + t58 * t35;
t23 = -t50 * qJD(6) - t110 * t31 + t115 * t97;
t36 = t58 * mrSges(7,1) - t50 * mrSges(7,3);
t15 = m(7) * (t110 * t18 + t115 * t17) + t23 * mrSges(7,3) - t29 * mrSges(7,2) + t49 * t32 - t58 * t36;
t56 = -t98 * mrSges(6,2) + t62 * mrSges(6,3);
t57 = t98 * mrSges(6,1) - t63 * mrSges(6,3);
t132 = m(6) * t122 - t30 * mrSges(6,1) + t31 * mrSges(6,2) + t110 * t15 + t115 * t14 - t62 * t56 + t63 * t57;
t65 = -t104 * mrSges(5,2) + t74 * mrSges(5,3);
t66 = t104 * mrSges(5,1) - t75 * mrSges(5,3);
t127 = m(5) * t123 - t54 * mrSges(5,1) + t55 * mrSges(5,2) - t74 * t65 + t75 * t66 + t132;
t125 = m(4) * ((pkin(2) * qJD(2) - (2 * qJD(3))) * t141 + t131) - t127 - t86 * mrSges(4,1);
t91 = mrSges(4,2) * t140 + qJD(2) * mrSges(4,3);
t145 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t140 + t91;
t88 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t141;
t89 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t141;
t155 = ((t88 - t89) * t113 - t118 * t145) * qJD(1) + m(3) * t76 - t86 * mrSges(3,1) + (mrSges(3,2) - mrSges(4,3)) * t85 + t125;
t43 = -t62 * mrSges(6,1) + t63 * mrSges(6,2);
t10 = m(6) * t148 - t97 * mrSges(6,2) + t30 * mrSges(6,3) - t110 * t14 + t115 * t15 + t62 * t43 - t98 * t57;
t133 = -t111 * t22 + t116 * t20;
t126 = m(7) * (-t97 * pkin(5) - t96 * pkin(10) + t63 * t45 - t133) - t23 * mrSges(7,1) + t24 * mrSges(7,2) - t49 * t35 + t50 * t36;
t11 = m(6) * t133 + t97 * mrSges(6,1) - t31 * mrSges(6,3) - t63 * t43 + t98 * t56 - t126;
t64 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t7 = m(5) * t135 + t103 * mrSges(5,1) - t55 * mrSges(5,3) + t111 * t10 + t104 * t65 + t116 * t11 - t75 * t64;
t8 = m(5) * t147 - t103 * mrSges(5,2) + t54 * mrSges(5,3) + t116 * t10 - t104 * t66 - t111 * t11 + t74 * t64;
t83 = (-mrSges(4,1) * t118 - mrSges(4,3) * t113) * qJD(1);
t130 = m(4) * t128 + qJDD(2) * mrSges(4,3) + qJD(2) * t89 - t112 * t7 + t117 * t8 + t83 * t140;
t149 = mrSges(3,3) + mrSges(4,2);
t84 = (-mrSges(3,1) * t118 + mrSges(3,2) * t113) * qJD(1);
t4 = m(3) * t138 - qJDD(2) * mrSges(3,2) - qJD(2) * t88 + t140 * t84 + t149 * t86 + t130;
t129 = -m(4) * t53 - t112 * t8 - t117 * t7;
t5 = m(3) * t144 - t149 * t85 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t145 * qJD(2) + (-t83 - t84) * t141 + t129;
t151 = t113 * t4 + t118 * t5;
t9 = m(2) * t142 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t155;
t1 = m(2) * t134 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t113 * t5 + t118 * t4;
t2 = [-m(1) * g(1) + t119 * t1 - t114 * t9, t1, t4, t86 * mrSges(4,2) + t130, t8, t10, t15; -m(1) * g(2) + t114 * t1 + t119 * t9, t9, t5, t125 + (-t113 * t89 - t118 * t91) * qJD(1) - t85 * mrSges(4,3), t7, t11, t14; (-m(1) - m(2)) * g(3) + t151, -m(2) * g(3) + t151, t155, -qJDD(2) * mrSges(4,1) + t85 * mrSges(4,2) - qJD(2) * t91 + t141 * t83 - t129, t127, t132, t126;];
f_new  = t2;
