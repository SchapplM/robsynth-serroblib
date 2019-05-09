% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR7
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
% Datum: 2019-05-06 22:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:16:20
% EndTime: 2019-05-06 22:16:28
% DurationCPUTime: 2.88s
% Computational Cost: add. (36586->204), mult. (74271->259), div. (0->0), fcn. (48510->10), ass. (0->101)
t116 = sin(qJ(2));
t121 = cos(qJ(2));
t141 = qJD(1) * qJD(2);
t138 = t121 * t141;
t124 = qJD(1) ^ 2;
t117 = sin(qJ(1));
t122 = cos(qJ(1));
t144 = t117 * g(1) - t122 * g(2);
t78 = -qJDD(1) * pkin(1) - t124 * pkin(7) - t144;
t88 = t116 * qJDD(1) + t138;
t139 = t116 * t141;
t89 = t121 * qJDD(1) - t139;
t133 = -t89 * pkin(2) + t78 + (-t138 - t88) * qJ(3);
t113 = sin(qJ(6));
t118 = cos(qJ(6));
t114 = sin(qJ(5));
t119 = cos(qJ(5));
t107 = -qJD(2) + qJD(4);
t143 = qJD(1) * t116;
t145 = t121 ^ 2 * t124;
t154 = 2 * qJD(3);
t95 = -qJD(2) * pkin(3) - pkin(8) * t143;
t126 = -pkin(2) * t139 + t89 * pkin(3) - pkin(8) * t145 - t133 + (t154 + t95) * t143;
t115 = sin(qJ(4));
t120 = cos(qJ(4));
t77 = (-t115 * t121 + t116 * t120) * qJD(1);
t54 = -t77 * qJD(4) - t115 * t88 - t120 * t89;
t142 = qJD(1) * t121;
t76 = -t115 * t143 - t120 * t142;
t55 = t76 * qJD(4) - t115 * t89 + t120 * t88;
t25 = t126 + (t107 * t77 - t54) * pkin(4) + (-t107 * t76 - t55) * pkin(9);
t105 = t107 ^ 2;
t106 = -qJDD(2) + qJDD(4);
t123 = qJD(2) ^ 2;
t135 = -t122 * g(1) - t117 * g(2);
t79 = -t124 * pkin(1) + qJDD(1) * pkin(7) + t135;
t140 = -t116 * g(3) + t121 * t79;
t85 = (-pkin(2) * t121 - qJ(3) * t116) * qJD(1);
t128 = -t123 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t154 + t85 * t142 + t140;
t40 = -pkin(3) * t145 - t89 * pkin(8) + qJD(2) * t95 + t128;
t146 = -t121 * g(3) - t116 * t79;
t53 = -qJDD(2) * pkin(2) - t123 * qJ(3) + t85 * t143 + qJDD(3) - t146;
t41 = (-t88 + t138) * pkin(8) + (-t116 * t121 * t124 - qJDD(2)) * pkin(3) + t53;
t149 = t115 * t41 + t120 * t40;
t62 = -t76 * pkin(4) - t77 * pkin(9);
t28 = -t105 * pkin(4) + t106 * pkin(9) + t76 * t62 + t149;
t137 = -t114 * t28 + t119 * t25;
t64 = t119 * t107 - t114 * t77;
t34 = t64 * qJD(5) + t114 * t106 + t119 * t55;
t52 = qJDD(5) - t54;
t65 = t114 * t107 + t119 * t77;
t75 = qJD(5) - t76;
t16 = (t64 * t75 - t34) * pkin(10) + (t64 * t65 + t52) * pkin(5) + t137;
t150 = t114 * t25 + t119 * t28;
t33 = -t65 * qJD(5) + t119 * t106 - t114 * t55;
t58 = t75 * pkin(5) - t65 * pkin(10);
t63 = t64 ^ 2;
t17 = -t63 * pkin(5) + t33 * pkin(10) - t75 * t58 + t150;
t44 = -t113 * t65 + t118 * t64;
t22 = t44 * qJD(6) + t113 * t33 + t118 * t34;
t45 = t113 * t64 + t118 * t65;
t30 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t70 = qJD(6) + t75;
t38 = -t70 * mrSges(7,2) + t44 * mrSges(7,3);
t48 = qJDD(6) + t52;
t14 = m(7) * (-t113 * t17 + t118 * t16) - t22 * mrSges(7,3) + t48 * mrSges(7,1) - t45 * t30 + t70 * t38;
t21 = -t45 * qJD(6) - t113 * t34 + t118 * t33;
t39 = t70 * mrSges(7,1) - t45 * mrSges(7,3);
t15 = m(7) * (t113 * t16 + t118 * t17) + t21 * mrSges(7,3) - t48 * mrSges(7,2) + t44 * t30 - t70 * t39;
t46 = -t64 * mrSges(6,1) + t65 * mrSges(6,2);
t56 = -t75 * mrSges(6,2) + t64 * mrSges(6,3);
t11 = m(6) * t137 + t52 * mrSges(6,1) - t34 * mrSges(6,3) + t113 * t15 + t118 * t14 - t65 * t46 + t75 * t56;
t57 = t75 * mrSges(6,1) - t65 * mrSges(6,3);
t12 = m(6) * t150 - t52 * mrSges(6,2) + t33 * mrSges(6,3) - t113 * t14 + t118 * t15 + t64 * t46 - t75 * t57;
t66 = -t107 * mrSges(5,2) + t76 * mrSges(5,3);
t67 = t107 * mrSges(5,1) - t77 * mrSges(5,3);
t134 = m(5) * t126 - t54 * mrSges(5,1) + t55 * mrSges(5,2) + t119 * t11 + t114 * t12 - t76 * t66 + t77 * t67;
t129 = m(4) * ((pkin(2) * qJD(2) - (2 * qJD(3))) * t143 + t133) - t89 * mrSges(4,1) - t134;
t94 = mrSges(4,2) * t142 + qJD(2) * mrSges(4,3);
t147 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t142 + t94;
t91 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t143;
t92 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t143;
t157 = ((t91 - t92) * t116 - t147 * t121) * qJD(1) + m(3) * t78 - t89 * mrSges(3,1) + (mrSges(3,2) - mrSges(4,3)) * t88 + t129;
t136 = -t115 * t40 + t120 * t41;
t27 = -t106 * pkin(4) - t105 * pkin(9) + t77 * t62 - t136;
t130 = t21 * mrSges(7,1) + t44 * t38 - m(7) * (-t33 * pkin(5) - t63 * pkin(10) + t65 * t58 + t27) - t22 * mrSges(7,2) - t45 * t39;
t125 = m(6) * t27 - t33 * mrSges(6,1) + t34 * mrSges(6,2) - t64 * t56 + t65 * t57 - t130;
t61 = -t76 * mrSges(5,1) + t77 * mrSges(5,2);
t13 = m(5) * t136 + t106 * mrSges(5,1) - t55 * mrSges(5,3) + t107 * t66 - t77 * t61 - t125;
t8 = m(5) * t149 - t106 * mrSges(5,2) + t54 * mrSges(5,3) - t107 * t67 - t114 * t11 + t119 * t12 + t76 * t61;
t86 = (-mrSges(4,1) * t121 - mrSges(4,3) * t116) * qJD(1);
t132 = m(4) * t128 + qJDD(2) * mrSges(4,3) + qJD(2) * t92 - t115 * t13 + t120 * t8 + t86 * t142;
t151 = mrSges(3,3) + mrSges(4,2);
t87 = (-mrSges(3,1) * t121 + mrSges(3,2) * t116) * qJD(1);
t4 = m(3) * t140 - qJDD(2) * mrSges(3,2) - qJD(2) * t91 + t87 * t142 + t151 * t89 + t132;
t131 = -m(4) * t53 - t115 * t8 - t120 * t13;
t5 = m(3) * t146 - t151 * t88 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t147 * qJD(2) + (-t86 - t87) * t143 + t131;
t153 = t116 * t4 + t121 * t5;
t6 = m(2) * t144 + qJDD(1) * mrSges(2,1) - t124 * mrSges(2,2) - t157;
t1 = m(2) * t135 - t124 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t116 * t5 + t121 * t4;
t2 = [-m(1) * g(1) + t122 * t1 - t117 * t6, t1, t4, t89 * mrSges(4,2) + t132, t8, t12, t15; -m(1) * g(2) + t117 * t1 + t122 * t6, t6, t5, -t88 * mrSges(4,3) + (-t116 * t92 - t121 * t94) * qJD(1) + t129, t13, t11, t14; (-m(1) - m(2)) * g(3) + t153, -m(2) * g(3) + t153, t157, -qJDD(2) * mrSges(4,1) + t88 * mrSges(4,2) - qJD(2) * t94 + t86 * t143 - t131, t134, t125, -t130;];
f_new  = t2;
