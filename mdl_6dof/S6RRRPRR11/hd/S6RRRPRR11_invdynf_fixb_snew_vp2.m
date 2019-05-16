% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 14:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 14:21:55
% EndTime: 2019-05-07 14:22:08
% DurationCPUTime: 4.60s
% Computational Cost: add. (62835->213), mult. (134767->278), div. (0->0), fcn. (105131->12), ass. (0->111)
t113 = sin(pkin(6));
t122 = cos(qJ(2));
t143 = qJD(1) * t122;
t138 = t113 * t143;
t104 = qJD(3) - t138;
t114 = cos(pkin(6));
t110 = t114 * qJD(1) + qJD(2);
t117 = sin(qJ(3));
t118 = sin(qJ(2));
t144 = qJD(1) * t113;
t139 = t118 * t144;
t158 = cos(qJ(3));
t82 = -t158 * t110 + t117 * t139;
t148 = t104 * t82;
t108 = t110 ^ 2;
t109 = t114 * qJDD(1) + qJDD(2);
t145 = t113 * t122;
t124 = qJD(1) ^ 2;
t119 = sin(qJ(1));
t123 = cos(qJ(1));
t137 = t119 * g(1) - t123 * g(2);
t92 = t124 * t113 * pkin(8) + qJDD(1) * pkin(1) + t137;
t147 = t114 * t92;
t135 = -t123 * g(1) - t119 * g(2);
t141 = qJDD(1) * t113;
t93 = -t124 * pkin(1) + pkin(8) * t141 + t135;
t140 = -g(3) * t145 - t118 * t93 + t122 * t147;
t95 = (-pkin(2) * t122 - pkin(9) * t118) * t144;
t48 = -t109 * pkin(2) - t108 * pkin(9) + t95 * t139 - t140;
t83 = t117 * t110 + t158 * t139;
t96 = (qJD(2) * t143 + qJDD(1) * t118) * t113;
t62 = t83 * qJD(3) - t158 * t109 + t117 * t96;
t63 = -t82 * qJD(3) + t117 * t109 + t158 * t96;
t129 = t62 * pkin(3) + t48 + (t148 - t63) * qJ(4);
t115 = sin(qJ(6));
t120 = cos(qJ(6));
t157 = pkin(3) * t104;
t159 = 2 * qJD(4);
t75 = -t104 * pkin(4) - t83 * pkin(10);
t81 = t82 ^ 2;
t125 = -t62 * pkin(4) - t81 * pkin(10) - t129 + (-t157 + t159 + t75) * t83;
t116 = sin(qJ(5));
t121 = cos(qJ(5));
t103 = t104 ^ 2;
t149 = t118 * t147 + t122 * t93;
t49 = -t108 * pkin(2) + t109 * pkin(9) + (-g(3) * t118 + t95 * t143) * t113 + t149;
t156 = t114 * g(3);
t97 = -qJD(2) * t139 + t122 * t141;
t50 = -t97 * pkin(2) - t96 * pkin(9) - t156 + (-t92 + (pkin(2) * t118 - pkin(9) * t122) * t110 * qJD(1)) * t113;
t136 = -t117 * t49 + t158 * t50;
t67 = t82 * pkin(3) - t83 * qJ(4);
t89 = qJDD(3) - t97;
t28 = -t89 * pkin(3) - t103 * qJ(4) + t83 * t67 + qJDD(4) - t136;
t22 = (-t63 - t148) * pkin(10) + (t82 * t83 - t89) * pkin(4) + t28;
t152 = t117 * t50 + t158 * t49;
t131 = -t103 * pkin(3) + t89 * qJ(4) + t104 * t159 - t82 * t67 + t152;
t25 = -t81 * pkin(4) + t62 * pkin(10) + t104 * t75 + t131;
t153 = t116 * t22 + t121 * t25;
t65 = -t116 * t83 + t121 * t82;
t66 = t116 * t82 + t121 * t83;
t44 = -t65 * pkin(5) - t66 * pkin(11);
t88 = qJDD(5) - t89;
t100 = qJD(5) - t104;
t99 = t100 ^ 2;
t19 = -t99 * pkin(5) + t88 * pkin(11) + t65 * t44 + t153;
t36 = -t66 * qJD(5) - t116 * t63 + t121 * t62;
t37 = t65 * qJD(5) + t116 * t62 + t121 * t63;
t20 = t125 + (t100 * t66 - t36) * pkin(5) + (-t100 * t65 - t37) * pkin(11);
t53 = t120 * t100 - t115 * t66;
t32 = t53 * qJD(6) + t115 * t88 + t120 * t37;
t35 = qJDD(6) - t36;
t54 = t115 * t100 + t120 * t66;
t38 = -t53 * mrSges(7,1) + t54 * mrSges(7,2);
t64 = qJD(6) - t65;
t39 = -t64 * mrSges(7,2) + t53 * mrSges(7,3);
t16 = m(7) * (-t115 * t19 + t120 * t20) - t32 * mrSges(7,3) + t35 * mrSges(7,1) - t54 * t38 + t64 * t39;
t31 = -t54 * qJD(6) - t115 * t37 + t120 * t88;
t40 = t64 * mrSges(7,1) - t54 * mrSges(7,3);
t17 = m(7) * (t115 * t20 + t120 * t19) + t31 * mrSges(7,3) - t35 * mrSges(7,2) + t53 * t38 - t64 * t40;
t55 = -t100 * mrSges(6,2) + t65 * mrSges(6,3);
t56 = t100 * mrSges(6,1) - t66 * mrSges(6,3);
t132 = m(6) * t125 - t36 * mrSges(6,1) + t37 * mrSges(6,2) + t115 * t17 + t120 * t16 - t65 * t55 + t66 * t56;
t74 = -t82 * mrSges(5,2) + t104 * mrSges(5,3);
t128 = m(5) * ((-(2 * qJD(4)) + t157) * t83 + t129) + t62 * mrSges(5,1) + t82 * t74 - t132;
t71 = -t104 * mrSges(4,2) - t82 * mrSges(4,3);
t72 = t104 * mrSges(4,1) - t83 * mrSges(4,3);
t73 = -t104 * mrSges(5,1) + t83 * mrSges(5,2);
t160 = m(4) * t48 + t62 * mrSges(4,1) + (t72 - t73) * t83 + (mrSges(4,2) - mrSges(5,3)) * t63 + t82 * t71 + t128;
t154 = -mrSges(4,3) - mrSges(5,2);
t68 = t82 * mrSges(5,1) - t83 * mrSges(5,3);
t151 = -t82 * mrSges(4,1) - t83 * mrSges(4,2) - t68;
t146 = t113 * t118;
t91 = -t110 * mrSges(3,2) + mrSges(3,3) * t138;
t94 = (-mrSges(3,1) * t122 + mrSges(3,2) * t118) * t144;
t10 = m(3) * t140 + t109 * mrSges(3,1) - t96 * mrSges(3,3) + t110 * t91 - t94 * t139 - t160;
t43 = -t65 * mrSges(6,1) + t66 * mrSges(6,2);
t12 = m(6) * t153 - t88 * mrSges(6,2) + t36 * mrSges(6,3) - t100 * t56 - t115 * t16 + t120 * t17 + t65 * t43;
t134 = -t116 * t25 + t121 * t22;
t127 = m(7) * (-t88 * pkin(5) - t99 * pkin(11) + t66 * t44 - t134) - t31 * mrSges(7,1) + t32 * mrSges(7,2) - t53 * t39 + t54 * t40;
t13 = m(6) * t134 + t88 * mrSges(6,1) - t37 * mrSges(6,3) + t100 * t55 - t66 * t43 - t127;
t133 = m(5) * t131 + t89 * mrSges(5,3) + t104 * t73 - t116 * t13 + t121 * t12;
t7 = m(4) * t152 - t89 * mrSges(4,2) - t104 * t72 + t151 * t82 + t154 * t62 + t133;
t130 = -m(5) * t28 - t116 * t12 - t121 * t13;
t8 = m(4) * t136 + (mrSges(4,1) + mrSges(5,1)) * t89 + t151 * t83 + t154 * t63 + (t71 + t74) * t104 + t130;
t90 = t110 * mrSges(3,1) - mrSges(3,3) * t139;
t4 = m(3) * (-g(3) * t146 + t149) + t97 * mrSges(3,3) - t109 * mrSges(3,2) + t94 * t138 - t110 * t90 + t158 * t7 - t117 * t8;
t6 = m(3) * (-t113 * t92 - t156) + t96 * mrSges(3,2) - t97 * mrSges(3,1) + t117 * t7 + t158 * t8 + (t118 * t90 - t122 * t91) * t144;
t142 = t10 * t145 + t114 * t6 + t4 * t146;
t2 = m(2) * t135 - t124 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t118 * t10 + t122 * t4;
t1 = m(2) * t137 + qJDD(1) * mrSges(2,1) - t124 * mrSges(2,2) - t113 * t6 + (t122 * t10 + t118 * t4) * t114;
t3 = [-m(1) * g(1) - t119 * t1 + t123 * t2, t2, t4, t7, -t62 * mrSges(5,2) - t82 * t68 + t133, t12, t17; -m(1) * g(2) + t123 * t1 + t119 * t2, t1, t10, t8, -t63 * mrSges(5,3) - t83 * t73 + t128, t13, t16; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t6, t160, -t89 * mrSges(5,1) + t63 * mrSges(5,2) - t104 * t74 + t83 * t68 - t130, t132, t127;];
f_new  = t3;
