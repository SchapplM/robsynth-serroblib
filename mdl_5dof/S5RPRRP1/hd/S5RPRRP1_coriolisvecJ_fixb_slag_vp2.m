% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:14
% EndTime: 2019-12-05 17:59:23
% DurationCPUTime: 2.89s
% Computational Cost: add. (1913->242), mult. (4215->338), div. (0->0), fcn. (2450->4), ass. (0->119)
t155 = Ifges(5,4) + Ifges(6,4);
t182 = Ifges(5,1) + Ifges(6,1);
t181 = Ifges(5,5) + Ifges(6,5);
t180 = Ifges(5,2) + Ifges(6,2);
t179 = Ifges(6,6) + Ifges(5,6);
t118 = sin(qJ(3));
t119 = cos(qJ(4));
t120 = cos(qJ(3));
t164 = sin(qJ(4));
t91 = -t119 * t118 - t120 * t164;
t83 = t91 * qJD(1);
t186 = t155 * t83;
t136 = t164 * t118;
t143 = qJD(1) * t120;
t84 = -qJD(1) * t136 + t119 * t143;
t185 = t155 * t84;
t115 = qJD(3) + qJD(4);
t184 = t179 * t115 + t180 * t83 + t185;
t183 = t181 * t115 + t182 * t84 + t186;
t162 = mrSges(6,3) * t83;
t62 = -mrSges(6,2) * t115 + t162;
t163 = mrSges(5,3) * t83;
t63 = -mrSges(5,2) * t115 + t163;
t178 = t62 + t63;
t158 = t84 * mrSges(6,3);
t64 = mrSges(6,1) * t115 - t158;
t159 = t84 * mrSges(5,3);
t65 = mrSges(5,1) * t115 - t159;
t177 = -t65 - t64;
t144 = qJD(1) * t118;
t121 = -pkin(1) - pkin(6);
t101 = qJD(1) * t121 + qJD(2);
t147 = t101 * t118;
t73 = -pkin(7) * t144 + t147;
t66 = t164 * t73;
t94 = t120 * t101;
t74 = -pkin(7) * t143 + t94;
t69 = qJD(3) * pkin(3) + t74;
t27 = t119 * t69 - t66;
t75 = t84 * qJ(5);
t13 = -t75 + t27;
t176 = (qJ(2) * (m(3) + m(4)));
t139 = qJD(1) * qJD(3);
t134 = t120 * t139;
t175 = qJD(1) ^ 2;
t172 = -t83 / 0.2e1;
t169 = t84 / 0.2e1;
t167 = pkin(4) * t84;
t161 = pkin(3) * t119;
t52 = t115 * t91;
t45 = t52 * qJD(1);
t160 = t45 * mrSges(6,3);
t154 = pkin(7) - t121;
t32 = t119 * t74 - t66;
t96 = t154 * t118;
t97 = t154 * t120;
t54 = -t119 * t96 - t164 * t97;
t153 = Ifges(4,4) * t118;
t152 = Ifges(4,4) * t120;
t67 = t119 * t73;
t150 = t83 * qJ(5);
t149 = Ifges(4,5) * qJD(3);
t148 = Ifges(4,6) * qJD(3);
t145 = t119 * t120;
t107 = t118 * pkin(3) + qJ(2);
t95 = pkin(3) * t134 + qJD(1) * qJD(2);
t100 = pkin(3) * t144 + qJD(1) * qJ(2);
t142 = qJD(3) * t118;
t141 = qJD(3) * t120;
t140 = qJD(4) * t119;
t102 = pkin(3) * t141 + qJD(2);
t138 = pkin(3) * t143;
t135 = qJD(4) * t164;
t133 = pkin(7) * qJD(1) - t101;
t127 = t136 - t145;
t46 = t115 * t127 * qJD(1);
t132 = pkin(3) * t46 * t164;
t31 = -t164 * t74 - t67;
t53 = -t119 * t97 + t164 * t96;
t129 = mrSges(4,1) * t118 + mrSges(4,2) * t120;
t128 = -t164 * t69 - t67;
t126 = qJ(2) * (mrSges(4,1) * t120 - mrSges(4,2) * t118);
t70 = t133 * t142;
t71 = t133 * t141;
t7 = -t119 * t71 - t135 * t73 + t69 * t140 + t164 * t70;
t89 = t154 * t142;
t90 = qJD(3) * t97;
t11 = -t119 * t90 + t135 * t96 - t97 * t140 + t164 * t89;
t10 = t115 * pkin(4) + t13;
t14 = -t128 + t150;
t2 = t46 * qJ(5) + t83 * qJD(5) + t7;
t8 = qJD(4) * t128 + t119 * t70 + t164 * t71;
t3 = -t45 * qJ(5) - t84 * qJD(5) + t8;
t51 = -qJD(3) * t136 + t115 * t145 - t118 * t135;
t125 = t10 * t52 - t127 * t3 + t14 * t51 - t2 * t91;
t124 = -t127 * t8 - t128 * t51 + t27 * t52 - t7 * t91;
t12 = -qJD(4) * t54 + t119 * t89 + t164 * t90;
t55 = -pkin(4) * t83 + qJD(5) + t100;
t123 = t8 * mrSges(5,1) + t3 * mrSges(6,1) - t7 * mrSges(5,2) - t2 * mrSges(6,2) + t10 * t162 - t100 * (mrSges(5,1) * t84 + mrSges(5,2) * t83) + t27 * t163 - t55 * (mrSges(6,1) * t84 + mrSges(6,2) * t83) + t179 * t46 + t181 * t45 - (t182 * t83 - t185) * t84 / 0.2e1 + t184 * t169 - (-t179 * t84 + t181 * t83) * t115 / 0.2e1 + (-t180 * t84 + t183 + t186) * t172;
t110 = pkin(4) + t161;
t99 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t143;
t98 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t144;
t93 = qJD(1) * t129;
t82 = t149 + (t120 * Ifges(4,1) - t153) * qJD(1);
t81 = t148 + (-t118 * Ifges(4,2) + t152) * qJD(1);
t72 = -pkin(4) * t91 + t107;
t61 = t138 + t167;
t50 = -mrSges(5,1) * t83 + mrSges(5,2) * t84;
t49 = -mrSges(6,1) * t83 + mrSges(6,2) * t84;
t40 = t45 * mrSges(6,2);
t39 = pkin(4) * t51 + t102;
t30 = qJ(5) * t91 + t54;
t29 = qJ(5) * t127 + t53;
t26 = -pkin(4) * t46 + t95;
t17 = -t75 + t32;
t16 = t31 - t150;
t5 = -t52 * qJ(5) + qJD(5) * t127 + t12;
t4 = -t51 * qJ(5) + t91 * qJD(5) + t11;
t1 = [(((2 * mrSges(3,3)) + t129 + (2 * t176)) * qJD(2) + (0.2e1 * t126 + (-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2)) * t118 * t120 + (-0.3e1 / 0.2e1 * t120 ^ 2 + 0.3e1 / 0.2e1 * t118 ^ 2) * Ifges(4,4)) * qJD(3)) * qJD(1) + m(6) * (t10 * t5 + t14 * t4 + t2 * t30 + t26 * t72 + t29 * t3 + t39 * t55) + t72 * t40 - t125 * mrSges(6,3) - t124 * mrSges(5,3) + t183 * t52 / 0.2e1 - t184 * t51 / 0.2e1 + (-t179 * t51 + t181 * t52) * t115 / 0.2e1 + (t155 * t52 - t180 * t51) * t83 / 0.2e1 + (-t155 * t51 + t182 * t52) * t169 + qJD(2) * t93 + t100 * (mrSges(5,1) * t51 + mrSges(5,2) * t52) + t102 * t50 + t55 * (mrSges(6,1) * t51 + mrSges(6,2) * t52) + t4 * t62 + t11 * t63 + t5 * t64 + t12 * t65 + t39 * t49 + (-t107 * mrSges(5,1) - t72 * mrSges(6,1) + t54 * mrSges(5,3) + t30 * mrSges(6,3) - t155 * t127 + t180 * t91) * t46 + ((-t81 / 0.2e1 + t121 * t98 - t148 / 0.2e1) * t120 + (-t82 / 0.2e1 - t121 * t99 - t149 / 0.2e1) * t118) * qJD(3) + (mrSges(5,2) * t107 - mrSges(5,3) * t53 - mrSges(6,3) * t29 - t127 * t182 + t155 * t91) * t45 + m(5) * (t100 * t102 + t107 * t95 - t11 * t128 + t12 * t27 + t53 * t8 + t54 * t7) + t26 * (-mrSges(6,1) * t91 - mrSges(6,2) * t127) + t95 * (-mrSges(5,1) * t91 - mrSges(5,2) * t127); -t177 * t52 + t178 * t51 + (-t118 * t99 + t120 * t98) * qJD(3) + m(5) * t124 + m(6) * t125 + (-mrSges(3,3) - t176) * t175 + (-m(5) * t100 - m(6) * t55 - t49 - t50 - t93) * qJD(1) + (mrSges(6,3) + mrSges(5,3)) * (t127 * t45 - t46 * t91); -t61 * t49 - t17 * t62 - t32 * t63 - t16 * t64 - t31 * t65 + t123 - t110 * t160 - t139 * Ifges(4,5) * t118 / 0.2e1 - t98 * t94 + t81 * t143 / 0.2e1 + t82 * t144 / 0.2e1 - t50 * t138 + t99 * t147 + t14 * t158 - t128 * t159 + mrSges(6,3) * t132 - Ifges(4,6) * t134 / 0.2e1 + (-mrSges(4,1) * t142 - mrSges(4,2) * t141) * t101 + (-t45 * t161 + t132) * mrSges(5,3) + (-t10 * t16 + t3 * t110 - t14 * t17 - t55 * t61) * m(6) + (-t100 * t138 + t128 * t32 - t27 * t31) * m(5) + (-t126 - t120 * (-Ifges(4,1) * t118 - t152) / 0.2e1 + t118 * (-Ifges(4,2) * t120 - t153) / 0.2e1) * t175 + (t178 * t140 + t177 * t135 + m(6) * (t164 * t2 + (-t10 * t164 + t119 * t14) * qJD(4)) + (t164 * t7 + t119 * t8 + (-t119 * t128 - t164 * t27) * qJD(4)) * m(5)) * pkin(3); (-t55 * t167 - (-t10 + t13) * t14 + t3 * pkin(4)) * m(6) + (-mrSges(5,3) * t128 + t14 * mrSges(6,3) - pkin(4) * t49) * t84 - pkin(4) * t160 - t13 * t62 - t27 * t63 + t14 * t64 - t128 * t65 + t123; -t46 * mrSges(6,1) - t83 * t62 + t84 * t64 + t40 + 0.2e1 * (t26 / 0.2e1 + t10 * t169 + t14 * t172) * m(6);];
tauc = t1(:);
