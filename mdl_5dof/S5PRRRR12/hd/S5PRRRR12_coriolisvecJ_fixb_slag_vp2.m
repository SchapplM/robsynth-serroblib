% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:21
% EndTime: 2024-09-28 18:07:23
% DurationCPUTime: 2.08s
% Computational Cost: add. (3370->232), mult. (8389->366), div. (0->0), fcn. (6513->12), ass. (0->131)
t95 = qJD(2) + qJD(3);
t92 = qJD(4) + t95;
t96 = sin(pkin(6));
t160 = t92 * t96;
t101 = sin(qJ(4));
t102 = sin(qJ(3));
t105 = cos(qJ(4));
t106 = cos(qJ(3));
t136 = qJD(4) * t105;
t137 = qJD(4) * t101;
t140 = t102 * t105;
t103 = sin(qJ(2));
t107 = cos(qJ(2));
t121 = t102 * t107 + t103 * t106;
t97 = sin(pkin(5));
t144 = qJD(1) * t97;
t70 = t121 * t144;
t120 = -t102 * t103 + t106 * t107;
t73 = t120 * t97;
t71 = qJD(1) * t73;
t91 = pkin(2) * t106 + pkin(3);
t159 = -t101 * t71 - t105 * t70 + t91 * t137 - (-t102 * t136 + (-t101 * t106 - t140) * qJD(3)) * pkin(2);
t141 = t101 * t102;
t169 = -t101 * t70 + t105 * t71 - t91 * t136 - (-t102 * t137 + (t105 * t106 - t141) * qJD(3)) * pkin(2);
t104 = cos(qJ(5));
t151 = t104 * Ifges(6,2);
t100 = sin(qJ(5));
t157 = Ifges(6,4) * t100;
t174 = (t151 + t157) * t160;
t99 = cos(pkin(5));
t143 = qJD(1) * t99;
t131 = t103 * t144;
t84 = qJD(2) * pkin(2) + t107 * t144;
t64 = t102 * t84 + t106 * t131;
t152 = t101 * t64;
t63 = -t102 * t131 + t106 * t84;
t59 = pkin(3) * t95 + t63;
t28 = t105 * t59 - t152;
t27 = pkin(4) * t92 + t28;
t98 = cos(pkin(6));
t119 = t143 * t96 + t27 * t98;
t146 = t105 * t64;
t29 = t101 * t59 + t146;
t26 = pkin(10) * t160 + t29;
t10 = -t100 * t26 + t104 * t119;
t153 = t100 * t98;
t114 = t120 * qJD(2);
t138 = qJD(3) * t106;
t139 = qJD(3) * t102;
t34 = t84 * t138 + (-t103 * t139 + t114) * t144;
t115 = t121 * qJD(2);
t35 = -t84 * t139 + (-t103 * t138 - t115) * t144;
t8 = t28 * qJD(4) + t101 * t35 + t105 * t34;
t9 = -qJD(4) * t29 - t101 * t34 + t105 * t35;
t3 = qJD(5) * t10 + t104 * t8 + t153 * t9;
t11 = t100 * t119 + t104 * t26;
t147 = t104 * t98;
t4 = -t11 * qJD(5) - t100 * t8 + t9 * t147;
t173 = t4 * mrSges(6,1) - t3 * mrSges(6,2);
t158 = pkin(2) * t140 + t101 * t91;
t93 = t96 * pkin(10);
t72 = t93 + t158;
t166 = -pkin(2) * t141 + t105 * t91;
t77 = pkin(4) + t166;
t39 = t104 * t72 + t153 * t77;
t171 = -qJD(5) * t39 + t169 * t100 - t159 * t147;
t38 = -t100 * t72 + t147 * t77;
t170 = qJD(5) * t38 - t169 * t104 - t159 * t153;
t155 = pkin(3) * qJD(4);
t30 = -t101 * t63 - t146;
t31 = t105 * t63 - t152;
t86 = pkin(3) * t101 + t93;
t90 = pkin(3) * t105 + pkin(4);
t65 = -t100 * t86 + t147 * t90;
t168 = -t104 * t31 - t153 * t30 + t65 * qJD(5) + (-t101 * t153 + t104 * t105) * t155;
t66 = t104 * t86 + t153 * t90;
t167 = t100 * t31 - t147 * t30 - t66 * qJD(5) + (-t100 * t105 - t101 * t147) * t155;
t125 = pkin(3) * t137 + t30;
t74 = t121 * t97;
t44 = -t101 * t74 + t105 * t73;
t124 = t44 * t98 + t96 * t99;
t45 = t101 * t73 + t105 * t74;
t23 = t100 * t124 + t104 * t45;
t94 = t96 ^ 2;
t78 = (-mrSges(6,1) * t104 + mrSges(6,2) * t100) * t96;
t165 = t9 * t78;
t164 = t9 * t94;
t163 = t100 / 0.2e1;
t25 = t143 * t98 - t27 * t96;
t161 = t25 * t96;
t156 = Ifges(6,4) * t104;
t154 = t100 * t96;
t148 = t104 * t96;
t145 = mrSges(6,3) * qJD(5);
t142 = qJD(5) * t96;
t135 = qJD(5) * t100;
t134 = t92 * t154;
t133 = t92 * t148;
t132 = mrSges(6,3) * t142;
t130 = t96 * t135;
t123 = mrSges(6,1) * t100 + mrSges(6,2) * t104;
t122 = t10 * t104 + t100 * t11;
t85 = t92 * t98 + qJD(5);
t118 = t85 * (Ifges(6,5) * t104 - Ifges(6,6) * t100);
t79 = pkin(4) * t147 - pkin(10) * t154;
t80 = pkin(4) * t153 + pkin(10) * t148;
t116 = t100 * (Ifges(6,1) * t104 - t157);
t113 = t123 * t142;
t22 = -t100 * t45 + t104 * t124;
t48 = Ifges(6,6) * t85 + t174;
t83 = Ifges(6,4) * t133;
t49 = Ifges(6,1) * t134 + Ifges(6,5) * t85 + t83;
t82 = Ifges(6,5) * qJD(5) * t133;
t110 = t25 * t113 + t9 * mrSges(5,1) + (t116 + t104 * (-Ifges(6,2) * t100 + t156)) * qJD(5) * t92 * t94 - (t48 + t174) * t130 / 0.2e1 + (t3 * t148 - t4 * t154) * mrSges(6,3) + (-Ifges(6,6) * t130 * t92 + t82 / 0.2e1 + t173) * t98 + (t118 + (t49 + t92 * (Ifges(6,5) * t98 + (t100 * Ifges(6,1) + t156) * t96)) * t104) * t142 / 0.2e1;
t109 = t35 * mrSges(4,1) - t8 * mrSges(5,2) + t110;
t76 = t80 * qJD(5);
t75 = t79 * qJD(5);
t69 = t92 * t78;
t68 = -mrSges(6,2) * t85 + mrSges(6,3) * t133;
t67 = mrSges(6,1) * t85 - mrSges(6,3) * t134;
t62 = t92 * t113;
t47 = (-qJD(3) * t121 - t115) * t97;
t46 = (qJD(3) * t120 + t114) * t97;
t32 = -t44 * t96 + t98 * t99;
t15 = -qJD(4) * t45 - t101 * t46 + t105 * t47;
t14 = qJD(4) * t44 + t101 * t47 + t105 * t46;
t13 = t104 * t28 - t153 * t29;
t12 = -t100 * t28 - t147 * t29;
t6 = -t23 * qJD(5) - t100 * t14 + t15 * t147;
t5 = qJD(5) * t22 + t104 * t14 + t15 * t153;
t1 = [-t96 * t15 * t69 + t32 * t62 + t5 * t68 + t6 * t67 + (mrSges(4,1) * t47 - mrSges(4,2) * t46) * t95 + (-mrSges(3,1) * t103 - mrSges(3,2) * t107) * t97 * qJD(2) ^ 2 + m(4) * (t34 * t74 + t35 * t73 + t46 * t64 + t47 * t63) + m(5) * (t14 * t29 + t15 * t28 + t44 * t9 + t45 * t8) + m(6) * (t10 * t6 + t11 * t5 + t22 * t4 + t23 * t3 + (-t15 * t25 - t32 * t9) * t96) + (mrSges(5,1) * t15 - mrSges(5,2) * t14 + (-t100 * t23 - t104 * t22) * t132) * t92; t170 * t68 + t171 * t67 - t34 * mrSges(4,2) + (-t77 * t62 - t165 + t159 * t69 + ((-t38 * t92 - t10) * t104 + (-t39 * t92 - t11) * t100) * t145) * t96 + (-t159 * mrSges(5,1) + t169 * mrSges(5,2)) * t92 + ((-pkin(2) * t138 + t71) * mrSges(4,2) + (-pkin(2) * t139 + t70) * mrSges(4,1)) * t95 + t109 + (t171 * t10 + t170 * t11 + t159 * t161 + t164 * t77 + t3 * t39 + t38 * t4) * m(6) + (t63 * t70 - t64 * t71 + (t102 * t34 + t106 * t35 + t64 * t138 - t63 * t139) * pkin(2)) * m(4) + (t158 * t8 - t159 * t28 + t166 * t9 - t169 * t29) * m(5); t168 * t68 + t167 * t67 + t64 * t95 * mrSges(4,1) + ((-t100 * t66 - t104 * t65) * t132 + (-pkin(3) * t136 + t31) * mrSges(5,2) - t125 * mrSges(5,1)) * t92 + (-t122 * t145 + t125 * t69 - t90 * t62 - t165) * t96 + t109 + (t63 * t95 - t34) * mrSges(4,2) + ((t101 * t8 + t105 * t9 + t29 * t136 - t28 * t137) * pkin(3) - t28 * t30 - t29 * t31) * m(5) + (t167 * t10 + t168 * t11 + t125 * t161 + t90 * t164 + t3 * t66 + t4 * t65) * m(6); (t28 * t92 - t8) * mrSges(5,2) + (-t76 - t12) * t67 + (t75 - t13) * t68 + m(6) * (pkin(4) * t164 - t10 * t76 + t11 * t75 + t3 * t80 + t4 * t79) + t29 * t92 * mrSges(5,1) + (-pkin(4) * t62 - t165 + (-m(6) * t25 - t69) * t29 + ((-t79 * t92 - t10) * t104 + (-t80 * t92 - t11) * t100) * t145) * t96 + t110 - m(6) * (t10 * t12 + t11 * t13); -t10 * t68 + t11 * t67 + t82 + (-Ifges(6,6) * t135 - t25 * t123 + t48 * t163 - t118 / 0.2e1 + (-t116 / 0.2e1 + t151 * t163) * t160 + t122 * mrSges(6,3) - (t49 + t83) * t104 / 0.2e1) * t160 + t173;];
tauc = t1(:);
