% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPPRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:29
% EndTime: 2019-03-08 19:21:38
% DurationCPUTime: 4.32s
% Computational Cost: add. (2726->358), mult. (5948->529), div. (0->0), fcn. (3835->10), ass. (0->179)
t103 = cos(qJ(2));
t94 = sin(pkin(6));
t171 = qJD(1) * t94;
t146 = t103 * t171;
t100 = sin(qJ(2));
t147 = t100 * t171;
t93 = sin(pkin(11));
t95 = cos(pkin(11));
t52 = t146 * t93 - t147 * t95;
t222 = qJD(3) * t93 - t52;
t102 = cos(qJ(5));
t99 = sin(qJ(5));
t133 = -pkin(5) * t99 + pkin(9) * t102;
t115 = t133 * qJD(5);
t104 = -pkin(2) - pkin(3);
t181 = t95 * qJ(3) + t93 * t104;
t77 = -pkin(8) + t181;
t221 = -qJD(6) * t102 * t77 + t115 + t222;
t170 = qJD(2) * t99;
t152 = mrSges(6,3) * t170;
t101 = cos(qJ(6));
t159 = qJD(5) * t101;
t98 = sin(qJ(6));
t74 = t170 * t98 + t159;
t162 = t98 * qJD(5);
t75 = t101 * t170 - t162;
t182 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t74 - mrSges(7,2) * t75 - t152;
t158 = qJD(5) * t102;
t164 = qJD(6) * t99;
t109 = t101 * t158 - t164 * t98;
t153 = qJD(5) * qJD(6);
t48 = -qJD(2) * t109 + t101 * t153;
t210 = t48 / 0.2e1;
t143 = t98 * t158;
t157 = qJD(6) * t101;
t110 = t157 * t99 + t143;
t49 = qJD(2) * t110 - t153 * t98;
t209 = t49 / 0.2e1;
t165 = qJD(5) * t99;
t168 = qJD(3) * t95;
t215 = t102 * t168 - t77 * t165;
t132 = pkin(5) * t102 + pkin(9) * t99;
t136 = -t93 * qJ(3) + t104 * t95;
t76 = pkin(4) - t136;
t56 = t132 + t76;
t108 = qJD(6) * t56 + t215;
t174 = t102 * t98;
t58 = (t100 * t93 + t103 * t95) * t94;
t55 = qJD(1) * t58;
t220 = t221 * t101 - t108 * t98 + t174 * t55;
t163 = t101 * t102;
t219 = t101 * t108 - t163 * t55 + t221 * t98;
t149 = t93 * t165;
t64 = t163 * t93 - t95 * t98;
t218 = -qJD(6) * t64 + t149 * t98 - (t101 * t93 - t174 * t95) * qJD(2);
t63 = -t101 * t95 - t174 * t93;
t217 = qJD(6) * t63 - t101 * t149 - (t163 * t95 + t93 * t98) * qJD(2);
t142 = qJD(2) * t171;
t72 = (qJD(3) + t146) * qJD(2);
t216 = -t103 * t142 + t72;
t119 = qJD(3) - t146;
t65 = qJD(2) * t104 + t119;
t81 = qJD(2) * qJ(3) + t147;
t36 = t93 * t65 + t95 * t81;
t31 = -qJD(2) * pkin(8) + t36;
t96 = cos(pkin(6));
t87 = -qJD(1) * t96 + qJD(4);
t24 = t102 * t31 + t87 * t99;
t166 = qJD(5) * t24;
t135 = t100 * t142;
t40 = t135 * t93 + t72 * t95;
t10 = t40 * t99 + t166;
t196 = t10 * t99;
t23 = t102 * t87 - t31 * t99;
t167 = qJD(5) * t23;
t9 = t102 * t40 + t167;
t214 = -t102 * t9 - t196;
t39 = -t95 * t135 + t72 * t93;
t27 = qJD(2) * t115 + t39;
t20 = qJD(5) * pkin(9) + t24;
t35 = t65 * t95 - t93 * t81;
t30 = qJD(2) * pkin(4) - t35;
t26 = qJD(2) * t132 + t30;
t5 = t101 * t26 - t20 * t98;
t1 = qJD(6) * t5 + t101 * t9 + t27 * t98;
t6 = t101 * t20 + t26 * t98;
t2 = -qJD(6) * t6 + t101 * t27 - t9 * t98;
t130 = t1 * t101 - t2 * t98;
t120 = Ifges(7,5) * t101 - Ifges(7,6) * t98;
t179 = Ifges(7,4) * t101;
t122 = -Ifges(7,2) * t98 + t179;
t199 = Ifges(7,4) * t98;
t124 = Ifges(7,1) * t101 - t199;
t126 = mrSges(7,1) * t98 + mrSges(7,2) * t101;
t203 = t6 * t98;
t129 = -t101 * t5 - t203;
t71 = Ifges(7,4) * t74;
t156 = t102 * qJD(2);
t88 = qJD(6) + t156;
t34 = -Ifges(7,1) * t75 + Ifges(7,5) * t88 + t71;
t178 = t101 * t34;
t19 = -qJD(5) * pkin(5) - t23;
t205 = t88 / 0.2e1;
t206 = -t75 / 0.2e1;
t207 = t74 / 0.2e1;
t200 = Ifges(7,4) * t75;
t33 = Ifges(7,2) * t74 + Ifges(7,6) * t88 - t200;
t212 = t120 * t205 + t122 * t207 + t124 * t206 + t19 * t126 + t129 * mrSges(7,3) + t178 / 0.2e1 - t98 * t33 / 0.2e1;
t154 = qJD(2) * qJD(5);
t141 = t99 * t154;
t211 = Ifges(7,4) * t210 + Ifges(7,2) * t209 - Ifges(7,6) * t141 / 0.2e1;
t202 = mrSges(5,2) * t95;
t201 = Ifges(6,4) * t99;
t59 = (t100 * t95 - t103 * t93) * t94;
t42 = t102 * t96 + t59 * t99;
t197 = t10 * t42;
t194 = t23 * mrSges(6,3);
t193 = t24 * mrSges(6,3);
t192 = t30 * t93;
t191 = t39 * t58;
t190 = t39 * t95;
t188 = t74 * Ifges(7,6);
t187 = t75 * Ifges(7,5);
t186 = t88 * Ifges(7,3);
t185 = t95 * t99;
t184 = t98 * t99;
t25 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t183 = t99 * t25;
t180 = Ifges(6,4) * t102;
t177 = t101 * t99;
t176 = t102 * t24;
t148 = mrSges(6,3) * t156;
t85 = -qJD(5) * mrSges(6,2) - t148;
t175 = t102 * t85;
t173 = Ifges(6,5) * qJD(5);
t172 = Ifges(6,6) * qJD(5);
t161 = qJD(2) * t100;
t160 = qJD(2) * t103;
t155 = qJD(2) * qJD(3);
t128 = -t35 * t93 + t36 * t95;
t127 = mrSges(6,1) * t102 - mrSges(6,2) * t99;
t125 = Ifges(7,1) * t98 + t179;
t123 = Ifges(7,2) * t101 + t199;
t121 = Ifges(7,5) * t98 + Ifges(7,6) * t101;
t37 = -mrSges(7,1) * t141 - mrSges(7,3) * t48;
t38 = mrSges(7,2) * t141 + mrSges(7,3) * t49;
t118 = t101 * t38 - t98 * t37;
t43 = t102 * t59 - t96 * t99;
t17 = t101 * t58 - t43 * t98;
t18 = t101 * t43 + t58 * t98;
t50 = -mrSges(7,2) * t88 + mrSges(7,3) * t74;
t51 = mrSges(7,1) * t88 + mrSges(7,3) * t75;
t117 = -t101 * t51 - t98 * t50;
t116 = -t23 * t99 + t176;
t114 = Ifges(7,5) * t48 + Ifges(7,6) * t49 - Ifges(7,3) * t141;
t113 = (-mrSges(6,1) * t99 - mrSges(6,2) * t102) * qJD(5);
t112 = (-t99 * Ifges(6,1) - t180) * qJD(2);
t111 = (-t102 * Ifges(6,2) - t201) * qJD(2);
t107 = (-t102 * t23 - t24 * t99) * qJD(5) - t214;
t105 = qJD(2) ^ 2;
t80 = t133 * qJD(2);
t78 = t127 * qJD(2);
t73 = -qJD(2) * pkin(2) + t119;
t70 = qJD(2) * t113;
t67 = t112 + t173;
t66 = t111 + t172;
t54 = qJD(2) * t58;
t53 = (t160 * t93 - t161 * t95) * t94;
t32 = t186 - t187 + t188;
t29 = t163 * t77 + t56 * t98;
t28 = t101 * t56 - t174 * t77;
t16 = -qJD(5) * t42 + t102 * t54;
t15 = t158 * t59 - t165 * t96 + t54 * t99;
t14 = t48 * Ifges(7,1) + t49 * Ifges(7,4) - Ifges(7,5) * t141;
t12 = t101 * t23 + t80 * t98;
t11 = t101 * t80 - t23 * t98;
t4 = qJD(6) * t17 + t101 * t16 + t53 * t98;
t3 = -qJD(6) * t18 + t101 * t53 - t16 * t98;
t7 = [t16 * t85 + t17 * t37 + t18 * t38 + t42 * t25 + t3 * t51 + t4 * t50 + t53 * t78 + t58 * t70 + t182 * t15 + (t53 * mrSges(5,1) + t54 * mrSges(5,2) + (-t102 * t42 + t43 * t99) * qJD(5) * mrSges(6,3)) * qJD(2) + m(6) * (-t15 * t23 + t16 * t24 + t30 * t53 + t43 * t9 + t191 + t197) + m(5) * (-t35 * t53 + t36 * t54 + t40 * t59 + t191) + m(7) * (t1 * t18 + t15 * t19 + t17 * t2 + t3 * t5 + t4 * t6 + t197) + (((-mrSges(3,2) + mrSges(4,3)) * t103 + (-mrSges(3,1) - mrSges(4,1)) * t100) * t105 + (t100 * t216 + t81 * t160 + t73 * t161) * m(4)) * t94; t102 * t114 / 0.2e1 + t19 * (-t110 * mrSges(7,1) - t109 * mrSges(7,2)) + (t39 * t76 + (t116 * t95 + t192) * qJD(3) + t107 * t77 - t116 * t55 - t30 * t52) * m(6) + t219 * t50 + t220 * t51 + (t155 + t216) * mrSges(4,3) + (qJD(3) * t128 - t136 * t39 + t181 * t40 + t35 * t52 - t36 * t55) * m(5) + (-qJD(2) * t52 + t155 * t93 + t39) * mrSges(5,1) + (m(7) * t19 + t182) * (t77 * t158 + (t168 - t55) * t99) + t214 * mrSges(6,3) + t215 * t85 + (-pkin(2) * t135 + qJ(3) * t72 + qJD(3) * t81 - (t100 * t73 + t103 * t81) * t171) * m(4) + t76 * t70 + t29 * t38 + t28 * t37 - t126 * t196 - (-Ifges(6,1) * t102 + t201) * t141 + t1 * (-mrSges(7,2) * t102 + mrSges(7,3) * t184) - t14 * t177 / 0.2e1 + t2 * (mrSges(7,1) * t102 + mrSges(7,3) * t177) - t55 * t175 + t5 * (-mrSges(7,1) * t165 + mrSges(7,3) * t109) + t6 * (mrSges(7,2) * t165 + mrSges(7,3) * t110) + (t1 * t29 + t77 * t196 + t2 * t28 + t219 * t6 + t220 * t5) * m(7) + t33 * t143 / 0.2e1 + t222 * t78 + (-qJD(2) * t55 + t40) * mrSges(5,2) + t39 * t127 + (t101 * t33 + t98 * t34) * t164 / 0.2e1 + (t111 + t66) * t165 / 0.2e1 - (qJD(2) * (t102 * Ifges(7,3) - t120 * t99) + t32) * t165 / 0.2e1 - (t178 + t112 + t67) * t158 / 0.2e1 + t155 * t202 + (t121 * t164 + (-Ifges(7,3) * t99 - t102 * t120) * qJD(5)) * t205 + (t125 * t164 + (-Ifges(7,5) * t99 - t102 * t124) * qJD(5)) * t206 + (t123 * t164 + (-Ifges(7,6) * t99 - t102 * t122) * qJD(5)) * t207 + (Ifges(7,6) * t102 - t122 * t99) * t209 + (Ifges(7,5) * t102 - t124 * t99) * t210 + t184 * t211 + t165 * t193 + t158 * t194 + t77 * t183 + t30 * t113 + qJD(5) ^ 2 * (-Ifges(6,5) * t102 + Ifges(6,6) * t99) / 0.2e1 - t102 * (Ifges(6,2) * t99 - t180) * t154; t63 * t37 + t64 * t38 - t95 * t70 + t218 * t51 + t217 * t50 + (-mrSges(4,3) - t202) * t105 + (-t105 * mrSges(5,1) + t183 + (t102 * t182 - t85 * t99) * qJD(5)) * t93 + m(6) * (t107 * t93 - t190) + m(5) * (t40 * t93 - t190) + (-t93 * t78 + (-t81 + t147) * m(4) + (-t182 * t99 - t175) * t95 - m(6) * (t176 * t95 - t185 * t23 + t192) - m(5) * t128) * qJD(2) + (t1 * t64 + t2 * t63 + (t158 * t19 + t196) * t93 - t19 * t185 * qJD(2) + t217 * t6 + t218 * t5) * m(7); (-t25 + (t101 * t50 - t51 * t98 + t148 + t85) * qJD(5) + m(6) * (-t10 + t166) + m(7) * (t159 * t6 - t162 * t5 - t10)) * t102 + (t117 * qJD(6) + (t152 + t182) * qJD(5) + m(6) * (t9 - t167) + m(7) * (qJD(5) * t19 - qJD(6) * t203 - t157 * t5 + t130) + t118) * t99; t125 * t210 + t123 * t209 + t101 * t211 + t98 * t14 / 0.2e1 - t23 * t85 - t12 * t50 - t11 * t51 - pkin(5) * t25 - t9 * mrSges(6,2) - t182 * t24 + (-mrSges(7,1) * t101 + mrSges(7,2) * t98 - mrSges(6,1)) * t10 + t130 * mrSges(7,3) + t212 * qJD(6) + ((-t66 / 0.2e1 + t32 / 0.2e1 + t186 / 0.2e1 - t187 / 0.2e1 + t188 / 0.2e1 + t5 * mrSges(7,1) - t6 * mrSges(7,2) + t30 * mrSges(6,1) + t172 / 0.2e1 - qJD(5) * t121 / 0.2e1 - t193 + Ifges(6,4) * t170 / 0.2e1) * t99 + (t67 / 0.2e1 + t30 * mrSges(6,2) - t173 / 0.2e1 - t194 + (-t180 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t99) * qJD(2) + t212) * t102) * qJD(2) + (-pkin(5) * t10 - t11 * t5 - t12 * t6 - t19 * t24) * m(7) + (t118 + m(7) * t130 + (m(7) * t129 + t117) * qJD(6)) * pkin(9); -t1 * mrSges(7,2) + t2 * mrSges(7,1) - t19 * (-mrSges(7,1) * t75 + mrSges(7,2) * t74) + t75 * (Ifges(7,1) * t74 + t200) / 0.2e1 + t33 * t206 - t88 * (Ifges(7,5) * t74 + Ifges(7,6) * t75) / 0.2e1 - t5 * t50 + t6 * t51 + (t5 * t74 - t6 * t75) * mrSges(7,3) + t114 - (Ifges(7,2) * t75 + t34 + t71) * t74 / 0.2e1;];
tauc  = t7(:);
