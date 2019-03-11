% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:14:05
% EndTime: 2019-03-08 19:14:14
% DurationCPUTime: 4.36s
% Computational Cost: add. (3925->337), mult. (10090->470), div. (0->0), fcn. (7839->12), ass. (0->166)
t114 = sin(pkin(11));
t121 = sin(qJ(2));
t115 = sin(pkin(6));
t154 = qJD(1) * t115;
t146 = t121 * t154;
t107 = t114 * t146;
t117 = cos(pkin(11));
t124 = cos(qJ(2));
t145 = t124 * t154;
t142 = t117 * t145;
t80 = -t107 + t142;
t219 = t80 - qJD(4);
t119 = sin(qJ(6));
t122 = cos(qJ(6));
t113 = sin(pkin(12));
t116 = cos(pkin(12));
t120 = sin(qJ(5));
t123 = cos(qJ(5));
t100 = t113 * t123 + t116 * t120;
t93 = t100 * qJD(2);
t75 = qJD(5) * t122 - t119 * t93;
t217 = -t75 / 0.2e1;
t76 = qJD(5) * t119 + t122 * t93;
t216 = -t76 / 0.2e1;
t215 = t93 / 0.2e1;
t99 = t113 * t120 - t123 * t116;
t92 = t99 * qJD(2);
t212 = qJD(6) + t92;
t214 = -t212 / 0.2e1;
t109 = pkin(2) * t114 + qJ(4);
t169 = pkin(8) + t109;
t96 = t169 * t113;
t97 = t169 * t116;
t57 = t120 * t97 + t123 * t96;
t205 = -qJD(5) * t57 + t219 * t99;
t86 = (t114 * t124 + t117 * t121) * t115;
t77 = qJD(1) * t86;
t94 = t99 * qJD(5);
t95 = t100 * qJD(5);
t213 = pkin(5) * t95 + pkin(9) * t94 - t77;
t155 = t113 ^ 2 + t116 ^ 2;
t203 = mrSges(5,3) * t155;
t58 = -t120 * t96 + t123 * t97;
t204 = qJD(5) * t58 - t100 * t219;
t161 = Ifges(6,5) * qJD(5);
t211 = t161 / 0.2e1 + Ifges(6,1) * t215;
t160 = Ifges(6,6) * qJD(5);
t210 = t160 / 0.2e1 + Ifges(7,3) * t214 + Ifges(7,6) * t217 + Ifges(7,5) * t216 + Ifges(6,4) * t215;
t87 = qJD(2) * t94;
t46 = qJD(6) * t75 - t122 * t87;
t193 = t46 / 0.2e1;
t47 = -qJD(6) * t76 + t119 * t87;
t192 = t47 / 0.2e1;
t88 = qJD(2) * t95;
t188 = t88 / 0.2e1;
t147 = -pkin(4) * t116 - pkin(3);
t180 = pkin(2) * t117;
t106 = t147 - t180;
t54 = pkin(5) * t99 - pkin(9) * t100 + t106;
t27 = t119 * t54 + t122 * t58;
t209 = -qJD(6) * t27 - t205 * t119 + t213 * t122;
t26 = -t119 * t58 + t122 * t54;
t208 = qJD(6) * t26 + t213 * t119 + t205 * t122;
t206 = Ifges(6,2) * t92;
t168 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t75 + mrSges(7,2) * t76 + t93 * mrSges(6,3);
t137 = -mrSges(5,1) * t116 + mrSges(5,2) * t113;
t202 = mrSges(6,1) * t92 + mrSges(6,2) * t93 + t137 * qJD(2);
t21 = -mrSges(7,1) * t47 + mrSges(7,2) * t46;
t201 = -t87 * mrSges(6,3) + t21;
t85 = (t114 * t121 - t117 * t124) * t115;
t186 = t122 / 0.2e1;
t184 = Ifges(7,4) * t76;
t29 = Ifges(7,2) * t75 + Ifges(7,6) * t212 + t184;
t70 = Ifges(7,4) * t75;
t30 = t76 * Ifges(7,1) + Ifges(7,5) * t212 + t70;
t200 = t30 * t186 - t119 * t29 / 0.2e1;
t78 = qJD(2) * t86;
t73 = qJD(1) * t78;
t36 = pkin(5) * t88 + pkin(9) * t87 + t73;
t118 = cos(pkin(6));
t108 = qJD(1) * t118 + qJD(3);
t105 = t116 * t108;
t165 = pkin(8) * qJD(2);
t103 = qJD(2) * pkin(2) + t145;
t72 = t114 * t103 + t117 * t146;
t68 = qJD(2) * qJ(4) + t72;
t44 = t105 + (-t68 - t165) * t113;
t51 = t113 * t108 + t116 * t68;
t45 = t116 * t165 + t51;
t20 = t120 * t44 + t123 * t45;
t18 = qJD(5) * pkin(9) + t20;
t71 = t103 * t117 - t107;
t138 = qJD(4) - t71;
t59 = qJD(2) * t147 + t138;
t33 = pkin(5) * t92 - pkin(9) * t93 + t59;
t5 = -t119 * t18 + t122 * t33;
t19 = -t120 * t45 + t123 * t44;
t102 = qJD(2) * t142;
t69 = t102 + (qJD(4) - t107) * qJD(2);
t9 = qJD(5) * t19 - t69 * t99;
t1 = qJD(6) * t5 + t119 * t36 + t122 * t9;
t6 = t119 * t33 + t122 * t18;
t2 = -qJD(6) * t6 - t119 * t9 + t122 * t36;
t141 = t1 * t122 - t119 * t2;
t133 = Ifges(7,5) * t122 - Ifges(7,6) * t119;
t166 = Ifges(7,4) * t122;
t134 = -Ifges(7,2) * t119 + t166;
t167 = Ifges(7,4) * t119;
t135 = Ifges(7,1) * t122 - t167;
t136 = mrSges(7,1) * t119 + mrSges(7,2) * t122;
t140 = t119 * t6 + t122 * t5;
t17 = -qJD(5) * pkin(5) - t19;
t187 = t212 / 0.2e1;
t189 = t76 / 0.2e1;
t190 = t75 / 0.2e1;
t198 = -t140 * mrSges(7,3) + t133 * t187 + t134 * t190 + t135 * t189 + t136 * t17 + t200;
t197 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t90 = Ifges(6,4) * t92;
t196 = t59 * mrSges(6,2) - t19 * mrSges(6,3) - t90 / 0.2e1 + t211;
t195 = t20 * mrSges(6,3) + t6 * mrSges(7,2) - t206 / 0.2e1 - t5 * mrSges(7,1) - t59 * mrSges(6,1) + t210;
t194 = Ifges(7,1) * t193 + Ifges(7,4) * t192 + Ifges(7,5) * t188;
t10 = qJD(5) * t20 + t100 * t69;
t66 = -t113 * t86 + t116 * t118;
t67 = t113 * t118 + t116 * t86;
t34 = t120 * t67 - t123 * t66;
t178 = t10 * t34;
t177 = t10 * t57;
t176 = t10 * t99;
t52 = t73 * t85;
t170 = t88 * mrSges(6,3);
t163 = t119 * t94;
t162 = t122 * t94;
t159 = qJD(2) * t80;
t158 = t100 * t119;
t157 = t100 * t122;
t153 = qJD(6) * t119;
t152 = qJD(6) * t122;
t151 = Ifges(7,5) * t46 + Ifges(7,6) * t47 + Ifges(7,3) * t88;
t53 = t88 * mrSges(6,1) - t87 * mrSges(6,2);
t50 = -t113 * t68 + t105;
t132 = -t113 * t50 + t116 * t51;
t31 = mrSges(7,1) * t88 - mrSges(7,3) * t46;
t32 = -mrSges(7,2) * t88 + mrSges(7,3) * t47;
t131 = -t119 * t31 + t122 * t32;
t35 = t120 * t66 + t123 * t67;
t23 = t119 * t85 + t122 * t35;
t22 = -t119 * t35 + t122 * t85;
t48 = -mrSges(7,2) * t212 + mrSges(7,3) * t75;
t49 = mrSges(7,1) * t212 - mrSges(7,3) * t76;
t130 = -t119 * t48 - t122 * t49;
t128 = t100 * t152 - t163;
t127 = t100 * t153 + t162;
t125 = qJD(2) ^ 2;
t81 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t92;
t79 = qJD(2) * t85;
t74 = -qJD(2) * t107 + t102;
t65 = -qJD(2) * pkin(3) + t138;
t61 = pkin(5) * t93 + pkin(9) * t92;
t15 = Ifges(7,4) * t46 + Ifges(7,2) * t47 + Ifges(7,6) * t88;
t14 = qJD(5) * t35 - t100 * t79;
t13 = -qJD(5) * t34 + t79 * t99;
t12 = t119 * t61 + t122 * t19;
t11 = -t119 * t19 + t122 * t61;
t4 = -qJD(6) * t23 - t119 * t13 + t122 * t78;
t3 = qJD(6) * t22 + t119 * t78 + t122 * t13;
t7 = [t13 * t81 + t34 * t21 + t22 * t31 + t23 * t32 + t3 * t48 + t4 * t49 + t85 * t53 + t202 * t78 + t168 * t14 + (-mrSges(3,1) * t121 - mrSges(3,2) * t124) * t125 * t115 + (-t34 * t87 - t35 * t88) * mrSges(6,3) + (-t78 * mrSges(4,1) - (-mrSges(4,2) + t203) * t79) * qJD(2) + m(6) * (t13 * t20 - t14 * t19 + t35 * t9 + t59 * t78 + t178 + t52) + m(5) * (t65 * t78 + t52 + (-t51 * t79 + t67 * t69) * t116 + (t50 * t79 - t66 * t69) * t113) + m(7) * (t1 * t23 + t14 * t17 + t2 * t22 + t3 * t6 + t4 * t5 + t178) + m(4) * (-t71 * t78 - t72 * t79 + t74 * t86 + t52); ((t114 * t74 - t117 * t73) * pkin(2) - t72 * t80) * m(4) + (t159 - t74) * mrSges(4,2) + t27 * t32 + t26 * t31 + (t137 - mrSges(4,1)) * t73 + (-Ifges(7,5) * t127 - Ifges(7,6) * t128) * t187 + t17 * (mrSges(7,1) * t128 - mrSges(7,2) * t127) - (t196 + t200 + t211) * t94 + t204 * t168 + (-Ifges(7,1) * t127 - Ifges(7,4) * t128) * t189 + (t92 * t94 / 0.2e1 - t93 * t95 / 0.2e1 + t99 * t87) * Ifges(6,4) + (-(t119 * t30 + t122 * t29) * qJD(6) / 0.2e1 - Ifges(6,4) * t88 + t73 * mrSges(6,2) - Ifges(6,1) * t87 + t133 * t188 + t134 * t192 + t135 * t193 + (t136 + mrSges(6,3)) * t10) * t100 + (t73 * (-pkin(3) - t180) + t155 * t69 * t109 - t219 * t132) * m(5) + t208 * t48 + t209 * t49 + (t1 * t27 + t17 * t204 + t2 * t26 + t208 * t6 + t209 * t5 + t177) * m(7) + (t206 / 0.2e1 - t160 / 0.2e1 + Ifges(7,3) * t187 + Ifges(7,5) * t189 + Ifges(7,6) * t190 - t195) * t95 + t205 * t81 + (t106 * t73 - t19 * t204 + t20 * t205 + t58 * t9 + t177) * m(6) + (-t1 * t158 + t127 * t5 - t128 * t6 - t157 * t2) * mrSges(7,3) + t201 * t57 + (m(4) * t71 - m(5) * t65 - m(6) * t59 + mrSges(4,1) * qJD(2) - t202) * t77 + (-t9 * mrSges(6,3) + Ifges(6,2) * t88 + t73 * mrSges(6,1) + Ifges(7,3) * t188 + Ifges(7,6) * t192 + Ifges(7,5) * t193 + t151 / 0.2e1 + t197) * t99 + (-Ifges(7,4) * t127 - Ifges(7,2) * t128) * t190 + t106 * t53 + t157 * t194 + (qJD(2) * qJD(4) - t159 + t69) * t203 - t15 * t158 / 0.2e1 - t58 * t170; t201 * t99 + t168 * t95 - (-t119 * t49 + t122 * t48 + t81) * t94 + m(7) * (-t162 * t6 + t163 * t5 + t17 * t95 + t176) + m(6) * (-t19 * t95 - t20 * t94 + t176) + (-t170 + t130 * qJD(6) + m(7) * (-t152 * t5 - t153 * t6 + t141) + m(6) * t9 + t131) * t100; t92 * t81 - t168 * t93 - t125 * t203 + (t212 * t48 + t31) * t122 + (-t212 * t49 + t32) * t119 + t53 + (t1 * t119 + t122 * t2 - t17 * t93 + t212 * (-t119 * t5 + t122 * t6)) * m(7) + (t19 * t93 + t20 * t92 + t73) * m(6) + (-qJD(2) * t132 + t73) * m(5); -pkin(5) * t21 - t9 * mrSges(6,2) - t12 * t48 - t11 * t49 - t19 * t81 - Ifges(6,5) * t87 - Ifges(6,6) * t88 + t119 * t194 + t15 * t186 + (Ifges(7,1) * t119 + t166) * t193 + (Ifges(7,2) * t122 + t167) * t192 + (Ifges(7,5) * t119 + Ifges(7,6) * t122) * t188 - t168 * t20 + (-mrSges(7,1) * t122 + mrSges(7,2) * t119 - mrSges(6,1)) * t10 + t141 * mrSges(7,3) + (t195 + t210) * t93 - (-t161 / 0.2e1 + t90 / 0.2e1 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t93 - t196 - t198) * t92 + t198 * qJD(6) + (-pkin(5) * t10 - t11 * t5 - t12 * t6 - t17 * t20) * m(7) + (t131 + m(7) * t141 + (-m(7) * t140 + t130) * qJD(6)) * pkin(9); -t17 * (mrSges(7,1) * t76 + mrSges(7,2) * t75) + (Ifges(7,1) * t75 - t184) * t216 + t29 * t189 + (Ifges(7,5) * t75 - Ifges(7,6) * t76) * t214 - t5 * t48 + t6 * t49 + (t5 * t75 + t6 * t76) * mrSges(7,3) + t151 + (-Ifges(7,2) * t76 + t30 + t70) * t217 + t197;];
tauc  = t7(:);
