% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:42:04
% EndTime: 2019-03-08 18:42:11
% DurationCPUTime: 2.72s
% Computational Cost: add. (4088->343), mult. (11337->503), div. (0->0), fcn. (9720->14), ass. (0->177)
t210 = qJD(5) / 0.2e1;
t116 = sin(qJ(5));
t163 = qJD(3) * t116;
t209 = t163 / 0.2e1;
t147 = Ifges(6,5) * t210;
t119 = cos(qJ(5));
t162 = qJD(3) * t119;
t106 = Ifges(6,4) * t162;
t102 = qJD(6) - t162;
t115 = sin(qJ(6));
t118 = cos(qJ(6));
t179 = Ifges(7,4) * t118;
t136 = -Ifges(7,2) * t115 + t179;
t180 = Ifges(7,4) * t115;
t138 = Ifges(7,1) * t118 - t180;
t139 = mrSges(7,1) * t115 + mrSges(7,2) * t118;
t107 = sin(pkin(13));
t111 = cos(pkin(13));
t114 = cos(pkin(6));
t100 = qJD(1) * t114 + qJD(2);
t109 = sin(pkin(7));
t117 = sin(qJ(3));
t171 = t109 * t117;
t108 = sin(pkin(12));
t110 = sin(pkin(6));
t120 = cos(qJ(3));
t112 = cos(pkin(12));
t113 = cos(pkin(7));
t169 = t112 * t113;
t207 = (t108 * t120 + t117 * t169) * t110;
t58 = qJD(1) * t207 + t100 * t171;
t53 = t111 * t58;
t164 = qJD(1) * t110;
t149 = t112 * t164;
t172 = t108 * t117;
t170 = t109 * t120;
t92 = t100 * t170;
t57 = t113 * t120 * t149 - t164 * t172 + t92;
t55 = qJD(3) * pkin(3) + t57;
t32 = t107 * t55 + t53;
t30 = qJD(3) * pkin(9) + t32;
t77 = t100 * t113 - t109 * t149 + qJD(4);
t22 = t116 * t77 + t119 * t30;
t20 = qJD(5) * pkin(10) + t22;
t129 = -pkin(5) * t119 - pkin(10) * t116 - pkin(4);
t52 = t107 * t58;
t31 = t111 * t55 - t52;
t24 = qJD(3) * t129 - t31;
t5 = -t115 * t20 + t118 * t24;
t6 = t115 * t24 + t118 * t20;
t141 = t115 * t6 + t118 * t5;
t177 = Ifges(7,6) * t115;
t178 = Ifges(7,5) * t118;
t21 = -t116 * t30 + t119 * t77;
t19 = -qJD(5) * pkin(5) - t21;
t192 = t118 / 0.2e1;
t193 = -t115 / 0.2e1;
t161 = qJD(5) * t115;
t94 = t118 * t163 + t161;
t198 = t94 / 0.2e1;
t191 = Ifges(7,4) * t94;
t160 = qJD(5) * t118;
t93 = -t115 * t163 + t160;
t60 = Ifges(7,2) * t93 + Ifges(7,6) * t102 + t191;
t90 = Ifges(7,4) * t93;
t61 = Ifges(7,1) * t94 + Ifges(7,5) * t102 + t90;
t122 = -t141 * mrSges(7,3) + t60 * t193 + t61 * t192 + t19 * t139 + t93 * t136 / 0.2e1 + t138 * t198 + t102 * (-t177 + t178) / 0.2e1;
t29 = -qJD(3) * pkin(4) - t31;
t208 = t122 + t29 * mrSges(6,2) - t21 * mrSges(6,3) + Ifges(6,1) * t209 + t106 / 0.2e1 + t147;
t67 = t114 * t171 + t207;
t146 = -Ifges(6,6) * qJD(5) / 0.2e1;
t83 = (t107 * t117 - t111 * t120) * t109;
t127 = (t120 * t169 - t172) * t110;
t50 = (qJD(1) * t127 + t92) * qJD(3);
t51 = t58 * qJD(3);
t27 = t107 * t50 + t111 * t51;
t143 = pkin(5) * t116 - pkin(10) * t119;
t97 = t143 * qJD(5);
t23 = qJD(3) * t97 + t27;
t166 = t21 * qJD(5);
t28 = -t107 * t51 + t111 * t50;
t7 = t119 * t28 + t166;
t1 = qJD(6) * t5 + t115 * t23 + t118 * t7;
t2 = -qJD(6) * t6 - t115 * t7 + t118 * t23;
t142 = t1 * t118 - t115 * t2;
t155 = qJD(5) * qJD(6);
t158 = qJD(6) * t115;
t159 = qJD(5) * t119;
t75 = t118 * t155 + (-t116 * t158 + t118 * t159) * qJD(3);
t157 = qJD(6) * t118;
t76 = -t115 * t155 + (-t115 * t159 - t116 * t157) * qJD(3);
t206 = -t2 * mrSges(7,1) + t1 * mrSges(7,2) - Ifges(7,5) * t75 - Ifges(7,6) * t76;
t204 = 2 * m(6);
t34 = t111 * t57 - t52;
t203 = -t34 / 0.2e1;
t202 = t75 / 0.2e1;
t201 = t76 / 0.2e1;
t200 = -t93 / 0.2e1;
t199 = -t94 / 0.2e1;
t66 = t114 * t170 + t127;
t38 = t107 * t66 + t111 * t67;
t85 = -t109 * t110 * t112 + t113 * t114;
t25 = t116 * t38 - t119 * t85;
t165 = t22 * qJD(5);
t8 = t116 * t28 + t165;
t197 = t25 * t8;
t84 = (t107 * t120 + t111 * t117) * t109;
t131 = t113 * t119 - t116 * t84;
t196 = t131 * t8;
t195 = -t102 / 0.2e1;
t104 = pkin(3) * t107 + pkin(9);
t194 = t104 / 0.2e1;
t190 = pkin(3) * t111;
t37 = t107 * t67 - t111 * t66;
t186 = t27 * t37;
t185 = t27 * t83;
t151 = mrSges(6,3) * t163;
t181 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t93 + mrSges(7,2) * t94 + t151;
t173 = qJD(5) * mrSges(6,3);
t168 = t115 * t119;
t167 = t118 * t119;
t156 = qJD(5) * qJD(3);
t154 = t8 * t194;
t150 = mrSges(6,3) * t162;
t148 = qJD(5) * t104 * t116;
t145 = t116 * t156;
t140 = mrSges(7,1) * t118 - mrSges(7,2) * t115;
t137 = Ifges(7,1) * t115 + t179;
t135 = Ifges(7,2) * t118 + t180;
t134 = Ifges(7,5) * t115 + Ifges(7,6) * t118;
t26 = t116 * t85 + t119 * t38;
t14 = t115 * t37 + t118 * t26;
t13 = -t115 * t26 + t118 * t37;
t64 = mrSges(7,1) * t145 - mrSges(7,3) * t75;
t65 = -mrSges(7,2) * t145 + mrSges(7,3) * t76;
t133 = -t115 * t64 + t118 * t65;
t72 = t113 * t116 + t119 * t84;
t46 = t115 * t83 + t118 * t72;
t45 = -t115 * t72 + t118 * t83;
t78 = -mrSges(7,2) * t102 + mrSges(7,3) * t93;
t79 = mrSges(7,1) * t102 - mrSges(7,3) * t94;
t132 = -t115 * t78 - t118 * t79;
t89 = t129 - t190;
t69 = t104 * t167 + t115 * t89;
t68 = -t104 * t168 + t118 * t89;
t124 = t29 * mrSges(6,1) + t5 * mrSges(7,1) + t102 * Ifges(7,3) + t94 * Ifges(7,5) + t93 * Ifges(7,6) + t146 - (Ifges(6,4) * t116 + t119 * Ifges(6,2)) * qJD(3) / 0.2e1 - t22 * mrSges(6,3) - t6 * mrSges(7,2);
t105 = -pkin(4) - t190;
t101 = Ifges(7,3) * t145;
t99 = -qJD(5) * mrSges(6,2) + t150;
t96 = t143 * qJD(3);
t95 = (-mrSges(6,1) * t119 + mrSges(6,2) * t116) * qJD(3);
t88 = (mrSges(6,1) * t116 + mrSges(6,2) * t119) * t156;
t81 = qJD(3) * t83;
t80 = qJD(3) * t84;
t63 = t67 * qJD(3);
t62 = t66 * qJD(3);
t47 = -mrSges(7,1) * t76 + mrSges(7,2) * t75;
t44 = qJD(5) * t72 - t116 * t81;
t43 = qJD(5) * t131 - t119 * t81;
t42 = Ifges(7,1) * t75 + Ifges(7,4) * t76 + Ifges(7,5) * t145;
t41 = Ifges(7,4) * t75 + Ifges(7,2) * t76 + Ifges(7,6) * t145;
t40 = -qJD(6) * t69 + t115 * t148 + t118 * t97;
t39 = qJD(6) * t68 + t115 * t97 - t118 * t148;
t36 = -t107 * t63 + t111 * t62;
t35 = t107 * t62 + t111 * t63;
t33 = t107 * t57 + t53;
t18 = t115 * t96 + t118 * t21;
t17 = -t115 * t21 + t118 * t96;
t16 = -qJD(6) * t46 - t115 * t43 + t118 * t80;
t15 = qJD(6) * t45 + t115 * t80 + t118 * t43;
t12 = t115 * t33 + t167 * t34;
t11 = t118 * t33 - t168 * t34;
t10 = qJD(5) * t26 + t116 * t36;
t9 = -qJD(5) * t25 + t119 * t36;
t4 = -qJD(6) * t14 - t115 * t9 + t118 * t35;
t3 = qJD(6) * t13 + t115 * t35 + t118 * t9;
t48 = [t13 * t64 + t14 * t65 + t25 * t47 + t3 * t78 + t35 * t95 + t37 * t88 + t4 * t79 + t9 * t99 + t181 * t10 + m(4) * (t50 * t67 - t51 * t66 - t57 * t63 + t58 * t62) + m(6) * (-t10 * t21 + t22 * t9 + t26 * t7 + t29 * t35 + t186 + t197) + m(5) * (t28 * t38 - t31 * t35 + t32 * t36 + t186) + m(7) * (t1 * t14 + t10 * t19 + t13 * t2 + t3 * t6 + t4 * t5 + t197) + (-t63 * mrSges(4,1) - t35 * mrSges(5,1) - t62 * mrSges(4,2) - t36 * mrSges(5,2) + (-t116 * t26 + t119 * t25) * t173) * qJD(3); t15 * t78 + t16 * t79 + t43 * t99 + t45 * t64 + t46 * t65 - t131 * t47 + t80 * t95 + t83 * t88 + t181 * t44 + m(6) * (-t21 * t44 + t22 * t43 + t29 * t80 + t7 * t72 + t185 - t196) + m(5) * (t28 * t84 - t31 * t80 - t32 * t81 + t185) + m(7) * (t1 * t46 + t15 * t6 + t16 * t5 + t19 * t44 + t2 * t45 - t196) + (-t80 * mrSges(5,1) + t81 * mrSges(5,2) + (-t116 * t72 - t119 * t131) * t173) * qJD(3) + (m(4) * (t117 * t50 - t120 * t51 + (-t117 * t57 + t120 * t58) * qJD(3)) + (-mrSges(4,1) * t117 - mrSges(4,2) * t120) * qJD(3) ^ 2) * t109; -t51 * mrSges(4,1) - t27 * mrSges(5,1) - t50 * mrSges(4,2) - t28 * mrSges(5,2) + t105 * t88 - t33 * t95 + t68 * t64 + t69 * t65 + (-t11 + t40) * t79 + (-t12 + t39) * t78 + (t58 * mrSges(4,1) + t33 * mrSges(5,1) + t57 * mrSges(4,2) + t34 * mrSges(5,2)) * qJD(3) - m(7) * (t11 * t5 + t12 * t6) + m(7) * (t1 * t69 + t2 * t68 + t39 * t6 + t40 * t5) + (-t29 * t33 / 0.2e1 + t105 * t27 / 0.2e1) * t204 + (t7 * mrSges(6,3) - t34 * t99 - t27 * mrSges(6,1) - t101 / 0.2e1 + (t194 * t7 + t203 * t22) * t204 + (0.3e1 / 0.2e1 * t106 + t147 + (-m(6) * t21 + m(7) * t19 + t181) * t104 + t208) * qJD(5) + t206) * t119 + (t104 * t47 + t41 * t193 + t42 * t192 + t27 * mrSges(6,2) + t138 * t202 + t136 * t201 + (mrSges(6,3) + t139) * t8 - t181 * t34 + (-t1 * t115 - t2 * t118) * mrSges(7,3) + 0.2e1 * (t19 * t203 + t154) * m(7) + (t154 + t21 * t34 / 0.2e1) * t204 + (t61 * t193 - t118 * t60 / 0.2e1 + t19 * t140 + t135 * t200 + t137 * t199 + t134 * t195 + (t115 * t5 - t118 * t6) * mrSges(7,3)) * qJD(6) + (t146 + (-m(6) * t22 - t99) * t104 + ((-0.3e1 / 0.2e1 * Ifges(6,4) + t178 / 0.2e1 - t177 / 0.2e1) * t116 + (0.3e1 / 0.2e1 * Ifges(6,1) - 0.3e1 / 0.2e1 * Ifges(6,2) - Ifges(7,3) / 0.2e1) * t119) * qJD(3) + t124) * qJD(5)) * t116 + (t31 * t33 - t32 * t34 + (t107 * t28 - t111 * t27) * pkin(3)) * m(5); (-t47 + (-t115 * t79 + t118 * t78 - t150 + t99) * qJD(5) + m(6) * (-t8 + t165) + m(7) * (t160 * t6 - t161 * t5 - t8)) * t119 + (t132 * qJD(6) + (-t151 + t181) * qJD(5) + m(6) * (t7 - t166) + m(7) * (qJD(5) * t19 - t157 * t5 - t158 * t6 + t142) + t133) * t116; -t7 * mrSges(6,2) - pkin(5) * t47 - t18 * t78 - t17 * t79 - t21 * t99 + t115 * t42 / 0.2e1 + t41 * t192 + t137 * t202 + t135 * t201 + (-mrSges(6,1) - t140) * t8 - t181 * t22 + t142 * mrSges(7,3) + t122 * qJD(6) + ((Ifges(6,4) * t209 + t134 * t210 - t124 + t146) * t116 + (t147 - t106 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t163 - t208) * t119) * qJD(3) + (-pkin(5) * t8 - t17 * t5 - t18 * t6 - t19 * t22) * m(7) + (t133 + m(7) * t142 + (-m(7) * t141 + t132) * qJD(6)) * pkin(10); t101 - t19 * (mrSges(7,1) * t94 + mrSges(7,2) * t93) + (Ifges(7,1) * t93 - t191) * t199 + t60 * t198 + (Ifges(7,5) * t93 - Ifges(7,6) * t94) * t195 - t5 * t78 + t6 * t79 + (t5 * t93 + t6 * t94) * mrSges(7,3) + (-Ifges(7,2) * t94 + t61 + t90) * t200 - t206;];
tauc  = t48(:);
