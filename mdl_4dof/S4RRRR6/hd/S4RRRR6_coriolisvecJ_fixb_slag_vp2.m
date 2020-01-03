% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:22
% EndTime: 2019-12-31 17:29:34
% DurationCPUTime: 4.80s
% Computational Cost: add. (3641->436), mult. (10182->645), div. (0->0), fcn. (7334->8), ass. (0->201)
t137 = sin(qJ(4));
t140 = cos(qJ(4));
t158 = Ifges(5,5) * t140 - Ifges(5,6) * t137;
t194 = Ifges(5,4) * t140;
t160 = -Ifges(5,2) * t137 + t194;
t195 = Ifges(5,4) * t137;
t162 = Ifges(5,1) * t140 - t195;
t163 = mrSges(5,1) * t137 + mrSges(5,2) * t140;
t139 = sin(qJ(2));
t135 = sin(pkin(4));
t183 = qJD(1) * t135;
t174 = t139 * t183;
t136 = cos(pkin(4));
t142 = cos(qJ(2));
t209 = pkin(1) * t142;
t176 = t136 * t209;
t113 = -pkin(6) * t174 + qJD(1) * t176;
t130 = qJD(1) * t136 + qJD(2);
t85 = -pkin(2) * t130 - t113;
t138 = sin(qJ(3));
t141 = cos(qJ(3));
t97 = t130 * t141 - t138 * t174;
t98 = t130 * t138 + t141 * t174;
t42 = -pkin(3) * t97 - pkin(8) * t98 + t85;
t173 = t142 * t183;
t127 = qJD(3) - t173;
t133 = t136 * t139 * pkin(1);
t186 = t135 * t142;
t184 = pkin(6) * t186 + t133;
t116 = t184 * qJD(1);
t86 = t130 * pkin(7) + t116;
t110 = (-pkin(2) * t142 - pkin(7) * t139 - pkin(1)) * t135;
t92 = qJD(1) * t110;
t51 = t138 * t92 + t141 * t86;
t44 = pkin(8) * t127 + t51;
t10 = t137 * t42 + t140 * t44;
t9 = -t137 * t44 + t140 * t42;
t165 = t10 * t137 + t140 * t9;
t217 = t140 / 0.2e1;
t218 = -t137 / 0.2e1;
t70 = t127 * t137 + t140 * t98;
t214 = Ifges(5,4) * t70;
t69 = t127 * t140 - t137 * t98;
t94 = qJD(4) - t97;
t22 = Ifges(5,2) * t69 + Ifges(5,6) * t94 + t214;
t224 = t94 / 0.2e1;
t228 = t70 / 0.2e1;
t67 = Ifges(5,4) * t69;
t23 = Ifges(5,1) * t70 + Ifges(5,5) * t94 + t67;
t230 = t69 / 0.2e1;
t50 = -t138 * t86 + t141 * t92;
t43 = -pkin(3) * t127 - t50;
t239 = -t165 * mrSges(5,3) + t158 * t224 + t160 * t230 + t162 * t228 + t163 * t43 + t217 * t23 + t218 * t22;
t193 = Ifges(4,5) * t127;
t216 = Ifges(4,1) * t98;
t93 = Ifges(4,4) * t97;
t49 = t193 + t93 + t216;
t148 = t50 * mrSges(4,3) - t49 / 0.2e1 - t193 / 0.2e1 - t85 * mrSges(4,2);
t238 = t148 - t239;
t237 = -Ifges(3,6) * t130 / 0.2e1;
t187 = t135 * t139;
t131 = pkin(6) * t187;
t121 = -t131 + t176;
t182 = qJD(2) * t135;
t169 = qJD(1) * t182;
t168 = t139 * t169;
t153 = (pkin(2) * t139 - pkin(7) * t142) * t135;
t115 = qJD(2) * t153;
t105 = qJD(1) * t115;
t117 = t121 * qJD(2);
t106 = qJD(1) * t117;
t179 = qJD(3) * t141;
t180 = qJD(3) * t138;
t19 = t105 * t138 + t106 * t141 + t179 * t92 - t180 * t86;
t13 = pkin(8) * t168 + t19;
t118 = t184 * qJD(2);
t107 = qJD(1) * t118;
t181 = qJD(2) * t142;
t171 = t141 * t181;
t74 = t130 * t179 + (-t139 * t180 + t171) * t183;
t172 = t138 * t181;
t75 = t130 * t180 + (t139 * t179 + t172) * t183;
t28 = pkin(3) * t75 - pkin(8) * t74 + t107;
t1 = qJD(4) * t9 + t13 * t140 + t137 * t28;
t2 = -qJD(4) * t10 - t13 * t137 + t140 * t28;
t236 = t1 * t140 - t137 * t2;
t20 = -qJD(3) * t51 + t105 * t141 - t106 * t138;
t109 = pkin(7) * t136 + t184;
t197 = t109 * t141 + t110 * t138;
t34 = -qJD(3) * t197 + t115 * t141 - t117 * t138;
t31 = qJD(4) * t69 + t137 * t168 + t140 * t74;
t32 = -qJD(4) * t70 - t137 * t74 + t140 * t168;
t7 = Ifges(5,1) * t31 + Ifges(5,4) * t32 + Ifges(5,5) * t75;
t235 = t7 / 0.2e1;
t234 = -t22 / 0.2e1;
t233 = t31 / 0.2e1;
t232 = t32 / 0.2e1;
t231 = -t69 / 0.2e1;
t229 = -t70 / 0.2e1;
t227 = t75 / 0.2e1;
t226 = -t93 / 0.2e1;
t225 = -t94 / 0.2e1;
t223 = pkin(1) * mrSges(3,1);
t222 = pkin(1) * mrSges(3,2);
t119 = -t136 * t141 + t138 * t187;
t221 = -t119 / 0.2e1;
t120 = t136 * t138 + t141 * t187;
t220 = t120 / 0.2e1;
t219 = t136 / 0.2e1;
t215 = Ifges(4,4) * t98;
t30 = Ifges(5,5) * t31;
t213 = Ifges(5,5) * t70;
t212 = Ifges(4,2) * t97;
t29 = Ifges(5,6) * t32;
t211 = Ifges(5,6) * t69;
t210 = Ifges(5,3) * t94;
t208 = pkin(7) * t141;
t205 = t19 * mrSges(4,2);
t204 = t20 * mrSges(4,1);
t203 = t74 * Ifges(4,1);
t202 = t74 * Ifges(4,4);
t201 = t75 * Ifges(4,4);
t37 = -mrSges(5,1) * t69 + mrSges(5,2) * t70;
t77 = mrSges(4,1) * t127 - mrSges(4,3) * t98;
t198 = t77 - t37;
t196 = Ifges(3,4) * t139;
t191 = Ifges(4,6) * t127;
t114 = qJD(1) * t153;
t66 = t113 * t141 + t114 * t138;
t188 = -mrSges(3,1) * t130 - mrSges(4,1) * t97 + mrSges(4,2) * t98 + mrSges(3,3) * t174;
t185 = t141 * t142;
t126 = -pkin(3) * t141 - pkin(8) * t138 - pkin(2);
t178 = qJD(4) * t126;
t177 = qJD(4) * t141;
t5 = Ifges(5,3) * t75 + t29 + t30;
t175 = Ifges(4,5) * t74 - Ifges(4,6) * t75 + Ifges(4,3) * t168;
t170 = t139 * t182;
t167 = -t2 * mrSges(5,1) + t1 * mrSges(5,2);
t166 = pkin(3) * t138 - pkin(8) * t141;
t164 = mrSges(5,1) * t140 - mrSges(5,2) * t137;
t161 = Ifges(5,1) * t137 + t194;
t159 = Ifges(5,2) * t140 + t195;
t157 = Ifges(5,5) * t137 + Ifges(5,6) * t140;
t108 = t131 + (-pkin(2) - t209) * t136;
t52 = t119 * pkin(3) - t120 * pkin(8) + t108;
t54 = -pkin(8) * t186 + t197;
t15 = -t137 * t54 + t140 * t52;
t16 = t137 * t52 + t140 * t54;
t63 = -t109 * t138 + t110 * t141;
t65 = -t113 * t138 + t114 * t141;
t80 = -t120 * t137 - t140 * t186;
t151 = -t120 * t140 + t137 * t186;
t33 = -t109 * t180 + t110 * t179 + t115 * t138 + t117 * t141;
t128 = Ifges(3,4) * t173;
t149 = -t113 * mrSges(3,3) + t130 * Ifges(3,5) + Ifges(3,1) * t174 / 0.2e1 + t128 / 0.2e1;
t147 = t50 * mrSges(4,1) + t127 * Ifges(4,3) + t98 * Ifges(4,5) + t97 * Ifges(4,6) + t237 - (Ifges(3,2) * t142 + t196) * t183 / 0.2e1 - t116 * mrSges(3,3) - t51 * mrSges(4,2);
t21 = t210 + t211 + t213;
t48 = t191 + t212 + t215;
t146 = -t210 / 0.2e1 - t211 / 0.2e1 - t213 / 0.2e1 + t191 / 0.2e1 + t51 * mrSges(4,3) - t21 / 0.2e1 + t48 / 0.2e1 - t85 * mrSges(4,1) - t9 * mrSges(5,1) + t10 * mrSges(5,2) + t215 / 0.2e1;
t145 = -t212 / 0.2e1 - t146;
t125 = Ifges(3,5) * t142 * t169;
t123 = t166 * qJD(3);
t112 = -mrSges(3,2) * t130 + mrSges(3,3) * t173;
t104 = t126 * t137 + t140 * t208;
t103 = t126 * t140 - t137 * t208;
t89 = (t137 * t139 + t140 * t185) * t183;
t88 = (-t137 * t185 + t139 * t140) * t183;
t79 = -qJD(3) * t119 + t135 * t171;
t78 = qJD(3) * t120 + t135 * t172;
t76 = -mrSges(4,2) * t127 + mrSges(4,3) * t97;
t68 = (t133 + (pkin(6) + t166) * t186) * qJD(1);
t62 = pkin(3) * t98 - pkin(8) * t97;
t61 = -t137 * t178 + t123 * t140 + (t137 * t180 - t140 * t177) * pkin(7);
t60 = t140 * t178 + t123 * t137 + (-t137 * t177 - t140 * t180) * pkin(7);
t58 = -mrSges(4,2) * t168 - mrSges(4,3) * t75;
t57 = mrSges(4,1) * t168 - mrSges(4,3) * t74;
t56 = pkin(8) * t174 + t66;
t55 = -pkin(3) * t174 - t65;
t53 = pkin(3) * t186 - t63;
t46 = mrSges(5,1) * t94 - mrSges(5,3) * t70;
t45 = -mrSges(5,2) * t94 + mrSges(5,3) * t69;
t41 = qJD(4) * t80 + t137 * t170 + t140 * t79;
t40 = qJD(4) * t151 - t137 * t79 + t140 * t170;
t39 = mrSges(4,1) * t75 + mrSges(4,2) * t74;
t38 = pkin(3) * t78 - pkin(8) * t79 + t118;
t36 = Ifges(4,5) * t168 - t201 + t203;
t35 = -t75 * Ifges(4,2) + Ifges(4,6) * t168 + t202;
t27 = t137 * t68 + t140 * t56;
t26 = -t137 * t56 + t140 * t68;
t25 = -pkin(3) * t170 - t34;
t24 = pkin(8) * t170 + t33;
t18 = t137 * t62 + t140 * t50;
t17 = -t137 * t50 + t140 * t62;
t14 = -pkin(3) * t168 - t20;
t12 = -mrSges(5,2) * t75 + mrSges(5,3) * t32;
t11 = mrSges(5,1) * t75 - mrSges(5,3) * t31;
t8 = -mrSges(5,1) * t32 + mrSges(5,2) * t31;
t6 = Ifges(5,4) * t31 + Ifges(5,2) * t32 + Ifges(5,6) * t75;
t4 = -qJD(4) * t16 - t137 * t24 + t140 * t38;
t3 = qJD(4) * t15 + t137 * t38 + t140 * t24;
t47 = [-t175 * t186 / 0.2e1 - t75 * (Ifges(4,4) * t120 - Ifges(4,2) * t119 - Ifges(4,6) * t186) / 0.2e1 + t74 * (Ifges(4,1) * t120 - Ifges(4,4) * t119 - Ifges(4,5) * t186) / 0.2e1 + t98 * (Ifges(4,1) * t79 - Ifges(4,4) * t78) / 0.2e1 + t25 * t37 + t40 * t22 / 0.2e1 + t15 * t11 + t16 * t12 + (-t119 * t19 - t120 * t20 - t50 * t79 - t51 * t78) * mrSges(4,3) + t127 * (Ifges(4,5) * t79 - Ifges(4,6) * t78) / 0.2e1 + t97 * (Ifges(4,4) * t79 - Ifges(4,2) * t78) / 0.2e1 + (-Ifges(5,5) * t151 + Ifges(5,6) * t80 + Ifges(5,3) * t119) * t227 + (-Ifges(5,4) * t151 + Ifges(5,2) * t80 + Ifges(5,6) * t119) * t232 + (-Ifges(5,1) * t151 + Ifges(5,4) * t80 + Ifges(5,5) * t119) * t233 + t2 * (mrSges(5,1) * t119 + mrSges(5,3) * t151) + t14 * (-mrSges(5,1) * t80 - mrSges(5,2) * t151) - t151 * t235 + m(4) * (t107 * t108 + t118 * t85 + t19 * t197 + t20 * t63 + t33 * t51 + t34 * t50) + t197 * t58 + m(3) * (t106 * t184 - t107 * t121 - t113 * t118 + t116 * t117) + ((-t121 * mrSges(3,3) + Ifges(3,5) * t219 + (-0.2e1 * t222 + 0.3e1 / 0.2e1 * Ifges(3,4) * t142) * t135) * t142 + (-t184 * mrSges(3,3) + Ifges(4,5) * t220 + Ifges(4,6) * t221 - Ifges(3,6) * t136 + (-0.2e1 * t223 - 0.3e1 / 0.2e1 * t196 + (-Ifges(4,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t142) * t135) * t139) * t169 + (Ifges(5,5) * t41 + Ifges(5,6) * t40 + Ifges(5,3) * t78) * t224 + (Ifges(5,1) * t41 + Ifges(5,4) * t40 + Ifges(5,5) * t78) * t228 + (Ifges(5,4) * t41 + Ifges(5,2) * t40 + Ifges(5,6) * t78) * t230 + t125 * t219 + t36 * t220 + t35 * t221 + t186 * t205 + m(5) * (t1 * t16 + t10 * t3 + t14 * t53 + t15 * t2 + t25 * t43 + t4 * t9) - t186 * t204 + t106 * (-mrSges(3,2) * t136 + mrSges(3,3) * t186) + (-mrSges(3,1) * t136 + mrSges(4,1) * t119 + mrSges(4,2) * t120 + mrSges(3,3) * t187) * t107 + t188 * t118 + t119 * t5 / 0.2e1 + t1 * (-mrSges(5,2) * t119 + mrSges(5,3) * t80) + t108 * t39 + t117 * t112 + t85 * (mrSges(4,1) * t78 + mrSges(4,2) * t79) + t10 * (-mrSges(5,2) * t78 + mrSges(5,3) * t40) + t78 * t21 / 0.2e1 + t9 * (mrSges(5,1) * t78 - mrSges(5,3) * t41) - t78 * t48 / 0.2e1 + t79 * t49 / 0.2e1 + t80 * t6 / 0.2e1 + t33 * t76 + t34 * t77 + (t149 * t142 + (t237 + t147) * t139) * t182 + t63 * t57 + t53 * t8 + t41 * t23 / 0.2e1 + t43 * (-mrSges(5,1) * t40 + mrSges(5,2) * t41) + t3 * t45 + t4 * t46; t125 + ((qJD(2) * (Ifges(4,5) * t138 + Ifges(4,6) * t141) / 0.2e1 + (t223 + t196 / 0.2e1) * t183 + (t130 / 0.2e1 - qJD(2)) * Ifges(3,6) - t147) * t139 + (-t128 / 0.2e1 + (t222 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t139) * t183 + (t226 - t216 / 0.2e1 + t148) * t141 - t145 * t138 - t149) * t142) * t183 + ((t145 + (-m(4) * t51 - t76) * pkin(7)) * t138 + (t93 / 0.2e1 + t216 / 0.2e1 + (-m(4) * t50 + m(5) * t43 - t198) * pkin(7) - t238) * t141) * qJD(3) - m(5) * (t10 * t27 + t26 * t9 + t43 * t55) - m(4) * (t116 * t85 + t50 * t65 + t51 * t66) + (-t10 * t88 + t9 * t89) * mrSges(5,3) + m(5) * (t1 * t104 + t10 * t60 + t103 * t2 + t61 * t9) + (Ifges(5,5) * t89 + Ifges(5,6) * t88) * t225 + (Ifges(5,1) * t89 + Ifges(5,4) * t88) * t229 + (Ifges(5,4) * t89 + Ifges(5,2) * t88) * t231 + t88 * t234 + (-t20 * mrSges(4,3) + t203 / 0.2e1 - t201 / 0.2e1 + t107 * mrSges(4,2) + t36 / 0.2e1 + t14 * t163 + t162 * t233 + t160 * t232 + t158 * t227 + t6 * t218 + t7 * t217 + (-t1 * t137 - t140 * t2) * mrSges(5,3) + (-m(4) * t20 + m(5) * t14 - t57 + t8) * pkin(7) + (t23 * t218 + t140 * t234 + t157 * t225 + t159 * t231 + t161 * t229 + t43 * t164 + (-t10 * t140 + t137 * t9) * mrSges(5,3)) * qJD(4)) * t138 + (-m(4) * t107 - t39) * pkin(2) + (t35 / 0.2e1 - t5 / 0.2e1 - t30 / 0.2e1 - t29 / 0.2e1 - t107 * mrSges(4,1) + t202 / 0.2e1 + t19 * mrSges(4,3) + (-Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1) * t75 + (m(4) * t19 + t58) * pkin(7) + t167) * t141 - t188 * t116 - t106 * mrSges(3,2) - t107 * mrSges(3,1) - t113 * t112 + t103 * t11 + t104 * t12 - t89 * t23 / 0.2e1 - t43 * (-mrSges(5,1) * t88 + mrSges(5,2) * t89) - t66 * t76 - t65 * t77 + (t60 - t27) * t45 - t55 * t37 + (t61 - t26) * t46; t175 + t198 * t51 - t205 + t204 - pkin(3) * t8 + t236 * mrSges(5,3) - t14 * t164 + (t226 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t98 + t238) * t97 + t239 * qJD(4) + t157 * t227 + t159 * t232 + t161 * t233 + t146 * t98 + t6 * t217 + t137 * t235 - t50 * t76 - t18 * t45 - t17 * t46 + (-pkin(3) * t14 - t10 * t18 - t17 * t9 - t43 * t51) * m(5) + (-t11 * t137 + t12 * t140 + m(5) * t236 + (-m(5) * t165 - t137 * t45 - t140 * t46) * qJD(4)) * pkin(8); -t43 * (mrSges(5,1) * t70 + mrSges(5,2) * t69) + (Ifges(5,1) * t69 - t214) * t229 + t22 * t228 + (Ifges(5,5) * t69 - Ifges(5,6) * t70) * t225 - t9 * t45 + t10 * t46 + (t10 * t70 + t69 * t9) * mrSges(5,3) - t167 + t5 + (-Ifges(5,2) * t70 + t23 + t67) * t231;];
tauc = t47(:);
