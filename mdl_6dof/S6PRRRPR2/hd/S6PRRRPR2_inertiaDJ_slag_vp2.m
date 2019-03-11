% Calculate time derivative of joint inertia matrix for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:05:27
% EndTime: 2019-03-08 23:05:34
% DurationCPUTime: 3.02s
% Computational Cost: add. (5241->398), mult. (12766->596), div. (0->0), fcn. (12318->12), ass. (0->170)
t159 = sin(pkin(12));
t161 = cos(pkin(12));
t238 = -mrSges(6,1) * t161 + mrSges(6,2) * t159 - mrSges(5,1);
t160 = sin(pkin(6));
t166 = sin(qJ(2));
t197 = qJD(2) * t166;
t189 = t160 * t197;
t162 = cos(pkin(6));
t165 = sin(qJ(3));
t169 = cos(qJ(3));
t201 = t160 * t166;
t132 = t162 * t165 + t169 * t201;
t170 = cos(qJ(2));
t196 = qJD(2) * t170;
t188 = t160 * t196;
t119 = -t132 * qJD(3) - t165 * t188;
t131 = t162 * t169 - t165 * t201;
t120 = t131 * qJD(3) + t169 * t188;
t164 = sin(qJ(4));
t168 = cos(qJ(4));
t175 = t168 * t131 - t132 * t164;
t44 = t175 * qJD(4) + t119 * t164 + t120 * t168;
t29 = -t159 * t44 + t161 * t189;
t30 = t159 * t189 + t161 * t44;
t178 = -t159 * t29 + t161 * t30;
t163 = sin(qJ(6));
t167 = cos(qJ(6));
t137 = t159 * t167 + t161 * t163;
t174 = t159 * t163 - t161 * t167;
t112 = mrSges(7,1) * t174 + mrSges(7,2) * t137;
t237 = t112 + t238;
t228 = -pkin(9) - pkin(8);
t147 = t228 * t165;
t148 = t228 * t169;
t236 = t168 * t147 + t148 * t164;
t138 = t164 * t165 - t168 * t169;
t139 = t164 * t169 + t165 * t168;
t154 = -pkin(3) * t169 - pkin(2);
t105 = pkin(4) * t138 - qJ(5) * t139 + t154;
t124 = t147 * t164 - t148 * t168;
t65 = t161 * t105 - t124 * t159;
t66 = t159 * t105 + t161 * t124;
t235 = -t159 * t65 + t161 * t66;
t129 = t174 * qJD(6);
t234 = qJD(3) + qJD(4);
t218 = pkin(3) * qJD(4);
t233 = (-mrSges(5,2) * t168 + t237 * t164) * t218;
t232 = 2 * m(6);
t231 = 2 * m(7);
t190 = qJD(3) * t228;
t142 = t165 * t190;
t184 = t169 * t190;
t74 = t124 * qJD(4) + t142 * t164 - t168 * t184;
t230 = 0.2e1 * t74;
t130 = t137 * qJD(6);
t99 = mrSges(7,1) * t130 - t129 * mrSges(7,2);
t229 = 0.2e1 * t99;
t158 = t161 ^ 2;
t226 = pkin(3) * t168;
t115 = t234 * t138;
t116 = t234 * t139;
t192 = pkin(3) * qJD(3) * t165;
t55 = pkin(4) * t116 + qJ(5) * t115 - qJD(5) * t139 + t192;
t73 = qJD(4) * t236 + t168 * t142 + t164 * t184;
t18 = -t159 * t73 + t161 * t55;
t225 = t18 * mrSges(6,3);
t92 = t131 * t164 + t132 * t168;
t45 = t92 * qJD(4) - t168 * t119 + t120 * t164;
t24 = t175 * t45;
t35 = t174 * t115 - t139 * t130;
t36 = t137 * t115 + t129 * t139;
t11 = -t36 * mrSges(7,1) + t35 * mrSges(7,2);
t207 = t115 * t161;
t208 = t115 * t159;
t61 = -mrSges(6,1) * t208 - mrSges(6,2) * t207;
t224 = t11 + t61;
t19 = t159 * t55 + t161 * t73;
t222 = mrSges(7,3) * t174;
t221 = Ifges(6,4) * t159;
t220 = Ifges(6,4) * t161;
t219 = Ifges(6,2) * t159;
t217 = t236 * t74;
t216 = t130 * mrSges(7,3);
t215 = t159 * t18;
t212 = t161 * t19;
t102 = (mrSges(6,1) * t159 + mrSges(6,2) * t161) * t139;
t93 = t137 * t139;
t94 = t174 * t139;
t56 = mrSges(7,1) * t93 - mrSges(7,2) * t94;
t209 = t102 + t56;
t206 = t236 * t164;
t205 = t139 * t159;
t204 = t139 * t161;
t150 = pkin(3) * t164 + qJ(5);
t202 = t150 * t161;
t200 = t160 * t170;
t199 = -Ifges(7,5) * t129 - Ifges(7,6) * t130;
t198 = t159 ^ 2 + t158;
t195 = 0.2e1 * mrSges(7,3);
t194 = 0.2e1 * t165;
t193 = Ifges(7,5) * t35 + Ifges(7,6) * t36 + Ifges(7,3) * t116;
t191 = t164 * t218;
t151 = -pkin(5) * t161 - pkin(4);
t100 = -Ifges(7,4) * t129 - Ifges(7,2) * t130;
t101 = -Ifges(7,1) * t129 - Ifges(7,4) * t130;
t113 = Ifges(7,4) * t137 - Ifges(7,2) * t174;
t114 = Ifges(7,1) * t137 - Ifges(7,4) * t174;
t187 = -t100 * t174 + t137 * t101 - t130 * t113 - t129 * t114;
t149 = t168 * t218 + qJD(5);
t186 = t198 * t149;
t185 = t198 * qJD(5);
t183 = t160 ^ 2 * t166 * t196;
t181 = -t175 * t74 - t236 * t45;
t180 = -mrSges(4,1) * t169 + mrSges(4,2) * t165;
t179 = Ifges(6,5) * t161 - Ifges(6,6) * t159;
t77 = -t159 * t92 - t161 * t200;
t78 = -t159 * t200 + t161 * t92;
t177 = -t159 * t77 + t161 * t78;
t48 = pkin(5) * t138 - pkin(10) * t204 + t65;
t54 = -pkin(10) * t205 + t66;
t14 = -t163 * t54 + t167 * t48;
t15 = t163 * t48 + t167 * t54;
t39 = -t163 * t78 + t167 * t77;
t40 = t163 * t77 + t167 * t78;
t176 = 0.2e1 * t198 * mrSges(6,3);
t134 = (-pkin(10) - t150) * t159;
t155 = t161 * pkin(10);
t135 = t155 + t202;
t97 = t134 * t167 - t135 * t163;
t98 = t134 * t163 + t135 * t167;
t144 = (-pkin(10) - qJ(5)) * t159;
t146 = qJ(5) * t161 + t155;
t121 = t144 * t167 - t146 * t163;
t122 = t144 * t163 + t146 * t167;
t173 = -t119 * t165 + t120 * t169 + (-t131 * t169 - t132 * t165) * qJD(3);
t5 = t39 * qJD(6) + t163 * t29 + t167 * t30;
t6 = -t40 * qJD(6) - t163 * t30 + t167 * t29;
t172 = -t40 * t216 - t5 * t222 - t175 * t99 + (t129 * t39 - t137 * t6) * mrSges(7,3) - t44 * mrSges(5,2) + t178 * mrSges(6,3) + t237 * t45;
t10 = Ifges(7,1) * t35 + Ifges(7,4) * t36 + Ifges(7,5) * t116;
t12 = pkin(5) * t116 + pkin(10) * t207 + t18;
t16 = pkin(10) * t208 + t19;
t2 = t14 * qJD(6) + t12 * t163 + t16 * t167;
t3 = -t15 * qJD(6) + t12 * t167 - t16 * t163;
t46 = t116 * Ifges(6,6) - (-t219 + t220) * t115;
t47 = t116 * Ifges(6,5) - (Ifges(6,1) * t161 - t221) * t115;
t49 = -Ifges(7,4) * t94 - Ifges(7,2) * t93 + Ifges(7,6) * t138;
t50 = -Ifges(7,1) * t94 - Ifges(7,4) * t93 + Ifges(7,5) * t138;
t51 = -pkin(5) * t208 + t74;
t89 = pkin(5) * t205 - t236;
t9 = Ifges(7,4) * t35 + Ifges(7,2) * t36 + Ifges(7,6) * t116;
t171 = (t129 * t14 - t137 * t3) * mrSges(7,3) + mrSges(6,3) * t212 - t15 * t216 - t2 * t222 + t161 * t46 / 0.2e1 + t159 * t47 / 0.2e1 - t174 * t9 / 0.2e1 + t137 * t10 / 0.2e1 - t129 * t50 / 0.2e1 - t130 * t49 / 0.2e1 + t51 * t112 + t36 * t113 / 0.2e1 + t35 * t114 / 0.2e1 - Ifges(5,5) * t115 - Ifges(5,6) * t116 + t89 * t99 - t93 * t100 / 0.2e1 - t94 * t101 / 0.2e1 - t73 * mrSges(5,2) - (Ifges(6,1) * t159 + t220) * t207 / 0.2e1 + (Ifges(6,2) * t161 + t221) * t208 / 0.2e1 + t138 * t199 / 0.2e1 + t238 * t74 + (Ifges(6,5) * t159 + Ifges(7,5) * t137 + Ifges(6,6) * t161 - Ifges(7,6) * t174) * t116 / 0.2e1;
t153 = -pkin(4) - t226;
t143 = t151 - t226;
t141 = (mrSges(4,1) * t165 + mrSges(4,2) * t169) * qJD(3);
t117 = mrSges(5,1) * t138 + mrSges(5,2) * t139;
t107 = mrSges(6,1) * t138 - mrSges(6,3) * t204;
t106 = -mrSges(6,2) * t138 - mrSges(6,3) * t205;
t88 = t175 * t191;
t85 = -t137 * qJD(5) - t122 * qJD(6);
t84 = -t174 * qJD(5) + t121 * qJD(6);
t76 = mrSges(7,1) * t138 + mrSges(7,3) * t94;
t75 = -mrSges(7,2) * t138 - mrSges(7,3) * t93;
t69 = mrSges(5,1) * t116 - mrSges(5,2) * t115;
t68 = mrSges(6,1) * t116 + mrSges(6,3) * t207;
t67 = -mrSges(6,2) * t116 + mrSges(6,3) * t208;
t64 = -t98 * qJD(6) - t137 * t149;
t63 = t97 * qJD(6) - t174 * t149;
t21 = -mrSges(7,2) * t116 + mrSges(7,3) * t36;
t20 = mrSges(7,1) * t116 - mrSges(7,3) * t35;
t1 = [0.2e1 * m(6) * (t29 * t77 + t30 * t78 - t24) + 0.2e1 * m(7) * (t39 * t6 + t40 * t5 - t24) + 0.2e1 * m(5) * (t92 * t44 - t183 - t24) + 0.2e1 * m(4) * (t131 * t119 + t132 * t120 - t183); t30 * t106 + t29 * t107 + t39 * t20 + t40 * t21 + t5 * t75 + t6 * t76 + t78 * t67 + t77 * t68 - t224 * t175 + t209 * t45 + (t115 * t175 - t116 * t92 - t138 * t44 + t139 * t45) * mrSges(5,3) + t173 * mrSges(4,3) + ((-t141 - t69) * t170 + (-t170 * mrSges(3,2) + (-mrSges(3,1) + t117 + t180) * t166) * qJD(2)) * t160 + m(5) * (t124 * t44 + t73 * t92 + (t154 * t197 - t170 * t192) * t160 + t181) + m(6) * (t18 * t77 + t19 * t78 + t29 * t65 + t30 * t66 + t181) + m(7) * (t14 * t6 + t15 * t5 - t175 * t51 + t2 * t40 + t3 * t39 + t45 * t89) + (-pkin(2) * t189 + t173 * pkin(8)) * m(4); t116 * (-Ifges(7,5) * t94 - Ifges(7,6) * t93) + t102 * t230 + (t14 * t3 + t15 * t2 + t51 * t89) * t231 + 0.2e1 * t154 * t69 - 0.2e1 * pkin(2) * t141 - t94 * t10 + 0.2e1 * t19 * t106 + 0.2e1 * t18 * t107 + 0.2e1 * t89 * t11 - t93 * t9 + 0.2e1 * t2 * t75 + 0.2e1 * t3 * t76 + 0.2e1 * t65 * t68 + 0.2e1 * t66 * t67 + t36 * t49 + t35 * t50 + 0.2e1 * t51 * t56 + 0.2e1 * t14 * t20 + 0.2e1 * t15 * t21 + ((0.2e1 * Ifges(4,4) * t169 + (Ifges(4,1) - Ifges(4,2)) * t194) * t169 + (-Ifges(4,4) * t165 + pkin(3) * t117) * t194) * qJD(3) + (t18 * t65 + t19 * t66 - t217) * t232 + 0.2e1 * m(5) * (t124 * t73 + t154 * t192 - t217) + (-0.2e1 * (-Ifges(5,4) + t179) * t115 - 0.2e1 * t73 * mrSges(5,3) + ((2 * Ifges(5,2)) + (2 * Ifges(6,3)) + Ifges(7,3)) * t116 + t193) * t138 + (mrSges(5,3) * t230 - t159 * t46 + t161 * t47 + (-0.2e1 * Ifges(5,4) + t179) * t116 - (Ifges(6,1) * t158 + (2 * Ifges(5,1)) + (t219 - 0.2e1 * t220) * t159) * t115) * t139 - 0.2e1 * t236 * t61 + 0.2e1 * (t115 * t236 - t116 * t124) * mrSges(5,3); m(7) * (t143 * t45 + t39 * t64 + t40 * t63 + t5 * t98 + t6 * t97 - t88) + t119 * mrSges(4,1) - t120 * mrSges(4,2) + m(5) * (t164 * t44 - t168 * t45 + (-t164 * t175 + t168 * t92) * qJD(4)) * pkin(3) + t172 + m(6) * (t177 * t149 + t178 * t150 + t153 * t45 - t88); (Ifges(4,5) * t169 - Ifges(4,6) * t165 + t180 * pkin(8)) * qJD(3) + (t149 * t106 + t150 * t67) * t161 + (-t149 * t107 - t150 * t68 - t225) * t159 + t153 * t61 + t143 * t11 + m(6) * (t149 * t235 - t150 * t215 + t153 * t74 + t19 * t202) + m(7) * (t14 * t64 + t143 * t51 + t15 * t63 + t2 * t98 + t3 * t97) + t97 * t20 + t98 * t21 + t63 * t75 + t64 * t76 + t171 + (m(5) * (t164 * t73 - t168 * t74) + (t115 * t168 - t116 * t164) * mrSges(5,3) + (-t168 * t138 * mrSges(5,3) + m(5) * (t124 * t168 - t206) - m(6) * t206 + (m(7) * t89 + t139 * mrSges(5,3) + t209) * t164) * qJD(4)) * pkin(3); t143 * t229 + t149 * t176 + 0.2e1 * t233 + (t143 * t191 + t63 * t98 + t64 * t97) * t231 + (t150 * t186 + t153 * t191) * t232 + (t97 * t129 - t98 * t130 - t64 * t137 - t174 * t63) * t195 + t187; m(6) * (-pkin(4) * t45 + t178 * qJ(5) + t177 * qJD(5)) + m(7) * (t121 * t6 + t122 * t5 + t151 * t45 + t39 * t85 + t40 * t84) + t172; t151 * t11 + t121 * t20 + t122 * t21 + t84 * t75 + t85 * t76 - pkin(4) * t61 + t171 + (qJ(5) * t67 + qJD(5) * t106) * t161 + (-qJ(5) * t68 - qJD(5) * t107 - t225) * t159 + m(6) * (-pkin(4) * t74 + t235 * qJD(5) + (t212 - t215) * qJ(5)) + m(7) * (t121 * t3 + t122 * t2 + t14 * t85 + t15 * t84 + t151 * t51); (t143 + t151) * t99 + t233 + m(7) * (t121 * t64 + t122 * t63 + t151 * t191 + t84 * t98 + t85 * t97) + m(6) * (-pkin(4) * t191 + qJ(5) * t186 + t150 * t185) + (t186 + t185) * mrSges(6,3) + ((-t64 - t85) * t137 - (t63 + t84) * t174 - (t122 + t98) * t130 - (-t121 - t97) * t129) * mrSges(7,3) + t187; t151 * t229 + (t121 * t85 + t122 * t84) * t231 + (t198 * qJ(5) * t232 + t176) * qJD(5) + (t121 * t129 - t122 * t130 - t85 * t137 - t174 * t84) * t195 + t187; 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t45; m(6) * t74 + m(7) * t51 + t224; (m(6) + m(7)) * t191 + t99; t99; 0; mrSges(7,1) * t6 - mrSges(7,2) * t5; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t193; mrSges(7,1) * t64 - mrSges(7,2) * t63 + t199; mrSges(7,1) * t85 - mrSges(7,2) * t84 + t199; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
