% Calculate time derivative of joint inertia matrix for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:23
% EndTime: 2019-03-09 04:46:29
% DurationCPUTime: 3.13s
% Computational Cost: add. (2624->418), mult. (6162->595), div. (0->0), fcn. (4899->6), ass. (0->182)
t212 = 2 * qJD(2);
t227 = mrSges(6,3) + mrSges(7,2);
t226 = Ifges(7,4) + Ifges(6,5);
t225 = Ifges(7,6) - Ifges(6,6);
t140 = sin(qJ(4));
t141 = sin(qJ(3));
t187 = qJD(3) * t141;
t168 = t140 * t187;
t142 = cos(qJ(4));
t143 = cos(qJ(3));
t182 = qJD(4) * t143;
t169 = t142 * t182;
t150 = t168 - t169;
t224 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t223 = m(6) + m(7);
t139 = sin(pkin(9));
t193 = cos(pkin(9));
t161 = t193 * t142;
t159 = qJD(4) * t161;
t170 = t140 * t182;
t162 = t193 * t140;
t106 = t139 * t142 + t162;
t85 = t106 * t141;
t44 = qJD(3) * t85 + t139 * t170 - t143 * t159;
t192 = t139 * t140;
t154 = t161 - t192;
t218 = t106 * qJD(4);
t46 = -t143 * t218 - t154 * t187;
t10 = -t44 * mrSges(7,1) - t46 * mrSges(7,3);
t11 = -t44 * mrSges(6,1) + t46 * mrSges(6,2);
t222 = -t10 - t11;
t99 = t154 * qJD(4);
t54 = mrSges(7,1) * t218 - t99 * mrSges(7,3);
t55 = mrSges(6,1) * t218 + t99 * mrSges(6,2);
t221 = -t54 - t55;
t220 = m(7) * qJD(6);
t207 = pkin(8) * t143;
t209 = pkin(3) * t141;
t114 = qJ(2) - t207 + t209;
t144 = -pkin(1) - pkin(7);
t190 = t141 * t144;
t81 = t140 * t114 + t142 * t190;
t188 = t140 ^ 2 + t142 ^ 2;
t184 = qJD(4) * t140;
t167 = t141 * t184;
t186 = qJD(3) * t143;
t219 = -t142 * t186 + t167;
t174 = t142 * t187;
t151 = t170 + t174;
t208 = pkin(4) * t139;
t129 = qJ(6) + t208;
t217 = m(7) * t129 + mrSges(7,3);
t175 = t193 * pkin(4);
t131 = -t175 - pkin(5);
t211 = m(6) * pkin(4);
t216 = m(7) * t131 - t193 * t211 - mrSges(6,1) - mrSges(7,1);
t215 = t139 * t211 - mrSges(6,2) + t217;
t214 = 0.2e1 * m(6);
t213 = 0.2e1 * m(7);
t210 = -t140 / 0.2e1;
t206 = t99 * mrSges(7,2);
t205 = -qJ(5) - pkin(8);
t104 = qJD(2) + (pkin(3) * t143 + pkin(8) * t141) * qJD(3);
t147 = -t81 * qJD(4) + t142 * t104;
t164 = -t140 * t144 + pkin(4);
t181 = qJD(5) * t142;
t13 = qJ(5) * t174 + (qJ(5) * t184 + t164 * qJD(3) - t181) * t143 + t147;
t185 = qJD(3) * t144;
t171 = t143 * t185;
t183 = qJD(4) * t142;
t178 = t140 * t104 + t114 * t183 + t142 * t171;
t16 = -qJ(5) * t169 + (-qJD(5) * t143 + (qJ(5) * qJD(3) - qJD(4) * t144) * t141) * t140 + t178;
t4 = t139 * t13 + t193 * t16;
t28 = mrSges(7,2) * t44 + mrSges(7,3) * t186;
t29 = -mrSges(6,2) * t186 + mrSges(6,3) * t44;
t204 = t28 + t29;
t30 = mrSges(6,1) * t186 - mrSges(6,3) * t46;
t31 = -mrSges(7,1) * t186 + t46 * mrSges(7,2);
t203 = -t30 + t31;
t103 = t142 * t114;
t189 = t142 * t143;
t62 = -qJ(5) * t189 + t164 * t141 + t103;
t191 = t140 * t143;
t70 = -qJ(5) * t191 + t81;
t19 = t139 * t62 + t193 * t70;
t86 = t106 * t143;
t73 = -mrSges(7,2) * t86 + mrSges(7,3) * t141;
t74 = -mrSges(6,2) * t141 - mrSges(6,3) * t86;
t202 = t73 + t74;
t88 = t154 * t143;
t75 = mrSges(6,1) * t141 - mrSges(6,3) * t88;
t76 = -mrSges(7,1) * t141 + mrSges(7,2) * t88;
t201 = t75 - t76;
t200 = Ifges(5,4) * t140;
t199 = Ifges(5,4) * t142;
t198 = Ifges(5,5) * t140;
t197 = Ifges(5,6) * t142;
t27 = -t140 * t171 + t147;
t196 = t140 * t27;
t195 = t141 * Ifges(5,6);
t194 = -mrSges(5,1) * t142 + mrSges(5,2) * t140 - mrSges(4,1);
t180 = pkin(4) * t184;
t179 = pkin(8) * t184;
t177 = mrSges(5,1) * t186;
t176 = mrSges(5,1) * t183;
t132 = -pkin(4) * t142 - pkin(3);
t173 = t141 * t186;
t127 = t141 * t185;
t163 = qJD(4) * t205;
t149 = -qJD(5) * t140 + t142 * t163;
t95 = t140 * t163 + t181;
t47 = t139 * t95 - t193 * t149;
t48 = t139 * t149 + t193 * t95;
t121 = t205 * t142;
t71 = -t121 * t139 - t205 * t162;
t72 = -t193 * t121 + t205 * t192;
t166 = t47 * t71 + t72 * t48;
t165 = -Ifges(5,5) * t142 + (2 * Ifges(4,4));
t26 = -t144 * t167 + t178;
t80 = -t140 * t190 + t103;
t160 = -qJD(4) * t80 + t26;
t107 = pkin(4) * t191 - t143 * t144;
t158 = mrSges(5,1) * t140 + mrSges(5,2) * t142;
t157 = Ifges(5,1) * t142 - t200;
t123 = Ifges(5,1) * t140 + t199;
t156 = -Ifges(5,2) * t140 + t199;
t122 = Ifges(5,2) * t142 + t200;
t3 = t193 * t13 - t139 * t16;
t18 = -t139 * t70 + t193 * t62;
t43 = t219 * t139 - t141 * t159 - t162 * t186;
t45 = qJD(3) * t88 - t141 * t218;
t87 = t154 * t141;
t155 = -t43 * t71 + t72 * t45 + t47 * t85 + t48 * t87;
t77 = -pkin(4) * t150 + t127;
t152 = Ifges(5,6) * t168 + t186 * t224 - t225 * t44 + t226 * t46;
t136 = Ifges(5,5) * t183;
t113 = mrSges(5,1) * t141 - mrSges(5,3) * t189;
t112 = -mrSges(5,2) * t141 - mrSges(5,3) * t191;
t111 = t157 * qJD(4);
t110 = t156 * qJD(4);
t109 = t158 * qJD(4);
t101 = t158 * t143;
t94 = Ifges(7,4) * t99;
t93 = Ifges(6,5) * t99;
t92 = Ifges(6,6) * t218;
t91 = Ifges(7,6) * t218;
t84 = Ifges(5,5) * t141 + t157 * t143;
t83 = t156 * t143 + t195;
t79 = -mrSges(5,2) * t186 + t150 * mrSges(5,3);
t78 = t151 * mrSges(5,3) + t177;
t69 = Ifges(6,1) * t106 + Ifges(6,4) * t154;
t68 = Ifges(7,1) * t106 - Ifges(7,5) * t154;
t67 = Ifges(6,4) * t106 + Ifges(6,2) * t154;
t66 = Ifges(7,5) * t106 - Ifges(7,3) * t154;
t65 = -mrSges(6,1) * t154 + mrSges(6,2) * t106;
t64 = -mrSges(7,1) * t154 - mrSges(7,3) * t106;
t61 = -t150 * mrSges(5,1) - t151 * mrSges(5,2);
t60 = -pkin(5) * t154 - qJ(6) * t106 + t132;
t59 = Ifges(6,1) * t99 - Ifges(6,4) * t218;
t58 = Ifges(7,1) * t99 + Ifges(7,5) * t218;
t57 = Ifges(6,4) * t99 - Ifges(6,2) * t218;
t56 = Ifges(7,5) * t99 + Ifges(7,3) * t218;
t52 = -t123 * t182 + (Ifges(5,5) * t143 - t157 * t141) * qJD(3);
t51 = -t122 * t182 + (Ifges(5,6) * t143 - t156 * t141) * qJD(3);
t50 = mrSges(6,1) * t86 + mrSges(6,2) * t88;
t49 = mrSges(7,1) * t86 - mrSges(7,3) * t88;
t35 = Ifges(6,1) * t88 - Ifges(6,4) * t86 + Ifges(6,5) * t141;
t34 = Ifges(7,1) * t88 + Ifges(7,4) * t141 + Ifges(7,5) * t86;
t33 = Ifges(6,4) * t88 - Ifges(6,2) * t86 + Ifges(6,6) * t141;
t32 = Ifges(7,5) * t88 + Ifges(7,6) * t141 + Ifges(7,3) * t86;
t25 = t86 * pkin(5) - t88 * qJ(6) + t107;
t24 = pkin(5) * t218 - qJ(6) * t99 - qJD(6) * t106 + t180;
t17 = -t141 * pkin(5) - t18;
t15 = qJ(6) * t141 + t19;
t9 = Ifges(6,1) * t46 + Ifges(6,4) * t44 + Ifges(6,5) * t186;
t8 = Ifges(7,1) * t46 + Ifges(7,4) * t186 - Ifges(7,5) * t44;
t7 = Ifges(6,4) * t46 + Ifges(6,2) * t44 + Ifges(6,6) * t186;
t6 = Ifges(7,5) * t46 + Ifges(7,6) * t186 - Ifges(7,3) * t44;
t5 = -pkin(5) * t44 - qJ(6) * t46 - qJD(6) * t88 + t77;
t2 = -pkin(5) * t186 - t3;
t1 = qJ(6) * t186 + qJD(6) * t141 + t4;
t12 = [((mrSges(4,2) * t212) - t140 * t51 + t142 * t52 - 0.2e1 * t144 * t61 + (-t142 * t83 - t140 * t84 + t141 * (-t197 - t198)) * qJD(4) + (0.2e1 * qJ(2) * mrSges(4,1) + (-Ifges(5,6) * t140 - t165) * t143 + t226 * t88 + t225 * t86 + (-0.2e1 * m(5) * t144 ^ 2 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + t224) * t141) * qJD(3)) * t143 + 0.2e1 * t5 * t49 + 0.2e1 * t15 * t28 + 0.2e1 * t19 * t29 + 0.2e1 * t18 * t30 + 0.2e1 * t17 * t31 + 0.2e1 * t25 * t10 + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t212 + 0.2e1 * m(5) * (t81 * t26 + t80 * t27) + (t8 + t9) * t88 + (t6 - t7) * t86 + 0.2e1 * t1 * t73 + 0.2e1 * t4 * t74 + 0.2e1 * t3 * t75 + 0.2e1 * t2 * t76 + 0.2e1 * t77 * t50 + 0.2e1 * t80 * t78 + 0.2e1 * t81 * t79 + 0.2e1 * t107 * t11 + 0.2e1 * t26 * t112 + 0.2e1 * t27 * t113 + (mrSges(4,1) * t212 + (-0.2e1 * qJ(2) * mrSges(4,2) + 0.2e1 * t144 * t101 + t140 * t83 + t165 * t141 - t142 * t84) * qJD(3) + t152) * t141 + (t1 * t15 + t17 * t2 + t25 * t5) * t213 + (t107 * t77 + t18 * t3 + t19 * t4) * t214 + (-t32 + t33) * t44 + (t34 + t35) * t46; t204 * t87 + t203 * t85 + t202 * t45 + t201 * t43 + (-t61 + (t112 * t142 - t113 * t140) * qJD(3) + t222) * t143 + (-t140 * t78 + t142 * t79 + (-t112 * t140 - t113 * t142) * qJD(4) + (t101 + t49 + t50) * qJD(3)) * t141 + m(7) * (t1 * t87 - t143 * t5 + t15 * t45 - t17 * t43 + t25 * t187 + t2 * t85) + m(6) * (t107 * t187 - t143 * t77 + t18 * t43 + t19 * t45 - t3 * t85 + t4 * t87) + m(5) * ((-t140 * t80 + t142 * t81) * t186 + (-0.2e1 * t171 - t196 + t142 * t26 + (-t140 * t81 - t142 * t80) * qJD(4)) * t141); 0.2e1 * m(5) * (-0.1e1 + t188) * t173 + 0.2e1 * t223 * (-t43 * t85 + t87 * t45 - t173); (t32 / 0.2e1 - t33 / 0.2e1) * t218 + (-t106 * t3 + t154 * t4 - t18 * t99 - t19 * t218) * mrSges(6,3) + (t1 * t154 + t106 * t2 - t15 * t218 + t17 * t99) * mrSges(7,2) + (t110 * t210 + t142 * t111 / 0.2e1 - t144 * t109 + (-t142 * t122 / 0.2e1 + t123 * t210) * qJD(4) + (-t144 * mrSges(4,2) + t198 / 0.2e1 + t197 / 0.2e1 - Ifges(4,6) + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t106 - (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t154) * qJD(3)) * t143 - (t6 / 0.2e1 - t7 / 0.2e1) * t154 + (t68 / 0.2e1 + t69 / 0.2e1) * t46 + t25 * t54 + t60 * t10 - pkin(3) * t61 + t5 * t64 + t24 * t49 + m(7) * (t1 * t72 + t15 * t48 + t17 * t47 + t2 * t71 + t24 * t25 + t5 * t60) + (t34 / 0.2e1 + t35 / 0.2e1) * t99 + (-t66 / 0.2e1 + t67 / 0.2e1) * t44 + (t58 / 0.2e1 + t59 / 0.2e1) * t88 + (t56 / 0.2e1 - t57 / 0.2e1) * t86 + (t136 / 0.2e1 + t93 / 0.2e1 - t92 / 0.2e1 + t94 / 0.2e1 + t91 / 0.2e1 + (t194 * t144 - Ifges(4,5)) * qJD(3)) * t141 + (t122 * t187 / 0.2e1 - t27 * mrSges(5,3) - pkin(8) * t78 + t52 / 0.2e1 + (-pkin(8) * t112 - t81 * mrSges(5,3) + pkin(4) * t50 - t83 / 0.2e1 - t195 / 0.2e1) * qJD(4)) * t140 + m(5) * (-pkin(3) * t127 - pkin(8) * t196 - t81 * t179) + (-t123 * t187 / 0.2e1 + qJD(4) * t84 / 0.2e1 + t51 / 0.2e1 + t160 * mrSges(5,3) + (m(5) * t160 - qJD(4) * t113 + t79) * pkin(8)) * t142 + m(6) * (t107 * t180 + t132 * t77 - t18 * t47 + t19 * t48 - t3 * t71 + t4 * t72) + (t8 / 0.2e1 + t9 / 0.2e1) * t106 + t77 * t65 - t201 * t47 + t202 * t48 + t203 * t71 + t204 * t72 + t107 * t55 + t132 * t11; (-t109 + t221) * t143 + m(7) * (-t143 * t24 + t155) + m(6) * (-pkin(4) * t170 + t155) + ((t188 * mrSges(5,3) - mrSges(4,2)) * t143 + m(5) * (t188 * t207 - t209) + (m(6) * t132 + m(7) * t60 + t194 + t64 + t65) * t141) * qJD(3) + t227 * (-t106 * t43 + t154 * t45 - t218 * t87 + t85 * t99); -0.2e1 * pkin(3) * t109 + t142 * t110 + t140 * t111 + 0.2e1 * t132 * t55 + 0.2e1 * t24 * t64 + 0.2e1 * t60 * t54 + (t68 + t69) * t99 + (t66 - t67) * t218 + (t58 + t59) * t106 - (t56 - t57) * t154 + (t142 * t123 + (0.2e1 * pkin(4) * t65 - t122) * t140) * qJD(4) + (t24 * t60 + t166) * t213 + (t132 * t180 + t166) * t214 + 0.2e1 * t227 * (t106 * t47 + t154 * t48 - t218 * t72 + t71 * t99); -t26 * mrSges(5,2) + t27 * mrSges(5,1) + t1 * mrSges(7,3) - t2 * mrSges(7,1) + t3 * mrSges(6,1) - t4 * mrSges(6,2) + t152 - Ifges(5,6) * t169 + t30 * t175 + m(7) * (qJD(6) * t15 + t1 * t129 + t131 * t2) + (t139 * t4 + t193 * t3) * t211 + t29 * t208 + qJD(6) * t73 + t129 * t28 + t131 * t31 - t151 * Ifges(5,5); t219 * mrSges(5,2) - t140 * t177 - t141 * t176 + t215 * t45 - t216 * t43 + t87 * t220; t72 * t220 + mrSges(5,2) * t179 - Ifges(5,6) * t184 - pkin(8) * t176 + t131 * t206 + t136 + t91 - t92 + t93 + t94 + t215 * t48 + t216 * t47 + (-t99 * t175 - t208 * t218) * mrSges(6,3) + (qJD(6) * t154 - t129 * t218) * mrSges(7,2); 0.2e1 * t217 * qJD(6); m(6) * t77 + m(7) * t5 - t222; t223 * t187; m(6) * t180 + m(7) * t24 - t221; 0; 0; m(7) * t2 + t31; -m(7) * t43; m(7) * t47 + t206; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
