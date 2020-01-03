% Calculate time derivative of joint inertia matrix for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:45
% EndTime: 2019-12-31 22:27:55
% DurationCPUTime: 3.51s
% Computational Cost: add. (5478->415), mult. (13107->635), div. (0->0), fcn. (11449->8), ass. (0->182)
t160 = sin(qJ(3));
t161 = sin(qJ(2));
t164 = cos(qJ(3));
t193 = qJD(3) * t164;
t165 = cos(qJ(2));
t196 = qJD(2) * t165;
t168 = t160 * t196 + t161 * t193;
t159 = sin(qJ(4));
t163 = cos(qJ(4));
t172 = t159 * t160 - t163 * t164;
t112 = t172 * t161;
t218 = -pkin(8) - pkin(7);
t142 = t218 * t160;
t143 = t218 * t164;
t100 = t159 * t142 - t163 * t143;
t182 = t164 * t196;
t197 = qJD(2) * t161;
t226 = -Ifges(4,5) * t182 - Ifges(4,3) * t197;
t139 = -pkin(2) * t165 - t161 * pkin(7) - pkin(1);
t199 = t164 * t165;
t149 = pkin(6) * t199;
t107 = t160 * t139 + t149;
t225 = qJD(3) + qJD(4);
t126 = t159 * t164 + t160 * t163;
t94 = t225 * t126;
t60 = -t161 * t94 - t172 * t196;
t61 = t112 * t225 - t126 * t196;
t224 = -Ifges(5,5) * t60 - Ifges(5,6) * t61 - Ifges(5,3) * t197;
t223 = 2 * m(4);
t222 = 2 * m(5);
t221 = 2 * m(6);
t220 = -0.2e1 * pkin(1);
t219 = 0.2e1 * pkin(6);
t217 = -t160 / 0.2e1;
t216 = t164 / 0.2e1;
t214 = pkin(6) * t160;
t150 = pkin(3) * t163 + pkin(4);
t162 = cos(qJ(5));
t189 = qJD(5) * t162;
t158 = sin(qJ(5));
t190 = qJD(5) * t158;
t205 = t158 * t159;
t81 = t150 * t189 + (-t159 * t190 + (t162 * t163 - t205) * qJD(4)) * pkin(3);
t213 = t81 * mrSges(6,2);
t86 = -t126 * t158 - t162 * t172;
t93 = t225 * t172;
t31 = qJD(5) * t86 - t158 * t94 - t162 * t93;
t87 = t126 * t162 - t158 * t172;
t32 = -qJD(5) * t87 + t158 * t93 - t162 * t94;
t211 = Ifges(6,5) * t31 + Ifges(6,6) * t32;
t210 = -Ifges(5,5) * t93 - Ifges(5,6) * t94;
t124 = t164 * t139;
t202 = t161 * t164;
t88 = -pkin(8) * t202 + t124 + (-pkin(3) - t214) * t165;
t203 = t160 * t161;
t98 = -pkin(8) * t203 + t107;
t52 = t159 * t88 + t163 * t98;
t209 = Ifges(4,4) * t160;
t208 = Ifges(4,4) * t164;
t207 = Ifges(4,6) * t160;
t206 = Ifges(4,6) * t165;
t204 = t159 * t162;
t176 = Ifges(4,1) * t164 - t209;
t110 = -Ifges(4,5) * t165 + t161 * t176;
t201 = t164 * t110;
t141 = Ifges(4,1) * t160 + t208;
t200 = t164 * t141;
t137 = (pkin(2) * t161 - pkin(7) * t165) * qJD(2);
t198 = t164 * t137 + t197 * t214;
t138 = pkin(3) * t203 + t161 * pkin(6);
t195 = qJD(3) * t160;
t194 = qJD(3) * t161;
t192 = qJD(4) * t159;
t191 = qJD(4) * t163;
t111 = t126 * t161;
t71 = -t111 * t162 + t112 * t158;
t24 = qJD(5) * t71 + t158 * t61 + t162 * t60;
t72 = -t111 * t158 - t112 * t162;
t25 = -qJD(5) * t72 - t158 * t60 + t162 * t61;
t188 = -Ifges(6,5) * t24 - Ifges(6,6) * t25 - Ifges(6,3) * t197;
t187 = pkin(3) * t195;
t105 = t168 * pkin(3) + pkin(6) * t196;
t151 = -pkin(3) * t164 - pkin(2);
t186 = qJD(3) * t218;
t185 = t160 * t194;
t82 = -t150 * t190 + (-t159 * t189 + (-t158 * t163 - t204) * qJD(4)) * pkin(3);
t79 = t82 * mrSges(6,1);
t181 = t79 - t213;
t180 = (2 * Ifges(3,4)) + t207;
t51 = -t159 * t98 + t163 * t88;
t99 = t163 * t142 + t143 * t159;
t135 = t160 * t186;
t136 = t164 * t186;
t63 = t163 * t135 + t159 * t136 + t142 * t191 + t143 * t192;
t35 = -pkin(9) * t94 + t63;
t64 = -qJD(4) * t100 - t135 * t159 + t163 * t136;
t36 = pkin(9) * t93 + t64;
t76 = -pkin(9) * t126 + t99;
t77 = -pkin(9) * t172 + t100;
t39 = -t158 * t77 + t162 * t76;
t8 = qJD(5) * t39 + t158 * t36 + t162 * t35;
t40 = t158 * t76 + t162 * t77;
t9 = -qJD(5) * t40 - t158 * t35 + t162 * t36;
t179 = t9 * mrSges(6,1) - t8 * mrSges(6,2) + t211;
t178 = -mrSges(4,1) * t164 + mrSges(4,2) * t160;
t177 = mrSges(4,1) * t160 + mrSges(4,2) * t164;
t175 = -Ifges(4,2) * t160 + t208;
t140 = Ifges(4,2) * t164 + t209;
t174 = Ifges(4,5) * t160 + Ifges(4,6) * t164;
t37 = -pkin(4) * t165 + t112 * pkin(9) + t51;
t41 = -pkin(9) * t111 + t52;
t18 = -t158 * t41 + t162 * t37;
t19 = t158 * t37 + t162 * t41;
t67 = t160 * t137 + t139 * t193 + (-t164 * t197 - t165 * t195) * pkin(6);
t68 = -t107 * qJD(3) + t198;
t173 = -t160 * t68 + t164 * t67;
t48 = (pkin(3) * t161 - pkin(8) * t199) * qJD(2) + (-t149 + (pkin(8) * t161 - t139) * t160) * qJD(3) + t198;
t57 = -pkin(8) * t168 + t67;
t14 = -qJD(4) * t52 - t159 * t57 + t163 * t48;
t10 = pkin(4) * t197 - pkin(9) * t60 + t14;
t13 = t159 * t48 + t163 * t57 + t88 * t191 - t192 * t98;
t11 = pkin(9) * t61 + t13;
t2 = qJD(5) * t18 + t10 * t158 + t11 * t162;
t3 = -qJD(5) * t19 + t10 * t162 - t11 * t158;
t171 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t188;
t170 = (-mrSges(5,1) * t159 - mrSges(5,2) * t163) * qJD(4) * pkin(3);
t169 = t182 - t185;
t167 = t64 * mrSges(5,1) - t63 * mrSges(5,2) + t179 + t210;
t166 = t14 * mrSges(5,1) - t13 * mrSges(5,2) + t171 - t224;
t155 = Ifges(4,5) * t193;
t134 = -mrSges(4,1) * t165 - mrSges(4,3) * t202;
t133 = mrSges(4,2) * t165 - mrSges(4,3) * t203;
t132 = t176 * qJD(3);
t131 = t175 * qJD(3);
t130 = t177 * qJD(3);
t121 = (-mrSges(6,1) * t158 - mrSges(6,2) * t162) * qJD(5) * pkin(4);
t114 = pkin(3) * t204 + t150 * t158;
t113 = -pkin(3) * t205 + t150 * t162;
t109 = t161 * t175 - t206;
t108 = pkin(4) * t172 + t151;
t106 = -t165 * t214 + t124;
t104 = -mrSges(4,2) * t197 - mrSges(4,3) * t168;
t103 = mrSges(4,1) * t197 - mrSges(4,3) * t169;
t102 = -mrSges(5,1) * t165 + t112 * mrSges(5,3);
t101 = mrSges(5,2) * t165 - t111 * mrSges(5,3);
t97 = Ifges(5,1) * t126 - Ifges(5,4) * t172;
t96 = Ifges(5,4) * t126 - Ifges(5,2) * t172;
t95 = mrSges(5,1) * t172 + mrSges(5,2) * t126;
t91 = pkin(4) * t111 + t138;
t85 = mrSges(4,1) * t168 + mrSges(4,2) * t169;
t78 = pkin(4) * t94 + t187;
t75 = mrSges(5,1) * t111 - mrSges(5,2) * t112;
t74 = -t141 * t194 + (Ifges(4,5) * t161 + t165 * t176) * qJD(2);
t73 = -t140 * t194 + (Ifges(4,6) * t161 + t165 * t175) * qJD(2);
t70 = -Ifges(5,1) * t112 - Ifges(5,4) * t111 - Ifges(5,5) * t165;
t69 = -Ifges(5,4) * t112 - Ifges(5,2) * t111 - Ifges(5,6) * t165;
t66 = -mrSges(6,1) * t165 - t72 * mrSges(6,3);
t65 = mrSges(6,2) * t165 + t71 * mrSges(6,3);
t56 = -Ifges(5,1) * t93 - Ifges(5,4) * t94;
t55 = -Ifges(5,4) * t93 - Ifges(5,2) * t94;
t54 = mrSges(5,1) * t94 - mrSges(5,2) * t93;
t50 = -mrSges(5,2) * t197 + mrSges(5,3) * t61;
t49 = mrSges(5,1) * t197 - mrSges(5,3) * t60;
t47 = Ifges(6,1) * t87 + Ifges(6,4) * t86;
t46 = Ifges(6,4) * t87 + Ifges(6,2) * t86;
t45 = -mrSges(6,1) * t86 + mrSges(6,2) * t87;
t42 = -pkin(4) * t61 + t105;
t38 = -mrSges(6,1) * t71 + mrSges(6,2) * t72;
t34 = Ifges(6,1) * t72 + Ifges(6,4) * t71 - Ifges(6,5) * t165;
t33 = Ifges(6,4) * t72 + Ifges(6,2) * t71 - Ifges(6,6) * t165;
t28 = -mrSges(5,1) * t61 + mrSges(5,2) * t60;
t27 = Ifges(5,1) * t60 + Ifges(5,4) * t61 + Ifges(5,5) * t197;
t26 = Ifges(5,4) * t60 + Ifges(5,2) * t61 + Ifges(5,6) * t197;
t21 = -mrSges(6,2) * t197 + mrSges(6,3) * t25;
t20 = mrSges(6,1) * t197 - mrSges(6,3) * t24;
t17 = Ifges(6,1) * t31 + Ifges(6,4) * t32;
t16 = Ifges(6,4) * t31 + Ifges(6,2) * t32;
t15 = -mrSges(6,1) * t32 + mrSges(6,2) * t31;
t6 = -mrSges(6,1) * t25 + mrSges(6,2) * t24;
t5 = Ifges(6,1) * t24 + Ifges(6,4) * t25 + Ifges(6,5) * t197;
t4 = Ifges(6,4) * t24 + Ifges(6,2) * t25 + Ifges(6,6) * t197;
t1 = [(t85 * t219 - t160 * t73 + t164 * t74 + (-t164 * t109 - t160 * t110 + t165 * t174) * qJD(3) + (mrSges(3,1) * t220 + Ifges(6,5) * t72 + Ifges(6,6) * t71 - Ifges(5,5) * t112 - Ifges(5,6) * t111 + (Ifges(4,5) * t164 - t180) * t161 + (pkin(6) ^ 2 * t223 + t177 * t219 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(5,3) - Ifges(6,3)) * t165) * qJD(2)) * t161 + 0.2e1 * t138 * t28 + 0.2e1 * t67 * t133 + 0.2e1 * t68 * t134 - t111 * t26 - t112 * t27 + 0.2e1 * t13 * t101 + 0.2e1 * t14 * t102 + 0.2e1 * t105 * t75 + 0.2e1 * t106 * t103 + 0.2e1 * t107 * t104 + 0.2e1 * t91 * t6 + 0.2e1 * t2 * t65 + 0.2e1 * t3 * t66 + t61 * t69 + t60 * t70 + t71 * t4 + t72 * t5 + 0.2e1 * t51 * t49 + 0.2e1 * t52 * t50 + 0.2e1 * t42 * t38 + t25 * t33 + t24 * t34 + 0.2e1 * t18 * t20 + 0.2e1 * t19 * t21 + (t18 * t3 + t19 * t2 + t42 * t91) * t221 + (t105 * t138 + t13 * t52 + t14 * t51) * t222 + (t106 * t68 + t107 * t67) * t223 + ((mrSges(3,2) * t220 - t160 * t109 + t165 * t180 + t201) * qJD(2) + t188 + t224 + t226) * t165; (-t18 * t31 + t19 * t32 + t2 * t86 - t3 * t87) * mrSges(6,3) + m(5) * (t100 * t13 + t105 * t151 + t138 * t187 + t14 * t99 + t51 * t64 + t52 * t63) + ((-t106 * t164 - t107 * t160) * qJD(3) + t173) * mrSges(4,3) + (t140 * t217 + t200 / 0.2e1 + Ifges(3,5) + (-m(4) * pkin(2) - mrSges(3,1) + t178) * pkin(6)) * t196 + (t201 / 0.2e1 + (pkin(3) * t75 + t206 / 0.2e1 - t109 / 0.2e1) * t160) * qJD(3) + (-t134 * t193 - t133 * t195 + t164 * t104 + m(4) * (-t106 * t193 - t107 * t195 + t173) - t160 * t103) * pkin(7) + t160 * t74 / 0.2e1 + t151 * t28 + t138 * t54 + (-t126 * t14 - t13 * t172 + t51 * t93 - t52 * t94) * mrSges(5,3) + t126 * t27 / 0.2e1 - t111 * t55 / 0.2e1 - (t211 + t210 + t155) * t165 / 0.2e1 - t112 * t56 / 0.2e1 + (t132 * t216 + t131 * t217 - Ifges(3,6) * qJD(2) + (-t164 * t140 / 0.2e1 + t141 * t217) * qJD(3) + (qJD(2) * mrSges(3,2) + t130) * pkin(6) + (Ifges(5,5) * t126 + Ifges(6,5) * t87 - Ifges(5,6) * t172 + Ifges(6,6) * t86 + t174) * qJD(2) / 0.2e1) * t161 - t172 * t26 / 0.2e1 - t94 * t69 / 0.2e1 + t61 * t96 / 0.2e1 + t60 * t97 / 0.2e1 + t99 * t49 + t100 * t50 + t63 * t101 + t64 * t102 + t105 * t95 + t108 * t6 + t86 * t4 / 0.2e1 + t87 * t5 / 0.2e1 + t91 * t15 - t93 * t70 / 0.2e1 + t78 * t38 - pkin(2) * t85 + t8 * t65 + t9 * t66 + t71 * t16 / 0.2e1 + t72 * t17 / 0.2e1 + t42 * t45 + t25 * t46 / 0.2e1 + t24 * t47 / 0.2e1 + t39 * t20 + t40 * t21 + t32 * t33 / 0.2e1 + t31 * t34 / 0.2e1 + m(6) * (t108 * t42 + t18 * t9 + t19 * t8 + t2 * t40 + t3 * t39 + t78 * t91) + t73 * t216; -0.2e1 * pkin(2) * t130 + 0.2e1 * t108 * t15 - t172 * t55 + t126 * t56 + t164 * t131 + t160 * t132 + 0.2e1 * t151 * t54 + t86 * t16 + t87 * t17 + t31 * t47 + t32 * t46 + 0.2e1 * t78 * t45 - t93 * t97 - t94 * t96 + (t200 + (0.2e1 * pkin(3) * t95 - t140) * t160) * qJD(3) + (t100 * t63 + t151 * t187 + t64 * t99) * t222 + (t108 * t78 + t39 * t9 + t40 * t8) * t221 + 0.2e1 * (-t31 * t39 + t32 * t40 + t8 * t86 - t87 * t9) * mrSges(6,3) + 0.2e1 * (-t100 * t94 - t126 * t64 - t172 * t63 + t93 * t99) * mrSges(5,3); m(6) * (t113 * t3 + t114 * t2 + t18 * t82 + t19 * t81) - t168 * Ifges(4,6) + (t101 * t191 - t102 * t192 + t163 * t49 + m(5) * (t13 * t159 + t14 * t163 + t191 * t52 - t192 * t51) + t159 * t50) * pkin(3) - Ifges(4,5) * t185 + t113 * t20 + t114 * t21 + t81 * t65 + t82 * t66 - t67 * mrSges(4,2) + t68 * mrSges(4,1) + t166 - t226; m(6) * (t113 * t9 + t114 * t8 + t39 * t82 + t40 * t81) + t155 + (pkin(7) * t178 - t207) * qJD(3) + (-t113 * t31 + t114 * t32 + t81 * t86 - t82 * t87) * mrSges(6,3) + (m(5) * (t159 * t63 + t163 * t64 + (t100 * t163 - t159 * t99) * qJD(4)) + (-t159 * t94 + t163 * t93 + (t126 * t159 - t163 * t172) * qJD(4)) * mrSges(5,3)) * pkin(3) + t167; (t113 * t82 + t114 * t81) * t221 - 0.2e1 * t213 + 0.2e1 * t79 + 0.2e1 * t170; (m(6) * (t158 * t2 + t162 * t3 - t18 * t190 + t189 * t19) + t65 * t189 + t158 * t21 - t66 * t190 + t162 * t20) * pkin(4) + t166; (m(6) * (t158 * t8 + t162 * t9 + (-t158 * t39 + t162 * t40) * qJD(5)) + (t158 * t32 - t162 * t31 + (t158 * t87 + t162 * t86) * qJD(5)) * mrSges(6,3)) * pkin(4) + t167; t170 + (m(6) * (-t113 * t190 + t114 * t189 + t158 * t81 + t162 * t82) - mrSges(6,2) * t189 - mrSges(6,1) * t190) * pkin(4) + t181; 0.2e1 * t121; t171; t179; t181; t121; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
