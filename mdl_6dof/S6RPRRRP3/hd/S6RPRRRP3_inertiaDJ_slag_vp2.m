% Calculate time derivative of joint inertia matrix for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:03:11
% EndTime: 2019-03-09 06:03:18
% DurationCPUTime: 3.44s
% Computational Cost: add. (3924->448), mult. (9134->630), div. (0->0), fcn. (7526->8), ass. (0->181)
t228 = -mrSges(6,1) - mrSges(7,1);
t227 = mrSges(6,3) + mrSges(7,2);
t226 = Ifges(7,4) + Ifges(6,5);
t225 = -Ifges(7,2) - Ifges(6,3);
t224 = Ifges(7,6) - Ifges(6,6);
t142 = sin(qJ(5));
t143 = sin(qJ(4));
t145 = cos(qJ(5));
t146 = cos(qJ(4));
t105 = t142 * t143 - t145 * t146;
t147 = cos(qJ(3));
t189 = qJD(3) * t147;
t220 = qJD(4) + qJD(5);
t106 = t142 * t146 + t143 * t145;
t144 = sin(qJ(3));
t94 = t106 * t144;
t46 = -t105 * t189 - t220 * t94;
t187 = qJD(4) * t144;
t171 = t143 * t187;
t172 = t146 * t189;
t154 = -t171 + t172;
t174 = t143 * t189;
t185 = qJD(5) * t142;
t194 = t144 * t146;
t195 = t143 * t144;
t47 = -t185 * t195 + (t220 * t194 + t174) * t145 + t154 * t142;
t12 = t47 * mrSges(7,1) - t46 * mrSges(7,3);
t13 = t47 * mrSges(6,1) + t46 * mrSges(6,2);
t158 = -t12 - t13;
t186 = qJD(4) * t146;
t153 = t144 * t186 + t174;
t67 = t153 * mrSges(5,1) + t154 * mrSges(5,2);
t223 = t158 - t67;
t177 = -cos(pkin(10)) * pkin(1) - pkin(2);
t222 = 0.2e1 * t177;
t190 = qJD(3) * t144;
t221 = -Ifges(5,5) * t172 - Ifges(5,3) * t190;
t191 = t143 ^ 2 + t146 ^ 2;
t129 = sin(pkin(10)) * pkin(1) + pkin(7);
t103 = -pkin(3) * t147 - t144 * pkin(8) + t177;
t97 = t146 * t103;
t59 = -pkin(9) * t194 + t97 + (-t129 * t143 - pkin(4)) * t147;
t193 = t146 * t147;
t113 = t129 * t193;
t73 = t143 * t103 + t113;
t65 = -pkin(9) * t195 + t73;
t207 = t142 * t59 + t145 * t65;
t210 = pkin(8) * t147;
t211 = pkin(3) * t144;
t115 = (-t210 + t211) * qJD(3);
t175 = t129 * t190;
t192 = t146 * t115 + t143 * t175;
t22 = (pkin(4) * t144 - pkin(9) * t193) * qJD(3) + (-t113 + (pkin(9) * t144 - t103) * t143) * qJD(4) + t192;
t188 = qJD(4) * t143;
t169 = t147 * t188;
t44 = t103 * t186 + t143 * t115 + (-t146 * t190 - t169) * t129;
t25 = -pkin(9) * t153 + t44;
t6 = -qJD(5) * t207 - t142 * t25 + t145 * t22;
t219 = 2 * m(5);
t218 = 2 * m(6);
t217 = 2 * m(7);
t216 = 0.2e1 * pkin(4);
t215 = 0.2e1 * t129;
t214 = -pkin(9) - pkin(8);
t204 = Ifges(5,4) * t143;
t120 = Ifges(5,2) * t146 + t204;
t213 = -t120 / 0.2e1;
t212 = -t143 / 0.2e1;
t27 = -mrSges(7,2) * t47 + mrSges(7,3) * t190;
t30 = -mrSges(6,2) * t190 - mrSges(6,3) * t47;
t209 = t27 + t30;
t28 = mrSges(6,1) * t190 - mrSges(6,3) * t46;
t29 = -mrSges(7,1) * t190 + t46 * mrSges(7,2);
t208 = -t28 + t29;
t85 = -t94 * mrSges(7,2) - mrSges(7,3) * t147;
t86 = mrSges(6,2) * t147 - t94 * mrSges(6,3);
t206 = t85 + t86;
t95 = t105 * t144;
t87 = -mrSges(6,1) * t147 + t95 * mrSges(6,3);
t88 = mrSges(7,1) * t147 - t95 * mrSges(7,2);
t205 = -t87 + t88;
t203 = Ifges(5,4) * t146;
t202 = Ifges(5,5) * t143;
t201 = Ifges(5,6) * t143;
t200 = Ifges(5,6) * t146;
t45 = -t73 * qJD(4) + t192;
t199 = t143 * t45;
t198 = t147 * Ifges(5,6);
t119 = -mrSges(5,1) * t146 + mrSges(5,2) * t143;
t197 = -mrSges(4,1) + t119;
t196 = t129 * t147;
t98 = pkin(4) * t195 + t144 * t129;
t184 = qJD(5) * t145;
t183 = pkin(4) * t188;
t182 = pkin(4) * t185;
t181 = pkin(4) * t184;
t180 = t214 * t143;
t123 = t214 * t146;
t83 = -t142 * t123 - t145 * t180;
t179 = t83 * t185;
t178 = t94 * t185;
t116 = t129 * t189;
t82 = t153 * pkin(4) + t116;
t133 = -pkin(4) * t146 - pkin(3);
t176 = qJD(4) * t214;
t173 = t144 * t189;
t114 = t143 * t176;
t165 = t146 * t176;
t50 = -t83 * qJD(5) + t145 * t114 + t142 * t165;
t84 = -t145 * t123 + t142 * t180;
t51 = t84 * qJD(5) + t142 * t114 - t145 * t165;
t168 = t84 * t50 + t51 * t83;
t167 = (2 * Ifges(4,4)) + t201;
t72 = -t143 * t196 + t97;
t166 = -t72 * qJD(4) + t44;
t164 = mrSges(5,1) * t143 + mrSges(5,2) * t146;
t163 = Ifges(5,1) * t146 - t204;
t121 = Ifges(5,1) * t143 + t203;
t162 = -Ifges(5,2) * t143 + t203;
t18 = -t142 * t65 + t145 * t59;
t159 = t225 * t190 - t224 * t47 - t226 * t46;
t5 = t142 * t22 + t145 * t25 + t59 * t184 - t65 * t185;
t157 = t84 * t46 + t83 * t47 - t50 * t95 + t51 * t94;
t155 = -pkin(5) * t47 + qJ(6) * t46 - qJD(6) * t95;
t75 = t220 * t106;
t68 = Ifges(7,6) * t75;
t69 = Ifges(6,6) * t75;
t74 = t220 * t105;
t70 = Ifges(6,5) * t74;
t71 = Ifges(7,4) * t74;
t151 = t68 - t69 - t70 - t71 + t228 * t51 + (-mrSges(6,2) + mrSges(7,3)) * t50;
t2 = qJ(6) * t190 - qJD(6) * t147 + t5;
t3 = -pkin(5) * t190 - t6;
t150 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) - t159;
t127 = qJD(6) + t181;
t148 = -mrSges(6,2) * t181 + t127 * mrSges(7,3) + t228 * t182;
t140 = qJD(6) * mrSges(7,3);
t137 = Ifges(5,5) * t186;
t132 = -pkin(4) * t145 - pkin(5);
t130 = pkin(4) * t142 + qJ(6);
t112 = -mrSges(5,1) * t147 - mrSges(5,3) * t194;
t111 = mrSges(5,2) * t147 - mrSges(5,3) * t195;
t110 = t163 * qJD(4);
t109 = t162 * qJD(4);
t108 = t164 * qJD(4);
t102 = t164 * t144;
t93 = -Ifges(5,5) * t147 + t163 * t144;
t92 = t162 * t144 - t198;
t90 = -mrSges(5,2) * t190 - mrSges(5,3) * t153;
t89 = mrSges(5,1) * t190 - mrSges(5,3) * t154;
t81 = Ifges(6,1) * t106 - Ifges(6,4) * t105;
t80 = Ifges(7,1) * t106 + Ifges(7,5) * t105;
t79 = Ifges(6,4) * t106 - Ifges(6,2) * t105;
t78 = Ifges(7,5) * t106 + Ifges(7,3) * t105;
t77 = mrSges(6,1) * t105 + mrSges(6,2) * t106;
t76 = mrSges(7,1) * t105 - mrSges(7,3) * t106;
t66 = pkin(5) * t105 - qJ(6) * t106 + t133;
t64 = mrSges(6,1) * t94 - mrSges(6,2) * t95;
t63 = mrSges(7,1) * t94 + mrSges(7,3) * t95;
t62 = -t121 * t187 + (Ifges(5,5) * t144 + t163 * t147) * qJD(3);
t61 = -t120 * t187 + (Ifges(5,6) * t144 + t162 * t147) * qJD(3);
t58 = -Ifges(6,1) * t95 - Ifges(6,4) * t94 - Ifges(6,5) * t147;
t57 = -Ifges(7,1) * t95 - Ifges(7,4) * t147 + Ifges(7,5) * t94;
t56 = -Ifges(6,4) * t95 - Ifges(6,2) * t94 - Ifges(6,6) * t147;
t55 = -Ifges(7,5) * t95 - Ifges(7,6) * t147 + Ifges(7,3) * t94;
t52 = pkin(5) * t94 + qJ(6) * t95 + t98;
t36 = -Ifges(6,1) * t74 - Ifges(6,4) * t75;
t35 = -Ifges(7,1) * t74 + Ifges(7,5) * t75;
t34 = -Ifges(6,4) * t74 - Ifges(6,2) * t75;
t33 = -Ifges(7,5) * t74 + Ifges(7,3) * t75;
t32 = mrSges(6,1) * t75 - mrSges(6,2) * t74;
t31 = mrSges(7,1) * t75 + mrSges(7,3) * t74;
t21 = pkin(5) * t75 + qJ(6) * t74 - qJD(6) * t106 + t183;
t15 = pkin(5) * t147 - t18;
t14 = -qJ(6) * t147 + t207;
t11 = Ifges(6,1) * t46 - Ifges(6,4) * t47 + Ifges(6,5) * t190;
t10 = Ifges(7,1) * t46 + Ifges(7,4) * t190 + Ifges(7,5) * t47;
t9 = Ifges(6,4) * t46 - Ifges(6,2) * t47 + Ifges(6,6) * t190;
t8 = Ifges(7,5) * t46 + Ifges(7,6) * t190 + Ifges(7,3) * t47;
t7 = -t155 + t82;
t1 = [(t67 * t215 - t143 * t61 + t146 * t62 + (-t143 * t93 - t146 * t92 - t147 * (-t200 - t202)) * qJD(4) + (mrSges(4,1) * t222 + (Ifges(5,5) * t146 - t167) * t144 - t226 * t95 + t224 * t94 + (t129 ^ 2 * t219 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3) + t225) * t147) * qJD(3)) * t144 - (t10 + t11) * t95 + (t8 - t9) * t94 + ((mrSges(4,2) * t222 + t102 * t215 - t143 * t92 + t146 * t93 + t167 * t147) * qJD(3) + t159 + t221) * t147 + (t18 * t6 + t207 * t5 + t82 * t98) * t218 + 0.2e1 * t207 * t30 + (t55 - t56) * t47 + (t57 + t58) * t46 + (t14 * t2 + t15 * t3 + t52 * t7) * t217 + (t73 * t44 + t72 * t45) * t219 + 0.2e1 * t14 * t27 + 0.2e1 * t18 * t28 + 0.2e1 * t15 * t29 + 0.2e1 * t52 * t12 + 0.2e1 * t7 * t63 + 0.2e1 * t82 * t64 + 0.2e1 * t2 * t85 + 0.2e1 * t5 * t86 + 0.2e1 * t6 * t87 + 0.2e1 * t3 * t88 + 0.2e1 * t72 * t89 + 0.2e1 * t73 * t90 + 0.2e1 * t98 * t13 + 0.2e1 * t44 * t111 + 0.2e1 * t45 * t112; -t209 * t95 + t208 * t94 + t205 * t47 + t206 * t46 + ((t111 * t146 - t112 * t143) * qJD(3) + t223) * t147 + (-t143 * t89 + t146 * t90 + (-t111 * t143 - t112 * t146) * qJD(4) + (t102 + t63 + t64) * qJD(3)) * t144 + m(7) * (t14 * t46 - t147 * t7 + t15 * t47 + t52 * t190 - t2 * t95 + t3 * t94) + m(6) * (-t147 * t82 - t18 * t47 + t98 * t190 + t207 * t46 - t5 * t95 - t6 * t94) + m(5) * ((-t143 * t72 + t146 * t73 - t196) * t189 + (t175 - t199 + t146 * t44 + (-t143 * t73 - t146 * t72) * qJD(4)) * t144); 0.2e1 * m(5) * (-0.1e1 + t191) * t173 + 0.2e1 * (m(7) + m(6)) * (-t95 * t46 + t94 * t47 - t173); (t55 / 0.2e1 - t56 / 0.2e1) * t75 + m(5) * (-pkin(3) * t116 + (-t188 * t73 - t199) * pkin(8)) - (t57 / 0.2e1 + t58 / 0.2e1) * t74 + (t78 / 0.2e1 - t79 / 0.2e1) * t47 + m(6) * (t133 * t82 - t18 * t51 + t98 * t183 + t207 * t50 + t5 * t84 - t6 * t83) + (t80 / 0.2e1 + t81 / 0.2e1) * t46 + (-t105 * t5 - t106 * t6 + t18 * t74 - t207 * t75) * mrSges(6,3) + (-t105 * t2 + t106 * t3 - t14 * t75 - t15 * t74) * mrSges(7,2) + (t10 / 0.2e1 + t11 / 0.2e1) * t106 + (t8 / 0.2e1 - t9 / 0.2e1) * t105 + (-t137 / 0.2e1 + t70 / 0.2e1 + t69 / 0.2e1 + t71 / 0.2e1 - t68 / 0.2e1 + (t197 * t129 + Ifges(4,5)) * qJD(3)) * t147 + (t33 / 0.2e1 - t34 / 0.2e1) * t94 - (t35 / 0.2e1 + t36 / 0.2e1) * t95 + m(7) * (t14 * t50 + t15 * t51 + t2 * t84 + t21 * t52 + t3 * t83 + t66 * t7) + (t121 * t189 / 0.2e1 + qJD(4) * t93 / 0.2e1 + t61 / 0.2e1 + t166 * mrSges(5,3) + (m(5) * t166 - qJD(4) * t112 + t90) * pkin(8)) * t146 + t205 * t51 + t206 * t50 + t208 * t83 + t209 * t84 + (t129 * t108 + t109 * t212 + t146 * t110 / 0.2e1 + (t121 * t212 + t146 * t213) * qJD(4) + (t129 * mrSges(4,2) - Ifges(4,6) + t202 / 0.2e1 + t200 / 0.2e1 + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t106 + (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t105) * qJD(3)) * t144 + (t189 * t213 - t45 * mrSges(5,3) - pkin(8) * t89 + t62 / 0.2e1 + (-pkin(8) * t111 - t73 * mrSges(5,3) + pkin(4) * t64 - t92 / 0.2e1 + t198 / 0.2e1) * qJD(4)) * t143 + t52 * t31 + t21 * t63 + t66 * t12 - pkin(3) * t67 + t7 * t76 + t82 * t77 + t98 * t32 + t133 * t13; (-t108 - t31 - t32) * t147 + m(7) * (-t147 * t21 + t157) + m(6) * (-pkin(4) * t169 + t157) + ((t191 * mrSges(5,3) - mrSges(4,2)) * t147 + m(5) * (t191 * t210 - t211) + (m(6) * t133 + m(7) * t66 + t197 + t76 + t77) * t144) * qJD(3) + t227 * (-t105 * t46 + t106 * t47 - t74 * t94 + t75 * t95); -0.2e1 * pkin(3) * t108 + t146 * t109 + t143 * t110 + 0.2e1 * t133 * t32 + 0.2e1 * t21 * t76 + 0.2e1 * t66 * t31 + (t78 - t79) * t75 - (t80 + t81) * t74 + (t35 + t36) * t106 + (t33 - t34) * t105 + (t146 * t121 + (t77 * t216 - t120) * t143) * qJD(4) + (t21 * t66 + t168) * t217 + (t133 * t183 + t168) * t218 + 0.2e1 * t227 * (-t105 * t50 + t106 * t51 - t74 * t83 - t75 * t84); m(7) * (t127 * t14 + t130 * t2 + t132 * t3) + ((t28 + m(6) * t6 + (m(6) * t207 + t86) * qJD(5)) * t145 + (t30 + m(6) * t5 + (-m(6) * t18 + m(7) * t15 + t205) * qJD(5)) * t142) * pkin(4) + t150 - t153 * Ifges(5,6) - Ifges(5,5) * t171 - t44 * mrSges(5,2) + t45 * mrSges(5,1) + t127 * t85 + t130 * t27 + t132 * t29 - t221; m(7) * (-t127 * t95 + t130 * t46 + t132 * t47) + (m(7) * t178 / 0.2e1 + m(6) * (t142 * t46 - t145 * t47 - t184 * t95 + t178) / 0.2e1) * t216 + t223; m(7) * (t127 * t84 + t130 * t50 + t132 * t51) + t137 + (pkin(8) * t119 - t201) * qJD(4) + (-t105 * t127 - t130 * t75 - t132 * t74) * mrSges(7,2) + (t106 * mrSges(7,2) * t185 + m(7) * t179 + m(6) * (t142 * t50 - t145 * t51 + t184 * t84 + t179) + (-t142 * t75 + t145 * t74 + (-t105 * t145 + t106 * t142) * qJD(5)) * mrSges(6,3)) * pkin(4) + t151; 0.2e1 * m(7) * (t127 * t130 + t132 * t182) + 0.2e1 * t148; m(7) * (-pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t14) + t150 + qJ(6) * t27 - pkin(5) * t29 + qJD(6) * t85; m(7) * t155 + t158; m(7) * (-pkin(5) * t51 + qJ(6) * t50 + qJD(6) * t84) + (pkin(5) * t74 - qJ(6) * t75 - qJD(6) * t105) * mrSges(7,2) + t151; t140 + m(7) * (-pkin(5) * t182 + qJ(6) * t127 + qJD(6) * t130) + t148; 0.2e1 * m(7) * qJ(6) * qJD(6) + 0.2e1 * t140; m(7) * t3 + t29; m(7) * t47; m(7) * t51 - t74 * mrSges(7,2); m(7) * t182; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
