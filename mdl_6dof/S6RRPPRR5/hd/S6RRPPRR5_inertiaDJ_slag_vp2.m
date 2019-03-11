% Calculate time derivative of joint inertia matrix for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:20
% EndTime: 2019-03-09 09:07:27
% DurationCPUTime: 3.68s
% Computational Cost: add. (3543->463), mult. (8815->674), div. (0->0), fcn. (7562->8), ass. (0->191)
t240 = Ifges(4,4) + Ifges(3,5);
t144 = sin(qJ(6));
t147 = cos(qJ(6));
t141 = sin(pkin(6));
t146 = sin(qJ(2));
t184 = qJD(2) * t146;
t171 = t141 * t184;
t148 = cos(qJ(5));
t149 = cos(qJ(2));
t183 = qJD(2) * t149;
t170 = t141 * t183;
t142 = cos(pkin(6));
t145 = sin(qJ(5));
t192 = t141 * t146;
t82 = t142 * t148 + t145 * t192;
t61 = -qJD(5) * t82 + t148 * t170;
t191 = t141 * t149;
t83 = -t142 * t145 + t148 * t192;
t63 = t144 * t191 + t147 * t83;
t25 = -qJD(6) * t63 - t144 * t61 - t147 * t171;
t62 = -t144 * t83 + t147 * t191;
t26 = qJD(6) * t62 - t144 * t171 + t147 * t61;
t12 = -mrSges(7,1) * t25 + mrSges(7,2) * t26;
t75 = -t141 * pkin(1) - pkin(2) * t191 - qJ(3) * t192;
t65 = pkin(3) * t191 - t75;
t43 = (pkin(4) * t149 - pkin(9) * t146) * t141 + t65;
t216 = pkin(1) * t142;
t186 = pkin(8) * t191 + t146 * t216;
t74 = t142 * qJ(3) + t186;
t64 = -qJ(4) * t191 + t74;
t55 = -pkin(9) * t142 + t64;
t208 = t145 * t43 + t148 * t55;
t195 = qJD(5) * t208;
t150 = -pkin(2) - pkin(3);
t136 = pkin(4) - t150;
t187 = qJ(3) * t170 + qJD(3) * t192;
t35 = (-pkin(9) * t149 - t136 * t146) * t141 * qJD(2) + t187;
t127 = t149 * t216;
t112 = qJD(2) * t127;
t130 = t142 * qJD(3);
t49 = t112 + t130 + (-qJD(4) * t149 + (-pkin(8) + qJ(4)) * t184) * t141;
t9 = -t145 * t49 + t148 * t35 - t195;
t4 = pkin(5) * t171 - t9;
t239 = -m(7) * t4 - t12;
t177 = qJD(6) * t147;
t180 = qJD(5) * t148;
t153 = t144 * t180 + t145 * t177;
t45 = -mrSges(6,1) * t171 - mrSges(6,3) * t61;
t238 = m(6) * (t9 + t195) + t239 + t45;
t81 = t186 * qJD(2);
t53 = -qJ(4) * t170 - qJD(4) * t192 + t81;
t185 = t144 ^ 2 + t147 ^ 2;
t17 = pkin(10) * t191 + t208;
t84 = -pkin(8) * t192 + t127;
t76 = -t142 * pkin(2) - t84;
t56 = -t142 * pkin(3) - qJ(4) * t192 + t76;
t50 = t142 * pkin(4) - t56;
t27 = pkin(5) * t82 - pkin(10) * t83 + t50;
t10 = -t144 * t17 + t147 * t27;
t11 = t144 * t27 + t147 * t17;
t235 = -t10 * t144 + t11 * t147;
t197 = t148 * mrSges(6,2);
t164 = t145 * mrSges(6,1) + t197;
t90 = t164 * qJD(5);
t18 = -t145 * t55 + t148 * t43;
t16 = -pkin(5) * t191 - t18;
t31 = -mrSges(7,1) * t62 + mrSges(7,2) * t63;
t68 = mrSges(6,1) * t191 - mrSges(6,3) * t83;
t209 = t31 - t68;
t233 = -m(6) * t18 + m(7) * t16 + t209;
t232 = 0.2e1 * m(6);
t231 = 0.2e1 * m(7);
t230 = -0.2e1 * pkin(1);
t229 = 2 * mrSges(4,2);
t228 = -2 * mrSges(3,3);
t227 = -0.2e1 * t53;
t143 = qJ(3) - pkin(9);
t226 = 0.2e1 * t143;
t225 = t62 / 0.2e1;
t224 = t63 / 0.2e1;
t203 = Ifges(7,4) * t144;
t162 = Ifges(7,1) * t147 - t203;
t79 = Ifges(7,5) * t148 - t145 * t162;
t223 = t79 / 0.2e1;
t102 = Ifges(7,2) * t147 + t203;
t222 = t102 / 0.2e1;
t221 = t144 / 0.2e1;
t220 = -t145 / 0.2e1;
t219 = t145 / 0.2e1;
t218 = t147 / 0.2e1;
t217 = t148 / 0.2e1;
t215 = pkin(5) * t145;
t214 = pkin(10) * t148;
t213 = t53 * mrSges(5,1);
t80 = -pkin(8) * t171 + t112;
t212 = t80 * mrSges(3,2);
t211 = mrSges(5,3) - mrSges(4,2);
t181 = qJD(5) * t145;
t60 = -t142 * t181 + (t145 * t183 + t146 * t180) * t141;
t207 = Ifges(6,5) * t61 - Ifges(6,6) * t60;
t206 = mrSges(7,3) * t145;
t205 = Ifges(6,4) * t145;
t204 = Ifges(6,4) * t148;
t202 = Ifges(7,4) * t147;
t201 = Ifges(7,6) * t144;
t198 = t146 * t81;
t100 = -mrSges(7,1) * t147 + mrSges(7,2) * t144;
t196 = -mrSges(6,1) + t100;
t138 = t145 ^ 2;
t194 = t138 * t143;
t140 = t148 ^ 2;
t193 = t140 * t143;
t190 = t144 * t148;
t189 = t147 * t148;
t182 = qJD(3) * t148;
t179 = qJD(6) * t144;
t178 = qJD(6) * t145;
t176 = t148 * t231;
t175 = 0.2e1 * t149;
t174 = Ifges(5,4) + Ifges(4,5) - Ifges(3,4);
t5 = Ifges(7,5) * t26 + Ifges(7,6) * t25 + Ifges(7,3) * t60;
t173 = Ifges(4,6) * t171 + t170 * t240;
t172 = (-mrSges(4,1) - mrSges(3,1)) * t81;
t168 = t144 * t178;
t30 = t60 * mrSges(6,1) + t61 * mrSges(6,2);
t166 = -qJD(6) * t143 * t148 + (t214 - t215) * qJD(5);
t15 = pkin(5) * t60 - pkin(10) * t61 - t53;
t8 = t145 * t35 + t148 * t49 + t43 * t180 - t181 * t55;
t3 = -pkin(10) * t171 + t8;
t1 = qJD(6) * t10 + t144 * t15 + t147 * t3;
t2 = -qJD(6) * t11 - t144 * t3 + t147 * t15;
t165 = t1 * t147 - t144 * t2;
t163 = mrSges(7,1) * t144 + mrSges(7,2) * t147;
t104 = Ifges(7,1) * t144 + t202;
t161 = -Ifges(7,2) * t144 + t202;
t13 = -mrSges(7,2) * t60 + mrSges(7,3) * t25;
t14 = mrSges(7,1) * t60 - mrSges(7,3) * t26;
t160 = t147 * t13 - t144 * t14;
t36 = -mrSges(7,2) * t82 + mrSges(7,3) * t62;
t37 = mrSges(7,1) * t82 - mrSges(7,3) * t63;
t159 = t144 * t37 - t147 * t36;
t21 = Ifges(7,4) * t63 + Ifges(7,2) * t62 + Ifges(7,6) * t82;
t22 = Ifges(7,1) * t63 + Ifges(7,4) * t62 + Ifges(7,5) * t82;
t156 = -t144 * t21 / 0.2e1 + t22 * t218;
t154 = t147 * t180 - t168;
t99 = pkin(5) * t148 + pkin(10) * t145 + t136;
t152 = -qJD(6) * t99 + t143 * t181 - t182;
t46 = Ifges(7,5) * t168 + (-Ifges(7,5) * t189 - Ifges(7,3) * t145) * qJD(5) + t153 * Ifges(7,6);
t151 = -t10 * t177 - t11 * t179 + t165;
t129 = Ifges(7,5) * t177;
t128 = Ifges(6,6) * t181;
t119 = qJD(3) * t194;
t108 = mrSges(5,2) * t170;
t105 = -Ifges(6,1) * t145 - t204;
t103 = -Ifges(6,2) * t148 - t205;
t101 = Ifges(7,5) * t144 + Ifges(7,6) * t147;
t97 = mrSges(7,1) * t148 + t147 * t206;
t96 = -mrSges(7,2) * t148 + t144 * t206;
t95 = (-Ifges(6,1) * t148 + t205) * qJD(5);
t94 = t162 * qJD(6);
t93 = (Ifges(6,2) * t145 - t204) * qJD(5);
t92 = t161 * qJD(6);
t91 = -Ifges(7,6) * t179 + t129;
t89 = t163 * qJD(6);
t88 = mrSges(4,2) * t191 + mrSges(4,3) * t142;
t87 = mrSges(5,2) * t142 - mrSges(5,3) * t191;
t86 = t163 * t145;
t78 = Ifges(7,6) * t148 - t145 * t161;
t77 = Ifges(7,3) * t148 + (-Ifges(7,5) * t147 + t201) * t145;
t73 = t130 + t80;
t72 = mrSges(7,2) * t181 + mrSges(7,3) * t153;
t71 = -mrSges(7,1) * t181 + mrSges(7,3) * t154;
t70 = t143 * t189 + t144 * t99;
t69 = -t143 * t190 + t147 * t99;
t67 = -mrSges(6,2) * t191 - mrSges(6,3) * t82;
t66 = pkin(2) * t171 - t187;
t54 = -mrSges(7,1) * t153 - mrSges(7,2) * t154;
t52 = t150 * t171 + t187;
t48 = t104 * t178 + (-Ifges(7,5) * t145 - t148 * t162) * qJD(5);
t47 = t102 * t178 + (-Ifges(7,6) * t145 - t148 * t161) * qJD(5);
t44 = mrSges(6,2) * t171 - mrSges(6,3) * t60;
t41 = Ifges(6,1) * t83 - Ifges(6,4) * t82 + Ifges(6,5) * t191;
t40 = Ifges(6,4) * t83 - Ifges(6,2) * t82 + Ifges(6,6) * t191;
t33 = t144 * t152 + t147 * t166;
t32 = t144 * t166 - t147 * t152;
t29 = Ifges(6,1) * t61 - Ifges(6,4) * t60 - Ifges(6,5) * t171;
t28 = Ifges(6,4) * t61 - Ifges(6,2) * t60 - Ifges(6,6) * t171;
t20 = Ifges(7,5) * t63 + Ifges(7,6) * t62 + Ifges(7,3) * t82;
t7 = Ifges(7,1) * t26 + Ifges(7,4) * t25 + Ifges(7,5) * t60;
t6 = Ifges(7,4) * t26 + Ifges(7,2) * t25 + Ifges(7,6) * t60;
t19 = [(t20 - t40) * t60 + (mrSges(6,1) * t82 + mrSges(6,2) * t83) * t227 + (t1 * t11 + t10 * t2 + t16 * t4) * t231 + (t149 * t207 + mrSges(5,3) * t146 * t227 + t198 * t229 + 0.2e1 * t52 * (mrSges(5,1) * t149 + mrSges(5,2) * t146) + 0.2e1 * t66 * (-mrSges(4,1) * t149 - mrSges(4,3) * t146) + 0.2e1 * (t149 * t80 + t198) * mrSges(3,3) + ((t76 * t229 + t84 * t228 - 0.2e1 * t75 * mrSges(4,3) - 0.2e1 * t56 * mrSges(5,3) + (mrSges(3,2) * t230 - t174 * t175) * t141 + (-(2 * Ifges(5,5)) + t240) * t142) * t149 + (0.2e1 * t75 * mrSges(4,1) - 0.2e1 * t65 * mrSges(5,1) - 0.2e1 * t74 * mrSges(4,2) + t186 * t228 + 0.2e1 * t64 * mrSges(5,3) - Ifges(6,5) * t83 + Ifges(6,6) * t82 + (-(2 * Ifges(5,6)) + Ifges(4,6) - (2 * Ifges(3,6))) * t142 + (0.2e1 * t174 * t146 + mrSges(3,1) * t230 + (Ifges(3,1) + Ifges(4,1) + Ifges(5,1) - Ifges(3,2) - Ifges(5,2) - Ifges(4,3) - Ifges(6,3)) * t175) * t141) * t146) * qJD(2)) * t141 + 0.2e1 * t65 * t108 + (0.2e1 * t172 + t173 - 0.2e1 * t212 - 0.2e1 * t213) * t142 + (t18 * t9 + t208 * t8 - t50 * t53) * t232 + 0.2e1 * t208 * t44 + 0.2e1 * m(3) * (t186 * t80 - t81 * t84) + 0.2e1 * m(5) * (t49 * t64 + t52 * t65 + t53 * t56) + 0.2e1 * m(4) * (t66 * t75 + t73 * t74 + t76 * t81) + 0.2e1 * t11 * t13 + 0.2e1 * t10 * t14 + 0.2e1 * t16 * t12 + t25 * t21 + t26 * t22 + 0.2e1 * t4 * t31 + 0.2e1 * t1 * t36 + 0.2e1 * t2 * t37 + 0.2e1 * t18 * t45 + 0.2e1 * t50 * t30 + t61 * t41 + t62 * t6 + t63 * t7 + 0.2e1 * t8 * t67 + 0.2e1 * t9 * t68 + t82 * t5 - t82 * t28 + t83 * t29 + 0.2e1 * t49 * t87 + 0.2e1 * t73 * t88; (-t8 * mrSges(6,3) + qJD(3) * t67 + t143 * t44 + m(6) * (qJD(3) * t208 + t143 * t8) + t5 / 0.2e1 - t28 / 0.2e1 - t53 * mrSges(6,1) + (t18 * mrSges(6,3) - t41 / 0.2e1 + t233 * t143 - t156) * qJD(5)) * t148 - t213 - t212 + t173 + t172 + t47 * t225 + t26 * t223 + t48 * t224 + (t87 + t88) * qJD(3) + (t9 * mrSges(6,3) + t6 * t221 - t147 * t7 / 0.2e1 - t29 / 0.2e1 + t53 * mrSges(6,2) + (t21 * t218 + t22 * t221) * qJD(6) + (t208 * mrSges(6,3) - t20 / 0.2e1 + t40 / 0.2e1) * qJD(5) + t233 * qJD(3) + (-qJD(5) * t67 - t238) * t143) * t145 + m(7) * (t1 * t70 + t10 * t33 + t11 * t32 + t2 * t69) + m(4) * (-pkin(2) * t81 + qJ(3) * t73 + qJD(3) * t74) + m(5) * (qJ(3) * t49 + qJD(3) * t64 + t150 * t53) + (t149 * (-Ifges(6,5) * t180 + t128) / 0.2e1 + ((-mrSges(4,2) * pkin(2) - mrSges(5,3) * t150 - Ifges(5,5)) * t149 + (Ifges(6,5) * t219 + Ifges(6,6) * t217 + qJ(3) * t211 - Ifges(3,6) - Ifges(5,6)) * t146) * qJD(2)) * t141 + (t46 / 0.2e1 - t93 / 0.2e1) * t82 + (t77 / 0.2e1 - t103 / 0.2e1) * t60 + (-m(6) * t53 + t30) * t136 + t32 * t36 + t33 * t37 + t49 * mrSges(5,2) + t16 * t54 + t69 * t14 + t70 * t13 + t10 * t71 + t11 * t72 + t73 * mrSges(4,3) + t25 * t78 / 0.2e1 - t4 * t86 - t50 * t90 + t83 * t95 / 0.2e1 + t1 * t96 + t2 * t97 + t61 * t105 / 0.2e1; -0.2e1 * t136 * t90 + 0.2e1 * t32 * t96 + 0.2e1 * t33 * t97 + 0.2e1 * t69 * t71 + 0.2e1 * t70 * t72 + t119 * t232 + (t32 * t70 + t33 * t69 + t119) * t231 + (t193 * t232 + 0.2e1 * mrSges(5,2) + 0.2e1 * mrSges(4,3) + 0.2e1 * (m(5) + m(4)) * qJ(3) + 0.2e1 * (-t138 - t140) * mrSges(6,3)) * qJD(3) + (t46 - t93 + (t144 * t78 - t147 * t79 - t226 * t86 - t105) * qJD(5)) * t148 + (-0.2e1 * qJD(3) * t86 + t54 * t226 + t144 * t47 - t147 * t48 - t95 + (t144 * t79 + t147 * t78) * qJD(6) + (t143 ^ 2 * t176 + t103 - t77) * qJD(5)) * t145; -t144 * t13 - t147 * t14 + t159 * qJD(6) - t211 * t170 + m(7) * (-t235 * qJD(6) - t1 * t144 - t147 * t2) + m(4) * t81 + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t53 - t30; m(7) * (-t144 * t32 - t147 * t33 + (t144 * t69 - t147 * t70) * qJD(6)) - t96 * t177 - t144 * t72 + t97 * t179 - t147 * t71 + t90; 0; m(5) * t52 + t108 - mrSges(5,1) * t171 + ((m(7) * t235 - t159 + t67) * qJD(5) + t238) * t148 + (t44 + (-t144 * t36 - t147 * t37) * qJD(6) + t209 * qJD(5) + m(7) * (qJD(5) * t16 + t151) + m(6) * (-qJD(5) * t18 + t8) + t160) * t145; -t148 * t54 + (m(7) * (t189 * t70 - t190 * t69 - t193 + t194) + t96 * t189 - t97 * t190) * qJD(5) + (-qJD(5) * t86 + m(7) * (-t144 * t33 + t147 * t32 - t177 * t69 - t179 * t70 - t182) - t96 * t179 + t147 * t72 - t97 * t177 - t144 * t71) * t145; 0; (-0.1e1 + t185) * t176 * t181; -Ifges(6,3) * t171 - t8 * mrSges(6,2) + t9 * mrSges(6,1) + t16 * t89 + t82 * t91 / 0.2e1 + t92 * t225 + t94 * t224 + t4 * t100 + t60 * t101 / 0.2e1 + t25 * t222 + t26 * t104 / 0.2e1 + t7 * t221 + t6 * t218 + t156 * qJD(6) + t239 * pkin(5) + ((-t10 * t147 - t11 * t144) * qJD(6) + t165) * mrSges(7,3) + (m(7) * t151 - t177 * t37 - t179 * t36 + t160) * pkin(10) + t207; t128 + t145 * t143 * t89 + t91 * t217 + (t145 * t196 - t197) * qJD(3) + (t101 * t220 - Ifges(6,5) * t148 + (t145 * mrSges(6,2) + t148 * t196) * t143) * qJD(5) + (m(7) * (-qJD(3) * t145 - t143 * t180) - t54) * pkin(5) + (-t104 * t180 / 0.2e1 + t32 * mrSges(7,3) + t94 * t220 + t47 / 0.2e1 + (-t69 * mrSges(7,3) + t102 * t219 + t223) * qJD(6) + (m(7) * (-qJD(6) * t69 + t32) + t72 - qJD(6) * t97) * pkin(10)) * t147 + (t180 * t222 - t33 * mrSges(7,3) + t92 * t219 + t48 / 0.2e1 + (t104 * t219 - t70 * mrSges(7,3) - t78 / 0.2e1) * qJD(6) + (m(7) * (-qJD(6) * t70 - t33) - qJD(6) * t96 - t71) * pkin(10)) * t144; 0; -t148 * t89 + (t145 * t100 + m(7) * (t185 * t214 - t215) + t185 * t148 * mrSges(7,3) - t164) * qJD(5); -0.2e1 * pkin(5) * t89 + t144 * t94 + t147 * t92 + (-t102 * t144 + t104 * t147) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t33 - mrSges(7,2) * t32 + t46; t89; t54; t129 + (pkin(10) * t100 - t201) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
