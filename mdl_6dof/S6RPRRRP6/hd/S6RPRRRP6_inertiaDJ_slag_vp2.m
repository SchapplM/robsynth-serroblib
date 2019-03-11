% Calculate time derivative of joint inertia matrix for
% S6RPRRRP6
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:17
% EndTime: 2019-03-09 06:14:28
% DurationCPUTime: 4.84s
% Computational Cost: add. (6368->411), mult. (14415->580), div. (0->0), fcn. (14131->8), ass. (0->171)
t156 = sin(pkin(10));
t157 = cos(pkin(10));
t160 = sin(qJ(3));
t163 = cos(qJ(3));
t133 = t156 * t163 + t160 * t157;
t126 = t133 * qJD(3);
t234 = Ifges(6,6) + Ifges(7,6);
t240 = t234 * t126;
t235 = Ifges(6,5) + Ifges(7,5);
t239 = t235 * t126;
t238 = Ifges(6,4) + Ifges(7,4);
t158 = sin(qJ(5));
t159 = sin(qJ(4));
t161 = cos(qJ(5));
t162 = cos(qJ(4));
t135 = t158 * t162 + t159 * t161;
t225 = qJD(4) + qJD(5);
t104 = t225 * t135;
t132 = t156 * t160 - t163 * t157;
t125 = t132 * qJD(3);
t170 = t158 * t159 - t161 * t162;
t34 = -t104 * t133 + t125 * t170;
t86 = t170 * t133;
t35 = t135 * t125 + t225 * t86;
t237 = (Ifges(6,2) + Ifges(7,2)) * t35 + t238 * t34 + t240;
t236 = t238 * t35 + (Ifges(6,1) + Ifges(7,1)) * t34 + t239;
t233 = Ifges(6,3) + Ifges(7,3);
t190 = qJD(4) * t162;
t168 = -t125 * t159 + t133 * t190;
t103 = t225 * t170;
t232 = -t235 * t103 - t234 * t104;
t230 = pkin(4) * t161;
t229 = -mrSges(6,1) - mrSges(7,1);
t228 = pkin(4) * qJD(5);
t211 = pkin(7) + qJ(2);
t144 = t211 * t156;
t145 = t211 * t157;
t112 = -t160 * t144 + t145 * t163;
t102 = t162 * t112;
t185 = -pkin(2) * t157 - pkin(1);
t91 = pkin(3) * t132 - pkin(8) * t133 + t185;
t58 = t159 * t91 + t102;
t195 = t125 * t162;
t227 = -Ifges(5,5) * t195 + Ifges(5,3) * t126;
t226 = -t163 * t144 - t145 * t160;
t219 = -pkin(9) - pkin(8);
t148 = t219 * t159;
t149 = t219 * t162;
t114 = t158 * t148 - t161 * t149;
t174 = mrSges(5,1) * t159 + mrSges(5,2) * t162;
t139 = t174 * qJD(4);
t224 = 2 * m(6);
t223 = 2 * m(7);
t222 = -2 * mrSges(4,3);
t220 = -0.2e1 * t226;
t216 = -t133 / 0.2e1;
t206 = Ifges(5,4) * t159;
t146 = Ifges(5,2) * t162 + t206;
t213 = -t146 / 0.2e1;
t212 = pkin(5) * t104;
t27 = -mrSges(7,2) * t126 + mrSges(7,3) * t35;
t28 = -mrSges(6,2) * t126 + mrSges(6,3) * t35;
t210 = t27 + t28;
t193 = t133 * t162;
t57 = -t112 * t159 + t162 * t91;
t42 = pkin(4) * t132 - pkin(9) * t193 + t57;
t194 = t133 * t159;
t50 = -pkin(9) * t194 + t58;
t22 = t158 * t42 + t161 * t50;
t85 = t135 * t133;
t70 = -mrSges(7,2) * t132 - mrSges(7,3) * t85;
t71 = -mrSges(6,2) * t132 - mrSges(6,3) * t85;
t209 = t70 + t71;
t72 = mrSges(7,1) * t132 + mrSges(7,3) * t86;
t73 = mrSges(6,1) * t132 + mrSges(6,3) * t86;
t208 = t72 + t73;
t207 = -t104 * mrSges(7,1) + t103 * mrSges(7,2);
t205 = Ifges(5,4) * t162;
t204 = Ifges(5,6) * t159;
t78 = qJD(2) * t133 + qJD(3) * t112;
t203 = t226 * t78;
t202 = t126 * Ifges(5,5);
t201 = t126 * Ifges(5,6);
t200 = t132 * Ifges(5,6);
t188 = qJD(5) * t161;
t197 = (-t103 * t158 + t135 * t188) * pkin(4);
t191 = qJD(4) * t159;
t189 = qJD(5) * t158;
t187 = pkin(4) * t191;
t77 = -t132 * qJD(2) + qJD(3) * t226;
t90 = pkin(3) * t126 + pkin(8) * t125;
t178 = -t159 * t77 + t162 * t90;
t17 = pkin(9) * t195 + pkin(4) * t126 + (-t102 + (pkin(9) * t133 - t91) * t159) * qJD(4) + t178;
t23 = -t112 * t191 + t159 * t90 + t162 * t77 + t91 * t190;
t20 = -pkin(9) * t168 + t23;
t6 = -qJD(5) * t22 - t158 * t20 + t161 * t17;
t2 = pkin(5) * t126 - qJ(6) * t34 + qJD(6) * t86 + t6;
t25 = mrSges(7,1) * t126 - mrSges(7,3) * t34;
t186 = m(7) * t2 + t25;
t152 = -pkin(4) * t162 - pkin(3);
t184 = qJD(4) * t219;
t183 = t133 * t191;
t181 = t170 * t189;
t11 = -t35 * mrSges(7,1) + t34 * mrSges(7,2);
t180 = (-mrSges(6,2) - mrSges(7,2)) * t161;
t60 = t104 * mrSges(6,1) - t103 * mrSges(6,2);
t179 = -(2 * Ifges(4,4)) - t204;
t21 = -t158 * t50 + t161 * t42;
t177 = t126 * mrSges(4,1) - t125 * mrSges(4,2);
t113 = t161 * t148 + t149 * t158;
t142 = t159 * t184;
t143 = t162 * t184;
t67 = -qJD(5) * t114 - t142 * t158 + t161 * t143;
t41 = qJ(6) * t103 - qJD(6) * t135 + t67;
t175 = m(7) * t41 + t103 * mrSges(7,3);
t79 = pkin(4) * t194 - t226;
t173 = Ifges(5,1) * t162 - t206;
t172 = -Ifges(5,2) * t159 + t205;
t171 = -t60 + t207;
t169 = t126 * t233 + t234 * t35 + t235 * t34;
t5 = t158 * t17 + t161 * t20 + t42 * t188 - t189 * t50;
t167 = t183 + t195;
t66 = t161 * t142 + t158 * t143 + t148 * t188 + t149 * t189;
t166 = -t158 * t104 + (t135 * t158 - t161 * t170) * qJD(5);
t40 = -qJ(6) * t104 - qJD(6) * t170 + t66;
t165 = t67 * mrSges(6,1) + t41 * mrSges(7,1) - t66 * mrSges(6,2) - t40 * mrSges(7,2) + t232;
t3 = qJ(6) * t35 - qJD(6) * t85 + t5;
t164 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2) + t169;
t53 = pkin(4) * t168 + t78;
t153 = Ifges(5,5) * t190;
t151 = pkin(5) + t230;
t147 = Ifges(5,1) * t159 + t205;
t141 = t173 * qJD(4);
t140 = t172 * qJD(4);
t117 = pkin(5) * t170 + t152;
t110 = Ifges(6,1) * t135 - Ifges(6,4) * t170;
t109 = Ifges(7,1) * t135 - Ifges(7,4) * t170;
t108 = Ifges(6,4) * t135 - Ifges(6,2) * t170;
t107 = Ifges(7,4) * t135 - Ifges(7,2) * t170;
t106 = mrSges(6,1) * t170 + mrSges(6,2) * t135;
t105 = mrSges(7,1) * t170 + mrSges(7,2) * t135;
t94 = mrSges(5,1) * t132 - mrSges(5,3) * t193;
t93 = -mrSges(5,2) * t132 - mrSges(5,3) * t194;
t89 = t187 + t212;
t82 = -qJ(6) * t170 + t114;
t81 = -qJ(6) * t135 + t113;
t75 = Ifges(5,5) * t132 + t133 * t173;
t74 = t133 * t172 + t200;
t69 = -mrSges(5,2) * t126 - mrSges(5,3) * t168;
t68 = mrSges(5,1) * t126 + mrSges(5,3) * t167;
t64 = -Ifges(6,1) * t103 - Ifges(6,4) * t104;
t63 = -Ifges(7,1) * t103 - Ifges(7,4) * t104;
t62 = -Ifges(6,4) * t103 - Ifges(6,2) * t104;
t61 = -Ifges(7,4) * t103 - Ifges(7,2) * t104;
t56 = mrSges(5,1) * t168 - mrSges(5,2) * t167;
t55 = mrSges(6,1) * t85 - mrSges(6,2) * t86;
t54 = mrSges(7,1) * t85 - mrSges(7,2) * t86;
t52 = pkin(5) * t85 + t79;
t48 = -Ifges(6,1) * t86 - Ifges(6,4) * t85 + Ifges(6,5) * t132;
t47 = -Ifges(7,1) * t86 - Ifges(7,4) * t85 + Ifges(7,5) * t132;
t46 = -Ifges(6,4) * t86 - Ifges(6,2) * t85 + Ifges(6,6) * t132;
t45 = -Ifges(7,4) * t86 - Ifges(7,2) * t85 + Ifges(7,6) * t132;
t44 = -Ifges(5,1) * t167 - Ifges(5,4) * t168 + t202;
t43 = -Ifges(5,4) * t167 - Ifges(5,2) * t168 + t201;
t26 = mrSges(6,1) * t126 - mrSges(6,3) * t34;
t24 = -t58 * qJD(4) + t178;
t18 = -t35 * pkin(5) + t53;
t14 = -qJ(6) * t85 + t22;
t13 = pkin(5) * t132 + qJ(6) * t86 + t21;
t12 = -mrSges(6,1) * t35 + mrSges(6,2) * t34;
t1 = [(t77 * t222 - t179 * t125 + ((2 * Ifges(4,2)) + Ifges(5,3) + t233) * t126 + t169 + t227) * t132 + 0.2e1 * t185 * t177 + (t162 * t44 - t159 * t43 - 0.2e1 * Ifges(4,1) * t125 + (Ifges(5,5) * t162 + t179) * t126 + (-t159 * t75 - t162 * t74 + t132 * (-Ifges(5,5) * t159 - Ifges(5,6) * t162)) * qJD(4) + 0.2e1 * (mrSges(4,3) + t174) * t78) * t133 + (t45 + t46) * t35 + (t47 + t48) * t34 + t56 * t220 + (t13 * t2 + t14 * t3 + t18 * t52) * t223 + (t21 * t6 + t22 * t5 + t53 * t79) * t224 + t112 * t126 * t222 + 0.2e1 * m(5) * (t23 * t58 + t24 * t57 - t203) + 0.2e1 * m(4) * (t112 * t77 - t203) + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t156 ^ 2 + t157 ^ 2) * qJD(2) - (t236 + t239) * t86 - (t237 + t240) * t85 + 0.2e1 * t13 * t25 + 0.2e1 * t21 * t26 + 0.2e1 * t14 * t27 + 0.2e1 * t22 * t28 + 0.2e1 * t52 * t11 + 0.2e1 * t18 * t54 + 0.2e1 * t53 * t55 + 0.2e1 * t57 * t68 + 0.2e1 * t58 * t69 + 0.2e1 * t3 * t70 + 0.2e1 * t5 * t71 + 0.2e1 * t2 * t72 + 0.2e1 * t6 * t73 - (mrSges(4,3) * t220 - t159 * t74 + t162 * t75) * t125 + 0.2e1 * t79 * t12 + 0.2e1 * t23 * t93 + 0.2e1 * t24 * t94; t159 * t69 + t162 * t68 + t210 * t135 - (t25 + t26) * t170 - t208 * t104 - t209 * t103 + (-t159 * t94 + t162 * t93) * qJD(4) + m(7) * (-t103 * t14 - t104 * t13 + t135 * t3 - t170 * t2) + m(6) * (-t103 * t22 - t104 * t21 + t135 * t5 - t170 * t6) + m(5) * (t159 * t23 + t162 * t24 + (-t159 * t57 + t162 * t58) * qJD(4)) + t177; 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (-t103 * t135 + t104 * t170); t236 * t135 / 0.2e1 - t237 * t170 / 0.2e1 + (t235 * t135 - t234 * t170) * t126 / 0.2e1 + (t133 * t141 / 0.2e1 - t125 * t147 / 0.2e1 + t23 * mrSges(5,3) + t43 / 0.2e1 - t78 * mrSges(5,1) + t201 / 0.2e1 + (t75 / 0.2e1 + t133 * t213 - t57 * mrSges(5,3)) * qJD(4) + (m(5) * (-qJD(4) * t57 + t23) + t69 - qJD(4) * t94) * pkin(8)) * t162 + (t140 * t216 - t125 * t213 - t24 * mrSges(5,3) + t44 / 0.2e1 + t78 * mrSges(5,2) + t202 / 0.2e1 + (-m(5) * t24 - t68) * pkin(8) + (-t74 / 0.2e1 - t58 * mrSges(5,3) + t147 * t216 - t200 / 0.2e1 + (-m(5) * t58 - t93) * pkin(8) + (m(6) * t79 + t55) * pkin(4)) * qJD(4)) * t159 - t226 * t139 + (t107 / 0.2e1 + t108 / 0.2e1) * t35 + (t109 / 0.2e1 + t110 / 0.2e1) * t34 + (-m(5) * t78 - t56) * pkin(3) + (t153 + t232) * t132 / 0.2e1 + (t103 * t13 - t104 * t14 - t135 * t2 - t170 * t3) * mrSges(7,3) + (t103 * t21 - t104 * t22 - t135 * t6 - t170 * t5) * mrSges(6,3) - t52 * t207 - (t45 / 0.2e1 + t46 / 0.2e1) * t104 - (t47 / 0.2e1 + t48 / 0.2e1) * t103 + m(7) * (t117 * t18 + t13 * t41 + t14 * t40 + t2 * t81 + t3 * t82 + t52 * t89) - (t63 / 0.2e1 + t64 / 0.2e1) * t86 - (t61 / 0.2e1 + t62 / 0.2e1) * t85 + m(6) * (t113 * t6 + t114 * t5 + t152 * t53 + t21 * t67 + t22 * t66) + t40 * t70 + t66 * t71 + t41 * t72 + t67 * t73 - t77 * mrSges(4,2) - t78 * mrSges(4,1) + t79 * t60 + t81 * t25 + t82 * t27 + t89 * t54 + t18 * t105 + t53 * t106 + t113 * t26 + t114 * t28 + t117 * t11 - Ifges(4,5) * t125 - Ifges(4,6) * t126 + t152 * t12; m(6) * (-t103 * t114 - t104 * t113 + t135 * t66 - t170 * t67) + m(7) * (-t103 * t82 - t104 * t81 + t135 * t40 - t170 * t41); -0.2e1 * pkin(3) * t139 + 0.2e1 * t89 * t105 - 0.2e1 * t117 * t207 + t162 * t140 + t159 * t141 + 0.2e1 * t152 * t60 + (t63 + t64) * t135 - (t61 + t62) * t170 - (t107 + t108) * t104 - (t109 + t110) * t103 + (t162 * t147 + (0.2e1 * pkin(4) * t106 - t146) * t159) * qJD(4) + (t113 * t67 + t114 * t66 + t152 * t187) * t224 + (t117 * t89 + t40 * t82 + t41 * t81) * t223 + 0.2e1 * (t103 * t81 - t104 * t82 - t135 * t41 - t170 * t40) * mrSges(7,3) + 0.2e1 * (t103 * t113 - t104 * t114 - t135 * t67 - t170 * t66) * mrSges(6,3); t164 - t168 * Ifges(5,6) + t186 * t151 + (t161 * t26 + t210 * t158 + (-t158 * t208 + t161 * t209) * qJD(5) + m(7) * (-t13 * t189 + t14 * t188 + t158 * t3) + m(6) * (t158 * t5 + t161 * t6 + t188 * t22 - t189 * t21)) * pkin(4) - t23 * mrSges(5,2) + t24 * mrSges(5,1) - Ifges(5,5) * t183 + t227; -t139 + m(6) * ((-t104 * t161 + t181) * pkin(4) + t197) + m(7) * (pkin(4) * t181 - t104 * t151 + t197) + t171; t153 + t175 * t151 + (-t204 + (-mrSges(5,1) * t162 + mrSges(5,2) * t159) * pkin(8)) * qJD(4) + (m(7) * (t158 * t40 + t188 * t82 - t189 * t81) + m(6) * (-t113 * t189 + t114 * t188 + t158 * t66 + t161 * t67) + t166 * mrSges(7,3) + (t161 * t103 + t166) * mrSges(6,3)) * pkin(4) + t165; 0.2e1 * (t180 + ((-t151 + t230) * m(7) + t229) * t158) * t228; pkin(5) * t186 + t164; -m(7) * t212 + t171; pkin(5) * t175 + t165; (t180 + (-m(7) * pkin(5) + t229) * t158) * t228; 0; m(7) * t18 + t11; 0; m(7) * t89 - t207; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
