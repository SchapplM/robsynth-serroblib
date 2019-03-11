% Calculate time derivative of joint inertia matrix for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:40
% EndTime: 2019-03-09 03:40:45
% DurationCPUTime: 2.88s
% Computational Cost: add. (4950->405), mult. (11350->607), div. (0->0), fcn. (10497->10), ass. (0->165)
t144 = sin(pkin(11));
t152 = cos(qJ(3));
t171 = qJD(3) * t152;
t161 = t144 * t171;
t145 = cos(pkin(11));
t186 = mrSges(5,2) * t145;
t105 = mrSges(5,1) * t161 + t171 * t186;
t148 = sin(qJ(5));
t151 = cos(qJ(5));
t124 = t144 * t151 + t145 * t148;
t112 = t124 * qJD(5);
t149 = sin(qJ(3));
t155 = t144 * t148 - t145 * t151;
t74 = -t112 * t149 - t155 * t171;
t111 = t155 * qJD(5);
t75 = t111 * t149 - t124 * t171;
t37 = -t75 * mrSges(6,1) + t74 * mrSges(6,2);
t147 = sin(qJ(6));
t150 = cos(qJ(6));
t101 = t124 * t149;
t102 = t155 * t149;
t64 = -t101 * t150 + t102 * t147;
t24 = qJD(6) * t64 + t147 * t75 + t150 * t74;
t65 = -t101 * t147 - t102 * t150;
t25 = -qJD(6) * t65 - t147 * t74 + t150 * t75;
t6 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t153 = -t37 - t6;
t206 = -t105 + t153;
t163 = -cos(pkin(10)) * pkin(1) - pkin(2);
t205 = 0.2e1 * t163;
t188 = pkin(8) + qJ(4);
t129 = t188 * t144;
t131 = t188 * t145;
t92 = -t148 * t129 + t151 * t131;
t165 = t149 * qJD(4);
t191 = pkin(3) * t149;
t110 = -t165 + (-qJ(4) * t152 + t191) * qJD(3);
t138 = sin(pkin(10)) * pkin(1) + pkin(7);
t172 = qJD(3) * t149;
t162 = t138 * t172;
t82 = t145 * t110 + t144 * t162;
t98 = t144 * t110;
t83 = -t145 * t162 + t98;
t204 = -t82 * t144 + t83 * t145;
t203 = -Ifges(6,5) * t74 - Ifges(6,6) * t75 - Ifges(6,3) * t172;
t202 = 2 * m(5);
t201 = 2 * m(6);
t200 = 2 * m(7);
t143 = t145 ^ 2;
t199 = 0.2e1 * t138;
t198 = m(7) * pkin(5);
t84 = -t124 * t147 - t150 * t155;
t197 = t84 / 0.2e1;
t85 = t124 * t150 - t147 * t155;
t196 = t85 / 0.2e1;
t195 = -t155 / 0.2e1;
t194 = t124 / 0.2e1;
t193 = t145 / 0.2e1;
t190 = pkin(5) * t112;
t189 = pkin(5) * t152;
t41 = qJD(6) * t84 - t111 * t150 - t112 * t147;
t42 = -qJD(6) * t85 + t111 * t147 - t112 * t150;
t187 = Ifges(7,5) * t41 + Ifges(7,6) * t42;
t116 = -pkin(3) * t152 - t149 * qJ(4) + t163;
t104 = t145 * t116;
t176 = t145 * t149;
t70 = -pkin(8) * t176 + t104 + (-t138 * t144 - pkin(4)) * t152;
t178 = t144 * t149;
t175 = t145 * t152;
t87 = t144 * t116 + t138 * t175;
t77 = -pkin(8) * t178 + t87;
t36 = t148 * t70 + t151 * t77;
t185 = Ifges(5,4) * t144;
t184 = Ifges(5,4) * t145;
t183 = t144 * Ifges(5,2);
t182 = t148 * t77;
t179 = -mrSges(5,1) * t145 + mrSges(5,2) * t144 - mrSges(4,1);
t177 = t144 * t152;
t132 = t149 * t138;
t174 = -Ifges(6,5) * t111 - Ifges(6,6) * t112;
t128 = t138 * t171;
t100 = pkin(4) * t161 + t128;
t109 = pkin(4) * t178 + t132;
t173 = t144 ^ 2 + t143;
t170 = qJD(4) * t144;
t169 = qJD(4) * t145;
t168 = qJD(5) * t151;
t167 = qJD(6) * t147;
t166 = qJD(6) * t150;
t164 = -Ifges(7,5) * t24 - Ifges(7,6) * t25 - Ifges(7,3) * t172;
t139 = -pkin(4) * t145 - pkin(3);
t160 = t149 * t171;
t16 = -mrSges(7,1) * t42 + t41 * mrSges(7,2);
t35 = t151 * t70 - t182;
t159 = t173 * mrSges(5,3);
t158 = t173 * qJ(4);
t91 = -t151 * t129 - t131 * t148;
t71 = -pkin(9) * t124 + t91;
t72 = -pkin(9) * t155 + t92;
t33 = -t147 * t72 + t150 * t71;
t60 = -t129 * t168 + t151 * t169 + (-qJD(5) * t131 - t170) * t148;
t47 = -pkin(9) * t112 + t60;
t61 = -t124 * qJD(4) - qJD(5) * t92;
t48 = pkin(9) * t111 + t61;
t12 = qJD(6) * t33 + t147 * t48 + t150 * t47;
t34 = t147 * t71 + t150 * t72;
t13 = -qJD(6) * t34 - t147 * t47 + t150 * t48;
t157 = t13 * mrSges(7,1) - t12 * mrSges(7,2) + t187;
t156 = -Ifges(5,5) * t145 + Ifges(5,6) * t144;
t26 = t102 * pkin(9) - t189 + t35;
t27 = -pkin(9) * t101 + t36;
t10 = -t147 * t27 + t150 * t26;
t11 = t147 * t26 + t150 * t27;
t55 = (pkin(4) * t149 - pkin(8) * t175) * qJD(3) + t82;
t66 = t98 + (-pkin(8) * t177 - t132 * t145) * qJD(3);
t15 = -qJD(5) * t36 - t148 * t66 + t151 * t55;
t7 = pkin(5) * t172 - pkin(9) * t74 + t15;
t14 = -qJD(5) * t182 + t148 * t55 + t151 * t66 + t70 * t168;
t8 = pkin(9) * t75 + t14;
t2 = qJD(6) * t10 + t147 * t7 + t150 * t8;
t3 = -qJD(6) * t11 - t147 * t8 + t150 * t7;
t154 = t3 * mrSges(7,1) - t2 * mrSges(7,2) - t164;
t126 = -mrSges(5,1) * t152 - mrSges(5,3) * t176;
t125 = mrSges(5,2) * t152 - mrSges(5,3) * t178;
t118 = (-mrSges(7,1) * t147 - mrSges(7,2) * t150) * qJD(6) * pkin(5);
t115 = (mrSges(5,1) * t149 - mrSges(5,3) * t175) * qJD(3);
t114 = (-mrSges(5,2) * t149 - mrSges(5,3) * t177) * qJD(3);
t113 = (mrSges(5,1) * t144 + t186) * t149;
t106 = t111 * mrSges(6,2);
t97 = pkin(5) * t155 + t139;
t96 = (t149 * Ifges(5,5) + (t145 * Ifges(5,1) - t185) * t152) * qJD(3);
t95 = (t149 * Ifges(5,6) + (-t183 + t184) * t152) * qJD(3);
t94 = -mrSges(6,1) * t152 + t102 * mrSges(6,3);
t93 = mrSges(6,2) * t152 - t101 * mrSges(6,3);
t90 = Ifges(6,1) * t124 - Ifges(6,4) * t155;
t89 = Ifges(6,4) * t124 - Ifges(6,2) * t155;
t88 = mrSges(6,1) * t155 + mrSges(6,2) * t124;
t86 = -t138 * t177 + t104;
t81 = -Ifges(6,1) * t111 - Ifges(6,4) * t112;
t80 = -Ifges(6,4) * t111 - Ifges(6,2) * t112;
t79 = mrSges(6,1) * t112 - t106;
t78 = pkin(5) * t101 + t109;
t76 = mrSges(6,1) * t101 - mrSges(6,2) * t102;
t63 = -Ifges(6,1) * t102 - Ifges(6,4) * t101 - Ifges(6,5) * t152;
t62 = -Ifges(6,4) * t102 - Ifges(6,2) * t101 - Ifges(6,6) * t152;
t54 = -mrSges(6,2) * t172 + mrSges(6,3) * t75;
t53 = mrSges(6,1) * t172 - mrSges(6,3) * t74;
t52 = -mrSges(7,1) * t152 - t65 * mrSges(7,3);
t51 = mrSges(7,2) * t152 + t64 * mrSges(7,3);
t46 = -pkin(5) * t75 + t100;
t45 = Ifges(7,1) * t85 + Ifges(7,4) * t84;
t44 = Ifges(7,4) * t85 + Ifges(7,2) * t84;
t43 = -mrSges(7,1) * t84 + mrSges(7,2) * t85;
t32 = -mrSges(7,1) * t64 + mrSges(7,2) * t65;
t31 = Ifges(6,1) * t74 + Ifges(6,4) * t75 + Ifges(6,5) * t172;
t30 = Ifges(6,4) * t74 + Ifges(6,2) * t75 + Ifges(6,6) * t172;
t29 = Ifges(7,1) * t65 + Ifges(7,4) * t64 - Ifges(7,5) * t152;
t28 = Ifges(7,4) * t65 + Ifges(7,2) * t64 - Ifges(7,6) * t152;
t20 = -mrSges(7,2) * t172 + mrSges(7,3) * t25;
t19 = mrSges(7,1) * t172 - mrSges(7,3) * t24;
t18 = Ifges(7,1) * t41 + Ifges(7,4) * t42;
t17 = Ifges(7,4) * t41 + Ifges(7,2) * t42;
t5 = Ifges(7,1) * t24 + Ifges(7,4) * t25 + Ifges(7,5) * t172;
t4 = Ifges(7,4) * t24 + Ifges(7,2) * t25 + Ifges(7,6) * t172;
t1 = [(t100 * t109 + t14 * t36 + t15 * t35) * t201 + (t86 * t82 + t87 * t83) * t202 + (t10 * t3 + t11 * t2 + t46 * t78) * t200 + ((mrSges(4,1) * t205 + Ifges(7,5) * t65 + Ifges(7,6) * t64 - Ifges(6,5) * t102 - Ifges(6,6) * t101 + (-(2 * Ifges(4,4)) - t156) * t149) * t149 + (t113 * t199 + mrSges(4,2) * t205 + 0.2e1 * (Ifges(4,4) + t156) * t152 + (Ifges(5,1) * t143 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - (2 * Ifges(5,3)) - Ifges(7,3) - Ifges(6,3) + t138 ^ 2 * t202 + (t183 - 0.2e1 * t184) * t144) * t149) * t152) * qJD(3) + (t164 + t203) * t152 + 0.2e1 * t83 * t125 + 0.2e1 * t82 * t126 + 0.2e1 * t87 * t114 + 0.2e1 * t86 * t115 - t102 * t31 + 0.2e1 * t109 * t37 + 0.2e1 * t100 * t76 - t101 * t30 + 0.2e1 * t14 * t93 + 0.2e1 * t15 * t94 + 0.2e1 * t78 * t6 + t74 * t63 + t75 * t62 + t64 * t4 + t65 * t5 + 0.2e1 * t46 * t32 + 0.2e1 * t2 * t51 + 0.2e1 * t3 * t52 + 0.2e1 * t35 * t53 + 0.2e1 * t36 * t54 + t24 * t29 + t25 * t28 + 0.2e1 * t10 * t19 + 0.2e1 * t11 * t20 + (t105 * t199 - t144 * t95 + t145 * t96) * t149; -t101 * t53 - t102 * t54 + t64 * t19 + t65 * t20 + t24 * t51 + t25 * t52 + t74 * t93 + t75 * t94 + (t145 * t114 - t144 * t115) * t149 + t206 * t152 + ((t125 * t145 - t126 * t144) * t152 + (t113 + t32 + t76) * t149) * qJD(3) + m(7) * (t10 * t25 + t11 * t24 - t152 * t46 + t172 * t78 + t2 * t65 + t3 * t64) + m(6) * (-t100 * t152 - t15 * t101 - t14 * t102 + t109 * t172 + t35 * t75 + t36 * t74) + m(5) * (t204 * t149 + (t138 * t149 ^ 2 + (-t138 * t152 - t144 * t86 + t145 * t87) * t152) * qJD(3)); 0.2e1 * m(7) * (t65 * t24 + t64 * t25 - t160) + 0.2e1 * m(6) * (-t101 * t75 - t102 * t74 - t160) + 0.2e1 * m(5) * (-0.1e1 + t173) * t160; (-t10 * t41 + t11 * t42 + t2 * t84 - t3 * t85) * mrSges(7,3) + m(6) * (t100 * t139 + t14 * t92 + t15 * t91 + t35 * t61 + t36 * t60) + t31 * t194 + t30 * t195 + t5 * t196 + t4 * t197 + (t95 / 0.2e1 + t83 * mrSges(5,3) + qJ(4) * t114 + qJD(4) * t125) * t145 + (t96 / 0.2e1 - t82 * mrSges(5,3) - qJ(4) * t115 - qJD(4) * t126) * t144 - (t62 / 0.2e1 - pkin(5) * t32) * t112 + m(5) * (qJ(4) * t204 + t169 * t87 - t170 * t86) + (t111 * t35 - t112 * t36 - t124 * t15 - t14 * t155) * mrSges(6,3) + t139 * t37 + m(7) * (t10 * t13 + t11 * t12 + t190 * t78 + t2 * t34 + t3 * t33 + t46 * t97) - t111 * t63 / 0.2e1 - pkin(3) * t105 + t109 * t79 + t100 * t88 - t101 * t80 / 0.2e1 - t102 * t81 / 0.2e1 + t91 * t53 + t92 * t54 + t60 * t93 + t61 * t94 + t97 * t6 - (t187 + t174) * t152 / 0.2e1 + t75 * t89 / 0.2e1 + t74 * t90 / 0.2e1 + t78 * t16 + t64 * t17 / 0.2e1 + t65 * t18 / 0.2e1 + t46 * t43 + t12 * t51 + t13 * t52 + t41 * t29 / 0.2e1 + t42 * t28 / 0.2e1 + t25 * t44 / 0.2e1 + t24 * t45 / 0.2e1 + t33 * t19 + t34 * t20 + ((t138 * mrSges(4,2) - Ifges(4,6) + Ifges(7,5) * t196 + Ifges(7,6) * t197 + Ifges(6,5) * t194 + Ifges(6,6) * t195 + Ifges(5,5) * t144 / 0.2e1 + Ifges(5,6) * t193) * t149 + ((Ifges(5,1) * t144 + t184) * t193 - t144 * (Ifges(5,2) * t145 + t185) / 0.2e1 + Ifges(4,5) + (-m(5) * pkin(3) + t179) * t138) * t152) * qJD(3); (-t16 - t79) * t152 + ((-mrSges(4,2) + t159) * t152 + (t43 + t88 + t179) * t149) * qJD(3) + m(7) * (-t112 * t189 + t12 * t65 + t13 * t64 + t172 * t97 + t34 * t24 + t33 * t25) + m(6) * (-t101 * t61 - t102 * t60 + t139 * t172 + t74 * t92 + t75 * t91) + m(5) * (t173 * t165 + (t152 * t158 - t191) * qJD(3)) + (t24 * t84 - t25 * t85 - t41 * t64 + t42 * t65) * mrSges(7,3) + (-t101 * t111 + t102 * t112 - t124 * t75 - t155 * t74) * mrSges(6,3); -t111 * t90 - t155 * t80 + t124 * t81 + 0.2e1 * t139 * t79 + 0.2e1 * t97 * t16 + t84 * t17 + t85 * t18 + t41 * t45 + t42 * t44 - (-0.2e1 * pkin(5) * t43 + t89) * t112 + (t12 * t34 + t13 * t33 + t190 * t97) * t200 + (t60 * t92 + t61 * t91) * t201 + 0.2e1 * (t12 * t84 - t13 * t85 - t33 * t41 + t34 * t42) * mrSges(7,3) + 0.2e1 * (t111 * t91 - t112 * t92 - t124 * t61 - t155 * t60) * mrSges(6,3) + (t158 * t202 + 0.2e1 * t159) * qJD(4); m(5) * t128 + m(6) * t100 + m(7) * t46 - t206; (m(5) + m(6) + m(7)) * t172; -t106 - (-mrSges(6,1) - t198) * t112 + t16; 0; t15 * mrSges(6,1) - t14 * mrSges(6,2) + (m(7) * (-t10 * t167 + t11 * t166 + t147 * t2 + t150 * t3) - t52 * t167 + t150 * t19 + t51 * t166 + t147 * t20) * pkin(5) + t154 - t203; (t147 * t24 + t150 * t25 + (-t147 * t64 + t150 * t65) * qJD(6)) * t198 + t153; t61 * mrSges(6,1) - t60 * mrSges(6,2) + (m(7) * (t12 * t147 + t13 * t150 + (-t147 * t33 + t150 * t34) * qJD(6)) + (t147 * t42 - t150 * t41 + (t147 * t85 + t150 * t84) * qJD(6)) * mrSges(7,3)) * pkin(5) + t157 + t174; 0; 0.2e1 * t118; t154; -t6; t157; 0; t118; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
