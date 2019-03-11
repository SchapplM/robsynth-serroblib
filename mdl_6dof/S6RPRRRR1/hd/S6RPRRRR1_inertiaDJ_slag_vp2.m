% Calculate time derivative of joint inertia matrix for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:54:13
% EndTime: 2019-03-09 06:54:19
% DurationCPUTime: 3.14s
% Computational Cost: add. (7403->327), mult. (15361->499), div. (0->0), fcn. (14810->10), ass. (0->159)
t217 = qJD(3) + qJD(4);
t115 = sin(qJ(6));
t112 = t115 ^ 2;
t119 = cos(qJ(6));
t113 = t119 ^ 2;
t166 = t112 + t113;
t117 = sin(qJ(4));
t121 = cos(qJ(4));
t118 = sin(qJ(3));
t106 = sin(pkin(11)) * pkin(1) + pkin(7);
t193 = pkin(8) + t106;
t153 = t193 * t118;
t122 = cos(qJ(3));
t221 = t193 * t122;
t63 = -t117 * t221 - t121 * t153;
t162 = qJD(6) * t119;
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t93 = -t117 * t118 + t121 * t122;
t94 = t117 * t122 + t118 * t121;
t68 = t116 * t94 - t120 * t93;
t74 = t217 * t93;
t75 = t217 * t94;
t40 = -qJD(5) * t68 - t116 * t75 + t120 * t74;
t183 = t115 * t40;
t69 = t116 * t93 + t120 * t94;
t224 = t162 * t69 + t183;
t64 = -t117 * t153 + t121 * t221;
t130 = -t94 * pkin(9) + t63;
t50 = pkin(9) * t93 + t64;
t30 = t116 * t130 + t120 * t50;
t156 = -cos(pkin(11)) * pkin(1) - pkin(2);
t98 = -pkin(3) * t122 + t156;
t76 = -pkin(4) * t93 + t98;
t42 = pkin(5) * t68 - pkin(10) * t69 + t76;
t17 = t115 * t42 + t119 * t30;
t167 = t17 * qJD(6);
t45 = t217 * t64;
t125 = -t74 * pkin(9) - t45;
t165 = qJD(5) * t116;
t44 = t217 * t63;
t222 = -pkin(9) * t75 + qJD(5) * t130 + t44;
t10 = t116 * t125 + t120 * t222 - t165 * t50;
t41 = qJD(5) * t69 + t116 * t74 + t120 * t75;
t207 = pkin(5) * t41;
t61 = pkin(3) * qJD(3) * t118 + pkin(4) * t75;
t15 = -pkin(10) * t40 + t207 + t61;
t3 = -t10 * t115 + t119 * t15 - t167;
t223 = -t3 - t167;
t163 = qJD(6) * t115;
t159 = t69 * t163;
t177 = t119 * t40;
t135 = t159 - t177;
t220 = t166 * t120;
t107 = pkin(4) * t116 + pkin(10);
t219 = t166 * t107;
t16 = -t115 * t30 + t119 * t42;
t218 = -t115 * t16 + t119 * t17;
t216 = -mrSges(5,1) * t45 - t44 * mrSges(5,2) + Ifges(5,5) * t74 - Ifges(5,6) * t75;
t164 = qJD(5) * t120;
t160 = pkin(4) * t164;
t147 = mrSges(7,3) * t160;
t102 = t112 * t147;
t103 = t113 * t147;
t108 = -pkin(4) * t120 - pkin(5);
t142 = mrSges(7,1) * t115 + mrSges(7,2) * t119;
t95 = t142 * qJD(6);
t80 = t108 * t95;
t191 = mrSges(7,1) * t119;
t99 = mrSges(7,2) * t115 - t191;
t90 = pkin(4) * t99 * t165;
t215 = t102 + t103 + t80 + t90;
t214 = 2 * m(6);
t213 = 2 * m(7);
t11 = t116 * t222 - t120 * t125 + t164 * t50;
t212 = 0.2e1 * t11;
t211 = 0.2e1 * t61;
t210 = 0.2e1 * t98;
t209 = m(6) / 0.2e1;
t208 = m(5) * pkin(3);
t206 = pkin(5) * t95;
t29 = t116 * t50 - t120 * t130;
t204 = t11 * t29;
t203 = t115 * t3;
t175 = qJD(6) * t16;
t2 = t10 * t119 + t115 * t15 + t175;
t202 = t119 * t2;
t109 = pkin(3) * t121 + pkin(4);
t168 = t117 * t120;
t66 = t109 * t165 + (t117 * t164 + (t116 * t121 + t168) * qJD(4)) * pkin(3);
t201 = t29 * t66;
t200 = t41 * t68;
t169 = t116 * t117;
t65 = t109 * t164 + (-t117 * t165 + (t120 * t121 - t169) * qJD(4)) * pkin(3);
t198 = t65 * mrSges(6,2);
t197 = t65 * t69;
t196 = t66 * t68;
t195 = t74 * t94;
t194 = t75 * t93;
t192 = Ifges(7,5) * t177 + Ifges(7,3) * t41;
t190 = mrSges(7,3) * t112;
t189 = mrSges(7,3) * t113;
t188 = Ifges(7,4) * t115;
t187 = Ifges(7,4) * t119;
t186 = Ifges(7,6) * t115;
t185 = pkin(4) * qJD(5);
t182 = t115 * t69;
t180 = t116 * t68;
t176 = t119 * t69;
t86 = pkin(3) * t168 + t109 * t116;
t84 = pkin(10) + t86;
t174 = qJD(6) * t84;
t161 = 0.2e1 * t122;
t12 = -mrSges(7,1) * t224 + mrSges(7,2) * t135;
t157 = m(7) * t11 - t12;
t155 = mrSges(5,1) * t75 + t74 * mrSges(5,2);
t154 = mrSges(6,1) * t41 + t40 * mrSges(6,2);
t151 = -t163 / 0.2e1;
t150 = -(2 * Ifges(6,4)) - t186;
t149 = t166 * t40;
t148 = t166 * t65;
t145 = t11 * t68 + t29 * t41;
t143 = -t116 * mrSges(6,1) - t120 * mrSges(6,2);
t141 = Ifges(7,1) * t119 - t188;
t140 = -Ifges(7,2) * t115 + t187;
t139 = Ifges(7,5) * t115 + Ifges(7,6) * t119;
t13 = mrSges(7,1) * t41 + mrSges(7,3) * t135;
t14 = -mrSges(7,2) * t41 - mrSges(7,3) * t224;
t138 = -t115 * t13 + t119 * t14;
t47 = -mrSges(7,2) * t68 - mrSges(7,3) * t182;
t48 = mrSges(7,1) * t68 - mrSges(7,3) * t176;
t137 = -t115 * t48 + t119 * t47;
t85 = -pkin(3) * t169 + t109 * t120;
t100 = Ifges(7,2) * t119 + t188;
t101 = Ifges(7,1) * t115 + t187;
t96 = t140 * qJD(6);
t97 = t141 * qJD(6);
t134 = -t100 * t163 + t101 * t162 + t115 * t97 + t119 * t96;
t133 = t41 * t99 + t68 * t95 - t154 + (t189 + t190) * t40;
t132 = (-mrSges(5,1) * t117 - mrSges(5,2) * t121) * qJD(4) * pkin(3);
t131 = -t203 + (-t115 * t17 - t119 * t16) * qJD(6);
t129 = t133 - t155;
t56 = t66 * t99;
t58 = t65 * t190;
t59 = t65 * t189;
t62 = t66 * mrSges(6,1);
t83 = -pkin(5) - t85;
t73 = t83 * t95;
t127 = t134 + t56 + t58 + t59 - t62 + t73 - t198;
t110 = Ifges(7,5) * t162;
t27 = Ifges(7,6) * t68 + t140 * t69;
t28 = Ifges(7,5) * t68 + t141 * t69;
t6 = -Ifges(7,4) * t135 - Ifges(7,2) * t224 + Ifges(7,6) * t41;
t7 = -Ifges(7,1) * t135 - Ifges(7,4) * t224 + Ifges(7,5) * t41;
t126 = -t10 * mrSges(6,2) + mrSges(7,3) * t202 + t29 * t95 + t27 * t151 + t28 * t162 / 0.2e1 + Ifges(6,5) * t40 + t115 * t7 / 0.2e1 + t119 * t6 / 0.2e1 - t96 * t182 / 0.2e1 + t97 * t176 / 0.2e1 + t68 * (-Ifges(7,6) * t163 + t110) / 0.2e1 + (t139 / 0.2e1 - Ifges(6,6)) * t41 - t224 * t100 / 0.2e1 + (t99 - mrSges(6,1)) * t11 + (t177 / 0.2e1 + t69 * t151) * t101;
t124 = -t48 * t162 - t47 * t163 + m(7) * (-t16 * t162 - t163 * t17 + t202 - t203) + t138;
t123 = mrSges(7,3) * t131 + t126;
t46 = t142 * t69;
t1 = [t28 * t177 - t27 * t183 + 0.2e1 * t16 * t13 + 0.2e1 * t17 * t14 - 0.2e1 * t29 * t12 + t46 * t212 + 0.2e1 * t2 * t47 + 0.2e1 * t3 * t48 + 0.2e1 * t76 * t154 - 0.2e1 * Ifges(5,2) * t194 + 0.2e1 * Ifges(5,1) * t195 + t155 * t210 + 0.2e1 * m(5) * (t44 * t64 - t45 * t63) + (t10 * t30 + t61 * t76 + t204) * t214 + (t16 * t3 + t17 * t2 + t204) * t213 + ((mrSges(4,2) * t156 + Ifges(4,4) * t122) * t161 + (0.2e1 * pkin(3) * (-mrSges(5,1) * t93 + mrSges(5,2) * t94) + t208 * t210 + 0.2e1 * t156 * mrSges(4,1) - 0.2e1 * Ifges(4,4) * t118 + (Ifges(4,1) - Ifges(4,2)) * t161) * t118) * qJD(3) + (mrSges(6,1) * t211 - 0.2e1 * t10 * mrSges(6,3) + ((2 * Ifges(6,2)) + Ifges(7,3)) * t41 + t150 * t40 + t192) * t68 + (mrSges(6,2) * t211 + mrSges(6,3) * t212 + 0.2e1 * Ifges(6,1) * t40 - t115 * t6 + t119 * t7 + (Ifges(7,5) * t119 + t150) * t41 + (-t115 * t28 - t119 * t27 - t139 * t68) * qJD(6)) * t69 + 0.2e1 * (t29 * t40 - t30 * t41) * mrSges(6,3) + 0.2e1 * (t74 * t93 - t75 * t94) * Ifges(5,4) + 0.2e1 * (t44 * t93 + t45 * t94 - t63 * t74 - t64 * t75) * mrSges(5,3); -t68 * t12 + t41 * t46 + t137 * t40 + ((-t115 * t47 - t119 * t48) * qJD(6) + t138) * t69 + m(7) * (t218 * t40 + (t131 + t202) * t69 + t145) + m(6) * (t10 * t69 + t30 * t40 + t145) + m(5) * (t44 * t94 - t45 * t93 - t63 * t75 + t64 * t74); 0.2e1 * m(7) * (t149 * t69 + t200) + 0.2e1 * m(6) * (t40 * t69 + t200) + 0.2e1 * m(5) * (-t194 + t195); t126 + (-mrSges(7,3) * t175 - t48 * t174 + m(7) * (-t16 * t174 + t17 * t65 + t2 * t84) + t84 * t14 + t65 * t47) * t119 + (m(5) * (t117 * t44 - t121 * t45 + (-t117 * t63 + t121 * t64) * qJD(4)) + (-t117 * t75 - t121 * t74 + (t117 * t94 + t121 * t93) * qJD(4)) * mrSges(5,3)) * pkin(3) + (-t47 * t174 + t223 * mrSges(7,3) + (-m(7) * t16 - t48) * t65 + (m(7) * t223 - t13) * t84) * t115 + (-t40 * t85 - t41 * t86 - t65 * t68 + t66 * t69) * mrSges(6,3) + (Ifges(4,5) * t122 - Ifges(4,6) * t118 + (-mrSges(4,1) * t122 + mrSges(4,2) * t118) * t106) * qJD(3) + m(7) * (t11 * t83 + t201) + m(6) * (t10 * t86 - t11 * t85 + t30 * t65 + t201) + t66 * t46 - t83 * t12 + t216; (-mrSges(4,1) * t118 - mrSges(4,2) * t122) * qJD(3) + m(7) * (t41 * t83 + t196 + t166 * (t40 * t84 + t197)) + m(6) * (t40 * t86 - t41 * t85 + t196 + t197) + (t117 * t74 - t121 * t75 + (-t117 * t93 + t121 * t94) * qJD(4)) * t208 + t129; -0.2e1 * t198 + 0.2e1 * t56 + 0.2e1 * t58 + 0.2e1 * t59 - 0.2e1 * t62 + 0.2e1 * t73 + 0.2e1 * t132 + (t148 * t84 + t66 * t83) * t213 + (t65 * t86 - t66 * t85) * t214 + t134; t123 + t157 * t108 + t124 * t107 + (m(6) * (t10 * t116 - t11 * t120) + (-t116 * t41 - t120 * t40) * mrSges(6,3) + ((m(6) * t30 + m(7) * t218 - t68 * mrSges(6,3) + t137) * t120 + (t69 * mrSges(6,3) + t46 + (m(7) + m(6)) * t29) * t116) * qJD(5)) * pkin(4) + t216; m(7) * (t108 * t41 + t219 * t40) + 0.2e1 * ((t116 * t40 - t120 * t41) * t209 + (m(7) * (t220 * t69 + t180) / 0.2e1 + (t120 * t69 + t180) * t209) * qJD(5)) * pkin(4) + t129; t127 + t132 + m(7) * (t108 * t66 + t219 * t65) + (m(6) * (t116 * t65 - t120 * t66) + (m(7) * (t116 * t83 + t220 * t84) + m(6) * (-t116 * t85 + t120 * t86) + t143) * qJD(5)) * pkin(4) + t215; 0.2e1 * t102 + 0.2e1 * t103 + 0.2e1 * t80 + 0.2e1 * t90 + 0.2e1 * (m(7) * (t107 * t220 + t108 * t116) + t143) * t185 + t134; -pkin(5) * t157 + pkin(10) * t124 + t123; m(7) * (pkin(10) * t149 - t207) + t133; m(7) * (-pkin(5) * t66 + pkin(10) * t148) - t206 + t127; -t206 + (m(7) * (-pkin(5) * t116 + pkin(10) * t220) + t143) * t185 + t134 + t215; t134 - 0.2e1 * t206; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t159 - Ifges(7,6) * t224 + t192; t12; t110 - t142 * t65 + (-t84 * t191 + (mrSges(7,2) * t84 - Ifges(7,6)) * t115) * qJD(6); t110 - t142 * t160 + (t107 * t99 - t186) * qJD(6); t110 + (pkin(10) * t99 - t186) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
