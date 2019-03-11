% Calculate time derivative of joint inertia matrix for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:58
% EndTime: 2019-03-08 20:41:03
% DurationCPUTime: 2.27s
% Computational Cost: add. (2819->310), mult. (6442->468), div. (0->0), fcn. (5923->10), ass. (0->156)
t106 = sin(qJ(6));
t109 = cos(qJ(6));
t82 = -mrSges(7,1) * t109 + mrSges(7,2) * t106;
t233 = t82 - mrSges(6,1);
t102 = t106 ^ 2;
t103 = t109 ^ 2;
t229 = t102 + t103;
t107 = sin(qJ(5));
t216 = cos(qJ(5));
t217 = cos(qJ(4));
t141 = t216 * t217;
t215 = sin(qJ(4));
t149 = qJD(4) * t215;
t74 = t107 * t215 - t141;
t53 = qJD(4) * t141 - qJD(5) * t74 - t107 * t149;
t147 = t229 * t53;
t175 = qJD(6) * t106;
t140 = t216 * t215;
t150 = qJD(4) * t217;
t75 = t107 * t217 + t140;
t52 = -qJD(4) * t140 - qJD(5) * t75 - t107 * t150;
t184 = t109 * t52;
t126 = t74 * t175 + t184;
t203 = t74 * t52;
t232 = -Ifges(5,1) + Ifges(5,2);
t111 = -pkin(2) - pkin(8);
t231 = -pkin(9) + t111;
t230 = mrSges(5,1) * t215 + mrSges(5,2) * t217 + mrSges(4,3);
t228 = -mrSges(5,1) * t149 - mrSges(5,2) * t150;
t174 = qJD(6) * t109;
t191 = t106 * t52;
t127 = t74 * t174 - t191;
t190 = t106 * t74;
t46 = -t75 * mrSges(7,2) + mrSges(7,3) * t190;
t183 = t109 * t74;
t47 = t75 * mrSges(7,1) + mrSges(7,3) * t183;
t227 = -t47 * t174 - t46 * t175;
t95 = t215 * pkin(4) + qJ(3);
t41 = t75 * pkin(5) + t74 * pkin(10) + t95;
t123 = t231 * t217;
t80 = t231 * t215;
t60 = t107 * t123 + t216 * t80;
t20 = -t106 * t60 + t109 * t41;
t21 = t106 * t41 + t109 * t60;
t226 = -t20 * t174 - t21 * t175;
t104 = sin(pkin(6));
t108 = sin(qJ(2));
t179 = t104 * t108;
t152 = qJD(2) * t179;
t105 = cos(pkin(6));
t110 = cos(qJ(2));
t178 = t104 * t110;
t67 = -t105 * t215 - t178 * t217;
t57 = qJD(4) * t67 + t152 * t215;
t121 = -t105 * t217 + t178 * t215;
t58 = qJD(4) * t121 + t152 * t217;
t225 = t215 * t57 + t217 * t58;
t224 = -t106 * t47 + t109 * t46;
t223 = -t107 * t53 - t216 * t52;
t120 = qJD(4) * t80;
t73 = t123 * qJD(4);
t24 = qJD(5) * t60 + t107 * t73 + t120 * t216;
t117 = t216 * t123;
t59 = t107 * t80 - t117;
t137 = t74 * t24 - t52 * t59;
t128 = t107 * t121 + t216 * t67;
t39 = t107 * t67 - t121 * t216;
t14 = qJD(5) * t39 + t107 * t57 - t216 * t58;
t222 = -t128 * t52 - t74 * t14;
t221 = 2 * qJD(3);
t220 = m(6) * pkin(4);
t134 = mrSges(7,1) * t106 + mrSges(7,2) * t109;
t76 = t134 * qJD(6);
t218 = pkin(5) * t76;
t87 = pkin(4) * t150 + qJD(3);
t19 = t53 * pkin(5) - t52 * pkin(10) + t87;
t176 = qJD(5) * t107;
t23 = qJD(5) * t117 - t107 * t120 - t176 * t80 + t216 * t73;
t3 = -qJD(6) * t21 - t106 * t23 + t109 * t19;
t213 = t106 * t3;
t2 = qJD(6) * t20 + t106 * t19 + t109 * t23;
t212 = t109 * t2;
t211 = t128 * t14;
t208 = t53 * mrSges(6,3);
t207 = t59 * t24;
t206 = t74 * mrSges(6,3);
t202 = t75 * mrSges(6,3);
t170 = t216 * pkin(4);
t98 = -t170 - pkin(5);
t201 = t98 * t76;
t177 = qJD(2) * t110;
t151 = t104 * t177;
t200 = qJ(3) * t151 + qJD(3) * t179;
t199 = mrSges(7,3) * t109;
t198 = Ifges(7,4) * t106;
t197 = Ifges(7,4) * t109;
t196 = Ifges(7,6) * t106;
t195 = pkin(4) * qJD(5);
t194 = t106 * mrSges(7,3);
t16 = t53 * mrSges(7,1) - mrSges(7,3) * t126;
t193 = t106 * t16;
t189 = t107 * t128;
t187 = t107 * t59;
t186 = t107 * t74;
t97 = pkin(4) * t107 + pkin(10);
t182 = t109 * t97;
t181 = t53 * t106;
t180 = t53 * t109;
t173 = -0.2e1 * t208;
t172 = m(7) * t195;
t169 = pkin(4) * t176;
t156 = t106 * t216;
t153 = t109 * t216;
t146 = mrSges(6,1) * t169;
t145 = t82 * t169;
t144 = qJD(5) * t170;
t139 = -t128 * t24 + t59 * t14;
t136 = mrSges(7,3) * t144;
t135 = mrSges(6,2) * t144;
t133 = Ifges(7,1) * t109 - t198;
t132 = -Ifges(7,2) * t106 + t197;
t13 = qJD(5) * t128 + t107 * t58 + t216 * t57;
t28 = -t106 * t39 + t109 * t179;
t5 = qJD(6) * t28 + t106 * t151 + t109 * t13;
t131 = -t13 * mrSges(6,2) - t128 * t76 + t14 * t233 + t5 * t199;
t130 = t102 * t136;
t129 = t103 * t136;
t29 = t106 * t179 + t109 * t39;
t78 = t132 * qJD(6);
t79 = t133 * qJD(6);
t84 = Ifges(7,2) * t109 + t198;
t85 = Ifges(7,1) * t106 + t197;
t125 = t106 * t79 + t109 * t78 + t85 * t174 - t175 * t84;
t124 = -t53 * mrSges(6,2) + mrSges(7,3) * t147 - t233 * t52 + t74 * t76;
t122 = t229 * t216;
t119 = -t213 + (-t106 * t21 - t109 * t20) * qJD(6);
t6 = -qJD(6) * t29 - t106 * t13 + t109 * t151;
t118 = -t106 * t6 + (-t106 * t29 - t109 * t28) * qJD(6);
t116 = t109 * t5 + t118;
t115 = Ifges(7,5) * t126 + Ifges(7,6) * t127 + Ifges(7,3) * t53;
t114 = (-t121 * t217 - t215 * t67) * qJD(4) + t225;
t10 = Ifges(7,1) * t126 + Ifges(7,4) * t127 + Ifges(7,5) * t53;
t31 = Ifges(7,6) * t75 - t132 * t74;
t32 = Ifges(7,5) * t75 - t133 * t74;
t9 = Ifges(7,4) * t126 + Ifges(7,2) * t127 + Ifges(7,6) * t53;
t99 = Ifges(7,5) * t174;
t113 = -t23 * mrSges(6,2) + t2 * t199 - t31 * t175 / 0.2e1 - t84 * t191 / 0.2e1 + t85 * t184 / 0.2e1 + t53 * (Ifges(7,5) * t106 + Ifges(7,6) * t109) / 0.2e1 + t59 * t76 - Ifges(6,6) * t53 + Ifges(6,5) * t52 - t79 * t183 / 0.2e1 + t75 * (-Ifges(7,6) * t175 + t99) / 0.2e1 + t106 * t10 / 0.2e1 + t109 * t9 / 0.2e1 + t233 * t24 + (qJD(6) * t85 + t78) * t190 / 0.2e1 + (t74 * t84 + t32) * t174 / 0.2e1;
t17 = -t53 * mrSges(7,2) + mrSges(7,3) * t127;
t112 = m(7) * (t212 - t213 + t226) + t109 * t17 - t193;
t81 = t104 ^ 2 * t108 * t177;
t77 = (mrSges(5,1) * t217 - mrSges(5,2) * t215) * qJD(4);
t56 = mrSges(6,1) * t75 - mrSges(6,2) * t74;
t40 = t134 * t74;
t25 = mrSges(6,1) * t53 + mrSges(6,2) * t52;
t15 = -mrSges(7,1) * t127 + mrSges(7,2) * t126;
t1 = [0.2e1 * m(7) * (t28 * t6 + t29 * t5 - t211) + 0.2e1 * m(6) * (t39 * t13 - t211 + t81) + 0.2e1 * m(5) * (-t121 * t57 + t67 * t58 + t81); m(7) * (t2 * t29 + t20 * t6 + t21 * t5 + t28 * t3 + t139) + t28 * t16 + t29 * t17 - t128 * t15 - t14 * t40 + t5 * t46 + t6 * t47 + m(6) * (t60 * t13 + t23 * t39 + (t108 * t87 + t177 * t95) * t104 + t139) + m(5) * (t111 * t114 + t200) + m(4) * t200 - t13 * t202 - t39 * t208 + (t25 + t77) * t179 + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2)) * t152 + t222 * mrSges(6,3) + (t56 - mrSges(3,2) + t230) * t151 + (t121 * t150 + t149 * t67 - t225) * mrSges(5,3); t9 * t190 + t60 * t173 + 0.2e1 * (-t75 * t52 + t74 * t53) * Ifges(6,4) - 0.2e1 * t137 * mrSges(6,3) + ((m(5) + m(4)) * t221 + 0.2e1 * t77) * qJ(3) + t75 * t115 + ((-Ifges(7,5) * t109 + t196) * t74 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t75) * t53 + t230 * t221 + (0.2e1 * Ifges(5,4) * t215 + t217 * t232) * t149 + (-0.2e1 * Ifges(5,4) * t217 + t215 * t232) * t150 - 0.2e1 * Ifges(6,1) * t203 + t126 * t32 + t127 * t31 + 0.2e1 * t20 * t16 + 0.2e1 * t21 * t17 - 0.2e1 * t24 * t40 + 0.2e1 * t2 * t46 + 0.2e1 * t3 * t47 + 0.2e1 * t59 * t15 + 0.2e1 * t87 * t56 + 0.2e1 * t95 * t25 - t10 * t183 - 0.2e1 * t23 * t202 + 0.2e1 * m(7) * (t2 * t21 + t20 * t3 + t207) + 0.2e1 * m(6) * (t60 * t23 + t87 * t95 + t207); m(4) * t152 + m(7) * ((-t106 * t28 + t109 * t29) * t53 + t116 * t75 - t222) + m(6) * (t13 * t75 + t39 * t53 - t222) + m(5) * t114; t74 * t15 + t224 * t53 + (t40 + 0.2e1 * t206) * t52 + m(7) * (t180 * t21 - t181 * t20 + t137) + m(6) * (t60 * t53 + t137) + (t173 + (-t106 * t46 - t109 * t47) * qJD(6) + m(6) * t23 + t112) * t75; 0.2e1 * m(6) * (t75 * t53 - t203) + 0.2e1 * m(7) * (t147 * t75 - t203); m(7) * (t98 * t14 + (t153 * t29 - t156 * t28 - t189) * t195 + t116 * t97) - t6 * t194 + (-t216 * t14 + t107 * t13 + (t216 * t39 - t189) * qJD(5)) * t220 - t57 * mrSges(5,2) + t58 * mrSges(5,1) + t131 + (-t174 * t28 - t175 * t29) * mrSges(7,3); m(7) * (t98 * t24 + (t153 * t21 - t156 * t20 + t187) * t195) + (-t216 * t24 + t107 * t23 + (t216 * t60 + t187) * qJD(5)) * t220 + t113 - t3 * t194 - Ifges(5,5) * t149 - Ifges(5,6) * t150 + t17 * t182 + t98 * t15 + t223 * mrSges(6,3) * pkin(4) + (-t40 - t206) * t169 + t228 * t111 + t226 * mrSges(7,3) + (m(7) * (t119 + t212) - t193 + t227) * t97 + (-t202 + t224) * t144; ((t216 * t75 + t186) * qJD(5) - t223) * t220 + m(7) * (-t98 * t52 + t97 * t147 + (t122 * t75 + t186) * t195) + t124 + t228; 0.2e1 * t145 + 0.2e1 * t201 + 0.2e1 * (t107 * t98 + t122 * t97) * t172 - 0.2e1 * t135 - 0.2e1 * t146 + 0.2e1 * t129 + 0.2e1 * t130 + t125; m(7) * (-pkin(5) * t14 + pkin(10) * t116) + t118 * mrSges(7,3) + t131; (-m(7) * t24 - t15) * pkin(5) + t113 + (t112 + t227) * pkin(10) + t119 * mrSges(7,3); m(7) * (pkin(5) * t52 + pkin(10) * t147) + t124; -t218 + t145 + t201 + (-pkin(5) * t107 + pkin(10) * t122) * t172 + t129 + t130 - t146 - t135 + t125; t125 - 0.2e1 * t218; mrSges(7,1) * t6 - mrSges(7,2) * t5; t3 * mrSges(7,1) - t2 * mrSges(7,2) + t115; (t175 * t75 - t180) * mrSges(7,2) + (-t174 * t75 - t181) * mrSges(7,1); t99 + (-mrSges(7,1) * t156 - mrSges(7,2) * t153) * t195 + (-mrSges(7,1) * t182 + (mrSges(7,2) * t97 - Ifges(7,6)) * t106) * qJD(6); t99 + (pkin(10) * t82 - t196) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
