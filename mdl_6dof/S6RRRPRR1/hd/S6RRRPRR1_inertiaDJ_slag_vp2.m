% Calculate time derivative of joint inertia matrix for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:02:14
% EndTime: 2019-03-09 18:02:21
% DurationCPUTime: 3.46s
% Computational Cost: add. (13540->357), mult. (29124->535), div. (0->0), fcn. (30030->10), ass. (0->175)
t144 = sin(qJ(6));
t148 = cos(qJ(6));
t184 = qJD(6) * t148;
t145 = sin(qJ(5));
t222 = cos(qJ(5));
t146 = sin(qJ(3));
t147 = sin(qJ(2));
t149 = cos(qJ(3));
t150 = cos(qJ(2));
t117 = -t146 * t147 + t149 * t150;
t118 = t146 * t150 + t149 * t147;
t142 = sin(pkin(11));
t143 = cos(pkin(11));
t86 = t117 * t143 - t118 * t142;
t87 = t117 * t142 + t118 * t143;
t161 = -t145 * t87 + t222 * t86;
t236 = qJD(2) + qJD(3);
t92 = t236 * t117;
t93 = t236 * t118;
t68 = -t142 * t92 - t143 * t93;
t69 = -t142 * t93 + t143 * t92;
t36 = qJD(5) * t161 + t145 * t68 + t222 * t69;
t63 = t145 * t86 + t222 * t87;
t160 = t144 * t36 + t63 * t184;
t225 = -pkin(8) - pkin(7);
t132 = t225 * t147;
t133 = t225 * t150;
t241 = t149 * t132 + t133 * t146;
t82 = -qJ(4) * t118 + t241;
t95 = t146 * t132 - t149 * t133;
t83 = qJ(4) * t117 + t95;
t57 = -t142 * t83 + t143 * t82;
t162 = -pkin(9) * t87 + t57;
t58 = t142 * t82 + t143 * t83;
t45 = pkin(9) * t86 + t58;
t22 = t145 * t162 + t222 * t45;
t137 = -pkin(2) * t150 - pkin(1);
t104 = -t117 * pkin(3) + t137;
t74 = -t86 * pkin(4) + t104;
t40 = -pkin(5) * t161 - t63 * pkin(10) + t74;
t17 = t144 * t40 + t148 * t22;
t189 = t17 * qJD(6);
t37 = qJD(5) * t63 + t145 * t69 - t222 * t68;
t84 = qJD(2) * t147 * pkin(2) + pkin(3) * t93;
t48 = -pkin(4) * t68 + t84;
t13 = pkin(5) * t37 - pkin(10) * t36 + t48;
t179 = qJD(2) * t225;
t125 = t147 * t179;
t126 = t150 * t179;
t165 = -t125 * t146 + t126 * t149;
t166 = t125 * t149 + t126 * t146;
t238 = -qJ(4) * t93 + qJD(4) * t117;
t239 = -qJ(4) * t92 - qJD(4) * t118;
t38 = -t142 * (t166 + t238) + t143 * (t165 + t239) + (-t142 * t241 - t143 * t95) * qJD(3);
t151 = -t69 * pkin(9) + t38;
t71 = qJD(3) * t241 + t166;
t72 = -qJD(3) * t95 + t165;
t39 = t143 * (t71 + t238) + t142 * (t72 + t239);
t19 = pkin(9) * t68 + t39;
t240 = -t145 * t45 + t222 * t162;
t6 = qJD(5) * t240 + t145 * t151 + t222 * t19;
t3 = t13 * t148 - t144 * t6 - t189;
t242 = -t3 - t189;
t135 = pkin(3) * t143 + pkin(4);
t221 = pkin(3) * t142;
t113 = t145 * t135 + t222 * t221;
t16 = -t144 * t22 + t148 * t40;
t196 = t148 * t17;
t237 = -t144 * t16 + t196;
t127 = -mrSges(7,1) * t148 + mrSges(7,2) * t144;
t190 = t143 * t146;
t207 = pkin(2) * qJD(3);
t108 = (-t142 * t149 - t190) * t207;
t191 = t142 * t146;
t109 = (t143 * t149 - t191) * t207;
t136 = pkin(2) * t149 + pkin(3);
t110 = -pkin(2) * t191 + t143 * t136;
t105 = pkin(4) + t110;
t112 = pkin(2) * t190 + t136 * t142;
t81 = t145 * t105 + t222 * t112;
t61 = qJD(5) * t81 - t108 * t222 + t145 * t109;
t49 = t61 * t127;
t140 = t144 ^ 2;
t212 = mrSges(7,3) * t140;
t80 = t105 * t222 - t145 * t112;
t60 = qJD(5) * t80 + t145 * t108 + t109 * t222;
t50 = t60 * t212;
t141 = t148 ^ 2;
t211 = mrSges(7,3) * t141;
t51 = t60 * t211;
t54 = t61 * mrSges(6,1);
t171 = mrSges(7,1) * t144 + mrSges(7,2) * t148;
t122 = t171 * qJD(6);
t75 = -pkin(5) - t80;
t73 = t75 * t122;
t235 = t49 + t50 + t51 + t73 - t54;
t234 = t72 * mrSges(4,1) + t38 * mrSges(5,1) - t71 * mrSges(4,2) - t39 * mrSges(5,2) + Ifges(4,5) * t92 + Ifges(5,5) * t69 - Ifges(4,6) * t93 + Ifges(5,6) * t68;
t233 = 2 * m(5);
t232 = 2 * m(6);
t231 = 2 * m(7);
t7 = qJD(5) * t22 + t145 * t19 - t151 * t222;
t230 = 0.2e1 * t7;
t229 = -2 * mrSges(6,3);
t228 = -0.2e1 * t240;
t227 = 0.2e1 * t48;
t226 = 0.2e1 * t137;
t224 = t240 * t7;
t220 = pkin(5) * t122;
t2 = qJD(6) * t16 + t13 * t144 + t148 * t6;
t219 = t148 * t2;
t218 = t240 * t61;
t217 = t3 * t144;
t215 = t60 * mrSges(6,2);
t195 = t148 * t36;
t213 = Ifges(7,5) * t195 + Ifges(7,3) * t37;
t210 = Ifges(7,4) * t144;
t209 = Ifges(7,4) * t148;
t208 = Ifges(7,6) * t144;
t111 = t135 * t222 - t145 * t221;
t102 = t111 * qJD(5);
t204 = t102 * mrSges(6,2);
t103 = t113 * qJD(5);
t203 = t103 * t240;
t202 = t108 * mrSges(5,1);
t201 = t109 * mrSges(5,2);
t198 = t144 * t63;
t194 = t148 * t63;
t76 = pkin(10) + t81;
t193 = t148 * t76;
t188 = t140 + t141;
t185 = qJD(6) * t144;
t183 = 0.2e1 * t150;
t182 = t63 * t185;
t181 = t16 * t184;
t178 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t177 = t37 * mrSges(6,1) + t36 * mrSges(6,2);
t176 = -t185 / 0.2e1;
t175 = -(2 * Ifges(6,4)) - t208;
t174 = t188 * t60;
t173 = t188 * t102;
t107 = pkin(10) + t113;
t172 = t188 * t107;
t170 = Ifges(7,1) * t148 - t210;
t169 = -Ifges(7,2) * t144 + t209;
t168 = Ifges(7,5) * t144 + Ifges(7,6) * t148;
t42 = mrSges(7,2) * t161 - mrSges(7,3) * t198;
t43 = -mrSges(7,1) * t161 - mrSges(7,3) * t194;
t167 = -t144 * t43 + t148 * t42;
t159 = t182 - t195;
t123 = t169 * qJD(6);
t124 = t170 * qJD(6);
t128 = Ifges(7,2) * t148 + t210;
t129 = Ifges(7,1) * t144 + t209;
t158 = t148 * t123 + t144 * t124 - t128 * t185 + t129 * t184;
t157 = (-mrSges(4,1) * t146 - mrSges(4,2) * t149) * t207;
t106 = -pkin(5) - t111;
t90 = t106 * t122;
t91 = t103 * t127;
t96 = t102 * t212;
t97 = t102 * t211;
t98 = t103 * mrSges(6,1);
t155 = t158 + t90 + t91 + t96 + t97 - t98;
t10 = -Ifges(7,4) * t159 - Ifges(7,2) * t160 + Ifges(7,6) * t37;
t11 = -Ifges(7,1) * t159 - Ifges(7,4) * t160 + Ifges(7,5) * t37;
t138 = Ifges(7,5) * t184;
t27 = -Ifges(7,6) * t161 + t169 * t63;
t28 = -Ifges(7,5) * t161 + t170 * t63;
t154 = -t6 * mrSges(6,2) + mrSges(7,3) * t219 - t240 * t122 + t27 * t176 + t28 * t184 / 0.2e1 + Ifges(6,5) * t36 - t123 * t198 / 0.2e1 + t124 * t194 / 0.2e1 - t161 * (-Ifges(7,6) * t185 + t138) / 0.2e1 + t144 * t11 / 0.2e1 + t148 * t10 / 0.2e1 + (-mrSges(6,1) + t127) * t7 + (t168 / 0.2e1 - Ifges(6,6)) * t37 - t160 * t128 / 0.2e1 + (t195 / 0.2e1 + t63 * t176) * t129;
t14 = mrSges(7,1) * t37 + mrSges(7,3) * t159;
t15 = -mrSges(7,2) * t37 - mrSges(7,3) * t160;
t153 = -t43 * t184 - t42 * t185 - t144 * t14 + t148 * t15 + m(7) * (-t17 * t185 - t181 - t217 + t219);
t152 = (-t217 + (-t144 * t17 - t148 * t16) * qJD(6)) * mrSges(7,3) + t154;
t41 = t171 * t63;
t12 = mrSges(7,1) * t160 - mrSges(7,2) * t159;
t1 = [(mrSges(6,3) * t228 - t144 * t27 + t148 * t28) * t36 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t150) * t183 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t117 + mrSges(4,2) * t118) - 0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t226 - 0.2e1 * Ifges(3,4) * t147 + (Ifges(3,1) - Ifges(3,2)) * t183) * t147) * qJD(2) + 0.2e1 * (t68 * t87 + t69 * t86) * Ifges(5,4) + 0.2e1 * (-t38 * t87 + t39 * t86 - t57 * t69 + t58 * t68) * mrSges(5,3) + 0.2e1 * t74 * t177 + 0.2e1 * t104 * t178 + 0.2e1 * t86 * Ifges(5,2) * t68 - 0.2e1 * t117 * Ifges(4,2) * t93 + t22 * t37 * t229 + 0.2e1 * t69 * t87 * Ifges(5,1) + 0.2e1 * t92 * t118 * Ifges(4,1) + 0.2e1 * m(4) * (t241 * t72 + t71 * t95) + 0.2e1 * (t117 * t92 - t118 * t93) * Ifges(4,4) + 0.2e1 * (t117 * t71 - t118 * t72 - t241 * t92 - t93 * t95) * mrSges(4,3) + (t93 * mrSges(4,1) + t92 * mrSges(4,2)) * t226 + (t16 * t3 + t17 * t2 - t224) * t231 + (t22 * t6 + t48 * t74 - t224) * t232 - (mrSges(6,1) * t227 + t6 * t229 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t37 + t175 * t36 + t213) * t161 + (mrSges(6,2) * t227 + mrSges(6,3) * t230 + 0.2e1 * Ifges(6,1) * t36 - t144 * t10 + t148 * t11 + (Ifges(7,5) * t148 + t175) * t37 + (-t144 * t28 - t148 * t27 + t161 * t168) * qJD(6)) * t63 + 0.2e1 * t84 * (-mrSges(5,1) * t86 + mrSges(5,2) * t87) + 0.2e1 * t2 * t42 + 0.2e1 * t3 * t43 + 0.2e1 * t16 * t14 + 0.2e1 * t17 * t15 + t12 * t228 + t41 * t230 + (t104 * t84 + t38 * t57 + t39 * t58) * t233; (t242 * mrSges(7,3) + (-m(7) * t16 - t43) * t60 + (m(7) * t242 - qJD(6) * t42 - t14) * t76) * t144 + t234 + m(6) * (t22 * t60 + t6 * t81 - t7 * t80 - t218) + m(7) * (-t76 * t181 + t2 * t193 + t60 * t196 + t7 * t75 - t218) + m(5) * (t108 * t57 + t109 * t58 + t110 * t38 + t112 * t39) + (m(4) * (t146 * t71 + t149 * t72 + (-t146 * t241 + t149 * t95) * qJD(3)) + (-t146 * t93 - t149 * t92 + (t117 * t149 + t118 * t146) * qJD(3)) * mrSges(4,3)) * pkin(2) + (-t108 * t87 + t109 * t86 - t110 * t69 + t112 * t68) * mrSges(5,3) + t154 + t75 * t12 + t61 * t41 + (t76 * t15 + t60 * t42 + (-mrSges(7,3) * t16 - t43 * t76) * qJD(6)) * t148 + (Ifges(3,5) * t150 - Ifges(3,6) * t147 + (-mrSges(3,1) * t150 + mrSges(3,2) * t147) * pkin(7)) * qJD(2) + (t161 * t60 - t36 * t80 - t37 * t81 + t61 * t63) * mrSges(6,3); 0.2e1 * t202 - 0.2e1 * t201 - 0.2e1 * t215 + 0.2e1 * t49 + 0.2e1 * t50 + 0.2e1 * t51 - 0.2e1 * t54 + 0.2e1 * t73 + 0.2e1 * t157 + (t174 * t76 + t61 * t75) * t231 + (t60 * t81 - t61 * t80) * t232 + (t108 * t110 + t109 * t112) * t233 + t158; m(6) * (-t111 * t7 + t113 * t6 - t203) + (t103 * t63 - t111 * t36 - t113 * t37) * mrSges(6,3) + m(7) * (t106 * t7 - t203) + t234 + t153 * t107 + (m(5) * (t142 * t39 + t143 * t38) + (t142 * t68 - t143 * t69) * mrSges(5,3)) * pkin(3) + t152 + t103 * t41 + t106 * t12 + (m(6) * t22 + m(7) * t237 + mrSges(6,3) * t161 + t167) * t102; t155 + m(6) * (t102 * t81 - t103 * t80 - t111 * t61 + t113 * t60) + (-t102 - t60) * mrSges(6,2) + m(7) * (t103 * t75 + t106 * t61 + t172 * t60 + t173 * t76) + t157 + m(5) * (t108 * t143 + t109 * t142) * pkin(3) + t202 - t201 + t235; -0.2e1 * t204 + 0.2e1 * t90 + 0.2e1 * t91 + 0.2e1 * t96 + 0.2e1 * t97 - 0.2e1 * t98 + (t102 * t172 + t103 * t106) * t231 + (t102 * t113 - t103 * t111) * t232 + t158; t148 * t14 + t144 * t15 + t167 * qJD(6) + m(7) * (qJD(6) * t237 + t144 * t2 + t148 * t3) + m(6) * t48 + m(5) * t84 + t177 + t178; 0; 0; 0; (-m(7) * t7 - t12) * pkin(5) + t153 * pkin(10) + t152; -t220 + m(7) * (-pkin(5) * t61 + pkin(10) * t174) - t215 + t158 + t235; -t220 + m(7) * (-pkin(5) * t103 + pkin(10) * t173) - t204 + t155; 0; t158 - 0.2e1 * t220; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t182 - Ifges(7,6) * t160 + t213; t138 - t171 * t60 + (-mrSges(7,1) * t193 + (mrSges(7,2) * t76 - Ifges(7,6)) * t144) * qJD(6); t138 - t171 * t102 + (t107 * t127 - t208) * qJD(6); -t122; t138 + (pkin(10) * t127 - t208) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
