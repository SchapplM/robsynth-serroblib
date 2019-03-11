% Calculate time derivative of joint inertia matrix for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:19
% EndTime: 2019-03-09 06:30:27
% DurationCPUTime: 3.61s
% Computational Cost: add. (3794->462), mult. (8601->659), div. (0->0), fcn. (6987->6), ass. (0->185)
t222 = 2 * qJD(2);
t235 = -mrSges(6,1) - mrSges(7,1);
t234 = -mrSges(6,2) + mrSges(7,3);
t233 = mrSges(6,3) + mrSges(7,2);
t232 = Ifges(7,4) + Ifges(6,5);
t231 = Ifges(7,2) + Ifges(6,3);
t230 = Ifges(7,6) - Ifges(6,6);
t146 = sin(qJ(5));
t147 = sin(qJ(4));
t149 = cos(qJ(5));
t150 = cos(qJ(4));
t109 = t146 * t150 + t147 * t149;
t226 = qJD(4) + qJD(5);
t229 = t226 * t109;
t148 = sin(qJ(3));
t200 = qJD(3) * t148;
t177 = t147 * t200;
t151 = cos(qJ(3));
t195 = qJD(4) * t151;
t179 = t150 * t195;
t159 = t177 - t179;
t199 = qJD(3) * t151;
t228 = Ifges(5,6) * t177 + Ifges(5,3) * t199;
t218 = pkin(8) * t151;
t219 = pkin(3) * t148;
t117 = qJ(2) - t218 + t219;
t152 = -pkin(1) - pkin(7);
t203 = t148 * t152;
t130 = t150 * t203;
t92 = t147 * t117 + t130;
t201 = t147 ^ 2 + t150 ^ 2;
t180 = t147 * t195;
t160 = t150 * t200 + t180;
t197 = qJD(4) * t147;
t181 = t148 * t197;
t227 = t150 * t199 - t181;
t105 = t150 * t117;
t173 = -t147 * t152 + pkin(4);
t202 = t150 * t151;
t67 = -pkin(9) * t202 + t173 * t148 + t105;
t204 = t147 * t151;
t75 = -pkin(9) * t204 + t92;
t215 = t146 * t67 + t149 * t75;
t106 = qJD(2) + (pkin(3) * t151 + pkin(8) * t148) * qJD(3);
t100 = t150 * t106;
t22 = t100 + (-t130 + (pkin(9) * t151 - t117) * t147) * qJD(4) + (pkin(9) * t148 * t150 + t173 * t151) * qJD(3);
t198 = qJD(3) * t152;
t182 = t151 * t198;
t196 = qJD(4) * t150;
t53 = t147 * t106 + t117 * t196 + t150 * t182 - t152 * t181;
t30 = t159 * pkin(9) + t53;
t6 = -qJD(5) * t215 - t146 * t30 + t149 * t22;
t225 = 2 * m(6);
t224 = 2 * m(7);
t223 = 0.2e1 * pkin(4);
t221 = -pkin(9) - pkin(8);
t220 = -t147 / 0.2e1;
t171 = t226 * t150;
t194 = qJD(5) * t146;
t178 = t147 * t194;
t47 = -t151 * t178 + (t151 * t171 - t177) * t149 - t160 * t146;
t26 = -t47 * mrSges(7,2) + mrSges(7,3) * t199;
t29 = -mrSges(6,2) * t199 - t47 * mrSges(6,3);
t217 = t26 + t29;
t108 = t146 * t147 - t149 * t150;
t45 = t108 * t200 - t151 * t229;
t27 = mrSges(6,1) * t199 - t45 * mrSges(6,3);
t28 = -mrSges(7,1) * t199 + t45 * mrSges(7,2);
t216 = -t27 + t28;
t96 = t109 * t151;
t84 = -mrSges(7,2) * t96 + mrSges(7,3) * t148;
t85 = -mrSges(6,2) * t148 - mrSges(6,3) * t96;
t214 = t84 + t85;
t98 = t108 * t151;
t86 = mrSges(6,1) * t148 + mrSges(6,3) * t98;
t87 = -mrSges(7,1) * t148 - mrSges(7,2) * t98;
t213 = -t86 + t87;
t212 = Ifges(5,4) * t147;
t211 = Ifges(5,4) * t150;
t210 = Ifges(5,5) * t147;
t209 = Ifges(5,6) * t147;
t208 = Ifges(5,6) * t150;
t54 = -t92 * qJD(4) - t147 * t182 + t100;
t207 = t147 * t54;
t206 = t148 * Ifges(5,6);
t125 = -mrSges(5,1) * t150 + mrSges(5,2) * t147;
t205 = -mrSges(4,1) + t125;
t193 = qJD(5) * t149;
t192 = pkin(4) * t197;
t191 = pkin(4) * t194;
t190 = pkin(4) * t193;
t189 = t221 * t147;
t128 = t221 * t150;
t82 = -t146 * t128 - t149 * t189;
t188 = t82 * t194;
t95 = t109 * t148;
t187 = t95 * t194;
t138 = -pkin(4) * t150 - pkin(3);
t186 = qJD(4) * t221;
t185 = t147 * t199;
t184 = t148 * t199;
t133 = t148 * t198;
t116 = t147 * t186;
t170 = t150 * t186;
t50 = -t82 * qJD(5) + t149 * t116 + t146 * t170;
t83 = -t149 * t128 + t146 * t189;
t51 = t83 * qJD(5) + t146 * t116 - t149 * t170;
t175 = t83 * t50 + t51 * t82;
t174 = -Ifges(5,5) * t150 + (2 * Ifges(4,4));
t91 = -t147 * t203 + t105;
t172 = -t91 * qJD(4) + t53;
t107 = pkin(4) * t204 - t151 * t152;
t169 = mrSges(5,1) * t147 + mrSges(5,2) * t150;
t168 = Ifges(5,1) * t150 - t212;
t127 = Ifges(5,1) * t147 + t211;
t167 = -Ifges(5,2) * t147 + t211;
t126 = Ifges(5,2) * t150 + t212;
t24 = -t146 * t75 + t149 * t67;
t164 = t199 * t231 + t230 * t47 + t232 * t45;
t5 = t146 * t22 + t149 * t30 + t67 * t193 - t75 * t194;
t44 = -qJD(3) * t98 - t148 * t229;
t46 = -t148 * t178 + (t148 * t171 + t185) * t149 + t227 * t146;
t163 = t234 * t44 + t235 * t46;
t97 = t108 * t148;
t162 = t83 * t44 + t46 * t82 - t50 * t97 + t51 * t95;
t88 = -pkin(4) * t159 + t133;
t69 = Ifges(7,6) * t229;
t70 = Ifges(6,6) * t229;
t73 = t226 * t108;
t71 = Ifges(6,5) * t73;
t72 = Ifges(7,4) * t73;
t157 = t234 * t50 + t235 * t51 + t69 - t70 - t71 - t72;
t2 = qJ(6) * t199 + qJD(6) * t148 + t5;
t3 = -pkin(5) * t199 - t6;
t156 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) + t164;
t132 = qJD(6) + t190;
t153 = -mrSges(6,2) * t190 + t132 * mrSges(7,3) + t191 * t235;
t145 = qJD(6) * mrSges(7,3);
t142 = Ifges(5,5) * t196;
t137 = -pkin(4) * t149 - pkin(5);
t135 = pkin(4) * t146 + qJ(6);
t115 = mrSges(5,1) * t148 - mrSges(5,3) * t202;
t114 = -mrSges(5,2) * t148 - mrSges(5,3) * t204;
t113 = t168 * qJD(4);
t112 = t167 * qJD(4);
t111 = t169 * qJD(4);
t103 = t169 * t151;
t94 = Ifges(5,5) * t148 + t168 * t151;
t93 = t167 * t151 + t206;
t90 = -mrSges(5,2) * t199 + t159 * mrSges(5,3);
t89 = mrSges(5,1) * t199 + t160 * mrSges(5,3);
t81 = Ifges(6,1) * t109 - Ifges(6,4) * t108;
t80 = Ifges(7,1) * t109 + Ifges(7,5) * t108;
t79 = Ifges(6,4) * t109 - Ifges(6,2) * t108;
t78 = Ifges(7,5) * t109 + Ifges(7,3) * t108;
t77 = mrSges(6,1) * t108 + mrSges(6,2) * t109;
t76 = mrSges(7,1) * t108 - mrSges(7,3) * t109;
t66 = -t159 * mrSges(5,1) - t160 * mrSges(5,2);
t65 = pkin(5) * t108 - qJ(6) * t109 + t138;
t62 = mrSges(6,1) * t96 - mrSges(6,2) * t98;
t61 = mrSges(7,1) * t96 + mrSges(7,3) * t98;
t60 = -t127 * t195 + (Ifges(5,5) * t151 - t168 * t148) * qJD(3);
t59 = -t126 * t195 + (Ifges(5,6) * t151 - t167 * t148) * qJD(3);
t58 = -Ifges(6,1) * t98 - Ifges(6,4) * t96 + Ifges(6,5) * t148;
t57 = -Ifges(7,1) * t98 + Ifges(7,4) * t148 + Ifges(7,5) * t96;
t56 = -Ifges(6,4) * t98 - Ifges(6,2) * t96 + Ifges(6,6) * t148;
t55 = -Ifges(7,5) * t98 + Ifges(7,6) * t148 + Ifges(7,3) * t96;
t52 = t96 * pkin(5) + t98 * qJ(6) + t107;
t36 = -Ifges(6,1) * t73 - Ifges(6,4) * t229;
t35 = -Ifges(7,1) * t73 + Ifges(7,5) * t229;
t34 = -Ifges(6,4) * t73 - Ifges(6,2) * t229;
t33 = -Ifges(7,5) * t73 + Ifges(7,3) * t229;
t32 = mrSges(6,1) * t229 - mrSges(6,2) * t73;
t31 = mrSges(7,1) * t229 + mrSges(7,3) * t73;
t21 = -pkin(5) * t148 - t24;
t20 = qJ(6) * t148 + t215;
t16 = pkin(5) * t229 + qJ(6) * t73 - qJD(6) * t109 + t192;
t13 = mrSges(6,1) * t47 + mrSges(6,2) * t45;
t12 = mrSges(7,1) * t47 - mrSges(7,3) * t45;
t11 = Ifges(6,1) * t45 - Ifges(6,4) * t47 + Ifges(6,5) * t199;
t10 = Ifges(7,1) * t45 + Ifges(7,4) * t199 + Ifges(7,5) * t47;
t9 = Ifges(6,4) * t45 - Ifges(6,2) * t47 + Ifges(6,6) * t199;
t8 = Ifges(7,5) * t45 + Ifges(7,6) * t199 + Ifges(7,3) * t47;
t7 = t47 * pkin(5) - t45 * qJ(6) + qJD(6) * t98 + t88;
t1 = [(t107 * t88 + t215 * t5 + t24 * t6) * t225 + 0.2e1 * t215 * t29 + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t222 + (t55 - t56) * t47 + (t57 + t58) * t45 + ((mrSges(4,2) * t222) - t147 * t59 + t150 * t60 - 0.2e1 * t152 * t66 + (-t150 * t93 - t147 * t94 + t148 * (-t208 - t210)) * qJD(4) + (0.2e1 * qJ(2) * mrSges(4,1) + (-t174 - t209) * t151 - t232 * t98 + t230 * t96 + (-0.2e1 * m(5) * t152 ^ 2 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + Ifges(5,3) + t231) * t148) * qJD(3)) * t151 + (mrSges(4,1) * t222 + (-0.2e1 * qJ(2) * mrSges(4,2) + 0.2e1 * t152 * t103 + t147 * t93 + t174 * t148 - t150 * t94) * qJD(3) + t164 + t228) * t148 + (t2 * t20 + t21 * t3 + t52 * t7) * t224 + 0.2e1 * t20 * t26 + 0.2e1 * t24 * t27 + 0.2e1 * t21 * t28 + 0.2e1 * t52 * t12 + 0.2e1 * t7 * t61 + 0.2e1 * t2 * t84 + 0.2e1 * t5 * t85 + 0.2e1 * t6 * t86 + 0.2e1 * t3 * t87 + 0.2e1 * t88 * t62 + 0.2e1 * t91 * t89 + 0.2e1 * t92 * t90 + 0.2e1 * t107 * t13 + 0.2e1 * t53 * t114 + 0.2e1 * t54 * t115 + (t8 - t9) * t96 - (t10 + t11) * t98 + 0.2e1 * m(5) * (t92 * t53 + t91 * t54); -t217 * t97 + t216 * t95 + t213 * t46 + t214 * t44 + (-t12 - t13 - t66 + (t114 * t150 - t115 * t147) * qJD(3)) * t151 + (-t147 * t89 + t150 * t90 + (-t114 * t147 - t115 * t150) * qJD(4) + (t103 + t61 + t62) * qJD(3)) * t148 + m(7) * (-t151 * t7 - t2 * t97 + t44 * t20 + t52 * t200 + t46 * t21 + t3 * t95) + m(5) * ((-t147 * t91 + t150 * t92) * t199 + (-0.2e1 * t182 - t207 + t150 * t53 + (-t147 * t92 - t150 * t91) * qJD(4)) * t148) + m(6) * (t107 * t200 - t151 * t88 + t215 * t44 - t46 * t24 - t5 * t97 - t6 * t95); 0.2e1 * m(5) * (-0.1e1 + t201) * t184 + 0.2e1 * (m(6) + m(7)) * (-t97 * t44 + t46 * t95 - t184); (t55 / 0.2e1 - t56 / 0.2e1) * t229 + (-t108 * t5 - t109 * t6 - t215 * t229 + t24 * t73) * mrSges(6,3) + (-t108 * t2 + t109 * t3 - t20 * t229 - t21 * t73) * mrSges(7,2) + (t142 / 0.2e1 - t71 / 0.2e1 - t70 / 0.2e1 - t72 / 0.2e1 + t69 / 0.2e1 + (t205 * t152 - Ifges(4,5)) * qJD(3)) * t148 + m(6) * (t107 * t192 + t138 * t88 + t215 * t50 - t51 * t24 + t5 * t83 - t6 * t82) + m(5) * (-pkin(3) * t133 + (-t197 * t92 - t207) * pkin(8)) + (t78 / 0.2e1 - t79 / 0.2e1) * t47 - (t57 / 0.2e1 + t58 / 0.2e1) * t73 + (t10 / 0.2e1 + t11 / 0.2e1) * t109 + (t8 / 0.2e1 - t9 / 0.2e1) * t108 + (t80 / 0.2e1 + t81 / 0.2e1) * t45 + m(7) * (t16 * t52 + t2 * t83 + t50 * t20 + t51 * t21 + t3 * t82 + t65 * t7) - (t35 / 0.2e1 + t36 / 0.2e1) * t98 + (t112 * t220 + t150 * t113 / 0.2e1 - t152 * t111 + (-t150 * t126 / 0.2e1 + t127 * t220) * qJD(4) + (-t152 * mrSges(4,2) - Ifges(4,6) + t210 / 0.2e1 + t208 / 0.2e1 + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t109 + (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t108) * qJD(3)) * t151 + t216 * t82 + t217 * t83 + (t33 / 0.2e1 - t34 / 0.2e1) * t96 + t213 * t51 + t214 * t50 + (t126 * t200 / 0.2e1 - t54 * mrSges(5,3) - pkin(8) * t89 + t60 / 0.2e1 + (-pkin(8) * t114 - t92 * mrSges(5,3) + pkin(4) * t62 - t93 / 0.2e1 - t206 / 0.2e1) * qJD(4)) * t147 + (-t127 * t200 / 0.2e1 + qJD(4) * t94 / 0.2e1 + t59 / 0.2e1 + t172 * mrSges(5,3) + (m(5) * t172 - qJD(4) * t115 + t90) * pkin(8)) * t150 + t52 * t31 + t16 * t61 + t65 * t12 - pkin(3) * t66 + t7 * t76 + t88 * t77 + t107 * t32 + t138 * t13; (-t111 - t31 - t32) * t151 + m(6) * (-pkin(4) * t180 + t162) + m(7) * (-t151 * t16 + t162) + ((t201 * mrSges(5,3) - mrSges(4,2)) * t151 + m(5) * (t201 * t218 - t219) + (m(6) * t138 + m(7) * t65 + t205 + t76 + t77) * t148) * qJD(3) + t233 * (-t108 * t44 + t109 * t46 + t229 * t97 - t73 * t95); -0.2e1 * pkin(3) * t111 + t150 * t112 + t147 * t113 + 0.2e1 * t138 * t32 + 0.2e1 * t16 * t76 + 0.2e1 * t65 * t31 + (t78 - t79) * t229 - (t80 + t81) * t73 + (t35 + t36) * t109 + (t33 - t34) * t108 + (t150 * t127 + (t77 * t223 - t126) * t147) * qJD(4) + (t65 * t16 + t175) * t224 + (t138 * t192 + t175) * t225 + 0.2e1 * t233 * (-t108 * t50 + t109 * t51 - t229 * t83 - t73 * t82); m(7) * (t132 * t20 + t135 * t2 + t137 * t3) + ((t27 + m(6) * t6 + (m(6) * t215 + t85) * qJD(5)) * t149 + (t29 + m(6) * t5 + (-m(6) * t24 + m(7) * t21 + t213) * qJD(5)) * t146) * pkin(4) - Ifges(5,6) * t179 + t156 - t53 * mrSges(5,2) + t54 * mrSges(5,1) + t132 * t84 + t135 * t26 + t137 * t28 - t160 * Ifges(5,5) + t228; -t227 * mrSges(5,2) + (-t148 * t196 - t185) * mrSges(5,1) + m(7) * (-t132 * t97 + t135 * t44 + t137 * t46) + (m(6) * (t146 * t44 - t149 * t46 - t97 * t193 + t187) / 0.2e1 + m(7) * t187 / 0.2e1) * t223 + t163; m(7) * (t132 * t83 + t135 * t50 + t137 * t51) + t142 + (t125 * pkin(8) - t209) * qJD(4) + (-t108 * t132 - t135 * t229 - t137 * t73) * mrSges(7,2) + (t109 * mrSges(7,2) * t194 + m(7) * t188 + m(6) * (t146 * t50 - t149 * t51 + t83 * t193 + t188) + (-t146 * t229 + t149 * t73 + (-t108 * t149 + t109 * t146) * qJD(5)) * mrSges(6,3)) * pkin(4) + t157; 0.2e1 * m(7) * (t132 * t135 + t137 * t191) + 0.2e1 * t153; m(7) * (-pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t20) + t156 + qJ(6) * t26 - pkin(5) * t28 + qJD(6) * t84; m(7) * (-pkin(5) * t46 + qJ(6) * t44 - qJD(6) * t97) + t163; m(7) * (-pkin(5) * t51 + qJ(6) * t50 + qJD(6) * t83) + (pkin(5) * t73 - qJ(6) * t229 - qJD(6) * t108) * mrSges(7,2) + t157; m(7) * (-pkin(5) * t191 + qJ(6) * t132 + qJD(6) * t135) + t145 + t153; 0.2e1 * m(7) * qJ(6) * qJD(6) + 0.2e1 * t145; m(7) * t3 + t28; m(7) * t46; m(7) * t51 - t73 * mrSges(7,2); m(7) * t191; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
