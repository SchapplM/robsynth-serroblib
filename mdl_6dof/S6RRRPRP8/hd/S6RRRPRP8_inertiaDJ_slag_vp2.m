% Calculate time derivative of joint inertia matrix for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:18
% EndTime: 2019-03-09 17:14:32
% DurationCPUTime: 6.33s
% Computational Cost: add. (3889->516), mult. (9159->711), div. (0->0), fcn. (7262->6), ass. (0->206)
t244 = Ifges(5,4) + Ifges(4,5);
t265 = -Ifges(7,5) - Ifges(6,5);
t264 = -Ifges(7,6) - Ifges(6,6);
t268 = -Ifges(5,2) - Ifges(4,3);
t182 = cos(qJ(3));
t218 = qJD(3) * t182;
t179 = sin(qJ(3));
t220 = qJD(3) * t179;
t267 = Ifges(5,6) * t220 + t218 * t244;
t180 = sin(qJ(2));
t183 = cos(qJ(2));
t221 = qJD(2) * t183;
t211 = t179 * t221;
t189 = t180 * t218 + t211;
t178 = sin(qJ(5));
t181 = cos(qJ(5));
t191 = t178 * t179 + t181 * t182;
t260 = qJD(3) - qJD(5);
t76 = t260 * t191;
t192 = t178 * t182 - t179 * t181;
t77 = t260 * t192;
t266 = t264 * t77 - t265 * t76;
t110 = t192 * t180;
t184 = -pkin(3) - pkin(4);
t143 = -qJ(4) * t178 + t181 * t184;
t108 = t181 * qJD(4) + qJD(5) * t143;
t245 = mrSges(6,2) + mrSges(7,2);
t263 = t245 * t108;
t262 = t245 * t181;
t253 = pkin(8) - pkin(9);
t152 = t253 * t179;
t153 = t253 * t182;
t87 = t178 * t152 + t181 * t153;
t235 = Ifges(5,5) * t179;
t196 = Ifges(5,1) * t182 + t235;
t237 = Ifges(4,4) * t179;
t197 = Ifges(4,1) * t182 - t237;
t261 = (t196 + t197) * qJD(3);
t174 = t179 * qJ(4);
t193 = -t182 * pkin(3) - t174;
t208 = t182 * t221;
t222 = qJD(2) * t180;
t259 = -t189 * Ifges(5,6) - t244 * t208 + t268 * t222;
t258 = 2 * m(4);
t257 = 2 * m(6);
t256 = 2 * m(7);
t255 = -2 * pkin(1);
t254 = 2 * pkin(7);
t249 = pkin(7) * t179;
t247 = qJD(2) / 0.2e1;
t246 = (mrSges(6,1) + mrSges(7,1));
t243 = Ifges(6,3) + Ifges(7,3);
t39 = t180 * t76 - t192 * t221;
t23 = mrSges(7,2) * t222 + mrSges(7,3) * t39;
t24 = mrSges(6,2) * t222 + mrSges(6,3) * t39;
t242 = t23 + t24;
t146 = -pkin(2) * t183 - pkin(8) * t180 - pkin(1);
t163 = t183 * t249;
t177 = t183 * pkin(3);
t65 = pkin(4) * t183 + t163 + t177 + (-pkin(9) * t180 - t146) * t182;
t230 = t179 * t180;
t228 = t182 * t183;
t164 = pkin(7) * t228;
t102 = t179 * t146 + t164;
t88 = -qJ(4) * t183 + t102;
t69 = pkin(9) * t230 + t88;
t20 = t178 * t65 + t181 * t69;
t90 = -mrSges(7,2) * t183 - mrSges(7,3) * t110;
t91 = -mrSges(6,2) * t183 - mrSges(6,3) * t110;
t239 = t90 + t91;
t111 = t191 * t180;
t92 = mrSges(7,1) * t183 - mrSges(7,3) * t111;
t93 = mrSges(6,1) * t183 - mrSges(6,3) * t111;
t238 = t92 + t93;
t236 = Ifges(4,4) * t182;
t234 = Ifges(5,5) * t182;
t233 = Ifges(4,6) * t182;
t232 = t183 * Ifges(4,6);
t229 = t180 * t182;
t105 = -Ifges(5,4) * t183 + t196 * t180;
t106 = -Ifges(4,5) * t183 + t197 * t180;
t227 = t105 + t106;
t142 = (pkin(2) * t180 - pkin(8) * t183) * qJD(2);
t226 = t179 * t142 + t146 * t218;
t217 = qJD(4) * t182;
t225 = qJ(4) * t208 + t180 * t217;
t224 = qJ(4) * t218 + t179 * qJD(4);
t219 = qJD(3) * t180;
t216 = qJD(5) * t178;
t215 = qJD(5) * t181;
t145 = -pkin(2) + t193;
t144 = t181 * qJ(4) + t178 * t184;
t109 = -t178 * qJD(4) - qJD(5) * t144;
t214 = t178 * t108 + t181 * t109 + t144 * t215;
t213 = t184 * t179;
t212 = -pkin(3) - t249;
t210 = t179 * t219;
t148 = -Ifges(5,3) * t182 + t235;
t149 = Ifges(4,2) * t182 + t237;
t207 = t148 / 0.2e1 - t149 / 0.2e1;
t150 = Ifges(5,1) * t179 - t234;
t151 = Ifges(4,1) * t179 + t236;
t206 = -t151 / 0.2e1 - t150 / 0.2e1;
t40 = t110 * t260 + t191 * t221;
t10 = -t39 * mrSges(7,1) + t40 * mrSges(7,2);
t27 = t77 * mrSges(7,1) + t76 * mrSges(7,2);
t19 = -t178 * t69 + t181 * t65;
t205 = -t264 * t39 - t265 * t40;
t86 = t181 * t152 - t153 * t178;
t101 = t146 * t182 - t163;
t121 = t182 * pkin(4) - t145;
t204 = m(7) * pkin(5) + t246;
t203 = -pkin(7) + t213;
t194 = Ifges(5,3) * t179 + t234;
t103 = -Ifges(5,6) * t183 + t194 * t180;
t195 = -Ifges(4,2) * t179 + t236;
t104 = t195 * t180 - t232;
t202 = t103 - t104 + t232;
t201 = qJD(3) * t164 - t142 * t182 + t146 * t220;
t200 = -t182 * mrSges(4,1) + t179 * mrSges(4,2);
t199 = mrSges(4,1) * t179 + mrSges(4,2) * t182;
t147 = -t182 * mrSges(5,1) - t179 * mrSges(5,3);
t198 = mrSges(5,1) * t179 - mrSges(5,3) * t182;
t21 = pkin(9) * t210 + (-pkin(9) * t228 + (-pkin(4) + t212) * t180) * qJD(2) + t201;
t166 = qJ(4) * t222;
t22 = t166 + (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t229 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t179) * t183 + t226;
t3 = t178 * t21 + t181 * t22 + t65 * t215 - t69 * t216;
t98 = qJD(3) * t213 + t224;
t140 = t253 * t220;
t141 = qJD(3) * t153;
t41 = -t181 * t140 + t178 * t141 + t152 * t215 - t153 * t216;
t162 = qJ(4) * t229;
t85 = t203 * t180 + t162;
t190 = t208 - t210;
t188 = -t108 * t191 + t109 * t192 - t144 * t77;
t95 = mrSges(5,2) * t208 + (-mrSges(5,1) * qJD(2) - mrSges(5,2) * t220) * t180;
t4 = -t20 * qJD(5) - t178 * t22 + t181 * t21;
t42 = -qJD(5) * t87 + t140 * t178 + t181 * t141;
t47 = (-t182 * t222 - t183 * t220) * pkin(7) + t226;
t1 = -pkin(5) * t222 - qJ(6) * t40 - qJD(6) * t111 + t4;
t2 = qJ(6) * t39 - qJD(6) * t110 + t3;
t187 = t4 * mrSges(6,1) + t1 * mrSges(7,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2) + t205;
t12 = -qJ(6) * t77 - qJD(6) * t191 + t41;
t13 = -qJ(6) * t76 + qJD(6) * t192 + t42;
t185 = t42 * mrSges(6,1) + t13 * mrSges(7,1) - t41 * mrSges(6,2) - t12 * mrSges(7,2) + t266;
t33 = (t184 * t182 - t174) * t219 + t203 * t221 + t225;
t139 = -pkin(5) + t143;
t138 = -mrSges(5,2) * t230 - mrSges(5,3) * t183;
t137 = mrSges(5,1) * t183 + mrSges(5,2) * t229;
t136 = -mrSges(4,1) * t183 - mrSges(4,3) * t229;
t135 = mrSges(4,2) * t183 - mrSges(4,3) * t230;
t132 = t195 * qJD(3);
t131 = t194 * qJD(3);
t130 = t199 * qJD(3);
t129 = t198 * qJD(3);
t118 = t198 * t180;
t112 = pkin(3) * t220 - t224;
t107 = -t162 + (pkin(3) * t179 + pkin(7)) * t180;
t97 = -mrSges(5,2) * t189 + mrSges(5,3) * t222;
t96 = -mrSges(4,2) * t222 - mrSges(4,3) * t189;
t94 = mrSges(4,1) * t222 - mrSges(4,3) * t190;
t89 = -t101 + t177;
t84 = pkin(5) * t191 + t121;
t83 = -Ifges(6,1) * t192 - Ifges(6,4) * t191;
t82 = -Ifges(7,1) * t192 - Ifges(7,4) * t191;
t81 = -Ifges(6,4) * t192 - Ifges(6,2) * t191;
t80 = -Ifges(7,4) * t192 - Ifges(7,2) * t191;
t79 = mrSges(6,1) * t191 - mrSges(6,2) * t192;
t78 = mrSges(7,1) * t191 - mrSges(7,2) * t192;
t75 = t144 * t108;
t68 = mrSges(4,1) * t189 + mrSges(4,2) * t190;
t67 = mrSges(5,1) * t189 - mrSges(5,3) * t190;
t62 = mrSges(6,1) * t110 + mrSges(6,2) * t111;
t61 = mrSges(7,1) * t110 + mrSges(7,2) * t111;
t60 = -qJ(6) * t191 + t87;
t59 = qJ(6) * t192 + t86;
t58 = -t151 * t219 + (Ifges(4,5) * t180 + t197 * t183) * qJD(2);
t57 = -t150 * t219 + (Ifges(5,4) * t180 + t196 * t183) * qJD(2);
t56 = -t149 * t219 + (Ifges(4,6) * t180 + t195 * t183) * qJD(2);
t55 = -t148 * t219 + (Ifges(5,6) * t180 + t194 * t183) * qJD(2);
t53 = Ifges(6,1) * t111 - Ifges(6,4) * t110 + Ifges(6,5) * t183;
t52 = Ifges(7,1) * t111 - Ifges(7,4) * t110 + Ifges(7,5) * t183;
t51 = Ifges(6,4) * t111 - Ifges(6,2) * t110 + Ifges(6,6) * t183;
t50 = Ifges(7,4) * t111 - Ifges(7,2) * t110 + Ifges(7,6) * t183;
t49 = pkin(5) * t110 + t85;
t48 = t222 * t249 - t201;
t46 = pkin(3) * t189 + pkin(7) * t221 + qJ(4) * t210 - t225;
t45 = t212 * t222 + t201;
t44 = pkin(5) * t77 + t98;
t43 = -qJD(4) * t183 + t166 + t47;
t32 = Ifges(6,1) * t76 - Ifges(6,4) * t77;
t31 = Ifges(7,1) * t76 - Ifges(7,4) * t77;
t30 = Ifges(6,4) * t76 - Ifges(6,2) * t77;
t29 = Ifges(7,4) * t76 - Ifges(7,2) * t77;
t28 = mrSges(6,1) * t77 + mrSges(6,2) * t76;
t26 = -mrSges(6,1) * t222 - mrSges(6,3) * t40;
t25 = -mrSges(7,1) * t222 - mrSges(7,3) * t40;
t15 = -qJ(6) * t110 + t20;
t14 = pkin(5) * t183 - qJ(6) * t111 + t19;
t11 = -mrSges(6,1) * t39 + mrSges(6,2) * t40;
t9 = -pkin(5) * t39 + t33;
t8 = Ifges(6,1) * t40 + Ifges(6,4) * t39 - Ifges(6,5) * t222;
t7 = Ifges(7,1) * t40 + Ifges(7,4) * t39 - Ifges(7,5) * t222;
t6 = Ifges(6,4) * t40 + Ifges(6,2) * t39 - Ifges(6,6) * t222;
t5 = Ifges(7,4) * t40 + Ifges(7,2) * t39 - Ifges(7,6) * t222;
t16 = [(((mrSges(3,2) * t255) + 0.2e1 * Ifges(3,4) * t183 + t202 * t179 + t227 * t182) * qJD(2) + t205 + t259) * t183 + (t7 + t8) * t111 - (t5 + t6) * t110 + (t68 * t254 + (t57 + t58) * t182 + (t55 - t56) * t179 + (t202 * t182 + (t244 * t183 - t227) * t179) * qJD(3) + ((mrSges(3,1) * t255) + (-0.2e1 * Ifges(3,4) + t244 * t182 + (-Ifges(4,6) + Ifges(5,6)) * t179) * t180 + t265 * t111 - t264 * t110 + ((pkin(7) ^ 2 * t258) + t199 * t254 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(6,3)) - (2 * Ifges(7,3)) + t268) * t183) * qJD(2)) * t180 + 0.2e1 * m(5) * (t107 * t46 + t43 * t88 + t45 * t89) + (t1 * t14 + t15 * t2 + t49 * t9) * t256 + (t19 * t4 + t20 * t3 + t33 * t85) * t257 + (t101 * t48 + t102 * t47) * t258 + (t52 + t53) * t40 + 0.2e1 * t15 * t23 + 0.2e1 * t20 * t24 + 0.2e1 * t14 * t25 + 0.2e1 * t19 * t26 + 0.2e1 * t49 * t10 + 0.2e1 * t9 * t61 + 0.2e1 * t33 * t62 + 0.2e1 * t85 * t11 + 0.2e1 * t2 * t90 + 0.2e1 * t3 * t91 + 0.2e1 * t1 * t92 + 0.2e1 * t4 * t93 + 0.2e1 * t89 * t95 + 0.2e1 * t88 * t97 + 0.2e1 * t101 * t94 + 0.2e1 * t102 * t96 + 0.2e1 * t107 * t67 + 0.2e1 * t46 * t118 + 0.2e1 * t47 * t135 + 0.2e1 * t48 * t136 + 0.2e1 * t45 * t137 + 0.2e1 * t43 * t138 + (t50 + t51) * t39; m(5) * (t107 * t112 + t145 * t46) + (t82 / 0.2e1 + t83 / 0.2e1) * t40 + (t80 / 0.2e1 + t81 / 0.2e1) * t39 + (t1 * t192 - t14 * t76 - t15 * t77 - t191 * t2) * mrSges(7,3) + (-t19 * t76 - t191 * t3 + t192 * t4 - t20 * t77) * mrSges(6,3) - (t5 / 0.2e1 + t6 / 0.2e1) * t191 + (-Ifges(3,6) * qJD(2) + t233 * t247 + (mrSges(3,2) * qJD(2) + t130) * pkin(7) - (t191 * t264 + t192 * t265) * qJD(2) / 0.2e1 + (t131 / 0.2e1 - t132 / 0.2e1 + t206 * qJD(3) + t244 * t247) * t179 + (t261 / 0.2e1 - Ifges(5,6) * t247 + t207 * qJD(3)) * t182) * t180 - (t7 / 0.2e1 + t8 / 0.2e1) * t192 + (t31 / 0.2e1 + t32 / 0.2e1) * t111 - (t29 / 0.2e1 + t30 / 0.2e1) * t110 + (-t55 / 0.2e1 + t56 / 0.2e1 + t43 * mrSges(5,2) + t47 * mrSges(4,3) + (t89 * mrSges(5,2) - t101 * mrSges(4,3) + t105 / 0.2e1 + t106 / 0.2e1) * qJD(3)) * t182 - (t50 / 0.2e1 + t51 / 0.2e1) * t77 + (t52 / 0.2e1 + t53 / 0.2e1) * t76 + (t266 / 0.2e1 - t267 / 0.2e1) * t183 + (Ifges(3,5) - t206 * t182 + t207 * t179 + (-m(4) * pkin(2) - mrSges(3,1) + t200) * pkin(7)) * t221 + ((t96 + t97) * t182 + (-t94 + t95) * t179 + ((-t136 + t137) * t182 + (-t135 - t138) * t179) * qJD(3) + m(5) * (t45 * t179 + t43 * t182 + t89 * t218 - t88 * t220) + m(4) * (-t101 * t218 - t102 * t220 - t48 * t179 + t47 * t182)) * pkin(8) + (t57 / 0.2e1 + t58 / 0.2e1 + t45 * mrSges(5,2) - t48 * mrSges(4,3) + (t232 / 0.2e1 + t103 / 0.2e1 - t104 / 0.2e1 - t88 * mrSges(5,2) - t102 * mrSges(4,3)) * qJD(3)) * t179 + m(7) * (t1 * t59 + t12 * t15 + t13 * t14 + t2 * t60 + t44 * t49 + t84 * t9) + m(6) * (t121 * t33 + t19 * t42 + t20 * t41 + t3 * t87 + t4 * t86 + t85 * t98) + t49 * t27 + t59 * t25 + t60 * t23 + t44 * t61 - pkin(2) * t68 + t9 * t78 + t33 * t79 + t84 * t10 + t85 * t28 + t86 * t26 + t87 * t24 + t12 * t90 + t41 * t91 + t13 * t92 + t42 * t93 + t98 * t62 + t112 * t118 + t121 * t11 + t107 * t129 + t145 * t67 + t46 * t147; -0.2e1 * pkin(2) * t130 + 0.2e1 * t121 * t28 + 0.2e1 * t145 * t129 + 0.2e1 * t84 * t27 + 0.2e1 * t44 * t78 + 0.2e1 * t98 * t79 - (t80 + t81) * t77 + (t82 + t83) * t76 + (-t131 + t132) * t182 + t261 * t179 - (t31 + t32) * t192 - (t29 + t30) * t191 + (t12 * t60 + t13 * t59 + t44 * t84) * t256 + (t121 * t98 + t41 * t87 + t42 * t86) * t257 + ((t150 + t151) * t182 + (t148 - t149) * t179) * qJD(3) + 0.2e1 * (m(5) * t145 + t147) * t112 + 0.2e1 * (-t12 * t191 + t13 * t192 - t59 * t76 - t60 * t77) * mrSges(7,3) + 0.2e1 * (-t191 * t41 + t192 * t42 - t76 * t86 - t77 * t87) * mrSges(6,3); -t259 - t187 + m(5) * (-pkin(3) * t45 + qJ(4) * t43 + qJD(4) * t88) + m(7) * (t1 * t139 + t108 * t15 + t109 * t14 + t144 * t2) + m(6) * (t108 * t20 + t109 * t19 + t143 * t4 + t144 * t3) - Ifges(4,6) * t211 + t238 * t109 + t239 * t108 + t242 * t144 + (t243 * qJD(2) + (-t244 * t179 - t233) * qJD(3)) * t180 + t43 * mrSges(5,3) - t45 * mrSges(5,1) - t47 * mrSges(4,2) + t48 * mrSges(4,1) - pkin(3) * t95 + qJ(4) * t97 + qJD(4) * t138 + t139 * t25 + t143 * t26; -t185 + (m(5) * t217 + (m(5) * t193 + t147 + t200) * qJD(3)) * pkin(8) + (qJD(3) * t193 + t217) * mrSges(5,2) + m(7) * (t108 * t60 + t109 * t59 + t12 * t144 + t13 * t139) + m(6) * (t108 * t87 + t109 * t86 + t143 * t42 + t144 * t41) + (-t139 * t76 + t188) * mrSges(7,3) + (-t143 * t76 + t188) * mrSges(6,3) - Ifges(4,6) * t220 + t267; (t109 * t143 + t75) * t257 + (t109 * t139 + t75) * t256 - 0.2e1 * t246 * t109 + 0.2e1 * t263 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4); (t25 + t26) * t181 + t242 * t178 + (-t238 * t178 + t239 * t181) * qJD(5) + m(7) * (t1 * t181 + t178 * t2 + (-t14 * t178 + t15 * t181) * qJD(5)) + m(6) * (t178 * t3 + t181 * t4 + (-t178 * t19 + t181 * t20) * qJD(5)) + m(5) * t45 + t95; (m(5) * pkin(8) + mrSges(5,2)) * t218 + m(7) * (t12 * t178 + t13 * t181 + (-t178 * t59 + t181 * t60) * qJD(5)) + m(6) * (t178 * t41 + t181 * t42 + (-t178 * t86 + t181 * t87) * qJD(5)) + (mrSges(7,3) + mrSges(6,3)) * (-t178 * t77 - t181 * t76 + (-t178 * t192 - t181 * t191) * qJD(5)); (t246 * t178 + t262) * qJD(5) + m(6) * (-t143 * t216 + t214) + m(7) * (-t139 * t216 + t214); 0; -t243 * t222 + (m(7) * t1 + t25) * pkin(5) + t187; (m(7) * t13 - t76 * mrSges(7,3)) * pkin(5) + t185; t204 * t109 - t263; (-t204 * t178 - t262) * qJD(5); 0; m(7) * t9 + t10; m(7) * t44 + t27; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
