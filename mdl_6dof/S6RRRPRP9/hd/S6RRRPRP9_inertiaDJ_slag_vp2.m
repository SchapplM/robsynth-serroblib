% Calculate time derivative of joint inertia matrix for
% S6RRRPRP9
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:55
% EndTime: 2019-03-09 17:22:08
% DurationCPUTime: 6.37s
% Computational Cost: add. (3819->528), mult. (8873->724), div. (0->0), fcn. (6874->6), ass. (0->209)
t250 = Ifges(5,4) + Ifges(4,5);
t273 = Ifges(7,4) + Ifges(6,5);
t272 = Ifges(7,6) - Ifges(6,6);
t276 = -Ifges(5,2) - Ifges(4,3);
t177 = sin(qJ(2));
t176 = sin(qJ(3));
t180 = cos(qJ(2));
t221 = qJD(2) * t180;
t209 = t176 * t221;
t179 = cos(qJ(3));
t218 = qJD(3) * t179;
t186 = t177 * t218 + t209;
t220 = qJD(3) * t176;
t275 = Ifges(5,6) * t220 + t250 * t218;
t175 = sin(qJ(5));
t178 = cos(qJ(5));
t114 = t175 * t176 + t178 * t179;
t269 = qJD(3) - qJD(5);
t72 = t269 * t114;
t215 = qJD(5) * t178;
t216 = qJD(5) * t175;
t73 = t175 * t218 + t176 * t215 - t178 * t220 - t179 * t216;
t203 = t272 * t73 + t273 * t72;
t252 = mrSges(6,1) + mrSges(7,1);
t274 = mrSges(6,3) + mrSges(7,2);
t271 = m(7) * qJD(6);
t239 = Ifges(5,5) * t176;
t194 = Ifges(5,1) * t179 + t239;
t241 = Ifges(4,4) * t176;
t195 = Ifges(4,1) * t179 - t241;
t270 = (t194 + t195) * qJD(3);
t170 = t176 * qJ(4);
t191 = -t179 * pkin(3) - t170;
t181 = -pkin(3) - pkin(4);
t132 = t178 * qJ(4) + t175 * t181;
t208 = t179 * t221;
t219 = qJD(3) * t177;
t211 = t176 * t219;
t187 = t208 - t211;
t268 = m(7) * qJ(6) + mrSges(7,3);
t130 = (pkin(2) * t177 - pkin(8) * t180) * qJD(2);
t134 = -pkin(2) * t180 - pkin(8) * t177 - pkin(1);
t227 = t179 * t180;
t158 = pkin(7) * t227;
t199 = qJD(3) * t158 - t130 * t179 + t134 * t220;
t255 = pkin(7) * t176;
t212 = -pkin(3) - t255;
t20 = pkin(9) * t211 + (-pkin(9) * t227 + (-pkin(4) + t212) * t177) * qJD(2) + t199;
t222 = qJD(2) * t177;
t160 = qJ(4) * t222;
t226 = t176 * t130 + t134 * t218;
t228 = t177 * t179;
t21 = t160 + (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t228 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t176) * t180 + t226;
t157 = t180 * t255;
t174 = t180 * pkin(3);
t62 = pkin(4) * t180 + t157 + t174 + (-pkin(9) * t177 - t134) * t179;
t230 = t176 * t177;
t95 = t176 * t134 + t158;
t83 = -qJ(4) * t180 + t95;
t66 = pkin(9) * t230 + t83;
t246 = t175 * t62 + t178 * t66;
t4 = -qJD(5) * t246 - t175 * t21 + t178 * t20;
t267 = -m(7) * pkin(5) - t252;
t266 = -mrSges(6,2) + t268;
t265 = -t186 * Ifges(5,6) - t250 * t208 + t276 * t222;
t264 = 2 * m(4);
t263 = 2 * m(6);
t262 = 2 * m(7);
t261 = -0.2e1 * pkin(1);
t260 = 0.2e1 * pkin(7);
t259 = pkin(8) - pkin(9);
t253 = qJD(2) / 0.2e1;
t251 = mrSges(6,2) - mrSges(7,3);
t249 = Ifges(7,2) + Ifges(6,3);
t105 = t114 * t177;
t40 = qJD(5) * t105 + t175 * t187 - t178 * t186;
t22 = -mrSges(7,2) * t40 - mrSges(7,3) * t222;
t23 = mrSges(6,2) * t222 - mrSges(6,3) * t40;
t248 = t22 + t23;
t229 = t176 * t178;
t188 = t175 * t179 - t229;
t41 = t177 * t188 * t269 + t114 * t221;
t24 = -mrSges(6,1) * t222 - mrSges(6,3) * t41;
t25 = mrSges(7,1) * t222 + t41 * mrSges(7,2);
t247 = t25 - t24;
t104 = t175 * t228 - t177 * t229;
t85 = -mrSges(7,2) * t104 + mrSges(7,3) * t180;
t86 = -mrSges(6,2) * t180 - mrSges(6,3) * t104;
t243 = t85 + t86;
t87 = mrSges(6,1) * t180 - mrSges(6,3) * t105;
t88 = -mrSges(7,1) * t180 + mrSges(7,2) * t105;
t242 = -t87 + t88;
t240 = Ifges(4,4) * t179;
t238 = Ifges(5,5) * t179;
t237 = Ifges(4,6) * t179;
t102 = -qJ(4) * t216 + t178 * qJD(4) + t181 * t215;
t236 = t102 * mrSges(6,2);
t103 = t175 * qJD(4) + qJD(5) * t132;
t143 = t259 * t179;
t213 = t259 * t176;
t81 = t175 * t143 - t178 * t213;
t235 = t103 * t81;
t234 = t180 * Ifges(4,6);
t100 = -Ifges(4,5) * t180 + t177 * t195;
t99 = -Ifges(5,4) * t180 + t177 * t194;
t233 = t100 + t99;
t232 = t103 * t188;
t231 = t103 * t178;
t217 = qJD(4) * t179;
t225 = qJ(4) * t208 + t177 * t217;
t224 = qJ(4) * t218 + t176 * qJD(4);
t133 = -pkin(2) + t191;
t214 = t181 * t176;
t139 = -Ifges(5,3) * t179 + t239;
t140 = Ifges(4,2) * t179 + t241;
t207 = t139 / 0.2e1 - t140 / 0.2e1;
t141 = Ifges(5,1) * t176 - t238;
t142 = Ifges(4,1) * t176 + t240;
t206 = t141 / 0.2e1 + t142 / 0.2e1;
t129 = t259 * t220;
t202 = qJD(3) * t143;
t42 = -qJD(5) * t81 - t178 * t129 + t175 * t202;
t82 = t178 * t143 + t175 * t213;
t43 = qJD(5) * t82 - t175 * t129 - t178 * t202;
t205 = t82 * t42 + t43 * t81;
t204 = t272 * t40 + t273 * t41;
t94 = t134 * t179 - t157;
t112 = t179 * pkin(4) - t133;
t192 = Ifges(5,3) * t176 + t238;
t97 = -Ifges(5,6) * t180 + t177 * t192;
t193 = -Ifges(4,2) * t176 + t240;
t98 = t177 * t193 - t234;
t201 = t97 - t98 + t234;
t200 = -pkin(7) + t214;
t198 = -t179 * mrSges(4,1) + t176 * mrSges(4,2);
t197 = mrSges(4,1) * t176 + mrSges(4,2) * t179;
t138 = -t179 * mrSges(5,1) - t176 * mrSges(5,3);
t196 = mrSges(5,1) * t176 - mrSges(5,3) * t179;
t18 = -t175 * t66 + t178 * t62;
t131 = -t175 * qJ(4) + t178 * t181;
t3 = t175 * t20 + t178 * t21 + t62 * t215 - t216 * t66;
t93 = qJD(3) * t214 + t224;
t155 = qJ(4) * t228;
t80 = t177 * t200 + t155;
t90 = mrSges(5,2) * t208 + (-mrSges(5,1) * qJD(2) - mrSges(5,2) * t220) * t177;
t47 = (-t179 * t222 - t180 * t220) * pkin(7) + t226;
t1 = -qJ(6) * t222 + qJD(6) * t180 + t3;
t2 = pkin(5) * t222 - t4;
t183 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3) + t204;
t33 = (t179 * t181 - t170) * t219 + t200 * t221 + t225;
t128 = pkin(5) - t131;
t127 = -mrSges(5,2) * t230 - mrSges(5,3) * t180;
t126 = mrSges(5,1) * t180 + mrSges(5,2) * t228;
t125 = -mrSges(4,1) * t180 - mrSges(4,3) * t228;
t124 = mrSges(4,2) * t180 - mrSges(4,3) * t230;
t121 = t193 * qJD(3);
t120 = t192 * qJD(3);
t119 = t197 * qJD(3);
t118 = t196 * qJD(3);
t117 = -qJ(6) + t132;
t109 = t196 * t177;
t106 = pkin(3) * t220 - t224;
t101 = -t155 + (pkin(3) * t176 + pkin(7)) * t177;
t96 = -qJD(6) + t102;
t92 = -mrSges(5,2) * t186 + mrSges(5,3) * t222;
t91 = -mrSges(4,2) * t222 - mrSges(4,3) * t186;
t89 = mrSges(4,1) * t222 - mrSges(4,3) * t187;
t84 = t174 - t94;
t79 = -Ifges(6,1) * t188 - Ifges(6,4) * t114;
t78 = -Ifges(7,1) * t188 + Ifges(7,5) * t114;
t77 = -Ifges(6,4) * t188 - Ifges(6,2) * t114;
t76 = -Ifges(7,5) * t188 + Ifges(7,3) * t114;
t75 = mrSges(6,1) * t114 - mrSges(6,2) * t188;
t74 = mrSges(7,1) * t114 + mrSges(7,3) * t188;
t65 = mrSges(4,1) * t186 + mrSges(4,2) * t187;
t64 = mrSges(5,1) * t186 - mrSges(5,3) * t187;
t60 = mrSges(6,1) * t104 + mrSges(6,2) * t105;
t59 = mrSges(7,1) * t104 - mrSges(7,3) * t105;
t58 = -t142 * t219 + (Ifges(4,5) * t177 + t180 * t195) * qJD(2);
t57 = -t141 * t219 + (Ifges(5,4) * t177 + t180 * t194) * qJD(2);
t56 = -t140 * t219 + (Ifges(4,6) * t177 + t180 * t193) * qJD(2);
t55 = -t139 * t219 + (Ifges(5,6) * t177 + t180 * t192) * qJD(2);
t53 = pkin(5) * t114 + qJ(6) * t188 + t112;
t52 = Ifges(6,1) * t105 - Ifges(6,4) * t104 + Ifges(6,5) * t180;
t51 = Ifges(7,1) * t105 + Ifges(7,4) * t180 + Ifges(7,5) * t104;
t50 = Ifges(6,4) * t105 - Ifges(6,2) * t104 + Ifges(6,6) * t180;
t49 = Ifges(7,5) * t105 + Ifges(7,6) * t180 + Ifges(7,3) * t104;
t48 = t222 * t255 - t199;
t46 = pkin(3) * t186 + pkin(7) * t221 + qJ(4) * t211 - t225;
t45 = t212 * t222 + t199;
t44 = -qJD(4) * t180 + t160 + t47;
t32 = Ifges(6,1) * t72 - Ifges(6,4) * t73;
t31 = Ifges(7,1) * t72 + Ifges(7,5) * t73;
t30 = Ifges(6,4) * t72 - Ifges(6,2) * t73;
t29 = Ifges(7,5) * t72 + Ifges(7,3) * t73;
t28 = mrSges(6,1) * t73 + mrSges(6,2) * t72;
t27 = mrSges(7,1) * t73 - mrSges(7,3) * t72;
t26 = pkin(5) * t104 - qJ(6) * t105 + t80;
t15 = -pkin(5) * t180 - t18;
t14 = qJ(6) * t180 + t246;
t12 = pkin(5) * t73 - qJ(6) * t72 + qJD(6) * t188 + t93;
t11 = mrSges(6,1) * t40 + mrSges(6,2) * t41;
t10 = mrSges(7,1) * t40 - mrSges(7,3) * t41;
t9 = Ifges(6,1) * t41 - Ifges(6,4) * t40 - Ifges(6,5) * t222;
t8 = Ifges(7,1) * t41 - Ifges(7,4) * t222 + Ifges(7,5) * t40;
t7 = Ifges(6,4) * t41 - Ifges(6,2) * t40 - Ifges(6,6) * t222;
t6 = Ifges(7,5) * t41 - Ifges(7,6) * t222 + Ifges(7,3) * t40;
t5 = pkin(5) * t40 - qJ(6) * t41 - qJD(6) * t105 + t33;
t13 = [((mrSges(3,2) * t261 + 0.2e1 * Ifges(3,4) * t180 + t176 * t201 + t179 * t233) * qJD(2) + t204 + t265) * t180 + (t18 * t4 + t246 * t3 + t33 * t80) * t263 + 0.2e1 * t246 * t23 + (t8 + t9) * t105 + (t6 - t7) * t104 + 0.2e1 * m(5) * (t101 * t46 + t44 * t83 + t45 * t84) + (t1 * t14 + t15 * t2 + t26 * t5) * t262 + (t47 * t95 + t48 * t94) * t264 + (t51 + t52) * t41 + (t49 - t50) * t40 + 0.2e1 * t44 * t127 + 0.2e1 * t47 * t124 + 0.2e1 * t48 * t125 + 0.2e1 * t45 * t126 + 0.2e1 * t46 * t109 + 0.2e1 * t101 * t64 + 0.2e1 * t4 * t87 + 0.2e1 * t2 * t88 + 0.2e1 * t84 * t90 + 0.2e1 * t83 * t92 + 0.2e1 * t94 * t89 + 0.2e1 * t95 * t91 + 0.2e1 * t80 * t11 + 0.2e1 * t1 * t85 + 0.2e1 * t3 * t86 + 0.2e1 * t5 * t59 + 0.2e1 * t33 * t60 + 0.2e1 * t14 * t22 + 0.2e1 * t18 * t24 + 0.2e1 * t15 * t25 + 0.2e1 * t26 * t10 + (t65 * t260 + (t57 + t58) * t179 + (t55 - t56) * t176 + (t201 * t179 + (t180 * t250 - t233) * t176) * qJD(3) + (mrSges(3,1) * t261 + (-0.2e1 * Ifges(3,4) + t250 * t179 + (-Ifges(4,6) + Ifges(5,6)) * t176) * t177 - t273 * t105 - t272 * t104 + (pkin(7) ^ 2 * t264 + t197 * t260 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(7,2)) - (2 * Ifges(6,3)) + t276) * t180) * qJD(2)) * t177; (t78 / 0.2e1 + t79 / 0.2e1) * t41 + (t76 / 0.2e1 - t77 / 0.2e1) * t40 + m(5) * (t101 * t106 + t133 * t46) + (-Ifges(3,6) * qJD(2) + t237 * t253 + (qJD(2) * mrSges(3,2) + t119) * pkin(7) - (t272 * t114 - t188 * t273) * qJD(2) / 0.2e1 + (t120 / 0.2e1 - t121 / 0.2e1 - t206 * qJD(3) + t250 * t253) * t176 + (t270 / 0.2e1 - Ifges(5,6) * t253 + qJD(3) * t207) * t179) * t177 + (-t114 * t3 - t18 * t72 + t188 * t4 - t246 * t73) * mrSges(6,3) + (-t1 * t114 - t14 * t73 + t15 * t72 - t188 * t2) * mrSges(7,2) - (t8 / 0.2e1 + t9 / 0.2e1) * t188 + m(6) * (t112 * t33 - t18 * t43 + t246 * t42 + t3 * t82 - t4 * t81 + t80 * t93) + (t31 / 0.2e1 + t32 / 0.2e1) * t105 + (t29 / 0.2e1 - t30 / 0.2e1) * t104 + (t6 / 0.2e1 - t7 / 0.2e1) * t114 + (t203 / 0.2e1 - t275 / 0.2e1) * t180 + t247 * t81 + t248 * t82 + t242 * t43 + t243 * t42 + (t57 / 0.2e1 + t58 / 0.2e1 + t45 * mrSges(5,2) - t48 * mrSges(4,3) + (-t83 * mrSges(5,2) - t95 * mrSges(4,3) + t234 / 0.2e1 + t97 / 0.2e1 - t98 / 0.2e1) * qJD(3)) * t176 + ((t91 + t92) * t179 + (-t89 + t90) * t176 + ((-t125 + t126) * t179 + (-t124 - t127) * t176) * qJD(3) + m(5) * (t176 * t45 + t179 * t44 + t218 * t84 - t220 * t83) + m(4) * (-t176 * t48 + t179 * t47 - t218 * t94 - t220 * t95)) * pkin(8) + m(7) * (t1 * t82 + t12 * t26 + t14 * t42 + t15 * t43 + t2 * t81 + t5 * t53) + (Ifges(3,5) + t206 * t179 + t207 * t176 + (-m(4) * pkin(2) - mrSges(3,1) + t198) * pkin(7)) * t221 + (-t55 / 0.2e1 + t56 / 0.2e1 + t44 * mrSges(5,2) + t47 * mrSges(4,3) + (t84 * mrSges(5,2) - t94 * mrSges(4,3) + t99 / 0.2e1 + t100 / 0.2e1) * qJD(3)) * t179 + t133 * t64 + t46 * t138 + t101 * t118 + t106 * t109 + t112 * t11 + t93 * t60 + t5 * t74 + t33 * t75 + t80 * t28 + t53 * t10 + t12 * t59 - pkin(2) * t65 + t26 * t27 + (t51 / 0.2e1 + t52 / 0.2e1) * t72 + (t49 / 0.2e1 - t50 / 0.2e1) * t73; -0.2e1 * pkin(2) * t119 + 0.2e1 * t112 * t28 + 0.2e1 * t133 * t118 + 0.2e1 * t12 * t74 + 0.2e1 * t53 * t27 + 0.2e1 * t93 * t75 + (t76 - t77) * t73 + (t78 + t79) * t72 + (-t120 + t121) * t179 + t270 * t176 - (t31 + t32) * t188 + (-t30 + t29) * t114 + 0.2e1 * (m(5) * t133 + t138) * t106 + (t12 * t53 + t205) * t262 + (t112 * t93 + t205) * t263 + ((t141 + t142) * t179 + (t139 - t140) * t176) * qJD(3) + 0.2e1 * t274 * (-t114 * t42 - t188 * t43 + t72 * t81 - t73 * t82); -t183 - t265 + m(6) * (t102 * t246 - t103 * t18 + t131 * t4 + t132 * t3) + (t249 * qJD(2) + (-t176 * t250 - t237) * qJD(3)) * t177 + t242 * t103 - Ifges(4,6) * t209 + m(5) * (-pkin(3) * t45 + qJ(4) * t44 + qJD(4) * t83) + m(7) * (t1 * t117 + t103 * t15 + t128 * t2 + t14 * t96) + qJD(4) * t127 + t128 * t25 + t131 * t24 + t132 * t23 + t117 * t22 + t102 * t86 - pkin(3) * t90 + qJ(4) * t92 + t96 * t85 + t44 * mrSges(5,3) - t45 * mrSges(5,1) - t47 * mrSges(4,2) + t48 * mrSges(4,1); -Ifges(4,6) * t220 + t252 * t43 + t251 * t42 + m(7) * (t117 * t42 + t128 * t43 + t82 * t96 + t235) + m(6) * (t102 * t82 - t131 * t43 + t132 * t42 + t235) + (qJD(3) * t191 + t217) * mrSges(5,2) + (-t102 * t114 - t131 * t72 - t132 * t73 - t232) * mrSges(6,3) + (-t114 * t96 - t117 * t73 + t128 * t72 - t232) * mrSges(7,2) + (m(5) * t217 + (m(5) * t191 + t138 + t198) * qJD(3)) * pkin(8) - t203 + t275; 0.2e1 * t236 - 0.2e1 * t96 * mrSges(7,3) + (t102 * t132 - t103 * t131) * t263 + (t103 * t128 + t117 * t96) * t262 + 0.2e1 * t252 * t103 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4); -t247 * t178 + t248 * t175 + (t175 * t242 + t178 * t243) * qJD(5) + m(7) * (t1 * t175 - t178 * t2 + (t14 * t178 + t15 * t175) * qJD(5)) + m(6) * (t175 * t3 + t178 * t4 + (-t175 * t18 + t178 * t246) * qJD(5)) + m(5) * t45 + t90; (m(5) * pkin(8) + mrSges(5,2)) * t218 + (m(6) + m(7)) * (t175 * t42 - t178 * t43 + t82 * t215 + t216 * t81) + t274 * (-t175 * t73 - t178 * t72 + (-t114 * t178 - t175 * t188) * qJD(5)); (t252 * t175 + t251 * t178) * qJD(5) + m(6) * (t102 * t175 - t231 + (-t131 * t175 + t132 * t178) * qJD(5)) + m(7) * (-t231 + t175 * t96 + (t117 * t178 + t128 * t175) * qJD(5)); 0; m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t14) + qJD(6) * t85 + qJ(6) * t22 - pkin(5) * t25 - t249 * t222 + t183; t82 * t271 + (-pkin(5) * t72 - qJ(6) * t73 - qJD(6) * t114) * mrSges(7,2) + t203 + t267 * t43 + t266 * t42; -t236 + m(7) * (qJ(6) * t96 + qJD(6) * t117) + (t96 - qJD(6)) * mrSges(7,3) + t267 * t103; t266 * t215 + (qJD(5) * t267 + t271) * t175; 0.2e1 * t268 * qJD(6); m(7) * t2 + t25; m(7) * t43 + t72 * mrSges(7,2); m(7) * t103; m(7) * t216; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
