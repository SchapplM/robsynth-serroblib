% Calculate time derivative of joint inertia matrix for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:46
% EndTime: 2019-03-09 04:21:58
% DurationCPUTime: 5.86s
% Computational Cost: add. (9716->527), mult. (28125->780), div. (0->0), fcn. (29495->12), ass. (0->225)
t287 = Ifges(5,4) - Ifges(4,5);
t286 = Ifges(5,5) - Ifges(4,6);
t160 = sin(pkin(7));
t162 = cos(pkin(7));
t161 = sin(pkin(6));
t237 = cos(pkin(12));
t201 = t161 * t237;
t238 = cos(pkin(6));
t119 = -t160 * t201 + t162 * t238;
t159 = sin(pkin(12));
t111 = (-pkin(9) * t159 * t160 - pkin(2) * t237 - pkin(1)) * t161;
t165 = sin(qJ(3));
t259 = cos(qJ(3));
t181 = t259 * t201;
t202 = qJD(3) * t259;
t192 = t160 * t202;
t226 = qJD(2) * t161;
t207 = t159 * t226;
t193 = t162 * t207;
t210 = pkin(1) * t238;
t154 = t237 * t210;
t233 = t159 * t161;
t101 = t238 * pkin(2) + t154 + (-pkin(9) * t162 - qJ(2)) * t233;
t236 = t101 * t162;
t195 = t259 * t236;
t200 = t238 * t160;
t228 = qJ(2) * t201 + t159 * t210;
t98 = (t162 * t201 + t200) * pkin(9) + t228;
t50 = -t165 * (qJD(3) * t98 + t193) + qJD(2) * t181 + qJD(3) * t195 + t111 * t192;
t45 = -t119 * qJD(4) - t50;
t163 = sin(qJ(6));
t166 = cos(qJ(6));
t167 = cos(qJ(5));
t216 = qJD(6) * t167;
t164 = sin(qJ(5));
t223 = qJD(5) * t164;
t172 = t163 * t223 - t166 * t216;
t209 = t160 * t259;
t175 = -t164 * t162 - t167 * t209;
t225 = qJD(5) * t175;
t285 = -2 * mrSges(4,3);
t284 = t163 / 0.2e1;
t260 = t166 / 0.2e1;
t283 = Ifges(4,4) + Ifges(5,6);
t199 = t237 * t165;
t208 = t259 * t159;
t100 = t165 * t200 + (t162 * t199 + t208) * t161;
t99 = -t162 * t181 + t165 * t233 - t259 * t200;
t74 = t119 * t167 + t164 * t99;
t64 = t100 * t166 - t163 * t74;
t65 = t100 * t163 + t166 * t74;
t31 = -mrSges(7,1) * t64 + mrSges(7,2) * t65;
t67 = mrSges(6,1) * t100 - mrSges(6,3) * t74;
t251 = t31 - t67;
t73 = t119 * t164 - t99 * t167;
t56 = mrSges(6,1) * t73 + mrSges(6,2) * t74;
t77 = mrSges(5,1) * t99 - mrSges(5,3) * t119;
t282 = -t77 + t56;
t142 = t164 * mrSges(6,1) + t167 * mrSges(6,2);
t281 = mrSges(5,3) + t142;
t227 = t163 ^ 2 + t166 ^ 2;
t94 = t99 * qJD(3);
t95 = t100 * qJD(3);
t279 = t286 * t95 + t287 * t94;
t171 = -t111 * t209 + t165 * t98 - t195;
t267 = pkin(3) + pkin(10);
t35 = t100 * pkin(4) - t119 * t267 + t171;
t72 = -t101 * t160 + t162 * t111;
t177 = -qJ(4) * t100 + t72;
t39 = t267 * t99 + t177;
t250 = t164 * t35 + t167 * t39;
t278 = qJD(5) * t250;
t277 = 0.2e1 * t99;
t62 = -qJD(5) * t73 + t164 * t95;
t26 = -qJD(6) * t65 - t163 * t62 - t166 * t94;
t27 = qJD(6) * t64 - t163 * t94 + t166 * t62;
t12 = -mrSges(7,1) * t26 + mrSges(7,2) * t27;
t88 = t259 * t98;
t51 = (t162 * t208 + t199) * t226 + (t88 + (t111 * t160 + t236) * t165) * qJD(3);
t37 = -t94 * pkin(4) + t51;
t135 = t160 * t207;
t178 = qJ(4) * t94 - qJD(4) * t100 + t135;
t41 = t267 * t95 + t178;
t6 = -t164 * t41 + t167 * t37 - t278;
t4 = pkin(5) * t94 - t6;
t276 = -m(7) * t4 - t12;
t232 = t160 * t165;
t206 = qJD(3) * t232;
t104 = t164 * t206 + t225;
t194 = t164 * t209;
t121 = t167 * t162 - t194;
t221 = qJD(5) * t167;
t105 = -qJD(5) * t194 + t162 * t221 - t167 * t206;
t235 = t105 * t167;
t275 = qJD(5) * (t121 * t167 - t164 * t175) + t104 * t164 - t235;
t141 = -mrSges(7,1) * t166 + mrSges(7,2) * t163;
t274 = -m(7) * pkin(5) - mrSges(6,1) + t141;
t273 = 0.2e1 * m(7);
t272 = -0.2e1 * t267;
t271 = m(6) / 0.2e1;
t270 = m(7) / 0.2e1;
t243 = Ifges(7,4) * t163;
t144 = Ifges(7,2) * t166 + t243;
t242 = Ifges(7,4) * t166;
t188 = -Ifges(7,2) * t163 + t242;
t82 = -t144 * t216 + (Ifges(7,6) * t167 - t164 * t188) * qJD(5);
t269 = t82 / 0.2e1;
t146 = Ifges(7,1) * t163 + t242;
t189 = Ifges(7,1) * t166 - t243;
t83 = -t146 * t216 + (Ifges(7,5) * t167 - t164 * t189) * qJD(5);
t268 = t83 / 0.2e1;
t116 = Ifges(7,5) * t164 + t167 * t189;
t266 = t116 / 0.2e1;
t217 = qJD(6) * t166;
t156 = Ifges(7,5) * t217;
t218 = qJD(6) * t163;
t265 = -Ifges(7,6) * t218 / 0.2e1 + t156 / 0.2e1;
t264 = Ifges(7,5) * t284 + Ifges(7,6) * t260;
t263 = t144 / 0.2e1;
t262 = -t163 / 0.2e1;
t261 = -t166 / 0.2e1;
t258 = m(7) * t164;
t257 = pkin(5) * t164;
t256 = pkin(11) * t167;
t255 = t94 * mrSges(5,1);
t254 = t94 * Ifges(6,5);
t253 = t94 * Ifges(6,6);
t252 = -mrSges(4,1) + mrSges(5,2);
t48 = -mrSges(6,1) * t94 - mrSges(6,3) * t62;
t249 = t48 - t12;
t246 = mrSges(7,3) * t167;
t245 = Ifges(6,4) * t164;
t244 = Ifges(6,4) * t167;
t241 = Ifges(7,6) * t163;
t240 = t100 * Ifges(6,5);
t239 = t100 * Ifges(6,6);
t234 = t175 * t105;
t231 = t164 * t267;
t230 = t167 * t267;
t229 = qJ(4) * t192 + qJD(4) * t232;
t224 = qJD(5) * t163;
t222 = qJD(5) * t166;
t138 = qJ(4) - t256 + t257;
t112 = t166 * t138 + t163 * t231;
t220 = qJD(6) * t112;
t113 = t163 * t138 - t166 * t231;
t219 = qJD(6) * t113;
t63 = qJD(5) * t74 - t95 * t167;
t9 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t63;
t215 = Ifges(6,5) * t62 - Ifges(6,6) * t63 - Ifges(6,3) * t94;
t58 = t111 * t232 + t165 * t236 + t88;
t212 = t259 * t51;
t205 = t267 * t221;
t124 = qJD(4) + (pkin(5) * t167 + pkin(11) * t164) * qJD(5);
t79 = t163 * t124 - t166 * t205 + t220;
t198 = t79 - t220;
t80 = t166 * t124 + t163 * t205 - t219;
t197 = -t80 - t219;
t53 = -t119 * qJ(4) - t58;
t32 = -pkin(4) * t95 - t45;
t13 = pkin(5) * t63 - pkin(11) * t62 + t32;
t5 = t164 * t37 + t167 * t41 + t35 * t221 - t223 * t39;
t3 = -pkin(11) * t94 + t5;
t15 = pkin(11) * t100 + t250;
t42 = -pkin(4) * t99 - t53;
t20 = pkin(5) * t73 - pkin(11) * t74 + t42;
t7 = -t15 * t163 + t166 * t20;
t1 = qJD(6) * t7 + t13 * t163 + t166 * t3;
t8 = t15 * t166 + t163 * t20;
t2 = -qJD(6) * t8 + t13 * t166 - t163 * t3;
t191 = t1 * t166 - t163 * t2;
t190 = mrSges(7,1) * t163 + mrSges(7,2) * t166;
t18 = -mrSges(7,2) * t63 + mrSges(7,3) * t26;
t19 = mrSges(7,1) * t63 - mrSges(7,3) * t27;
t187 = -t163 * t19 + t166 * t18;
t106 = -t121 * t163 + t166 * t232;
t70 = qJD(6) * t106 + t166 * t104 + t163 * t192;
t107 = t121 * t166 + t163 * t232;
t71 = -qJD(6) * t107 - t163 * t104 + t166 * t192;
t186 = -t163 * t71 + t166 * t70;
t16 = -t164 * t39 + t167 * t35;
t22 = Ifges(7,4) * t65 + Ifges(7,2) * t64 + Ifges(7,6) * t73;
t23 = Ifges(7,1) * t65 + Ifges(7,4) * t64 + Ifges(7,5) * t73;
t179 = t22 * t262 + t23 * t260;
t173 = t163 * t216 + t164 * t222;
t170 = -t217 * t7 - t218 * t8 + t191;
t81 = -Ifges(7,5) * t173 + t172 * Ifges(7,6) + Ifges(7,3) * t221;
t147 = Ifges(6,1) * t167 - t245;
t145 = -Ifges(6,2) * t164 + t244;
t134 = mrSges(7,1) * t164 - t166 * t246;
t133 = -mrSges(7,2) * t164 - t163 * t246;
t131 = (-Ifges(6,1) * t164 - t244) * qJD(5);
t130 = t189 * qJD(6);
t129 = (-Ifges(6,2) * t167 - t245) * qJD(5);
t128 = t188 * qJD(6);
t126 = (mrSges(6,1) * t167 - mrSges(6,2) * t164) * qJD(5);
t125 = t190 * qJD(6);
t122 = t190 * t167;
t115 = Ifges(7,6) * t164 + t167 * t188;
t114 = Ifges(7,3) * t164 + (Ifges(7,5) * t166 - t241) * t167;
t110 = -mrSges(7,2) * t221 + mrSges(7,3) * t172;
t109 = mrSges(7,1) * t221 + mrSges(7,3) * t173;
t97 = -mrSges(7,1) * t172 - mrSges(7,2) * t173;
t78 = mrSges(5,1) * t100 + mrSges(5,2) * t119;
t76 = mrSges(4,1) * t119 - mrSges(4,3) * t100;
t75 = -mrSges(4,2) * t119 - mrSges(4,3) * t99;
t69 = mrSges(4,1) * t95 - mrSges(4,2) * t94;
t68 = -mrSges(5,2) * t95 + mrSges(5,3) * t94;
t66 = -mrSges(6,2) * t100 - mrSges(6,3) * t73;
t55 = pkin(3) * t95 + t178;
t54 = -t119 * pkin(3) + t171;
t52 = pkin(3) * t99 + t177;
t49 = mrSges(6,2) * t94 - mrSges(6,3) * t63;
t47 = Ifges(6,1) * t74 - Ifges(6,4) * t73 + t240;
t46 = Ifges(6,4) * t74 - Ifges(6,2) * t73 + t239;
t44 = mrSges(7,1) * t73 - mrSges(7,3) * t65;
t43 = -mrSges(7,2) * t73 + mrSges(7,3) * t64;
t30 = mrSges(6,1) * t63 + mrSges(6,2) * t62;
t29 = Ifges(6,1) * t62 - Ifges(6,4) * t63 - t254;
t28 = Ifges(6,4) * t62 - Ifges(6,2) * t63 - t253;
t21 = Ifges(7,5) * t65 + Ifges(7,6) * t64 + Ifges(7,3) * t73;
t14 = -pkin(5) * t100 - t16;
t11 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t63;
t10 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t63;
t17 = [(0.2e1 * mrSges(5,1) * t53 + t58 * t285 + (Ifges(4,2) + Ifges(5,3)) * t277 + t286 * t119 - 0.2e1 * t283 * t100) * t95 + (t171 * t285 - Ifges(6,5) * t74 + Ifges(6,6) * t73 + t283 * t277 + t287 * t119 + (-(2 * Ifges(4,1)) - (2 * Ifges(5,2)) - Ifges(6,3)) * t100) * t94 + 0.2e1 * m(4) * (t135 * t72 + t171 * t51 + t50 * t58) + (t21 - t46) * t63 + (t1 * t8 + t14 * t4 + t2 * t7) * t273 + 0.2e1 * t55 * (-mrSges(5,2) * t99 - mrSges(5,3) * t100) + t279 * t119 + 0.2e1 * t45 * t77 + 0.2e1 * t51 * t78 + 0.2e1 * t72 * t69 + t74 * t29 + 0.2e1 * t50 * t75 - 0.2e1 * t51 * t76 + t65 * t11 + 0.2e1 * t5 * t66 + 0.2e1 * t6 * t67 + 0.2e1 * t52 * t68 + t62 * t47 + t64 * t10 + 0.2e1 * t32 * t56 + 0.2e1 * t16 * t48 + 0.2e1 * t42 * t30 + 0.2e1 * t1 * t43 + 0.2e1 * t2 * t44 + t27 * t23 + 0.2e1 * t4 * t31 + t26 * t22 + 0.2e1 * t8 * t18 + 0.2e1 * t7 * t19 + 0.2e1 * t14 * t12 + 0.2e1 * m(5) * (t45 * t53 + t51 * t54 + t52 * t55) + 0.2e1 * t250 * t49 + 0.2e1 * m(6) * (t16 * t6 + t250 * t5 + t32 * t42) + (t9 - t28) * t73 + 0.2e1 * m(3) * (t228 * t237 + (qJ(2) * t233 - t154) * t159) * t226 + 0.2e1 * qJD(2) * (-mrSges(3,2) * t238 + mrSges(3,3) * t201) * t201 + 0.2e1 * (mrSges(4,1) * t99 + mrSges(4,2) * t100) * t135 + t100 * t215 - 0.2e1 * (mrSges(3,1) * t238 - mrSges(3,3) * t233) * t207 - 0.2e1 * t54 * t255; t30 * t232 + t121 * t49 + t104 * t66 + t106 * t19 + t107 * t18 + t70 * t43 + t71 * t44 + m(6) * (t104 * t250 + t121 * t5) + m(7) * (t1 * t107 + t106 * t2 + t7 * t71 + t70 * t8) + (-t76 + t78) * t206 + (m(4) * (t193 - t212 + t165 * t50 + (t165 * t171 + t259 * t58) * qJD(3)) + m(6) * (t165 * t32 + t202 * t42) + m(5) * (-t212 - t165 * t45 + (t165 * t54 - t259 * t53) * qJD(3))) * t160 + (m(5) * t55 + t68 + t69) * t162 - (-m(6) * t6 - t276 - t48) * t175 + (-m(6) * t16 + m(7) * t14 + t251) * t105 + (t75 + t282) * t192 + (mrSges(4,3) + mrSges(5,1)) * (t94 * t209 - t95 * t232); 0.2e1 * m(7) * (t106 * t71 + t107 * t70 - t234) + 0.2e1 * m(6) * (t160 ^ 2 * t165 * t202 + t121 * t104 - t234); (-t129 / 0.2e1 + t81 / 0.2e1) * t73 + m(5) * (-pkin(3) * t51 - qJ(4) * t45 - qJD(4) * t53) + t279 + ((-t46 / 0.2e1 + t21 / 0.2e1 - t250 * mrSges(6,3) - t239 / 0.2e1) * t167 + (-t47 / 0.2e1 + t16 * mrSges(6,3) - t240 / 0.2e1 - t179) * t164 - (t167 * t66 + t251 * t164 + t14 * t258 + m(6) * (-t16 * t164 + t167 * t250)) * t267) * qJD(5) + (-t267 * t49 - t5 * mrSges(6,3) - t28 / 0.2e1 + t9 / 0.2e1 + t253 / 0.2e1) * t164 + (t11 * t260 + t10 * t262 - t6 * mrSges(6,3) - t254 / 0.2e1 + t29 / 0.2e1 - t249 * t267 + (t22 * t261 + t23 * t262) * qJD(6)) * t167 + (-t145 / 0.2e1 + t114 / 0.2e1) * t63 + t27 * t266 + t65 * t268 + t64 * t269 + (pkin(3) * t94 - qJ(4) * t95) * mrSges(5,1) + t74 * t131 / 0.2e1 + t1 * t133 + t2 * t134 + t32 * t142 + t62 * t147 / 0.2e1 + t4 * t122 + t42 * t126 + t7 * t109 + t8 * t110 + t112 * t19 + t113 * t18 + t26 * t115 / 0.2e1 + t14 * t97 + t79 * t43 + t80 * t44 - t50 * mrSges(4,2) - t45 * mrSges(5,3) + qJ(4) * t30 + m(7) * (t113 * t1 + t112 * t2 + t230 * t4 + t80 * t7 + t79 * t8) + m(6) * (qJ(4) * t32 + qJD(4) * t42 - t230 * t6 - t231 * t5) + t282 * qJD(4) + t252 * t51; t105 * t122 + t106 * t109 + t107 * t110 - t175 * t97 + t70 * t133 + t71 * t134 + ((qJD(3) * t252 + t126) * t165 + (-mrSges(4,2) + t281) * t202) * t160 + m(7) * (t80 * t106 + t79 * t107 + t112 * t71 + t113 * t70) + m(6) * t229 + m(5) * (-pkin(3) * t206 + t229) + ((-t175 * t223 - t235) * t270 + t275 * t271) * t272 - t275 * mrSges(6,3); (t112 * t80 + t113 * t79) * t273 + 0.2e1 * t79 * t133 + 0.2e1 * t113 * t110 + 0.2e1 * t80 * t134 + 0.2e1 * t112 * t109 + 0.2e1 * qJ(4) * t126 + 0.2e1 * ((m(5) + m(6)) * qJ(4) + t281) * qJD(4) + (-t129 + t81 + (t115 * t163 - t116 * t166 + t122 * t272 - t147) * qJD(5)) * t164 + (-t163 * t82 + t166 * t83 + 0.2e1 * t267 * t97 + t131 + (-t115 * t166 - t116 * t163) * qJD(6) + (-0.2e1 * t258 * t267 ^ 2 + t114 - t145) * qJD(5)) * t167; m(5) * t51 - t255 + ((-t163 * t44 + t166 * t43 + t66) * qJD(5) + m(7) * (t222 * t8 - t224 * t7 - t4) + m(6) * (t6 + t278) + t249) * t167 + (t49 + (-t163 * t43 - t166 * t44) * qJD(6) + t251 * qJD(5) + m(7) * (qJD(5) * t14 + t170) + m(6) * (-qJD(5) * t16 + t5) + t187) * t164; m(5) * t206 + 0.2e1 * ((-t106 * t224 + t107 * t222 - t105) * t270 + (qJD(5) * t121 - t105) * t271) * t167 + 0.2e1 * ((-t106 * t217 - t107 * t218 + t186 - t225) * t270 + (t104 - t225) * t271) * t164; (-t97 + (m(7) * (-t112 * t163 + t113 * t166) + t166 * t133 - t163 * t134) * qJD(5)) * t167 + (m(7) * (-t112 * t217 - t113 * t218 - t163 * t80 + t166 * t79 + 0.2e1 * t205) - t133 * t218 + t166 * t110 - t134 * t217 - t163 * t109 + qJD(5) * t122) * t164; 0.2e1 * (-0.1e1 + t227) * t221 * t258; t10 * t260 + t11 * t284 + t65 * t130 / 0.2e1 + t4 * t141 + t63 * t264 + t26 * t263 + t27 * t146 / 0.2e1 + t14 * t125 + t73 * t265 + t64 * t128 / 0.2e1 - t5 * mrSges(6,2) + t6 * mrSges(6,1) + t179 * qJD(6) + t276 * pkin(5) + ((-t163 * t8 - t166 * t7) * qJD(6) + t191) * mrSges(7,3) + (m(7) * t170 - t217 * t44 - t218 * t43 + t187) * pkin(11) + t215; -t104 * mrSges(6,2) - t175 * t125 + (m(7) * pkin(11) + mrSges(7,3)) * ((-t106 * t166 - t107 * t163) * qJD(6) + t186) + t274 * t105; -pkin(5) * t97 + (t265 + (-t267 * t274 - Ifges(6,5)) * qJD(5)) * t164 + (qJD(6) * t266 - t146 * t223 / 0.2e1 + t269 + t198 * mrSges(7,3) + (m(7) * t198 - qJD(6) * t134 + t110) * pkin(11)) * t166 + (-qJD(6) * t115 / 0.2e1 + t223 * t263 + t268 + t197 * mrSges(7,3) + (m(7) * t197 - qJD(6) * t133 - t109) * pkin(11)) * t163 + (t130 * t260 + t128 * t262 + t267 * t125 + (t144 * t261 + t146 * t262) * qJD(6) + (mrSges(6,2) * t267 - Ifges(6,6) + t264) * qJD(5)) * t167; -t167 * t125 + (t164 * t141 + m(7) * (t227 * t256 - t257) + t227 * t246 - t142) * qJD(5); -0.2e1 * pkin(5) * t125 + t128 * t166 + t130 * t163 + (-t144 * t163 + t146 * t166) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t9; mrSges(7,1) * t71 - mrSges(7,2) * t70; mrSges(7,1) * t80 - mrSges(7,2) * t79 + t81; (t164 * t218 - t166 * t221) * mrSges(7,2) + (-t163 * t221 - t164 * t217) * mrSges(7,1); t156 + (pkin(11) * t141 - t241) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t17(1) t17(2) t17(4) t17(7) t17(11) t17(16); t17(2) t17(3) t17(5) t17(8) t17(12) t17(17); t17(4) t17(5) t17(6) t17(9) t17(13) t17(18); t17(7) t17(8) t17(9) t17(10) t17(14) t17(19); t17(11) t17(12) t17(13) t17(14) t17(15) t17(20); t17(16) t17(17) t17(18) t17(19) t17(20) t17(21);];
Mq  = res;
