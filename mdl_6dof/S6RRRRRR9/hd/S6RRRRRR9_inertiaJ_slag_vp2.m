% Calculate joint inertia matrix for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:45:05
% EndTime: 2018-11-23 18:45:08
% DurationCPUTime: 3.15s
% Computational Cost: add. (9190->599), mult. (23677->867), div. (0->0), fcn. (27070->14), ass. (0->222)
t293 = 2 * pkin(11);
t228 = sin(qJ(2));
t221 = sin(pkin(6));
t233 = cos(qJ(2));
t257 = t221 * t233;
t223 = cos(pkin(6));
t270 = pkin(1) * t223;
t168 = pkin(9) * t257 + t228 * t270;
t220 = sin(pkin(7));
t222 = cos(pkin(7));
t255 = t222 * t233;
t250 = t221 * t255;
t110 = (t220 * t223 + t250) * pkin(10) + t168;
t204 = t233 * t270;
t258 = t221 * t228;
t119 = pkin(2) * t223 + t204 + (-pkin(10) * t222 - pkin(9)) * t258;
t139 = (-pkin(10) * t220 * t228 - pkin(2) * t233 - pkin(1)) * t221;
t227 = sin(qJ(3));
t232 = cos(qJ(3));
t54 = -t227 * t110 + (t119 * t222 + t139 * t220) * t232;
t259 = t220 * t232;
t115 = -t223 * t259 + t227 * t258 - t232 * t250;
t225 = sin(qJ(5));
t230 = cos(qJ(5));
t260 = t220 * t227;
t116 = t223 * t260 + (t227 * t255 + t228 * t232) * t221;
t161 = -t220 * t257 + t222 * t223;
t226 = sin(qJ(4));
t231 = cos(qJ(4));
t89 = t116 * t231 + t161 * t226;
t58 = t115 * t230 - t225 * t89;
t59 = t115 * t225 + t230 * t89;
t88 = t116 * t226 - t231 * t161;
t24 = Ifges(6,4) * t59 + Ifges(6,2) * t58 + Ifges(6,6) * t88;
t292 = t24 / 0.2e1;
t25 = Ifges(6,1) * t59 + Ifges(6,4) * t58 + Ifges(6,5) * t88;
t291 = t25 / 0.2e1;
t42 = Ifges(5,1) * t89 - Ifges(5,4) * t88 + Ifges(5,5) * t115;
t290 = t42 / 0.2e1;
t163 = t222 * t226 + t231 * t260;
t124 = -t163 * t225 - t230 * t259;
t125 = t163 * t230 - t225 * t259;
t162 = -t231 * t222 + t226 * t260;
t69 = Ifges(6,4) * t125 + Ifges(6,2) * t124 + Ifges(6,6) * t162;
t289 = t69 / 0.2e1;
t70 = Ifges(6,1) * t125 + Ifges(6,4) * t124 + Ifges(6,5) * t162;
t288 = t70 / 0.2e1;
t287 = -pkin(13) - pkin(12);
t224 = sin(qJ(6));
t229 = cos(qJ(6));
t177 = t224 * t230 + t225 * t229;
t155 = t177 * t226;
t176 = -t224 * t225 + t229 * t230;
t156 = t176 * t226;
t100 = Ifges(7,4) * t156 - Ifges(7,2) * t155 - Ifges(7,6) * t231;
t286 = t100 / 0.2e1;
t101 = Ifges(7,1) * t156 - Ifges(7,4) * t155 - Ifges(7,5) * t231;
t285 = t101 / 0.2e1;
t104 = Ifges(5,1) * t163 - Ifges(5,4) * t162 - Ifges(5,5) * t259;
t284 = t104 / 0.2e1;
t122 = Ifges(7,4) * t177 + Ifges(7,2) * t176;
t283 = t122 / 0.2e1;
t123 = Ifges(7,1) * t177 + Ifges(7,4) * t176;
t282 = t123 / 0.2e1;
t264 = Ifges(6,4) * t230;
t152 = -Ifges(6,6) * t231 + (-Ifges(6,2) * t225 + t264) * t226;
t281 = t152 / 0.2e1;
t265 = Ifges(6,4) * t225;
t153 = -Ifges(6,5) * t231 + (Ifges(6,1) * t230 - t265) * t226;
t280 = t153 / 0.2e1;
t279 = -t155 / 0.2e1;
t278 = t156 / 0.2e1;
t277 = t176 / 0.2e1;
t276 = t177 / 0.2e1;
t187 = Ifges(6,2) * t230 + t265;
t275 = t187 / 0.2e1;
t189 = Ifges(6,1) * t225 + t264;
t274 = t189 / 0.2e1;
t190 = Ifges(5,1) * t226 + Ifges(5,4) * t231;
t273 = t190 / 0.2e1;
t272 = -t225 / 0.2e1;
t271 = t230 / 0.2e1;
t269 = pkin(2) * t232;
t268 = pkin(11) * t226;
t267 = pkin(11) * t231;
t266 = -Ifges(7,3) - Ifges(6,3);
t80 = -t119 * t220 + t222 * t139;
t46 = pkin(3) * t115 - pkin(11) * t116 + t80;
t256 = t222 * t227;
t55 = t232 * t110 + t119 * t256 + t139 * t260;
t50 = pkin(11) * t161 + t55;
t20 = t226 * t46 + t231 * t50;
t18 = pkin(12) * t115 + t20;
t49 = -pkin(3) * t161 - t54;
t28 = pkin(4) * t88 - pkin(12) * t89 + t49;
t7 = t230 * t18 + t225 * t28;
t198 = pkin(10) * t260;
t146 = t198 + (-pkin(3) - t269) * t222;
t90 = pkin(4) * t162 - pkin(12) * t163 + t146;
t167 = pkin(2) * t256 + pkin(10) * t259;
t147 = pkin(11) * t222 + t167;
t148 = (-pkin(3) * t232 - pkin(11) * t227 - pkin(2)) * t220;
t98 = t231 * t147 + t226 * t148;
t92 = -pkin(12) * t259 + t98;
t52 = t225 * t90 + t230 * t92;
t263 = Ifges(7,3) * t231;
t166 = -pkin(9) * t258 + t204;
t262 = t166 * mrSges(3,1);
t261 = t168 * mrSges(3,2);
t254 = t225 * t226;
t253 = t226 * t230;
t252 = Ifges(7,5) * t156 - Ifges(7,6) * t155;
t121 = Ifges(7,5) * t177 + Ifges(7,6) * t176;
t182 = -pkin(4) * t231 - pkin(12) * t226 - pkin(3);
t141 = t225 * t182 + t230 * t267;
t185 = Ifges(6,5) * t225 + Ifges(6,6) * t230;
t186 = Ifges(5,5) * t226 + Ifges(5,6) * t231;
t251 = t225 ^ 2 + t230 ^ 2;
t31 = -t224 * t59 + t229 * t58;
t32 = t224 * t58 + t229 * t59;
t8 = Ifges(7,5) * t32 + Ifges(7,6) * t31 + Ifges(7,3) * t88;
t23 = Ifges(6,5) * t59 + Ifges(6,6) * t58 + Ifges(6,3) * t88;
t77 = t124 * t229 - t125 * t224;
t78 = t124 * t224 + t125 * t229;
t35 = Ifges(7,5) * t78 + Ifges(7,6) * t77 + Ifges(7,3) * t162;
t40 = Ifges(5,5) * t89 - Ifges(5,6) * t88 + Ifges(5,3) * t115;
t65 = Ifges(4,5) * t116 - Ifges(4,6) * t115 + Ifges(4,3) * t161;
t68 = Ifges(6,5) * t125 + Ifges(6,6) * t124 + Ifges(6,3) * t162;
t142 = Ifges(4,5) * t260 + Ifges(4,6) * t259 + Ifges(4,3) * t222;
t249 = Ifges(3,5) * t258 + Ifges(3,6) * t257 + Ifges(3,3) * t223;
t248 = t185 / 0.2e1 + t121 / 0.2e1;
t6 = -t18 * t225 + t230 * t28;
t51 = -t225 * t92 + t230 * t90;
t19 = -t226 * t50 + t231 * t46;
t97 = -t226 * t147 + t148 * t231;
t41 = Ifges(5,4) * t89 - Ifges(5,2) * t88 + Ifges(5,6) * t115;
t247 = t8 / 0.2e1 + t23 / 0.2e1 - t41 / 0.2e1;
t103 = Ifges(5,4) * t163 - Ifges(5,2) * t162 - Ifges(5,6) * t259;
t246 = t35 / 0.2e1 + t68 / 0.2e1 - t103 / 0.2e1;
t245 = Ifges(6,5) * t253 - Ifges(6,6) * t254;
t151 = -Ifges(6,3) * t231 + t245;
t188 = Ifges(5,4) * t226 + Ifges(5,2) * t231;
t99 = t252 - t263;
t244 = -t188 / 0.2e1 + t99 / 0.2e1 + t151 / 0.2e1;
t91 = pkin(4) * t259 - t97;
t243 = mrSges(6,1) * t225 + mrSges(6,2) * t230;
t173 = t230 * t182;
t114 = -pkin(13) * t253 + t173 + (-pkin(11) * t225 - pkin(5)) * t231;
t127 = -pkin(13) * t254 + t141;
t72 = t114 * t229 - t127 * t224;
t73 = t114 * t224 + t127 * t229;
t241 = t72 * mrSges(7,1) - t73 * mrSges(7,2) + t252;
t102 = Ifges(5,5) * t163 - Ifges(5,6) * t162 - Ifges(5,3) * t259;
t192 = t287 * t225;
t193 = t287 * t230;
t132 = t192 * t229 + t193 * t224;
t133 = t192 * t224 - t193 * t229;
t240 = t132 * mrSges(7,1) - t133 * mrSges(7,2) + t121;
t4 = pkin(5) * t88 - pkin(13) * t59 + t6;
t5 = pkin(13) * t58 + t7;
t2 = -t224 * t5 + t229 * t4;
t3 = t224 * t4 + t229 * t5;
t239 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t8;
t17 = -pkin(4) * t115 - t19;
t34 = pkin(5) * t162 - pkin(13) * t125 + t51;
t43 = pkin(13) * t124 + t52;
t14 = -t224 * t43 + t229 * t34;
t15 = t224 * t34 + t229 * t43;
t238 = t14 * mrSges(7,1) - t15 * mrSges(7,2) + t35;
t237 = (mrSges(7,1) * t229 - mrSges(7,2) * t224) * pkin(5);
t235 = pkin(11) ^ 2;
t219 = t231 ^ 2;
t217 = t226 ^ 2;
t215 = t217 * t235;
t208 = -pkin(5) * t230 - pkin(4);
t184 = -mrSges(5,1) * t231 + mrSges(5,2) * t226;
t183 = -mrSges(6,1) * t230 + mrSges(6,2) * t225;
t181 = (pkin(5) * t225 + pkin(11)) * t226;
t179 = -mrSges(6,1) * t231 - mrSges(6,3) * t253;
t178 = mrSges(6,2) * t231 - mrSges(6,3) * t254;
t175 = -mrSges(4,2) * t222 + mrSges(4,3) * t259;
t174 = mrSges(4,1) * t222 - mrSges(4,3) * t260;
t169 = t243 * t226;
t165 = t222 * t269 - t198;
t164 = (-mrSges(4,1) * t232 + mrSges(4,2) * t227) * t220;
t144 = Ifges(4,5) * t222 + (Ifges(4,1) * t227 + Ifges(4,4) * t232) * t220;
t143 = Ifges(4,6) * t222 + (Ifges(4,4) * t227 + Ifges(4,2) * t232) * t220;
t140 = -t225 * t267 + t173;
t135 = -mrSges(7,1) * t231 - mrSges(7,3) * t156;
t134 = mrSges(7,2) * t231 - mrSges(7,3) * t155;
t131 = -mrSges(5,1) * t259 - mrSges(5,3) * t163;
t130 = mrSges(5,2) * t259 - mrSges(5,3) * t162;
t120 = -mrSges(7,1) * t176 + mrSges(7,2) * t177;
t106 = mrSges(5,1) * t162 + mrSges(5,2) * t163;
t105 = mrSges(7,1) * t155 + mrSges(7,2) * t156;
t96 = mrSges(6,1) * t162 - mrSges(6,3) * t125;
t95 = -mrSges(6,2) * t162 + mrSges(6,3) * t124;
t94 = mrSges(4,1) * t161 - mrSges(4,3) * t116;
t93 = -mrSges(4,2) * t161 - mrSges(4,3) * t115;
t79 = -mrSges(6,1) * t124 + mrSges(6,2) * t125;
t76 = mrSges(4,1) * t115 + mrSges(4,2) * t116;
t67 = Ifges(4,1) * t116 - Ifges(4,4) * t115 + Ifges(4,5) * t161;
t66 = Ifges(4,4) * t116 - Ifges(4,2) * t115 + Ifges(4,6) * t161;
t64 = -pkin(5) * t124 + t91;
t63 = mrSges(7,1) * t162 - mrSges(7,3) * t78;
t62 = -mrSges(7,2) * t162 + mrSges(7,3) * t77;
t61 = mrSges(5,1) * t115 - mrSges(5,3) * t89;
t60 = -mrSges(5,2) * t115 - mrSges(5,3) * t88;
t53 = mrSges(5,1) * t88 + mrSges(5,2) * t89;
t44 = -mrSges(7,1) * t77 + mrSges(7,2) * t78;
t39 = mrSges(6,1) * t88 - mrSges(6,3) * t59;
t38 = -mrSges(6,2) * t88 + mrSges(6,3) * t58;
t37 = Ifges(7,1) * t78 + Ifges(7,4) * t77 + Ifges(7,5) * t162;
t36 = Ifges(7,4) * t78 + Ifges(7,2) * t77 + Ifges(7,6) * t162;
t33 = -mrSges(6,1) * t58 + mrSges(6,2) * t59;
t22 = mrSges(7,1) * t88 - mrSges(7,3) * t32;
t21 = -mrSges(7,2) * t88 + mrSges(7,3) * t31;
t12 = -pkin(5) * t58 + t17;
t11 = -mrSges(7,1) * t31 + mrSges(7,2) * t32;
t10 = Ifges(7,1) * t32 + Ifges(7,4) * t31 + Ifges(7,5) * t88;
t9 = Ifges(7,4) * t32 + Ifges(7,2) * t31 + Ifges(7,6) * t88;
t1 = [((Ifges(3,5) * t228 + Ifges(3,6) * t233) * t223 + 0.2e1 * (-t166 * t228 + t168 * t233) * mrSges(3,3) + (m(3) * pkin(1) ^ 2 - 0.2e1 * pkin(1) * (-mrSges(3,1) * t233 + mrSges(3,2) * t228) + t228 * (Ifges(3,1) * t228 + Ifges(3,4) * t233) + t233 * (Ifges(3,4) * t228 + Ifges(3,2) * t233)) * t221) * t221 + (t40 - t66) * t115 + m(7) * (t12 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t17 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2 + t49 ^ 2) + m(4) * (t54 ^ 2 + t55 ^ 2 + t80 ^ 2) + (t8 + t23 - t41) * t88 + m(3) * (t166 ^ 2 + t168 ^ 2) + Ifges(2,3) + (t249 - 0.2e1 * t261 + 0.2e1 * t262) * t223 + 0.2e1 * t12 * t11 + 0.2e1 * t3 * t21 + 0.2e1 * t2 * t22 + t31 * t9 + t32 * t10 + 0.2e1 * t17 * t33 + 0.2e1 * t7 * t38 + 0.2e1 * t6 * t39 + 0.2e1 * t49 * t53 + t58 * t24 + t59 * t25 + 0.2e1 * t20 * t60 + 0.2e1 * t19 * t61 + 0.2e1 * t80 * t76 + t89 * t42 + 0.2e1 * t55 * t93 + 0.2e1 * t54 * t94 + t116 * t67 + t161 * t65; (t102 / 0.2e1 - t143 / 0.2e1) * t115 + m(7) * (t12 * t64 + t14 * t2 + t15 * t3) + m(6) * (t17 * t91 + t51 * t6 + t52 * t7) + m(5) * (t146 * t49 + t19 * t97 + t20 * t98) + m(4) * (-pkin(2) * t220 * t80 + t165 * t54 + t167 * t55) + t246 * t88 + t247 * t162 - t261 + t262 + t249 + t89 * t284 + t59 * t288 + t58 * t289 + t163 * t290 + t125 * t291 + t124 * t292 + (-pkin(2) * t76 + t227 * t67 / 0.2e1 + (-t40 / 0.2e1 + t66 / 0.2e1) * t232) * t220 + t222 * t65 / 0.2e1 + t15 * t21 + t14 * t22 + t31 * t36 / 0.2e1 + t32 * t37 / 0.2e1 + t12 * t44 + t51 * t39 + t52 * t38 + t3 * t62 + t2 * t63 + t64 * t11 + t77 * t9 / 0.2e1 + t78 * t10 / 0.2e1 + t17 * t79 + t91 * t33 + t7 * t95 + t6 * t96 + t97 * t61 + t98 * t60 + t49 * t106 + t20 * t130 + t19 * t131 + t116 * t144 / 0.2e1 + t146 * t53 + t161 * t142 / 0.2e1 + t80 * t164 + t165 * t94 + t167 * t93 + t54 * t174 + t55 * t175; t222 * t142 + Ifges(3,3) + 0.2e1 * t15 * t62 + 0.2e1 * t14 * t63 + 0.2e1 * t64 * t44 + t77 * t36 + t78 * t37 + 0.2e1 * t91 * t79 + 0.2e1 * t52 * t95 + 0.2e1 * t51 * t96 + t124 * t69 + t125 * t70 + 0.2e1 * t98 * t130 + 0.2e1 * t97 * t131 + 0.2e1 * t146 * t106 + t163 * t104 + 0.2e1 * t165 * t174 + 0.2e1 * t167 * t175 + (t35 + t68 - t103) * t162 + (-0.2e1 * pkin(2) * t164 + t227 * t144 + (-t102 + t143) * t232) * t220 + m(7) * (t14 ^ 2 + t15 ^ 2 + t64 ^ 2) + m(6) * (t51 ^ 2 + t52 ^ 2 + t91 ^ 2) + m(5) * (t146 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(4) * (pkin(2) ^ 2 * t220 ^ 2 + t165 ^ 2 + t167 ^ 2); t65 + t244 * t88 + (t20 * mrSges(5,3) + pkin(11) * t60 - t247) * t231 + t32 * t285 + t31 * t286 + t89 * t273 + t10 * t278 + t9 * t279 + t59 * t280 + t58 * t281 + t181 * t11 + t49 * t184 + t115 * t186 / 0.2e1 + m(7) * (t12 * t181 + t2 * t72 + t3 * t73) + m(6) * (t140 * t6 + t141 * t7 + t17 * t268) + (t290 - t19 * mrSges(5,3) + t24 * t272 + t25 * t271 + (t33 - t61) * pkin(11)) * t226 + m(5) * (-pkin(3) * t49 + (-t19 * t226 + t20 * t231) * pkin(11)) - pkin(3) * t53 + t54 * mrSges(4,1) - t55 * mrSges(4,2) + t72 * t22 + t73 * t21 + t12 * t105 + t3 * t134 + t2 * t135 + t140 * t39 + t141 * t38 + t17 * t169 + t7 * t178 + t6 * t179; t142 + t244 * t162 + (t98 * mrSges(5,3) + pkin(11) * t130 - t246) * t231 + m(7) * (t14 * t72 + t15 * t73 + t181 * t64) + t78 * t285 + t77 * t286 + t163 * t273 + t37 * t278 + t36 * t279 + t125 * t280 + t124 * t281 + t181 * t44 + t146 * t184 - t186 * t259 / 0.2e1 + m(6) * (t140 * t51 + t141 * t52 + t268 * t91) + (t284 - t97 * mrSges(5,3) + t69 * t272 + t70 * t271 + (t79 - t131) * pkin(11)) * t226 + m(5) * (-pkin(3) * t146 + (-t226 * t97 + t231 * t98) * pkin(11)) + t72 * t63 + t73 * t62 + t64 * t105 - pkin(3) * t106 + t15 * t134 + t14 * t135 + t140 * t96 + t141 * t95 + t165 * mrSges(4,1) - t167 * mrSges(4,2) + t91 * t169 + t52 * t178 + t51 * t179; -0.2e1 * pkin(3) * t184 - t155 * t100 + t156 * t101 + 0.2e1 * t181 * t105 + 0.2e1 * t73 * t134 + 0.2e1 * t72 * t135 + 0.2e1 * t140 * t179 + 0.2e1 * t141 * t178 + Ifges(4,3) + (t217 + t219) * mrSges(5,3) * t293 + (-t99 - t151 + t188) * t231 + m(7) * (t181 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(6) * (t140 ^ 2 + t141 ^ 2 + t215) + m(5) * (pkin(3) ^ 2 + t219 * t235 + t215) + (-t152 * t225 + t153 * t230 + t169 * t293 + t190) * t226; m(7) * (t12 * t208 + t132 * t2 + t133 * t3) + (t7 * mrSges(6,3) + pkin(12) * t38 + t292) * t230 + (t176 * t3 - t177 * t2) * mrSges(7,3) + t40 + (-t6 * mrSges(6,3) - pkin(12) * t39 + t291) * t225 + t208 * t11 + t17 * t183 + t58 * t275 + t59 * t274 + t248 * t88 + m(6) * (-pkin(4) * t17 + (-t225 * t6 + t230 * t7) * pkin(12)) + t19 * mrSges(5,1) - t20 * mrSges(5,2) - pkin(4) * t33 + t12 * t120 + t31 * t283 + t32 * t282 + t132 * t22 + t133 * t21 + t9 * t277 + t10 * t276; (t52 * mrSges(6,3) + pkin(12) * t95 + t289) * t230 + (-t51 * mrSges(6,3) - pkin(12) * t96 + t288) * t225 + t102 + m(7) * (t132 * t14 + t133 * t15 + t208 * t64) + t248 * t162 + (-t14 * t177 + t15 * t176) * mrSges(7,3) + t208 * t44 + t91 * t183 + t124 * t275 + t125 * t274 + m(6) * (-pkin(4) * t91 + (-t225 * t51 + t230 * t52) * pkin(12)) - pkin(4) * t79 + t97 * mrSges(5,1) - t98 * mrSges(5,2) + t64 * t120 + t77 * t283 + t78 * t282 + t132 * t63 + t133 * t62 + t36 * t277 + t37 * t276; t208 * t105 + t133 * t134 + t132 * t135 + t122 * t279 + t123 * t278 - pkin(4) * t169 + t100 * t277 + t101 * t276 + t181 * t120 + (-mrSges(5,1) + t183) * t268 + (t176 * t73 - t177 * t72) * mrSges(7,3) + (t141 * mrSges(6,3) + pkin(12) * t178 + t226 * t274 + t281) * t230 + (t280 - t140 * mrSges(6,3) - pkin(12) * t179 - t226 * t187 / 0.2e1) * t225 + m(7) * (t132 * t72 + t133 * t73 + t181 * t208) + m(6) * (-pkin(4) * t268 + (-t140 * t225 + t141 * t230) * pkin(12)) + (-pkin(11) * mrSges(5,2) - t248) * t231 + t186; -0.2e1 * pkin(4) * t183 + 0.2e1 * t208 * t120 + t176 * t122 + t177 * t123 + t230 * t187 + t225 * t189 + Ifges(5,3) + m(7) * (t132 ^ 2 + t133 ^ 2 + t208 ^ 2) + m(6) * (pkin(12) ^ 2 * t251 + pkin(4) ^ 2) + 0.2e1 * (-t132 * t177 + t133 * t176) * mrSges(7,3) + 0.2e1 * t251 * pkin(12) * mrSges(6,3); t6 * mrSges(6,1) - t7 * mrSges(6,2) + (t224 * t21 + t229 * t22 + m(7) * (t2 * t229 + t224 * t3)) * pkin(5) + t239 + t23; t51 * mrSges(6,1) - t52 * mrSges(6,2) + (m(7) * (t14 * t229 + t15 * t224) + t224 * t62 + t229 * t63) * pkin(5) + t238 + t68; t140 * mrSges(6,1) - t141 * mrSges(6,2) + t266 * t231 + (m(7) * (t224 * t73 + t229 * t72) + t224 * t134 + t229 * t135) * pkin(5) + t241 + t245; -t243 * pkin(12) + (m(7) * (t132 * t229 + t133 * t224) + (t176 * t224 - t177 * t229) * mrSges(7,3)) * pkin(5) + t240 + t185; m(7) * (t224 ^ 2 + t229 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t237 - t266; t239; t238; t241 - t263; t240; Ifges(7,3) + t237; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
