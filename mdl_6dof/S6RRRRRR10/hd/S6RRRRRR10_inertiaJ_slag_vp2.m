% Calculate joint inertia matrix for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_inertiaJ_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:49:49
% EndTime: 2019-03-10 05:49:59
% DurationCPUTime: 4.32s
% Computational Cost: add. (15302->672), mult. (42027->983), div. (0->0), fcn. (48834->16), ass. (0->248)
t327 = 2 * pkin(13);
t242 = sin(pkin(8));
t245 = cos(pkin(8));
t250 = sin(qJ(4));
t255 = cos(qJ(4));
t244 = sin(pkin(6));
t247 = cos(pkin(6));
t251 = sin(qJ(3));
t252 = sin(qJ(2));
t256 = cos(qJ(3));
t246 = cos(pkin(7));
t257 = cos(qJ(2));
t280 = t246 * t257;
t243 = sin(pkin(7));
t288 = t243 * t256;
t142 = t247 * t288 + (-t251 * t252 + t256 * t280) * t244;
t286 = t244 * t257;
t184 = -t243 * t286 + t246 * t247;
t262 = t142 * t245 + t184 * t242;
t303 = pkin(1) * t247;
t194 = pkin(10) * t286 + t252 * t303;
t135 = (t243 * t247 + t244 * t280) * pkin(11) + t194;
t227 = t257 * t303;
t287 = t244 * t252;
t148 = pkin(2) * t247 + t227 + (-pkin(11) * t246 - pkin(10)) * t287;
t162 = (-pkin(11) * t243 * t252 - pkin(2) * t257 - pkin(1)) * t244;
t282 = t246 * t251;
t289 = t243 * t251;
t78 = t256 * t135 + t148 * t282 + t162 * t289;
t53 = pkin(12) * t262 + t78;
t143 = t247 * t289 + (t251 * t280 + t252 * t256) * t244;
t300 = pkin(12) * t245;
t281 = t246 * t256;
t76 = -t135 * t251 + t148 * t281 + t162 * t288;
t56 = pkin(3) * t184 - t143 * t300 + t76;
t103 = -t148 * t243 + t246 * t162;
t301 = pkin(12) * t242;
t67 = -pkin(3) * t142 - t143 * t301 + t103;
t19 = -t250 * t53 + (t242 * t67 + t245 * t56) * t255;
t193 = pkin(2) * t282 + pkin(11) * t288;
t283 = t245 * t256;
t272 = t243 * t283;
t134 = (t242 * t246 + t272) * pkin(12) + t193;
t225 = pkin(2) * t281;
t147 = pkin(3) * t246 + t225 + (-pkin(11) - t300) * t289;
t161 = (-pkin(3) * t256 - t251 * t301 - pkin(2)) * t243;
t75 = -t250 * t134 + (t147 * t245 + t161 * t242) * t255;
t109 = -t142 * t242 + t184 * t245;
t249 = sin(qJ(5));
t254 = cos(qJ(5));
t83 = t143 * t255 + t250 * t262;
t54 = -t109 * t254 + t249 * t83;
t55 = t109 * t249 + t254 * t83;
t284 = t245 * t255;
t290 = t242 * t255;
t82 = -t142 * t284 + t143 * t250 - t184 * t290;
t25 = Ifges(6,1) * t55 - Ifges(6,4) * t54 + Ifges(6,5) * t82;
t326 = t25 / 0.2e1;
t291 = t242 * t250;
t144 = t246 * t291 + (t250 * t283 + t251 * t255) * t243;
t183 = -t242 * t288 + t245 * t246;
t110 = t144 * t249 - t254 * t183;
t111 = t144 * t254 + t183 * t249;
t141 = -t246 * t290 + t250 * t289 - t255 * t272;
t63 = Ifges(6,1) * t111 - Ifges(6,4) * t110 + Ifges(6,5) * t141;
t325 = t63 / 0.2e1;
t186 = t245 * t249 + t254 * t291;
t248 = sin(qJ(6));
t253 = cos(qJ(6));
t149 = -t186 * t248 - t253 * t290;
t150 = t186 * t253 - t248 * t290;
t185 = -t254 * t245 + t249 * t291;
t97 = Ifges(7,4) * t150 + Ifges(7,2) * t149 + Ifges(7,6) * t185;
t324 = t97 / 0.2e1;
t98 = Ifges(7,1) * t150 + Ifges(7,4) * t149 + Ifges(7,5) * t185;
t323 = t98 / 0.2e1;
t322 = t109 / 0.2e1;
t125 = Ifges(6,1) * t186 - Ifges(6,4) * t185 - Ifges(6,5) * t290;
t321 = t125 / 0.2e1;
t320 = t149 / 0.2e1;
t319 = t150 / 0.2e1;
t169 = Ifges(5,5) * t245 + (Ifges(5,1) * t250 + Ifges(5,4) * t255) * t242;
t318 = t169 / 0.2e1;
t294 = Ifges(7,4) * t253;
t176 = -Ifges(7,6) * t254 + (-Ifges(7,2) * t248 + t294) * t249;
t317 = t176 / 0.2e1;
t295 = Ifges(7,4) * t248;
t177 = -Ifges(7,5) * t254 + (Ifges(7,1) * t253 - t295) * t249;
t316 = t177 / 0.2e1;
t315 = t183 / 0.2e1;
t314 = t186 / 0.2e1;
t206 = Ifges(7,5) * t248 + Ifges(7,6) * t253;
t313 = t206 / 0.2e1;
t207 = Ifges(6,5) * t249 + Ifges(6,6) * t254;
t312 = t207 / 0.2e1;
t208 = Ifges(7,2) * t253 + t295;
t311 = t208 / 0.2e1;
t210 = Ifges(7,1) * t248 + t294;
t310 = t210 / 0.2e1;
t211 = Ifges(6,1) * t249 + Ifges(6,4) * t254;
t309 = t211 / 0.2e1;
t308 = t245 / 0.2e1;
t307 = -t248 / 0.2e1;
t306 = t248 / 0.2e1;
t305 = t250 / 0.2e1;
t304 = t253 / 0.2e1;
t302 = pkin(3) * t242;
t299 = pkin(13) * t249;
t298 = pkin(13) * t254;
t297 = pkin(14) * t248;
t296 = pkin(14) * t253;
t285 = t245 * t250;
t20 = t255 * t53 + t56 * t285 + t67 * t291;
t16 = pkin(13) * t109 + t20;
t27 = -t242 * t56 + t245 * t67;
t18 = pkin(4) * t82 - pkin(13) * t83 + t27;
t6 = t254 * t16 + t249 * t18;
t102 = -t147 * t242 + t245 * t161;
t68 = pkin(4) * t141 - pkin(13) * t144 + t102;
t77 = t255 * t134 + t147 * t285 + t161 * t291;
t71 = pkin(13) * t183 + t77;
t37 = t249 * t68 + t254 * t71;
t191 = -pkin(10) * t287 + t227;
t293 = t191 * mrSges(3,1);
t292 = t194 * mrSges(3,2);
t279 = t248 * t249;
t278 = t249 * t253;
t192 = pkin(3) * t285 + pkin(12) * t290;
t173 = pkin(13) * t245 + t192;
t174 = (-pkin(4) * t255 - pkin(13) * t250 - pkin(3)) * t242;
t122 = t254 * t173 + t249 * t174;
t277 = t248 ^ 2 + t253 ^ 2;
t32 = -t248 * t55 + t253 * t82;
t33 = t248 * t82 + t253 * t55;
t8 = Ifges(7,5) * t33 + Ifges(7,6) * t32 + Ifges(7,3) * t54;
t23 = Ifges(6,5) * t55 - Ifges(6,6) * t54 + Ifges(6,3) * t82;
t24 = Ifges(6,4) * t55 - Ifges(6,2) * t54 + Ifges(6,6) * t82;
t276 = -t24 / 0.2e1 + t8 / 0.2e1;
t38 = Ifges(5,5) * t83 - Ifges(5,6) * t82 + Ifges(5,3) * t109;
t86 = -t111 * t248 + t141 * t253;
t87 = t111 * t253 + t141 * t248;
t41 = Ifges(7,5) * t87 + Ifges(7,6) * t86 + Ifges(7,3) * t110;
t39 = Ifges(5,4) * t83 - Ifges(5,2) * t82 + Ifges(5,6) * t109;
t275 = t23 / 0.2e1 - t39 / 0.2e1;
t61 = Ifges(6,5) * t111 - Ifges(6,6) * t110 + Ifges(6,3) * t141;
t92 = Ifges(5,4) * t144 - Ifges(5,2) * t141 + Ifges(5,6) * t183;
t274 = t61 / 0.2e1 - t92 / 0.2e1;
t62 = Ifges(6,4) * t111 - Ifges(6,2) * t110 + Ifges(6,6) * t141;
t273 = -t62 / 0.2e1 + t41 / 0.2e1;
t124 = Ifges(6,4) * t186 - Ifges(6,2) * t185 - Ifges(6,6) * t290;
t96 = Ifges(7,5) * t150 + Ifges(7,6) * t149 + Ifges(7,3) * t185;
t271 = t96 / 0.2e1 - t124 / 0.2e1;
t90 = Ifges(5,5) * t144 - Ifges(5,6) * t141 + Ifges(5,3) * t183;
t91 = Ifges(4,5) * t143 + Ifges(4,6) * t142 + Ifges(4,3) * t184;
t165 = Ifges(5,5) * t291 + Ifges(5,6) * t290 + Ifges(5,3) * t245;
t166 = Ifges(4,5) * t289 + Ifges(4,6) * t288 + Ifges(4,3) * t246;
t270 = Ifges(3,5) * t287 + Ifges(3,6) * t286 + Ifges(3,3) * t247;
t123 = Ifges(6,5) * t186 - Ifges(6,6) * t185 - Ifges(6,3) * t290;
t167 = Ifges(5,6) * t245 + (Ifges(5,4) * t250 + Ifges(5,2) * t255) * t242;
t269 = t123 / 0.2e1 - t167 / 0.2e1;
t175 = Ifges(7,5) * t278 - Ifges(7,6) * t279 - Ifges(7,3) * t254;
t209 = Ifges(6,4) * t249 + Ifges(6,2) * t254;
t268 = -t209 / 0.2e1 + t175 / 0.2e1;
t4 = pkin(14) * t82 + t6;
t15 = -pkin(4) * t109 - t19;
t7 = pkin(5) * t54 - pkin(14) * t55 + t15;
t1 = -t248 * t4 + t253 * t7;
t2 = t248 * t7 + t253 * t4;
t267 = -t1 * t248 + t2 * t253;
t266 = mrSges(7,1) * t248 + mrSges(7,2) * t253;
t29 = pkin(14) * t141 + t37;
t70 = -pkin(4) * t183 - t75;
t44 = pkin(5) * t110 - pkin(14) * t111 + t70;
t11 = -t248 * t29 + t253 * t44;
t12 = t248 * t44 + t253 * t29;
t265 = -t11 * t248 + t12 * t253;
t5 = -t16 * t249 + t18 * t254;
t219 = pkin(12) * t291;
t172 = t219 + (-pkin(3) * t255 - pkin(4)) * t245;
t112 = pkin(5) * t185 - pkin(14) * t186 + t172;
t114 = -pkin(14) * t290 + t122;
t72 = t112 * t253 - t114 * t248;
t73 = t112 * t248 + t114 * t253;
t263 = -t248 * t72 + t253 * t73;
t36 = -t249 * t71 + t254 * t68;
t203 = -pkin(5) * t254 - pkin(14) * t249 - pkin(4);
t163 = t203 * t253 - t248 * t298;
t164 = t203 * t248 + t253 * t298;
t260 = -t163 * t248 + t164 * t253;
t121 = -t173 * t249 + t174 * t254;
t259 = pkin(13) ^ 2;
t241 = t254 ^ 2;
t239 = t249 ^ 2;
t237 = t239 * t259;
t205 = -mrSges(6,1) * t254 + mrSges(6,2) * t249;
t204 = -mrSges(7,1) * t253 + mrSges(7,2) * t248;
t201 = -mrSges(7,1) * t254 - mrSges(7,3) * t278;
t200 = mrSges(7,2) * t254 - mrSges(7,3) * t279;
t199 = -mrSges(4,2) * t246 + mrSges(4,3) * t288;
t198 = mrSges(4,1) * t246 - mrSges(4,3) * t289;
t197 = -mrSges(5,2) * t245 + mrSges(5,3) * t290;
t196 = mrSges(5,1) * t245 - mrSges(5,3) * t291;
t195 = t266 * t249;
t190 = -pkin(11) * t289 + t225;
t189 = pkin(3) * t284 - t219;
t188 = (-mrSges(4,1) * t256 + mrSges(4,2) * t251) * t243;
t187 = (-mrSges(5,1) * t255 + mrSges(5,2) * t250) * t242;
t170 = Ifges(4,5) * t246 + (Ifges(4,1) * t251 + Ifges(4,4) * t256) * t243;
t168 = Ifges(4,6) * t246 + (Ifges(4,4) * t251 + Ifges(4,2) * t256) * t243;
t157 = -mrSges(6,1) * t290 - mrSges(6,3) * t186;
t156 = mrSges(6,2) * t290 - mrSges(6,3) * t185;
t127 = mrSges(6,1) * t185 + mrSges(6,2) * t186;
t120 = mrSges(7,1) * t185 - mrSges(7,3) * t150;
t119 = -mrSges(7,2) * t185 + mrSges(7,3) * t149;
t118 = mrSges(5,1) * t183 - mrSges(5,3) * t144;
t117 = mrSges(4,1) * t184 - mrSges(4,3) * t143;
t116 = -mrSges(4,2) * t184 + mrSges(4,3) * t142;
t115 = -mrSges(5,2) * t183 - mrSges(5,3) * t141;
t113 = pkin(5) * t290 - t121;
t101 = -mrSges(7,1) * t149 + mrSges(7,2) * t150;
t100 = mrSges(5,1) * t141 + mrSges(5,2) * t144;
t99 = -mrSges(4,1) * t142 + mrSges(4,2) * t143;
t95 = Ifges(4,1) * t143 + Ifges(4,4) * t142 + Ifges(4,5) * t184;
t94 = Ifges(5,1) * t144 - Ifges(5,4) * t141 + Ifges(5,5) * t183;
t93 = Ifges(4,4) * t143 + Ifges(4,2) * t142 + Ifges(4,6) * t184;
t89 = mrSges(6,1) * t141 - mrSges(6,3) * t111;
t88 = -mrSges(6,2) * t141 - mrSges(6,3) * t110;
t74 = mrSges(6,1) * t110 + mrSges(6,2) * t111;
t60 = mrSges(7,1) * t110 - mrSges(7,3) * t87;
t59 = -mrSges(7,2) * t110 + mrSges(7,3) * t86;
t58 = mrSges(5,1) * t109 - mrSges(5,3) * t83;
t57 = -mrSges(5,2) * t109 - mrSges(5,3) * t82;
t46 = -mrSges(7,1) * t86 + mrSges(7,2) * t87;
t45 = mrSges(5,1) * t82 + mrSges(5,2) * t83;
t43 = Ifges(7,1) * t87 + Ifges(7,4) * t86 + Ifges(7,5) * t110;
t42 = Ifges(7,4) * t87 + Ifges(7,2) * t86 + Ifges(7,6) * t110;
t40 = Ifges(5,1) * t83 - Ifges(5,4) * t82 + Ifges(5,5) * t109;
t35 = mrSges(6,1) * t82 - mrSges(6,3) * t55;
t34 = -mrSges(6,2) * t82 - mrSges(6,3) * t54;
t28 = -pkin(5) * t141 - t36;
t26 = mrSges(6,1) * t54 + mrSges(6,2) * t55;
t22 = mrSges(7,1) * t54 - mrSges(7,3) * t33;
t21 = -mrSges(7,2) * t54 + mrSges(7,3) * t32;
t13 = -mrSges(7,1) * t32 + mrSges(7,2) * t33;
t10 = Ifges(7,1) * t33 + Ifges(7,4) * t32 + Ifges(7,5) * t54;
t9 = Ifges(7,4) * t33 + Ifges(7,2) * t32 + Ifges(7,6) * t54;
t3 = -pkin(5) * t82 - t5;
t14 = [(t23 - t39) * t82 + (t270 - 0.2e1 * t292 + 0.2e1 * t293) * t247 + (t8 - t24) * t54 + ((Ifges(3,5) * t252 + Ifges(3,6) * t257) * t247 + 0.2e1 * (-t191 * t252 + t194 * t257) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t257 + mrSges(3,2) * t252) + t252 * (Ifges(3,1) * t252 + Ifges(3,4) * t257) + t257 * (Ifges(3,4) * t252 + Ifges(3,2) * t257) + m(3) * pkin(1) ^ 2) * t244) * t244 + m(3) * (t191 ^ 2 + t194 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t15 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2 + t27 ^ 2) + m(4) * (t103 ^ 2 + t76 ^ 2 + t78 ^ 2) + t184 * t91 + Ifges(2,3) + 0.2e1 * t3 * t13 + 0.2e1 * t2 * t21 + 0.2e1 * t1 * t22 + 0.2e1 * t15 * t26 + t32 * t9 + t33 * t10 + 0.2e1 * t6 * t34 + 0.2e1 * t5 * t35 + 0.2e1 * t27 * t45 + t55 * t25 + 0.2e1 * t20 * t57 + 0.2e1 * t19 * t58 + t83 * t40 + 0.2e1 * t103 * t99 + t109 * t38 + 0.2e1 * t78 * t116 + 0.2e1 * t76 * t117 + t142 * t93 + t143 * t95; t276 * t110 + t270 + t273 * t54 + t274 * t82 + t275 * t141 + t38 * t315 + t90 * t322 + t55 * t325 + t111 * t326 + t246 * t91 / 0.2e1 + t193 * t116 + t76 * t198 + t78 * t199 + t103 * t188 + t190 * t117 + t184 * t166 / 0.2e1 - t292 + t293 + m(7) * (t1 * t11 + t12 * t2 + t28 * t3) + m(6) * (t15 * t70 + t36 * t5 + t37 * t6) + m(5) * (t102 * t27 + t19 * t75 + t20 * t77) + m(4) * (-pkin(2) * t103 * t243 + t190 * t76 + t193 * t78) + (-pkin(2) * t99 + t251 * t95 / 0.2e1 + t256 * t93 / 0.2e1) * t243 + t12 * t21 + t11 * t22 + t28 * t13 + t36 * t35 + t37 * t34 + t32 * t42 / 0.2e1 + t33 * t43 / 0.2e1 + t3 * t46 + t2 * t59 + t1 * t60 + t70 * t26 + t15 * t74 + t75 * t58 + t77 * t57 + t86 * t9 / 0.2e1 + t87 * t10 / 0.2e1 + t6 * t88 + t5 * t89 + t83 * t94 / 0.2e1 + t27 * t100 + t102 * t45 + t20 * t115 + t19 * t118 + t144 * t40 / 0.2e1 + t142 * t168 / 0.2e1 + t143 * t170 / 0.2e1; m(5) * (t102 ^ 2 + t75 ^ 2 + t77 ^ 2) + m(4) * (pkin(2) ^ 2 * t243 ^ 2 + t190 ^ 2 + t193 ^ 2) + (t61 - t92) * t141 + (-0.2e1 * pkin(2) * t188 + t256 * t168 + t251 * t170) * t243 + m(7) * (t11 ^ 2 + t12 ^ 2 + t28 ^ 2) + m(6) * (t36 ^ 2 + t37 ^ 2 + t70 ^ 2) + (t41 - t62) * t110 + t246 * t166 + 0.2e1 * t190 * t198 + 0.2e1 * t193 * t199 + t183 * t90 + Ifges(3,3) + 0.2e1 * t28 * t46 + 0.2e1 * t12 * t59 + 0.2e1 * t11 * t60 + 0.2e1 * t70 * t74 + t86 * t42 + t87 * t43 + 0.2e1 * t37 * t88 + 0.2e1 * t36 * t89 + 0.2e1 * t102 * t100 + t111 * t63 + 0.2e1 * t77 * t115 + 0.2e1 * t75 * t118 + t144 * t94; t91 + m(7) * (t1 * t72 + t113 * t3 + t2 * t73) + m(6) * (t121 * t5 + t122 * t6 + t15 * t172) + (-pkin(3) * t45 - t255 * t275 + t305 * t40) * t242 + t271 * t54 + t276 * t185 + t25 * t314 + t83 * t318 + t10 * t319 + t9 * t320 + t55 * t321 + t165 * t322 + t33 * t323 + t32 * t324 + t269 * t82 + t192 * t57 + t19 * t196 + t20 * t197 + t27 * t187 + t189 * t58 + t172 * t26 + t38 * t308 + m(5) * (t189 * t19 + t192 * t20 - t27 * t302) + t72 * t22 + t73 * t21 + t76 * mrSges(4,1) - t78 * mrSges(4,2) + t3 * t101 + t113 * t13 + t2 * t119 + t1 * t120 + t121 * t35 + t122 * t34 + t15 * t127 + t6 * t156 + t5 * t157; (-pkin(3) * t100 - t255 * t274 + t305 * t94) * t242 + t271 * t110 + t273 * t185 + t63 * t314 + t165 * t315 + t144 * t318 + t43 * t319 + t42 * t320 + t111 * t321 + t87 * t323 + t86 * t324 + m(7) * (t11 * t72 + t113 * t28 + t12 * t73) + m(6) * (t121 * t36 + t122 * t37 + t172 * t70) + t269 * t141 + t166 - t193 * mrSges(4,2) + t75 * t196 + t77 * t197 + t102 * t187 + t189 * t118 + t190 * mrSges(4,1) + t192 * t115 + t172 * t74 + t90 * t308 + m(5) * (-t102 * t302 + t189 * t75 + t192 * t77) + t72 * t60 + t73 * t59 + t28 * t101 + t113 * t46 + t12 * t119 + t11 * t120 + t121 * t89 + t122 * t88 + t70 * t127 + t37 * t156 + t36 * t157; 0.2e1 * t113 * t101 + 0.2e1 * t73 * t119 + 0.2e1 * t72 * t120 + 0.2e1 * t121 * t157 + 0.2e1 * t122 * t156 + t186 * t125 + 0.2e1 * t172 * t127 + t149 * t97 + t150 * t98 + t245 * t165 + 0.2e1 * t189 * t196 + 0.2e1 * t192 * t197 + Ifges(4,3) + (-t124 + t96) * t185 + (-0.2e1 * pkin(3) * t187 + t250 * t169 + (-t123 + t167) * t255) * t242 + m(6) * (t121 ^ 2 + t122 ^ 2 + t172 ^ 2) + m(5) * (pkin(3) ^ 2 * t242 ^ 2 + t189 ^ 2 + t192 ^ 2) + m(7) * (t113 ^ 2 + t72 ^ 2 + t73 ^ 2); t38 + m(6) * (-pkin(4) * t15 + (-t249 * t5 + t254 * t6) * pkin(13)) + t268 * t54 + t15 * t205 + t82 * t312 + t55 * t309 + t3 * t195 + t2 * t200 + t1 * t201 + t32 * t317 + t33 * t316 + m(7) * (t1 * t163 + t164 * t2 + t299 * t3) + (t6 * mrSges(6,3) + pkin(13) * t34 - t276) * t254 + (t326 - t5 * mrSges(6,3) + t9 * t307 + t10 * t304 + (t13 - t35) * pkin(13)) * t249 + t19 * mrSges(5,1) - t20 * mrSges(5,2) - pkin(4) * t26 + t163 * t22 + t164 * t21; t90 + m(6) * (-pkin(4) * t70 + (-t249 * t36 + t254 * t37) * pkin(13)) + t268 * t110 + (t325 - t36 * mrSges(6,3) + t42 * t307 + t43 * t304 + (t46 - t89) * pkin(13)) * t249 + (t37 * mrSges(6,3) + pkin(13) * t88 - t273) * t254 + t70 * t205 + t141 * t312 + t111 * t309 + t28 * t195 + t12 * t200 + t11 * t201 + t86 * t317 + t87 * t316 - pkin(4) * t74 + t75 * mrSges(5,1) - t77 * mrSges(5,2) + m(7) * (t11 * t163 + t12 * t164 + t28 * t299) + t163 * t60 + t164 * t59; m(6) * (-pkin(4) * t172 + (-t121 * t249 + t122 * t254) * pkin(13)) + m(7) * (t113 * t299 + t163 * t72 + t164 * t73) + (t321 - t121 * mrSges(6,3) + t97 * t307 + t98 * t304 + (t101 - t157) * pkin(13)) * t249 + (t122 * mrSges(6,3) + pkin(13) * t156 - t271) * t254 + t268 * t185 + t165 + t172 * t205 + t186 * t309 + t113 * t195 + t73 * t200 + t72 * t201 + t189 * mrSges(5,1) - t192 * mrSges(5,2) + t149 * t317 + t150 * t316 - t207 * t290 / 0.2e1 - pkin(4) * t127 + t163 * t120 + t164 * t119; -0.2e1 * pkin(4) * t205 + 0.2e1 * t163 * t201 + 0.2e1 * t164 * t200 + Ifges(5,3) + (t209 - t175) * t254 + (t239 + t241) * mrSges(6,3) * t327 + m(6) * (pkin(4) ^ 2 + t241 * t259 + t237) + m(7) * (t163 ^ 2 + t164 ^ 2 + t237) + (-t176 * t248 + t177 * t253 + t195 * t327 + t211) * t249; -pkin(5) * t13 + t3 * t204 + t54 * t313 + t32 * t311 + t33 * t310 - t22 * t297 + t21 * t296 + m(7) * (-pkin(5) * t3 + pkin(14) * t267) + t10 * t306 + t9 * t304 - t6 * mrSges(6,2) + t5 * mrSges(6,1) + t267 * mrSges(7,3) + t23; t59 * t296 + t28 * t204 + t86 * t311 + t110 * t313 + t87 * t310 + m(7) * (-pkin(5) * t28 + pkin(14) * t265) - pkin(5) * t46 + t43 * t306 + t42 * t304 - t60 * t297 - t37 * mrSges(6,2) + t36 * mrSges(6,1) + t265 * mrSges(7,3) + t61; t150 * t310 + m(7) * (-pkin(5) * t113 + pkin(14) * t263) + t119 * t296 - t120 * t297 + t98 * t306 + t97 * t304 + t113 * t204 - pkin(5) * t101 - t122 * mrSges(6,2) + t149 * t311 + t185 * t313 + t121 * mrSges(6,1) + t263 * mrSges(7,3) + t123; t177 * t306 + t176 * t304 - pkin(5) * t195 + (m(7) * t260 + t253 * t200 - t248 * t201) * pkin(14) + (t210 * t304 + t208 * t307 + (-m(7) * pkin(5) - mrSges(6,1) + t204) * pkin(13)) * t249 + (-pkin(13) * mrSges(6,2) - t206 / 0.2e1) * t254 + t260 * mrSges(7,3) + t207; Ifges(6,3) + t248 * t210 + t253 * t208 - 0.2e1 * pkin(5) * t204 + m(7) * (t277 * pkin(14) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t277 * pkin(14) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t8; mrSges(7,1) * t11 - mrSges(7,2) * t12 + t41; mrSges(7,1) * t72 - mrSges(7,2) * t73 + t96; mrSges(7,1) * t163 - mrSges(7,2) * t164 + t175; -pkin(14) * t266 + t206; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
