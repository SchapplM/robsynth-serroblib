% Calculate joint inertia matrix for
% S6RRRRRR8
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:48:24
% EndTime: 2019-03-10 04:48:32
% DurationCPUTime: 3.68s
% Computational Cost: add. (9322->549), mult. (23859->797), div. (0->0), fcn. (27433->14), ass. (0->212)
t229 = sin(qJ(6));
t234 = cos(qJ(6));
t230 = sin(qJ(5));
t231 = sin(qJ(4));
t235 = cos(qJ(5));
t236 = cos(qJ(4));
t186 = t230 * t231 - t235 * t236;
t187 = t230 * t236 + t231 * t235;
t214 = -pkin(4) * t236 - pkin(3);
t130 = pkin(5) * t186 - pkin(13) * t187 + t214;
t305 = -pkin(12) - pkin(11);
t199 = t305 * t236;
t257 = t305 * t231;
t154 = -t235 * t199 + t230 * t257;
t87 = t130 * t234 - t154 * t229;
t88 = t130 * t229 + t154 * t234;
t314 = -t229 * t87 + t234 * t88;
t225 = sin(pkin(7));
t237 = cos(qJ(3));
t272 = t225 * t237;
t227 = cos(pkin(7));
t232 = sin(qJ(3));
t269 = t227 * t232;
t177 = pkin(2) * t269 + pkin(10) * t272;
t164 = pkin(11) * t227 + t177;
t165 = (-pkin(3) * t237 - pkin(11) * t232 - pkin(2)) * t225;
t113 = -t164 * t231 + t236 * t165;
t273 = t225 * t232;
t173 = t227 * t231 + t236 * t273;
t85 = -pkin(4) * t272 - pkin(12) * t173 + t113;
t114 = t236 * t164 + t231 * t165;
t172 = t227 * t236 - t231 * t273;
t90 = pkin(12) * t172 + t114;
t54 = t230 * t85 + t235 * t90;
t45 = -pkin(13) * t272 + t54;
t121 = -t235 * t172 + t173 * t230;
t122 = t172 * t230 + t173 * t235;
t204 = pkin(10) * t273;
t290 = pkin(2) * t237;
t163 = t204 + (-pkin(3) - t290) * t227;
t123 = -pkin(4) * t172 + t163;
t62 = pkin(5) * t121 - pkin(13) * t122 + t123;
t24 = -t229 * t45 + t234 * t62;
t25 = t229 * t62 + t234 * t45;
t249 = -t229 * t24 + t234 * t25;
t228 = cos(pkin(6));
t226 = sin(pkin(6));
t238 = cos(qJ(2));
t270 = t226 * t238;
t171 = -t225 * t270 + t227 * t228;
t233 = sin(qJ(2));
t291 = pkin(1) * t228;
t178 = pkin(9) * t270 + t233 * t291;
t268 = t227 * t238;
t256 = t226 * t268;
t131 = (t225 * t228 + t256) * pkin(10) + t178;
t209 = t238 * t291;
t271 = t226 * t233;
t140 = pkin(2) * t228 + t209 + (-pkin(10) * t227 - pkin(9)) * t271;
t159 = (-pkin(10) * t225 * t233 - pkin(2) * t238 - pkin(1)) * t226;
t70 = -t232 * t131 + (t140 * t227 + t159 * t225) * t237;
t60 = -pkin(3) * t171 - t70;
t139 = t228 * t273 + (t232 * t268 + t233 * t237) * t226;
t97 = -t139 * t231 + t171 * t236;
t32 = -pkin(4) * t97 + t60;
t98 = t139 * t236 + t171 * t231;
t66 = t230 * t98 - t235 * t97;
t67 = t230 * t97 + t235 * t98;
t15 = pkin(5) * t66 - pkin(13) * t67 + t32;
t138 = -t228 * t272 + t232 * t271 - t237 * t256;
t93 = -t140 * t225 + t227 * t159;
t58 = pkin(3) * t138 - pkin(11) * t139 + t93;
t71 = t237 * t131 + t140 * t269 + t159 * t273;
t61 = pkin(11) * t171 + t71;
t29 = -t231 * t61 + t236 * t58;
t17 = pkin(4) * t138 - pkin(12) * t98 + t29;
t30 = t231 * t58 + t236 * t61;
t20 = pkin(12) * t97 + t30;
t9 = t230 * t17 + t235 * t20;
t6 = pkin(13) * t138 + t9;
t2 = t15 * t234 - t229 * t6;
t3 = t15 * t229 + t234 * t6;
t251 = -t2 * t229 + t234 * t3;
t250 = mrSges(7,1) * t229 + mrSges(7,2) * t234;
t129 = t250 * t187;
t152 = -t199 * t230 - t235 * t257;
t313 = m(7) * t152 + t129;
t277 = t187 * t229;
t136 = -mrSges(7,2) * t186 - mrSges(7,3) * t277;
t276 = t187 * t234;
t137 = mrSges(7,1) * t186 - mrSges(7,3) * t276;
t312 = m(7) * t314 + t234 * t136 - t229 * t137;
t311 = t152 ^ 2;
t310 = 0.2e1 * t152;
t28 = Ifges(6,1) * t67 - Ifges(6,4) * t66 + Ifges(6,5) * t138;
t309 = t28 / 0.2e1;
t50 = Ifges(5,4) * t98 + Ifges(5,2) * t97 + Ifges(5,6) * t138;
t308 = t50 / 0.2e1;
t51 = Ifges(5,1) * t98 + Ifges(5,4) * t97 + Ifges(5,5) * t138;
t307 = t51 / 0.2e1;
t76 = Ifges(6,1) * t122 - Ifges(6,4) * t121 - Ifges(6,5) * t272;
t306 = t76 / 0.2e1;
t284 = Ifges(7,4) * t234;
t109 = Ifges(7,6) * t186 + (-Ifges(7,2) * t229 + t284) * t187;
t304 = t109 / 0.2e1;
t285 = Ifges(7,4) * t229;
t110 = Ifges(7,5) * t186 + (Ifges(7,1) * t234 - t285) * t187;
t303 = t110 / 0.2e1;
t116 = Ifges(5,4) * t173 + Ifges(5,2) * t172 - Ifges(5,6) * t272;
t302 = t116 / 0.2e1;
t117 = Ifges(5,1) * t173 + Ifges(5,4) * t172 - Ifges(5,5) * t272;
t301 = t117 / 0.2e1;
t144 = Ifges(6,1) * t187 - Ifges(6,4) * t186;
t300 = t144 / 0.2e1;
t192 = Ifges(7,5) * t229 + Ifges(7,6) * t234;
t299 = t192 / 0.2e1;
t194 = Ifges(7,2) * t234 + t285;
t298 = t194 / 0.2e1;
t195 = Ifges(5,4) * t231 + Ifges(5,2) * t236;
t297 = t195 / 0.2e1;
t196 = Ifges(7,1) * t229 + t284;
t296 = t196 / 0.2e1;
t197 = Ifges(5,1) * t231 + Ifges(5,4) * t236;
t295 = t197 / 0.2e1;
t294 = -t229 / 0.2e1;
t293 = t229 / 0.2e1;
t292 = t234 / 0.2e1;
t289 = pkin(13) * t229;
t288 = pkin(13) * t234;
t176 = -pkin(9) * t271 + t209;
t283 = t176 * mrSges(3,1);
t282 = t178 * mrSges(3,2);
t212 = pkin(4) * t230 + pkin(13);
t275 = t212 * t229;
t274 = t212 * t234;
t265 = Ifges(6,5) * t122 - Ifges(6,6) * t121;
t264 = Ifges(5,5) * t173 + Ifges(5,6) * t172;
t142 = Ifges(6,5) * t187 - Ifges(6,6) * t186;
t193 = Ifges(5,5) * t231 + Ifges(5,6) * t236;
t263 = t229 ^ 2 + t234 ^ 2;
t262 = t231 ^ 2 + t236 ^ 2;
t38 = t138 * t234 - t229 * t67;
t39 = t138 * t229 + t234 * t67;
t12 = Ifges(7,5) * t39 + Ifges(7,6) * t38 + Ifges(7,3) * t66;
t26 = Ifges(6,5) * t67 - Ifges(6,6) * t66 + Ifges(6,3) * t138;
t49 = Ifges(5,5) * t98 + Ifges(5,6) * t97 + Ifges(5,3) * t138;
t261 = Ifges(6,3) * t272;
t27 = Ifges(6,4) * t67 - Ifges(6,2) * t66 + Ifges(6,6) * t138;
t260 = t12 / 0.2e1 - t27 / 0.2e1;
t103 = -t122 * t229 - t234 * t272;
t104 = t122 * t234 - t229 * t272;
t46 = Ifges(7,5) * t104 + Ifges(7,6) * t103 + Ifges(7,3) * t121;
t75 = Ifges(6,4) * t122 - Ifges(6,2) * t121 - Ifges(6,6) * t272;
t259 = t46 / 0.2e1 - t75 / 0.2e1;
t258 = t234 * t194 + t229 * t196 + Ifges(6,3);
t79 = Ifges(4,5) * t139 - Ifges(4,6) * t138 + Ifges(4,3) * t171;
t160 = Ifges(4,5) * t273 + Ifges(4,6) * t272 + Ifges(4,3) * t227;
t255 = Ifges(3,5) * t271 + Ifges(3,6) * t270 + Ifges(3,3) * t228;
t108 = Ifges(7,5) * t276 - Ifges(7,6) * t277 + Ifges(7,3) * t186;
t143 = Ifges(6,4) * t187 - Ifges(6,2) * t186;
t254 = -t143 / 0.2e1 + t108 / 0.2e1;
t253 = t193 / 0.2e1 + t142 / 0.2e1;
t252 = t263 * t212;
t8 = t17 * t235 - t20 * t230;
t53 = -t230 * t90 + t235 * t85;
t247 = 0.2e1 * mrSges(7,3) * t263;
t245 = (mrSges(6,1) * t235 - mrSges(6,2) * t230) * pkin(4);
t13 = Ifges(7,4) * t39 + Ifges(7,2) * t38 + Ifges(7,6) * t66;
t14 = Ifges(7,1) * t39 + Ifges(7,4) * t38 + Ifges(7,5) * t66;
t190 = -mrSges(7,1) * t234 + mrSges(7,2) * t229;
t5 = -pkin(5) * t138 - t8;
t244 = t8 * mrSges(6,1) - t9 * mrSges(6,2) + t251 * mrSges(7,3) + t13 * t292 + t14 * t293 + t5 * t190 + t39 * t296 + t38 * t298 + t66 * t299 + t26;
t44 = pkin(5) * t272 - t53;
t47 = Ifges(7,4) * t104 + Ifges(7,2) * t103 + Ifges(7,6) * t121;
t48 = Ifges(7,1) * t104 + Ifges(7,4) * t103 + Ifges(7,5) * t121;
t243 = t53 * mrSges(6,1) - t54 * mrSges(6,2) + t249 * mrSges(7,3) + t103 * t298 + t104 * t296 + t121 * t299 + t44 * t190 + t47 * t292 + t48 * t293 + t265;
t242 = -t154 * mrSges(6,2) + t109 * t292 + t110 * t293 + t142 - t194 * t277 / 0.2e1 + t276 * t296 + t186 * t299 + (-mrSges(6,1) + t190) * t152 + t314 * mrSges(7,3);
t213 = -pkin(4) * t235 - pkin(5);
t191 = -mrSges(5,1) * t236 + mrSges(5,2) * t231;
t183 = -mrSges(4,2) * t227 + mrSges(4,3) * t272;
t182 = mrSges(4,1) * t227 - mrSges(4,3) * t273;
t175 = t227 * t290 - t204;
t174 = (-mrSges(4,1) * t237 + mrSges(4,2) * t232) * t225;
t162 = Ifges(4,5) * t227 + (Ifges(4,1) * t232 + Ifges(4,4) * t237) * t225;
t161 = Ifges(4,6) * t227 + (Ifges(4,4) * t232 + Ifges(4,2) * t237) * t225;
t151 = -mrSges(5,1) * t272 - mrSges(5,3) * t173;
t150 = mrSges(5,2) * t272 + mrSges(5,3) * t172;
t141 = mrSges(6,1) * t186 + mrSges(6,2) * t187;
t124 = -mrSges(5,1) * t172 + mrSges(5,2) * t173;
t115 = -Ifges(5,3) * t272 + t264;
t112 = -mrSges(6,1) * t272 - mrSges(6,3) * t122;
t111 = mrSges(6,2) * t272 - mrSges(6,3) * t121;
t106 = mrSges(4,1) * t171 - mrSges(4,3) * t139;
t105 = -mrSges(4,2) * t171 - mrSges(4,3) * t138;
t89 = mrSges(4,1) * t138 + mrSges(4,2) * t139;
t82 = mrSges(6,1) * t121 + mrSges(6,2) * t122;
t81 = Ifges(4,1) * t139 - Ifges(4,4) * t138 + Ifges(4,5) * t171;
t80 = Ifges(4,4) * t139 - Ifges(4,2) * t138 + Ifges(4,6) * t171;
t78 = mrSges(5,1) * t138 - mrSges(5,3) * t98;
t77 = -mrSges(5,2) * t138 + mrSges(5,3) * t97;
t74 = -t261 + t265;
t73 = mrSges(7,1) * t121 - mrSges(7,3) * t104;
t72 = -mrSges(7,2) * t121 + mrSges(7,3) * t103;
t69 = -mrSges(7,1) * t103 + mrSges(7,2) * t104;
t68 = -mrSges(5,1) * t97 + mrSges(5,2) * t98;
t43 = mrSges(6,1) * t138 - mrSges(6,3) * t67;
t42 = -mrSges(6,2) * t138 - mrSges(6,3) * t66;
t31 = mrSges(6,1) * t66 + mrSges(6,2) * t67;
t23 = mrSges(7,1) * t66 - mrSges(7,3) * t39;
t22 = -mrSges(7,2) * t66 + mrSges(7,3) * t38;
t18 = -mrSges(7,1) * t38 + mrSges(7,2) * t39;
t1 = [m(3) * (t176 ^ 2 + t178 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t5 ^ 2) + m(6) * (t32 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2 + t60 ^ 2) + m(4) * (t70 ^ 2 + t71 ^ 2 + t93 ^ 2) + (t12 - t27) * t66 + (t26 + t49 - t80) * t138 + t171 * t79 + t139 * t81 + t98 * t51 + 0.2e1 * t71 * t105 + 0.2e1 * t70 * t106 + 0.2e1 * t93 * t89 + t97 * t50 + ((Ifges(3,5) * t233 + Ifges(3,6) * t238) * t228 + 0.2e1 * (-t176 * t233 + t178 * t238) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t238 + mrSges(3,2) * t233) + t233 * (Ifges(3,1) * t233 + Ifges(3,4) * t238) + t238 * (Ifges(3,4) * t233 + Ifges(3,2) * t238) + m(3) * pkin(1) ^ 2) * t226) * t226 + Ifges(2,3) + (t255 - 0.2e1 * t282 + 0.2e1 * t283) * t228 + 0.2e1 * t5 * t18 + 0.2e1 * t3 * t22 + 0.2e1 * t2 * t23 + 0.2e1 * t32 * t31 + t38 * t13 + t39 * t14 + 0.2e1 * t9 * t42 + 0.2e1 * t8 * t43 + t67 * t28 + 0.2e1 * t60 * t68 + 0.2e1 * t30 * t77 + 0.2e1 * t29 * t78; t255 + (-t161 / 0.2e1 + t115 / 0.2e1 + t74 / 0.2e1) * t138 - t282 + t283 + t227 * t79 / 0.2e1 + t70 * t182 + t71 * t183 + t171 * t160 / 0.2e1 + t93 * t174 + t175 * t106 + t177 * t105 + t139 * t162 / 0.2e1 + t163 * t68 + t30 * t150 + t29 * t151 + t123 * t31 + t60 * t124 + t9 * t111 + t8 * t112 + t113 * t78 + t114 * t77 + t103 * t13 / 0.2e1 + t104 * t14 / 0.2e1 + t98 * t301 + t97 * t302 + t67 * t306 + t173 * t307 + t172 * t308 + t122 * t309 + t259 * t66 + t260 * t121 + m(7) * (t2 * t24 + t25 * t3 + t44 * t5) + m(6) * (t123 * t32 + t53 * t8 + t54 * t9) + m(5) * (t113 * t29 + t114 * t30 + t163 * t60) + m(4) * (-pkin(2) * t225 * t93 + t175 * t70 + t177 * t71) + t24 * t23 + t25 * t22 + t44 * t18 + t38 * t47 / 0.2e1 + t39 * t48 / 0.2e1 + t53 * t43 + t54 * t42 + (-pkin(2) * t89 + t232 * t81 / 0.2e1 + (-t26 / 0.2e1 - t49 / 0.2e1 + t80 / 0.2e1) * t237) * t225 + t5 * t69 + t3 * t72 + t2 * t73 + t32 * t82; t227 * t160 + 0.2e1 * t175 * t182 + 0.2e1 * t177 * t183 + t172 * t116 + t173 * t117 + 0.2e1 * t163 * t124 + 0.2e1 * t114 * t150 + 0.2e1 * t113 * t151 + 0.2e1 * t123 * t82 + t122 * t76 + 0.2e1 * t54 * t111 + 0.2e1 * t53 * t112 + t103 * t47 + t104 * t48 + Ifges(3,3) + 0.2e1 * t44 * t69 + 0.2e1 * t25 * t72 + 0.2e1 * t24 * t73 + (t46 - t75) * t121 + (-0.2e1 * pkin(2) * t174 + t232 * t162 + (-t115 + t161 - t74) * t237) * t225 + m(7) * (t24 ^ 2 + t25 ^ 2 + t44 ^ 2) + m(6) * (t123 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(5) * (t113 ^ 2 + t114 ^ 2 + t163 ^ 2) + m(4) * (pkin(2) ^ 2 * t225 ^ 2 + t175 ^ 2 + t177 ^ 2); t79 + (t18 - t43) * t152 + m(5) * (-pkin(3) * t60 + (-t231 * t29 + t236 * t30) * pkin(11)) + m(7) * (t152 * t5 + t2 * t87 + t3 * t88) + m(6) * (-t152 * t8 + t154 * t9 + t214 * t32) + t60 * t191 + t214 * t31 + t154 * t42 + t3 * t136 + t2 * t137 + t32 * t141 + t5 * t129 + t253 * t138 + t254 * t66 + t98 * t295 + t97 * t297 + t67 * t300 + t39 * t303 + t38 * t304 + (-t9 * mrSges(6,3) + t260) * t186 + (-t29 * mrSges(5,3) - pkin(11) * t78 + t307) * t231 + (t30 * mrSges(5,3) + pkin(11) * t77 + t308) * t236 + (-t8 * mrSges(6,3) + t13 * t294 + t14 * t292 + t309) * t187 - pkin(3) * t68 + t70 * mrSges(4,1) - t71 * mrSges(4,2) + t87 * t23 + t88 * t22; m(7) * (t152 * t44 + t24 * t87 + t25 * t88) + m(6) * (t123 * t214 - t152 * t53 + t154 * t54) + t160 + m(5) * (-pkin(3) * t163 + (-t113 * t231 + t114 * t236) * pkin(11)) + (t69 - t112) * t152 + t163 * t191 + t214 * t82 + t175 * mrSges(4,1) - t177 * mrSges(4,2) + t154 * t111 + t25 * t136 + t24 * t137 + t123 * t141 - pkin(3) * t124 + t44 * t129 + t254 * t121 + t173 * t295 + t172 * t297 + t122 * t300 + t104 * t303 + t103 * t304 + (-t54 * mrSges(6,3) + t259) * t186 - t253 * t272 + (-t113 * mrSges(5,3) - pkin(11) * t151 + t301) * t231 + (t114 * mrSges(5,3) + pkin(11) * t150 + t302) * t236 + (-t53 * mrSges(6,3) + t292 * t48 + t294 * t47 + t306) * t187 + t87 * t73 + t88 * t72; -0.2e1 * pkin(3) * t191 + t129 * t310 + 0.2e1 * t88 * t136 + 0.2e1 * t87 * t137 + 0.2e1 * t214 * t141 + t236 * t195 + t231 * t197 + Ifges(4,3) + 0.2e1 * t262 * pkin(11) * mrSges(5,3) + (-0.2e1 * mrSges(6,3) * t154 + t108 - t143) * t186 + m(5) * (pkin(11) ^ 2 * t262 + pkin(3) ^ 2) + m(7) * (t87 ^ 2 + t88 ^ 2 + t311) + m(6) * (t154 ^ 2 + t214 ^ 2 + t311) + (mrSges(6,3) * t310 - t109 * t229 + t110 * t234 + t144) * t187; t244 + m(7) * (t212 * t251 + t213 * t5) + (t230 * t42 + t235 * t43 + m(6) * (t230 * t9 + t235 * t8)) * pkin(4) + t213 * t18 - t23 * t275 + t22 * t274 + t29 * mrSges(5,1) - t30 * mrSges(5,2) + t49; t243 + m(7) * (t212 * t249 + t213 * t44) + (t230 * t111 + t235 * t112 + m(6) * (t230 * t54 + t235 * t53)) * pkin(4) + t213 * t69 + t113 * mrSges(5,1) - t114 * mrSges(5,2) + (-Ifges(5,3) - Ifges(6,3)) * t272 - t73 * t275 + t72 * t274 + t264; t242 + (m(6) * (-t152 * t235 + t154 * t230) + (-t186 * t230 - t187 * t235) * mrSges(6,3)) * pkin(4) + (-mrSges(5,1) * t231 - mrSges(5,2) * t236) * pkin(11) + t193 + t313 * t213 + t312 * t212; 0.2e1 * t213 * t190 + Ifges(5,3) + 0.2e1 * t245 + t212 * t247 + m(7) * (t212 ^ 2 * t263 + t213 ^ 2) + m(6) * (t230 ^ 2 + t235 ^ 2) * pkin(4) ^ 2 + t258; t244 + m(7) * (-pkin(5) * t5 + pkin(13) * t251) - t23 * t289 + t22 * t288 - pkin(5) * t18; t243 + m(7) * (-pkin(5) * t44 + pkin(13) * t249) - t73 * t289 + t72 * t288 - t261 - pkin(5) * t69; -pkin(5) * t313 + t312 * pkin(13) + t242; m(7) * (-pkin(5) * t213 + pkin(13) * t252) + (t213 - pkin(5)) * t190 + t245 + (pkin(13) * t263 + t252) * mrSges(7,3) + t258; -0.2e1 * pkin(5) * t190 + m(7) * (pkin(13) ^ 2 * t263 + pkin(5) ^ 2) + pkin(13) * t247 + t258; mrSges(7,1) * t2 - mrSges(7,2) * t3 + t12; mrSges(7,1) * t24 - mrSges(7,2) * t25 + t46; mrSges(7,1) * t87 - mrSges(7,2) * t88 + t108; -t212 * t250 + t192; -pkin(13) * t250 + t192; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
