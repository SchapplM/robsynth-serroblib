% Calculate time derivative of joint inertia matrix for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:45
% EndTime: 2019-03-09 18:20:55
% DurationCPUTime: 4.87s
% Computational Cost: add. (7510->436), mult. (15769->615), div. (0->0), fcn. (14692->8), ass. (0->195)
t183 = sin(qJ(6));
t184 = sin(qJ(5));
t247 = qJD(6) * t183;
t249 = qJD(5) * t184;
t187 = cos(qJ(6));
t188 = cos(qJ(5));
t254 = t187 * t188;
t301 = qJD(5) + qJD(6);
t115 = -t183 * t249 - t184 * t247 + t254 * t301;
t207 = t183 * t188 + t187 * t184;
t116 = t301 * t207;
t149 = -t183 * t184 + t254;
t185 = sin(qJ(3));
t186 = sin(qJ(2));
t286 = cos(qJ(3));
t287 = cos(qJ(2));
t150 = t185 * t287 + t186 * t286;
t148 = t185 * t186 - t286 * t287;
t176 = -pkin(2) * t287 - pkin(1);
t198 = -t150 * qJ(4) + t176;
t290 = pkin(3) + pkin(9);
t83 = t148 * t290 + t198;
t225 = -pkin(10) * t148 - t83;
t309 = -pkin(8) - pkin(7);
t161 = t309 * t186;
t162 = t309 * t287;
t131 = -t286 * t161 - t162 * t185;
t93 = pkin(4) * t150 + t131;
t91 = t188 * t93;
t31 = pkin(5) * t150 + t184 * t225 + t91;
t260 = t148 * t188;
t90 = t184 * t93;
t47 = t188 * t83 + t90;
t37 = pkin(10) * t260 + t47;
t15 = -t183 * t37 + t187 * t31;
t16 = t183 * t31 + t187 * t37;
t302 = qJD(2) + qJD(3);
t117 = t148 * t302;
t118 = t302 * t150;
t251 = qJD(2) * t186;
t179 = pkin(2) * t251;
t199 = qJ(4) * t117 - qJD(4) * t150 + t179;
t36 = t118 * t290 + t199;
t132 = t185 * t161 - t162 * t286;
t197 = qJD(2) * t162;
t215 = qJD(2) * t161;
t73 = qJD(3) * t132 + t185 * t215 - t286 * t197;
t53 = -t117 * pkin(4) + t73;
t51 = t188 * t53;
t5 = -pkin(5) * t117 + t51 + (-pkin(10) * t118 - t36) * t184 + (t188 * t225 - t90) * qJD(5);
t248 = qJD(5) * t188;
t11 = t184 * t53 + t188 * t36 + t93 * t248 - t249 * t83;
t200 = -t118 * t188 + t148 * t249;
t6 = -pkin(10) * t200 + t11;
t3 = qJD(6) * t15 + t183 * t5 + t187 * t6;
t4 = -qJD(6) * t16 - t183 * t6 + t187 * t5;
t318 = -t115 * t16 + t15 * t116 - t4 * t149 - t207 * t3;
t317 = t149 * t183 - t207 * t187;
t244 = 2 * mrSges(7,3);
t316 = mrSges(5,2) - mrSges(4,1);
t95 = t149 * t148;
t24 = t207 * t118 + t301 * t95;
t25 = -t116 * t148 + t118 * t149;
t10 = Ifges(7,1) * t24 + Ifges(7,4) * t25 - t117 * Ifges(7,5);
t121 = mrSges(7,1) * t207 + mrSges(7,2) * t149;
t122 = Ifges(7,4) * t149 - Ifges(7,2) * t207;
t123 = Ifges(7,1) * t149 - Ifges(7,4) * t207;
t212 = mrSges(6,1) * t188 - mrSges(6,2) * t184;
t152 = t212 * qJD(5);
t279 = Ifges(6,4) * t184;
t210 = Ifges(6,2) * t188 + t279;
t153 = t210 * qJD(5);
t278 = Ifges(6,4) * t188;
t211 = Ifges(6,1) * t184 + t278;
t154 = t211 * qJD(5);
t158 = mrSges(6,1) * t184 + mrSges(6,2) * t188;
t159 = -Ifges(6,2) * t184 + t278;
t160 = Ifges(6,1) * t188 - t279;
t209 = Ifges(6,5) * t184 + Ifges(6,6) * t188;
t226 = t260 / 0.2e1;
t253 = -Ifges(7,5) * t116 - Ifges(7,6) * t115;
t223 = qJD(3) * t286;
t250 = qJD(3) * t185;
t72 = -t161 * t223 - t162 * t250 - t185 * t197 - t286 * t215;
t52 = -pkin(4) * t118 - t72;
t26 = pkin(5) * t200 + t52;
t261 = t148 * t184;
t201 = t118 * t184 + t148 * t248;
t29 = Ifges(6,4) * t201 - Ifges(6,2) * t200 - t117 * Ifges(6,6);
t30 = Ifges(6,1) * t201 - Ifges(6,4) * t200 - t117 * Ifges(6,5);
t313 = t159 * t188 + t160 * t184;
t46 = -t184 * t83 + t91;
t96 = t207 * t148;
t48 = Ifges(7,4) * t96 + Ifges(7,2) * t95 + Ifges(7,6) * t150;
t49 = Ifges(7,1) * t96 + Ifges(7,4) * t95 + Ifges(7,5) * t150;
t64 = mrSges(7,1) * t115 - mrSges(7,2) * t116;
t65 = -Ifges(7,4) * t116 - Ifges(7,2) * t115;
t66 = -Ifges(7,1) * t116 - Ifges(7,4) * t115;
t74 = (-pkin(5) * t188 - pkin(4)) * t148 + t132;
t81 = Ifges(6,6) * t150 + t148 * t210;
t82 = Ifges(6,5) * t150 + t148 * t211;
t9 = Ifges(7,4) * t24 + Ifges(7,2) * t25 - t117 * Ifges(7,6);
t94 = -t148 * pkin(4) + t132;
t315 = (Ifges(5,5) - Ifges(4,6) + t313 / 0.2e1) * t118 + t316 * t73 - t154 * t261 / 0.2e1 + t150 * t253 / 0.2e1 - t81 * t248 / 0.2e1 + (Ifges(5,4) - Ifges(4,5)) * t117 + (-t150 * t209 / 0.2e1 + t160 * t226) * qJD(5) - (t148 * t159 + t82) * t249 / 0.2e1 - (Ifges(6,5) * t188 + Ifges(7,5) * t149 - Ifges(6,6) * t184 - Ifges(7,6) * t207) * t117 / 0.2e1 - t207 * t9 / 0.2e1 + t46 * mrSges(6,3) * t249 + t318 * mrSges(7,3) - t153 * t226 + t74 * t64 + t95 * t65 / 0.2e1 + t96 * t66 / 0.2e1 - t115 * t48 / 0.2e1 - t116 * t49 / 0.2e1 + t26 * t121 + t25 * t122 / 0.2e1 + t24 * t123 / 0.2e1 + t149 * t10 / 0.2e1 + t94 * t152 + t52 * t158 - t184 * t29 / 0.2e1 + t188 * t30 / 0.2e1;
t105 = -mrSges(6,2) * t150 + mrSges(6,3) * t260;
t305 = qJD(5) * t47;
t12 = -t184 * t36 - t305 + t51;
t275 = t11 * t184;
t303 = t12 * t188 + t275;
t54 = mrSges(6,2) * t117 - mrSges(6,3) * t200;
t55 = -mrSges(6,1) * t117 - mrSges(6,3) * t201;
t191 = m(6) * ((-t184 * t46 + t188 * t47) * qJD(5) + t303) + t188 * t55 + t184 * t54;
t314 = t105 * t248 + t191;
t312 = -t115 * t207 + t116 * t149;
t252 = t184 ^ 2 + t188 ^ 2;
t310 = (-mrSges(4,2) * t286 + (-t252 * mrSges(6,3) + t316) * t185) * pkin(2) * qJD(3);
t308 = Ifges(3,1) - Ifges(3,2);
t306 = mrSges(5,3) + t158;
t236 = t286 * pkin(2);
t175 = -t236 - pkin(3);
t169 = -pkin(9) + t175;
t235 = pkin(2) * t250;
t304 = -t169 * t249 + t188 * t235;
t267 = t115 * t183;
t276 = pkin(5) * qJD(6);
t300 = (-pkin(5) * t267 + t276 * t317) * mrSges(7,3);
t296 = 2 * m(5);
t295 = 0.2e1 * m(6);
t294 = 2 * m(7);
t293 = 0.2e1 * t64;
t292 = 0.2e1 * t121;
t291 = 0.2e1 * t152;
t285 = pkin(2) * t185;
t282 = -pkin(10) + t169;
t281 = -pkin(10) - t290;
t273 = t117 * mrSges(5,1);
t265 = t116 * t187;
t216 = pkin(2) * t223;
t163 = t216 + qJD(4);
t170 = qJ(4) + t285;
t256 = t163 * t170;
t246 = qJD(6) * t187;
t239 = t312 * t244;
t237 = Ifges(7,5) * t24 + Ifges(7,6) * t25 - Ifges(7,3) * t117;
t233 = mrSges(7,3) * t265;
t224 = qJD(2) * t287;
t141 = t282 * t188;
t156 = t281 * t188;
t220 = -t116 * mrSges(7,1) - t115 * mrSges(7,2);
t219 = t150 * t235;
t218 = t184 * t235;
t208 = t131 * t73 - t132 * t72;
t140 = t282 * t184;
t98 = t140 * t187 + t141 * t183;
t97 = -t140 * t183 + t141 * t187;
t155 = t281 * t184;
t125 = t155 * t187 + t156 * t183;
t124 = -t155 * t183 + t156 * t187;
t177 = pkin(10) * t249;
t129 = t177 + t304;
t130 = qJD(5) * t141 + t218;
t44 = qJD(6) * t97 + t129 * t183 + t130 * t187;
t45 = -qJD(6) * t98 + t129 * t187 - t130 * t183;
t206 = t45 * mrSges(7,1) - t44 * mrSges(7,2) + t253;
t144 = t249 * t290 + t177;
t145 = qJD(5) * t156;
t62 = qJD(6) * t124 + t144 * t183 + t145 * t187;
t63 = -qJD(6) * t125 + t144 * t187 - t145 * t183;
t205 = t63 * mrSges(7,1) - t62 * mrSges(7,2) + t253;
t204 = qJ(4) * t163 + qJD(4) * t170;
t203 = t252 * t235;
t202 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t237;
t196 = Ifges(6,5) * t201 - Ifges(6,6) * t200 - Ifges(6,3) * t117;
t195 = t98 * t115 - t97 * t116 + t45 * t149 + t207 * t44;
t194 = t125 * t115 - t124 * t116 + t63 * t149 + t207 * t62;
t192 = -qJD(5) * t313 - t115 * t122 - t116 * t123 + t149 * t66 + t184 * t153 - t188 * t154 - t207 * t65;
t180 = t184 * pkin(5);
t178 = pkin(5) * t248;
t171 = qJ(4) + t180;
t164 = qJD(4) + t178;
t157 = t170 + t180;
t146 = t178 + t163;
t142 = (-mrSges(7,1) * t183 - mrSges(7,2) * t187) * t276;
t104 = mrSges(6,1) * t150 - mrSges(6,3) * t261;
t101 = t148 * pkin(3) + t198;
t100 = t212 * t148;
t78 = mrSges(7,1) * t150 - mrSges(7,3) * t96;
t77 = -mrSges(7,2) * t150 + mrSges(7,3) * t95;
t57 = -mrSges(7,1) * t95 + mrSges(7,2) * t96;
t56 = pkin(3) * t118 + t199;
t38 = mrSges(6,1) * t200 + mrSges(6,2) * t201;
t18 = mrSges(7,2) * t117 + mrSges(7,3) * t25;
t17 = -mrSges(7,1) * t117 - mrSges(7,3) * t24;
t13 = -mrSges(7,1) * t25 + mrSges(7,2) * t24;
t1 = [0.2e1 * m(4) * (t176 * t179 + t208) + t150 * t237 + (0.2e1 * Ifges(3,4) * t287 + t186 * t308) * t224 + (-0.2e1 * Ifges(3,4) * t186 + t287 * t308) * t251 - t200 * t81 + t201 * t82 + 0.2e1 * (mrSges(4,1) * t148 + mrSges(4,2) * t150) * t179 + t150 * t196 + 0.2e1 * (mrSges(4,1) * t176 - mrSges(5,2) * t101 + (Ifges(4,2) + Ifges(5,3)) * t148) * t118 + 0.2e1 * (mrSges(4,3) + mrSges(5,1)) * (-t117 * t131 - t118 * t132 + t148 * t72 + t150 * t73) - 0.2e1 * pkin(1) * (mrSges(3,1) * t186 + mrSges(3,2) * t287) * qJD(2) + (-0.2e1 * mrSges(4,2) * t176 + 0.2e1 * mrSges(5,3) * t101 - Ifges(7,5) * t96 - Ifges(7,6) * t95 - t148 * t209 + (-(2 * Ifges(4,1)) - (2 * Ifges(5,2)) - Ifges(6,3) - Ifges(7,3)) * t150) * t117 + 0.2e1 * (Ifges(5,6) + Ifges(4,4)) * (t117 * t148 - t118 * t150) + (t15 * t4 + t16 * t3 + t26 * t74) * t294 + (t11 * t47 + t12 * t46 + t52 * t94) * t295 + (t101 * t56 + t208) * t296 + t29 * t260 + t30 * t261 + 0.2e1 * t15 * t17 + 0.2e1 * t16 * t18 + t25 * t48 + t24 * t49 + 0.2e1 * t47 * t54 + 0.2e1 * t46 * t55 + 0.2e1 * t26 * t57 + 0.2e1 * t74 * t13 + 0.2e1 * t3 * t77 + 0.2e1 * t4 * t78 + 0.2e1 * t94 * t38 + t95 * t9 + t96 * t10 - 0.2e1 * t52 * t100 + 0.2e1 * t12 * t104 + 0.2e1 * t11 * t105 + 0.2e1 * t56 * (-mrSges(5,2) * t148 - mrSges(5,3) * t150); t315 - t175 * t273 - Ifges(3,6) * t251 + m(5) * (t131 * t235 + t132 * t163 - t170 * t72 + t175 * t73) + t304 * t104 + (-t248 * t47 - t303) * mrSges(6,3) + t314 * t169 + m(4) * (-t286 * t73 - t185 * t72 + (t131 * t185 + t132 * t286) * qJD(3)) * pkin(2) + m(7) * (t146 * t74 + t15 * t45 + t157 * t26 + t16 * t44 + t3 * t98 + t4 * t97) + t105 * t218 + m(6) * (t163 * t94 + t170 * t52 + (t184 * t47 + t188 * t46) * t235) + Ifges(3,5) * t224 - t72 * mrSges(5,3) + t72 * mrSges(4,2) + t44 * t77 + t45 * t78 + t97 * t17 + t98 * t18 + t146 * t57 + t157 * t13 - t163 * t100 + t170 * t38 + (-t118 * t170 - t148 * t163 + t219) * mrSges(5,1) + (t117 * t236 - t118 * t285 - t148 * t216 + t219) * mrSges(4,3) + (-mrSges(3,1) * t224 + mrSges(3,2) * t251) * pkin(7); t146 * t292 + t170 * t291 + t157 * t293 + 0.2e1 * t306 * t163 + 0.2e1 * t310 + (t169 * t203 + t256) * t295 + (t175 * t235 + t256) * t296 + (t146 * t157 + t44 * t98 + t45 * t97) * t294 - t195 * t244 + t192; -(-t104 * t249 + t314) * t290 + (-t275 + (-t12 - t305) * t188) * mrSges(6,3) + (-mrSges(5,3) + mrSges(4,2)) * t72 + m(6) * (qJ(4) * t52 + qJD(4) * t94) + m(7) * (t124 * t4 + t125 * t3 + t15 * t63 + t16 * t62 + t164 * t74 + t171 * t26) + m(5) * (-pkin(3) * t73 - qJ(4) * t72 + qJD(4) * t132) + (pkin(3) * t117 - qJ(4) * t118 - qJD(4) * t148) * mrSges(5,1) + qJ(4) * t38 + t62 * t77 + t63 * t78 - qJD(4) * t100 + t124 * t17 + t125 * t18 + t164 * t57 + t171 * t13 + t315; t192 + m(6) * (-t203 * t290 + t204) + m(5) * (-pkin(3) * t235 + t204) + m(7) * (t124 * t45 + t125 * t44 + t146 * t171 + t157 * t164 + t62 * t98 + t63 * t97) + t310 + (t157 + t171) * t64 + (qJ(4) + t170) * t152 + (t146 + t164) * t121 + ((-t45 - t63) * t149 - (t44 + t62) * t207 - (-t124 - t97) * t116 + (-t125 - t98) * t115) * mrSges(7,3) + t306 * (qJD(4) + t163); t164 * t292 + t171 * t293 + qJ(4) * t291 + (t124 * t63 + t125 * t62 + t164 * t171) * t294 + 0.2e1 * ((m(5) + m(6)) * qJ(4) + t306) * qJD(4) - t194 * t244 + t192; -t273 + t115 * t77 - t116 * t78 + t207 * t18 + t149 * t17 + (-t184 * t104 + t188 * t105) * qJD(5) - m(7) * t318 + m(5) * t73 + t191; m(7) * t195 + (m(6) * t252 + m(5)) * t235 + t239; m(7) * t194 + t239; -t312 * t294; t12 * mrSges(6,1) - t11 * mrSges(6,2) + (m(7) * (-t15 * t247 + t16 * t246 + t183 * t3 + t187 * t4) + t77 * t246 + t183 * t18 - t78 * t247 + t187 * t17) * pkin(5) + t196 + t202; t212 * t235 + (-t158 * t169 - t209) * qJD(5) + (m(7) * (t183 * t44 + t187 * t45 + t246 * t98 - t247 * t97) + t233) * pkin(5) + t206 + t300; ((mrSges(6,2) * t290 - Ifges(6,6)) * t188 + (mrSges(6,1) * t290 - Ifges(6,5)) * t184) * qJD(5) + (m(7) * (-t124 * t247 + t125 * t246 + t183 * t62 + t187 * t63) + t233) * pkin(5) + t205 + t300; -t158 * qJD(5) + m(7) * (-t317 * qJD(6) - t265 + t267) * pkin(5) + t220; 0.2e1 * t142; t202; t206; t205; t220; t142; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
