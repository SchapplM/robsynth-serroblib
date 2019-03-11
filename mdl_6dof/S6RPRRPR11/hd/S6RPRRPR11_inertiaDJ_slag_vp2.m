% Calculate time derivative of joint inertia matrix for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:39:25
% EndTime: 2019-03-09 05:39:40
% DurationCPUTime: 7.36s
% Computational Cost: add. (16660->646), mult. (48728->971), div. (0->0), fcn. (51978->14), ass. (0->266)
t233 = sin(pkin(7));
t237 = cos(pkin(7));
t238 = cos(pkin(6));
t232 = sin(pkin(12));
t234 = sin(pkin(6));
t236 = cos(pkin(12));
t277 = t234 * t236;
t302 = pkin(1) * t238;
t270 = qJ(2) * t277 + t232 * t302;
t147 = (t233 * t238 + t237 * t277) * pkin(9) + t270;
t241 = sin(qJ(3));
t244 = cos(qJ(3));
t223 = t236 * t302;
t281 = t232 * t234;
t151 = pkin(2) * t238 + t223 + (-pkin(9) * t237 - qJ(2)) * t281;
t172 = (-pkin(9) * t232 * t233 - pkin(2) * t236 - pkin(1)) * t234;
t248 = t151 * t237 + t172 * t233;
t90 = -t241 * t147 + t248 * t244;
t231 = sin(pkin(13));
t235 = cos(pkin(13));
t329 = mrSges(6,1) * t231 + mrSges(6,2) * t235;
t273 = t237 * t244;
t278 = t233 * t244;
t328 = t234 * (-t232 * t241 + t236 * t273) + t238 * t278;
t322 = 2 * pkin(10);
t327 = -2 * Ifges(4,4);
t239 = sin(qJ(6));
t242 = cos(qJ(6));
t247 = t231 * t239 - t235 * t242;
t306 = -t247 / 0.2e1;
t199 = t231 * t242 + t235 * t239;
t305 = t199 / 0.2e1;
t303 = t235 / 0.2e1;
t274 = t237 * t241;
t279 = t233 * t241;
t150 = t238 * t279 + (t232 * t244 + t236 * t274) * t234;
t146 = t150 * qJD(3);
t145 = t328 * qJD(3);
t268 = qJD(2) * t234;
t257 = t232 * t268;
t253 = t233 * t257;
t106 = pkin(3) * t146 - pkin(10) * t145 + t253;
t240 = sin(qJ(4));
t243 = cos(qJ(4));
t264 = qJD(4) * t243;
t265 = qJD(4) * t240;
t78 = (-t232 * t274 + t236 * t244) * t268 + t90 * qJD(3);
t112 = -t151 * t233 + t237 * t172;
t83 = -pkin(3) * t328 - pkin(10) * t150 + t112;
t187 = -t233 * t277 + t237 * t238;
t137 = t244 * t147;
t91 = t151 * t274 + t172 * t279 + t137;
t88 = pkin(10) * t187 + t91;
t29 = t106 * t243 - t240 * t78 - t88 * t264 - t83 * t265;
t22 = -pkin(4) * t146 - t29;
t113 = t150 * t240 - t187 * t243;
t96 = -qJD(4) * t113 + t145 * t243;
t72 = t146 * t235 - t231 * t96;
t73 = t146 * t231 + t235 * t96;
t44 = -t72 * mrSges(6,1) + t73 * mrSges(6,2);
t326 = -m(6) * t22 - t44;
t188 = t247 * qJD(6);
t324 = 0.2e1 * m(6);
t323 = 2 * m(7);
t321 = -2 * mrSges(4,3);
t320 = m(6) / 0.2e1;
t114 = t150 * t243 + t187 * t240;
t97 = -t114 * t231 - t235 * t328;
t98 = t114 * t235 - t231 * t328;
t60 = -t239 * t98 + t242 * t97;
t319 = t60 / 0.2e1;
t61 = t239 * t97 + t242 * t98;
t318 = t61 / 0.2e1;
t317 = t72 / 0.2e1;
t316 = t73 / 0.2e1;
t154 = Ifges(7,4) * t199 - Ifges(7,2) * t247;
t314 = t154 / 0.2e1;
t155 = Ifges(7,1) * t199 - Ifges(7,4) * t247;
t313 = t155 / 0.2e1;
t293 = Ifges(6,4) * t235;
t251 = -Ifges(6,2) * t231 + t293;
t170 = (Ifges(6,6) * t240 + t243 * t251) * qJD(4);
t312 = t170 / 0.2e1;
t294 = Ifges(6,4) * t231;
t252 = Ifges(6,1) * t235 - t294;
t171 = (Ifges(6,5) * t240 + t243 * t252) * qJD(4);
t311 = t171 / 0.2e1;
t181 = t199 * t240;
t310 = -t181 / 0.2e1;
t182 = t247 * t240;
t309 = -t182 / 0.2e1;
t308 = -t188 / 0.2e1;
t189 = t199 * qJD(6);
t307 = -t189 / 0.2e1;
t304 = -t231 / 0.2e1;
t301 = pkin(10) * t240;
t300 = pkin(11) + qJ(5);
t28 = t240 * t106 + t243 * t78 + t83 * t264 - t265 * t88;
t19 = qJ(5) * t146 - qJD(5) * t328 + t28;
t79 = (t232 * t273 + t236 * t241) * t268 + (t241 * t248 + t137) * qJD(3);
t95 = qJD(4) * t114 + t145 * t240;
t36 = t95 * pkin(4) - t96 * qJ(5) - t114 * qJD(5) + t79;
t12 = t235 * t19 + t231 * t36;
t46 = t240 * t83 + t243 * t88;
t42 = -qJ(5) * t328 + t46;
t87 = -t187 * pkin(3) - t90;
t56 = t113 * pkin(4) - t114 * qJ(5) + t87;
t21 = t231 * t56 + t235 * t42;
t75 = mrSges(5,1) * t146 - mrSges(5,3) * t96;
t299 = t44 - t75;
t296 = Ifges(5,4) * t240;
t295 = Ifges(5,4) * t243;
t292 = t146 * Ifges(5,5);
t291 = t146 * Ifges(5,6);
t290 = t328 * Ifges(5,6);
t289 = t244 * t79;
t288 = -t243 * mrSges(5,1) + t240 * mrSges(5,2) - mrSges(4,1);
t103 = -mrSges(5,1) * t328 - mrSges(5,3) * t114;
t62 = -mrSges(6,1) * t97 + mrSges(6,2) * t98;
t287 = -t103 + t62;
t209 = -mrSges(6,1) * t235 + mrSges(6,2) * t231;
t286 = t209 - mrSges(5,1);
t190 = -t243 * t237 + t240 * t279;
t266 = qJD(3) * t244;
t255 = t233 * t266;
t161 = -qJD(4) * t190 + t243 * t255;
t285 = t161 * t243;
t191 = t237 * t240 + t243 * t279;
t162 = qJD(4) * t191 + t240 * t255;
t284 = t162 * t240;
t124 = t190 * t162;
t283 = t231 * t240;
t282 = t231 * t243;
t276 = t235 * t240;
t275 = t235 * t243;
t131 = -t189 * t240 - t247 * t264;
t132 = t188 * t240 - t199 * t264;
t105 = -t132 * mrSges(7,1) + t131 * mrSges(7,2);
t183 = t329 * t264;
t272 = t105 + t183;
t271 = Ifges(4,5) * t145 - Ifges(4,6) * t146;
t142 = -Ifges(7,5) * t188 - Ifges(7,6) * t189;
t186 = -qJD(5) * t240 + (pkin(4) * t240 - qJ(5) * t243) * qJD(4);
t261 = pkin(10) * t265;
t156 = t235 * t186 + t231 * t261;
t207 = -pkin(4) * t243 - qJ(5) * t240 - pkin(3);
t174 = pkin(10) * t275 + t231 * t207;
t267 = qJD(3) * t241;
t26 = qJD(6) * t60 + t239 * t72 + t242 * t73;
t27 = -qJD(6) * t61 - t239 * t73 + t242 * t72;
t7 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t95;
t262 = Ifges(5,5) * t96 - Ifges(5,6) * t95 + Ifges(5,3) * t146;
t99 = Ifges(7,5) * t131 + Ifges(7,6) * t132 + Ifges(7,3) * t265;
t258 = pkin(5) * t231 + pkin(10);
t256 = t233 * t267;
t254 = Ifges(6,5) * t231 / 0.2e1 + Ifges(6,6) * t303 + Ifges(7,5) * t305 + Ifges(7,6) * t306;
t10 = -t27 * mrSges(7,1) + t26 * mrSges(7,2);
t11 = -t19 * t231 + t235 * t36;
t20 = -t231 * t42 + t235 * t56;
t45 = -t240 * t88 + t243 * t83;
t250 = Ifges(6,5) * t235 - Ifges(6,6) * t231;
t13 = pkin(5) * t113 - pkin(11) * t98 + t20;
t15 = pkin(11) * t97 + t21;
t3 = t13 * t242 - t15 * t239;
t4 = t13 * t239 + t15 * t242;
t126 = -t161 * t231 + t235 * t256;
t127 = t161 * t235 + t231 * t256;
t249 = -t126 * t231 + t127 * t235;
t197 = t235 * t207;
t148 = -pkin(11) * t276 + t197 + (-pkin(10) * t231 - pkin(5)) * t243;
t163 = -pkin(11) * t283 + t174;
t108 = t148 * t242 - t163 * t239;
t109 = t148 * t239 + t163 * t242;
t158 = -t231 * t191 - t235 * t278;
t159 = t235 * t191 - t231 * t278;
t110 = t158 * t242 - t159 * t239;
t111 = t158 * t239 + t159 * t242;
t208 = t300 * t231;
t210 = t300 * t235;
t164 = -t208 * t242 - t210 * t239;
t165 = -t208 * t239 + t210 * t242;
t43 = pkin(4) * t328 - t45;
t246 = t190 * t264 + t284;
t228 = Ifges(5,5) * t264;
t226 = -pkin(5) * t235 - pkin(4);
t216 = Ifges(5,1) * t240 + t295;
t215 = Ifges(5,2) * t243 + t296;
t213 = Ifges(6,1) * t231 + t293;
t212 = Ifges(6,2) * t235 + t294;
t206 = t258 * t240;
t204 = (Ifges(5,1) * t243 - t296) * qJD(4);
t203 = (-Ifges(5,2) * t240 + t295) * qJD(4);
t202 = (mrSges(5,1) * t240 + mrSges(5,2) * t243) * qJD(4);
t201 = -mrSges(6,1) * t243 - mrSges(6,3) * t276;
t200 = mrSges(6,2) * t243 - mrSges(6,3) * t283;
t195 = t258 * t264;
t194 = (mrSges(6,1) * t240 - mrSges(6,3) * t275) * qJD(4);
t193 = (-mrSges(6,2) * t240 - mrSges(6,3) * t282) * qJD(4);
t192 = t329 * t240;
t180 = -Ifges(6,5) * t243 + t240 * t252;
t179 = -Ifges(6,6) * t243 + t240 * t251;
t178 = -Ifges(6,3) * t243 + t240 * t250;
t175 = t231 * t186;
t173 = -pkin(10) * t282 + t197;
t169 = (Ifges(6,3) * t240 + t243 * t250) * qJD(4);
t168 = -mrSges(7,1) * t243 + mrSges(7,3) * t182;
t167 = mrSges(7,2) * t243 - mrSges(7,3) * t181;
t157 = -t235 * t261 + t175;
t152 = mrSges(7,1) * t247 + mrSges(7,2) * t199;
t144 = -Ifges(7,1) * t188 - Ifges(7,4) * t189;
t143 = -Ifges(7,4) * t188 - Ifges(7,2) * t189;
t141 = mrSges(7,1) * t189 - mrSges(7,2) * t188;
t134 = t175 + (-pkin(10) * t276 - pkin(11) * t282) * qJD(4);
t133 = mrSges(7,1) * t181 - mrSges(7,2) * t182;
t125 = (pkin(5) * t240 - pkin(11) * t275) * qJD(4) + t156;
t123 = -Ifges(7,1) * t182 - Ifges(7,4) * t181 - Ifges(7,5) * t243;
t122 = -Ifges(7,4) * t182 - Ifges(7,2) * t181 - Ifges(7,6) * t243;
t121 = -Ifges(7,5) * t182 - Ifges(7,6) * t181 - Ifges(7,3) * t243;
t120 = -qJD(5) * t199 - qJD(6) * t165;
t119 = -qJD(5) * t247 + qJD(6) * t164;
t118 = -mrSges(7,2) * t265 + mrSges(7,3) * t132;
t117 = mrSges(7,1) * t265 - mrSges(7,3) * t131;
t116 = mrSges(4,1) * t187 - mrSges(4,3) * t150;
t115 = -mrSges(4,2) * t187 + mrSges(4,3) * t328;
t107 = mrSges(4,1) * t146 + mrSges(4,2) * t145;
t102 = mrSges(5,2) * t328 - mrSges(5,3) * t113;
t101 = Ifges(7,1) * t131 + Ifges(7,4) * t132 + Ifges(7,5) * t265;
t100 = Ifges(7,4) * t131 + Ifges(7,2) * t132 + Ifges(7,6) * t265;
t89 = mrSges(5,1) * t113 + mrSges(5,2) * t114;
t74 = -mrSges(5,2) * t146 - mrSges(5,3) * t95;
t70 = Ifges(5,1) * t114 - Ifges(5,4) * t113 - Ifges(5,5) * t328;
t69 = Ifges(5,4) * t114 - Ifges(5,2) * t113 - t290;
t68 = mrSges(6,1) * t113 - mrSges(6,3) * t98;
t67 = -mrSges(6,2) * t113 + mrSges(6,3) * t97;
t66 = -qJD(6) * t109 + t125 * t242 - t134 * t239;
t65 = qJD(6) * t108 + t125 * t239 + t134 * t242;
t64 = -qJD(6) * t111 + t126 * t242 - t127 * t239;
t63 = qJD(6) * t110 + t126 * t239 + t127 * t242;
t59 = mrSges(5,1) * t95 + mrSges(5,2) * t96;
t58 = Ifges(5,1) * t96 - Ifges(5,4) * t95 + t292;
t57 = Ifges(5,4) * t96 - Ifges(5,2) * t95 + t291;
t53 = mrSges(6,1) * t95 - mrSges(6,3) * t73;
t52 = -mrSges(6,2) * t95 + mrSges(6,3) * t72;
t51 = Ifges(6,1) * t98 + Ifges(6,4) * t97 + Ifges(6,5) * t113;
t50 = Ifges(6,4) * t98 + Ifges(6,2) * t97 + Ifges(6,6) * t113;
t49 = Ifges(6,5) * t98 + Ifges(6,6) * t97 + Ifges(6,3) * t113;
t48 = mrSges(7,1) * t113 - mrSges(7,3) * t61;
t47 = -mrSges(7,2) * t113 + mrSges(7,3) * t60;
t40 = Ifges(6,1) * t73 + Ifges(6,4) * t72 + Ifges(6,5) * t95;
t39 = Ifges(6,4) * t73 + Ifges(6,2) * t72 + Ifges(6,6) * t95;
t38 = Ifges(6,5) * t73 + Ifges(6,6) * t72 + Ifges(6,3) * t95;
t37 = -pkin(5) * t97 + t43;
t33 = -mrSges(7,1) * t60 + mrSges(7,2) * t61;
t32 = Ifges(7,1) * t61 + Ifges(7,4) * t60 + Ifges(7,5) * t113;
t31 = Ifges(7,4) * t61 + Ifges(7,2) * t60 + Ifges(7,6) * t113;
t30 = Ifges(7,5) * t61 + Ifges(7,6) * t60 + Ifges(7,3) * t113;
t17 = -mrSges(7,2) * t95 + mrSges(7,3) * t27;
t16 = mrSges(7,1) * t95 - mrSges(7,3) * t26;
t14 = -pkin(5) * t72 + t22;
t9 = Ifges(7,1) * t26 + Ifges(7,4) * t27 + Ifges(7,5) * t95;
t8 = Ifges(7,4) * t26 + Ifges(7,2) * t27 + Ifges(7,6) * t95;
t6 = pkin(11) * t72 + t12;
t5 = pkin(5) * t95 - pkin(11) * t73 + t11;
t2 = -qJD(6) * t4 - t239 * t6 + t242 * t5;
t1 = qJD(6) * t3 + t239 * t5 + t242 * t6;
t18 = [0.2e1 * (t236 * (-mrSges(3,2) * t238 + mrSges(3,3) * t277) + m(3) * (t270 * t236 + (qJ(2) * t281 - t223) * t232)) * t268 + (t38 + t7 - t57) * t113 + (t1 * t4 + t14 * t37 + t2 * t3) * t323 + (t11 * t20 + t12 * t21 + t22 * t43) * t324 + 0.2e1 * m(5) * (t28 * t46 + t29 * t45 + t79 * t87) + (t49 + t30 - t69) * t95 + (t150 * t327 + Ifges(5,5) * t114 - Ifges(4,6) * t187 - Ifges(5,6) * t113 + t321 * t91 - ((2 * Ifges(4,2)) + Ifges(5,3)) * t328) * t146 + (0.2e1 * Ifges(4,1) * t150 + Ifges(4,5) * t187 + t321 * t90 - t327 * t328) * t145 + 0.2e1 * (-mrSges(4,1) * t328 + mrSges(4,2) * t150) * t253 - t328 * t262 + t114 * t58 + 0.2e1 * t78 * t115 - 0.2e1 * t79 * t116 + 0.2e1 * t112 * t107 + 0.2e1 * t28 * t102 + 0.2e1 * t29 * t103 + t96 * t70 + t97 * t39 + t98 * t40 + 0.2e1 * t87 * t59 + 0.2e1 * t79 * t89 + t72 * t50 + t73 * t51 + 0.2e1 * t46 * t74 + 0.2e1 * t45 * t75 + 0.2e1 * t12 * t67 + 0.2e1 * t11 * t68 + t60 * t8 + t61 * t9 + 0.2e1 * t22 * t62 + 0.2e1 * t2 * t48 + 0.2e1 * t21 * t52 + 0.2e1 * t20 * t53 + 0.2e1 * t43 * t44 + 0.2e1 * t1 * t47 + 0.2e1 * t37 * t10 + t27 * t31 + t26 * t32 + 0.2e1 * t14 * t33 + 0.2e1 * t4 * t17 + 0.2e1 * t3 * t16 + 0.2e1 * m(4) * (t112 * t253 + t78 * t91 - t79 * t90) + t187 * t271 - 0.2e1 * (mrSges(3,1) * t238 - mrSges(3,3) * t281) * t257; t161 * t102 + t237 * t107 + t110 * t16 + t111 * t17 + t126 * t68 + t127 * t67 + t158 * t53 + t159 * t52 + t191 * t74 + t63 * t47 + t64 * t48 + (t10 + t299) * t190 + (t33 + t287) * t162 + m(5) * (t161 * t46 - t162 * t45 - t190 * t29 + t191 * t28) + m(7) * (t1 * t111 + t110 * t2 + t14 * t190 + t162 * t37 + t3 * t64 + t4 * t63) + m(6) * (t11 * t158 + t12 * t159 + t126 * t20 + t127 * t21 + t162 * t43 + t190 * t22) + (-t244 * t59 + (-t145 * t244 - t146 * t241) * mrSges(4,3) + (t244 * t115 + (-t116 + t89) * t241) * qJD(3) + m(5) * (t267 * t87 - t289) + m(4) * (t237 * t257 + t241 * t78 + t266 * t91 - t267 * t90 - t289)) * t233; 0.2e1 * m(7) * (t110 * t64 + t111 * t63 + t124) + 0.2e1 * m(6) * (t126 * t158 + t127 * t159 + t124) + 0.2e1 * m(5) * (-t233 ^ 2 * t241 * t266 + t191 * t161 + t124); (-t215 / 0.2e1 + t178 / 0.2e1 + t121 / 0.2e1) * t95 + (-t203 / 0.2e1 + t169 / 0.2e1 + t99 / 0.2e1) * t113 + t271 + ((t290 / 0.2e1 + t30 / 0.2e1 + t49 / 0.2e1 - t69 / 0.2e1 - t46 * mrSges(5,3) + (-m(5) * t46 - t102) * pkin(10)) * t240 + (t70 / 0.2e1 - t45 * mrSges(5,3) + t50 * t304 + t51 * t303 + (-m(5) * t45 + m(6) * t43 + t287) * pkin(10)) * t243) * qJD(4) + t9 * t309 + t8 * t310 + t98 * t311 + t97 * t312 + t180 * t316 + t179 * t317 + t101 * t318 + t100 * t319 - t328 * t228 / 0.2e1 + t96 * t216 / 0.2e1 + t195 * t33 + t12 * t200 + t11 * t201 + t87 * t202 + t114 * t204 / 0.2e1 + t206 * t10 + t22 * t192 + t21 * t193 + t20 * t194 + t43 * t183 + t1 * t167 + t2 * t168 + t173 * t53 + t174 * t52 + t156 * t68 + t157 * t67 + t131 * t32 / 0.2e1 + t132 * t31 / 0.2e1 + t14 * t133 + t4 * t118 + t27 * t122 / 0.2e1 + t26 * t123 / 0.2e1 + t3 * t117 + t108 * t16 + t109 * t17 + t37 * t105 - t78 * mrSges(4,2) + t65 * t47 + t66 * t48 - pkin(3) * t59 + m(7) * (t1 * t109 + t108 * t2 + t14 * t206 + t195 * t37 + t3 * t66 + t4 * t65) + t288 * t79 + (t291 / 0.2e1 + t57 / 0.2e1 - t38 / 0.2e1 - t7 / 0.2e1 + t28 * mrSges(5,3) + pkin(10) * t74) * t243 + m(5) * (pkin(10) * t243 * t28 - pkin(3) * t79 - t29 * t301) + m(6) * (t11 * t173 + t12 * t174 + t156 * t20 + t157 * t21 + t22 * t301) + (t292 / 0.2e1 + t58 / 0.2e1 - t29 * mrSges(5,3) + t39 * t304 + t40 * t303 + t299 * pkin(10)) * t240; t110 * t117 + t111 * t118 + t126 * t201 + t127 * t200 + t158 * t194 + t159 * t193 + t63 * t167 + t64 * t168 + t272 * t190 + (t133 + t192) * t162 + (-t244 * t202 + (-mrSges(4,2) * t244 + t241 * t288) * qJD(3)) * t233 + m(7) * (t108 * t64 + t109 * t63 + t110 * t66 + t111 * t65 + t162 * t206 + t190 * t195) + m(6) * (t126 * t173 + t127 * t174 + t156 * t158 + t157 * t159) - m(5) * pkin(3) * t256 + (t246 * t320 + m(5) * (-t191 * t265 + t246 + t285) / 0.2e1) * t322 + (t285 + t284 + (t190 * t243 - t191 * t240) * qJD(4)) * mrSges(5,3); -0.2e1 * pkin(3) * t202 - t181 * t100 - t182 * t101 + 0.2e1 * t206 * t105 + 0.2e1 * t108 * t117 + 0.2e1 * t109 * t118 + t132 * t122 + t131 * t123 + 0.2e1 * t195 * t133 + 0.2e1 * t156 * t201 + 0.2e1 * t157 * t200 + 0.2e1 * t65 * t167 + 0.2e1 * t66 * t168 + 0.2e1 * t173 * t194 + 0.2e1 * t174 * t193 + (t156 * t173 + t157 * t174) * t324 + (t108 * t66 + t109 * t65 + t195 * t206) * t323 + (t203 - t169 - t99) * t243 + (-t170 * t231 + t171 * t235 + t183 * t322 + t204) * t240 + ((t121 + t178 - t215) * t240 + (-t231 * t179 + t235 * t180 + t216 + (m(6) * t301 + t192) * t322) * t243) * qJD(4); t326 * pkin(4) + t9 * t305 + t8 * t306 + t31 * t307 + t32 * t308 + t26 * t313 + t27 * t314 + t213 * t316 + t212 * t317 + t144 * t318 + t143 * t319 + (t40 / 0.2e1 - t11 * mrSges(6,3) - qJ(5) * t53 - qJD(5) * t68 + m(6) * (-qJ(5) * t11 - qJD(5) * t20)) * t231 + (t39 / 0.2e1 + t12 * mrSges(6,3) + qJ(5) * t52 + qJD(5) * t67 + m(6) * (qJ(5) * t12 + qJD(5) * t21)) * t235 + (-t1 * t247 + t188 * t3 - t189 * t4 - t199 * t2) * mrSges(7,3) + t262 + t226 * t10 + t22 * t209 + t164 * t16 + t165 * t17 + t14 * t152 + t37 * t141 + t113 * t142 / 0.2e1 + t119 * t47 + t120 * t48 - t28 * mrSges(5,2) + t29 * mrSges(5,1) + t254 * t95 + m(7) * (t1 * t165 + t119 * t4 + t120 * t3 + t14 * t226 + t164 * t2); -t161 * mrSges(5,2) + t190 * t141 + t249 * mrSges(6,3) + (t152 + t286) * t162 + m(7) * (t110 * t120 + t111 * t119 + t162 * t226 + t164 * t64 + t165 * t63) + m(6) * (-pkin(4) * t162 + (-t158 * t231 + t159 * t235) * qJD(5) + t249 * qJ(5)) + (t110 * t188 - t111 * t189 - t199 * t64 - t247 * t63) * mrSges(7,3); t228 - t243 * t142 / 0.2e1 + t226 * t105 + t195 * t152 + t100 * t306 + t101 * t305 + t206 * t141 + t123 * t308 + t122 * t307 + t143 * t310 + t144 * t309 - pkin(4) * t183 + t164 * t117 + t165 * t118 + t119 * t167 + t120 * t168 + t132 * t314 + t131 * t313 + m(7) * (t108 * t120 + t109 * t119 + t164 * t66 + t165 * t65 + t195 * t226) + (m(6) * (qJ(5) * t157 + qJD(5) * t174) + t312 + t157 * mrSges(6,3) + qJ(5) * t193 + qJD(5) * t200) * t235 + (m(6) * (-qJ(5) * t156 - qJD(5) * t173) + t311 - t156 * mrSges(6,3) - qJ(5) * t194 - qJD(5) * t201) * t231 + (t108 * t188 - t109 * t189 - t199 * t66 - t247 * t65) * mrSges(7,3) + ((pkin(10) * mrSges(5,2) - Ifges(5,6) + t254) * t240 + (t212 * t304 + t213 * t303 + (-m(6) * pkin(4) + t286) * pkin(10)) * t243) * qJD(4); 0.2e1 * t226 * t141 - t188 * t155 + t199 * t144 - t189 * t154 - t247 * t143 + (t119 * t165 + t120 * t164) * t323 + 0.2e1 * (-t119 * t247 - t120 * t199 + t164 * t188 - t165 * t189) * mrSges(7,3) + (qJ(5) * t324 + 0.2e1 * mrSges(6,3)) * qJD(5) * (t231 ^ 2 + t235 ^ 2); m(7) * t14 + t10 - t326; 0.2e1 * (m(7) / 0.2e1 + t320) * t162; m(6) * pkin(10) * t264 + m(7) * t195 + t272; t141; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t64 - mrSges(7,2) * t63; mrSges(7,1) * t66 - mrSges(7,2) * t65 + t99; mrSges(7,1) * t120 - t119 * mrSges(7,2) + t142; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t18(1) t18(2) t18(4) t18(7) t18(11) t18(16); t18(2) t18(3) t18(5) t18(8) t18(12) t18(17); t18(4) t18(5) t18(6) t18(9) t18(13) t18(18); t18(7) t18(8) t18(9) t18(10) t18(14) t18(19); t18(11) t18(12) t18(13) t18(14) t18(15) t18(20); t18(16) t18(17) t18(18) t18(19) t18(20) t18(21);];
Mq  = res;
