% Calculate time derivative of joint inertia matrix for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR11_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:50
% EndTime: 2019-12-31 22:39:03
% DurationCPUTime: 5.13s
% Computational Cost: add. (7114->570), mult. (18938->856), div. (0->0), fcn. (17731->10), ass. (0->242)
t210 = sin(qJ(4));
t214 = cos(qJ(4));
t211 = sin(qJ(3));
t244 = qJD(4) * t211;
t215 = cos(qJ(3));
t246 = qJD(3) * t215;
t218 = -t210 * t244 + t214 * t246;
t183 = -pkin(3) * t215 - pkin(9) * t211 - pkin(2);
t251 = t214 * t215;
t200 = pkin(8) * t251;
t145 = t210 * t183 + t200;
t294 = qJD(4) * t145;
t209 = sin(qJ(5));
t213 = cos(qJ(5));
t223 = t209 * t210 - t213 * t214;
t272 = -t223 / 0.2e1;
t169 = t209 * t214 + t210 * t213;
t271 = t169 / 0.2e1;
t293 = t210 / 0.2e1;
t267 = t214 / 0.2e1;
t156 = t223 * t211;
t207 = sin(pkin(5));
t212 = sin(qJ(2));
t255 = t207 * t212;
t196 = pkin(7) * t255;
t208 = cos(pkin(5));
t216 = cos(qJ(2));
t266 = pkin(1) * t216;
t161 = t208 * t266 - t196;
t184 = -mrSges(5,1) * t214 + mrSges(5,2) * t210;
t292 = -m(5) * pkin(3) + t184;
t291 = qJD(4) + qJD(5);
t290 = 0.2e1 * m(5);
t289 = 2 * m(6);
t288 = 0.2e1 * pkin(8);
t287 = -2 * mrSges(3,3);
t157 = -t208 * t215 + t211 * t255;
t249 = qJD(2) * t207;
t232 = t216 * t249;
t122 = -qJD(3) * t157 + t215 * t232;
t158 = t208 * t211 + t215 * t255;
t254 = t207 * t216;
t219 = -t158 * t214 + t210 * t254;
t248 = qJD(2) * t212;
t233 = t207 * t248;
t58 = qJD(4) * t219 - t122 * t210 + t214 * t233;
t285 = t58 / 0.2e1;
t123 = -t158 * t210 - t214 * t254;
t67 = t123 * t213 + t209 * t219;
t284 = t67 / 0.2e1;
t68 = t123 * t209 - t213 * t219;
t283 = t68 / 0.2e1;
t282 = -pkin(10) - pkin(9);
t115 = t291 * t223;
t281 = -t115 / 0.2e1;
t116 = t291 * t169;
t280 = -t116 / 0.2e1;
t119 = Ifges(6,4) * t169 - Ifges(6,2) * t223;
t279 = t119 / 0.2e1;
t120 = Ifges(6,1) * t169 - Ifges(6,4) * t223;
t278 = t120 / 0.2e1;
t277 = t123 / 0.2e1;
t276 = -t219 / 0.2e1;
t262 = Ifges(5,4) * t210;
t226 = Ifges(5,1) * t214 - t262;
t151 = -Ifges(5,5) * t215 + t211 * t226;
t275 = t151 / 0.2e1;
t155 = t169 * t211;
t274 = -t155 / 0.2e1;
t273 = -t156 / 0.2e1;
t261 = Ifges(5,4) * t214;
t188 = Ifges(5,1) * t210 + t261;
t270 = t188 / 0.2e1;
t269 = -t210 / 0.2e1;
t268 = -t214 / 0.2e1;
t265 = pkin(8) * t210;
t146 = t196 + (-pkin(2) - t266) * t208;
t87 = t157 * pkin(3) - t158 * pkin(9) + t146;
t162 = t208 * t212 * pkin(1) + pkin(7) * t254;
t147 = pkin(8) * t208 + t162;
t148 = (-pkin(2) * t216 - pkin(8) * t212 - pkin(1)) * t207;
t95 = t215 * t147 + t211 * t148;
t89 = -pkin(9) * t254 + t95;
t39 = t210 * t87 + t214 * t89;
t264 = Ifges(4,4) * t211;
t263 = Ifges(4,4) * t215;
t260 = Ifges(5,6) * t210;
t153 = t161 * qJD(2);
t259 = t153 * mrSges(3,2);
t154 = t162 * qJD(2);
t258 = t154 * mrSges(3,1);
t257 = t154 * mrSges(4,1);
t256 = t154 * mrSges(4,2);
t253 = t210 * t211;
t252 = t211 * t214;
t70 = -Ifges(6,5) * t115 - Ifges(6,6) * t116;
t181 = (pkin(3) * t211 - pkin(9) * t215) * qJD(3);
t247 = qJD(3) * t211;
t250 = t214 * t181 + t247 * t265;
t245 = qJD(4) * t210;
t243 = qJD(4) * t214;
t242 = qJD(5) * t209;
t241 = qJD(5) * t213;
t121 = qJD(3) * t158 + t211 * t232;
t59 = qJD(4) * t123 + t122 * t214 + t210 * t233;
t18 = qJD(5) * t67 + t209 * t58 + t213 * t59;
t19 = -qJD(5) * t68 - t209 * t59 + t213 * t58;
t4 = Ifges(6,5) * t18 + Ifges(6,6) * t19 + Ifges(6,3) * t121;
t23 = Ifges(5,5) * t59 + Ifges(5,6) * t58 + Ifges(5,3) * t121;
t78 = -t116 * t211 - t223 * t246;
t79 = t291 * t156 - t169 * t246;
t34 = Ifges(6,5) * t78 + Ifges(6,6) * t79 + Ifges(6,3) * t247;
t239 = pkin(4) * t245;
t238 = Ifges(4,6) * t254;
t54 = -Ifges(5,1) * t219 + Ifges(5,4) * t123 + Ifges(5,5) * t157;
t237 = t54 * t267;
t205 = Ifges(5,5) * t243;
t236 = t70 / 0.2e1 - Ifges(5,6) * t245 / 0.2e1 + t205 / 0.2e1;
t235 = Ifges(4,5) * t122 - Ifges(4,6) * t121 + Ifges(4,3) * t233;
t234 = qJD(4) * t282;
t229 = Ifges(5,5) * t293 + Ifges(6,5) * t271 + Ifges(5,6) * t267 + Ifges(6,6) * t272;
t38 = -t210 * t89 + t214 * t87;
t167 = t214 * t183;
t144 = -t215 * t265 + t167;
t92 = t210 * t181 + t183 * t243 + (-t214 * t247 - t215 * t245) * pkin(8);
t228 = -qJD(4) * t144 + t92;
t94 = -t211 * t147 + t148 * t215;
t88 = pkin(3) * t254 - t94;
t227 = mrSges(5,1) * t210 + mrSges(5,2) * t214;
t225 = -Ifges(5,2) * t210 + t261;
t186 = Ifges(5,2) * t214 + t262;
t152 = (pkin(2) * t212 - pkin(8) * t216) * t249;
t45 = -t147 * t247 + t148 * t246 + t211 * t152 + t215 * t153;
t43 = pkin(9) * t233 + t45;
t55 = t121 * pkin(3) - t122 * pkin(9) + t154;
t12 = t210 * t55 + t214 * t43 + t87 * t243 - t245 * t89;
t13 = -qJD(4) * t39 - t210 * t43 + t214 * t55;
t224 = t12 * t214 - t13 * t210;
t27 = pkin(4) * t157 + pkin(10) * t219 + t38;
t32 = pkin(10) * t123 + t39;
t10 = -t209 * t32 + t213 * t27;
t11 = t209 * t27 + t213 * t32;
t108 = -pkin(10) * t252 + t167 + (-pkin(4) - t265) * t215;
t125 = -pkin(10) * t253 + t145;
t65 = t108 * t213 - t125 * t209;
t66 = t108 * t209 + t125 * t213;
t190 = t282 * t210;
t191 = t282 * t214;
t128 = t190 * t213 + t191 * t209;
t129 = t190 * t209 - t191 * t213;
t179 = t210 * t234;
t180 = t214 * t234;
t82 = qJD(5) * t128 + t179 * t213 + t180 * t209;
t83 = -qJD(5) * t129 - t179 * t209 + t180 * t213;
t222 = t83 * mrSges(6,1) - t82 * mrSges(6,2) + t70;
t46 = -t147 * t246 - t148 * t247 + t152 * t215 - t211 * t153;
t8 = pkin(4) * t121 - pkin(10) * t59 + t13;
t9 = pkin(10) * t58 + t12;
t2 = qJD(5) * t10 + t209 * t8 + t213 * t9;
t3 = -qJD(5) * t11 - t209 * t9 + t213 * t8;
t221 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t4;
t62 = (pkin(4) * t211 - pkin(10) * t251) * qJD(3) + (-t200 + (pkin(10) * t211 - t183) * t210) * qJD(4) + t250;
t217 = t210 * t246 + t211 * t243;
t73 = -pkin(10) * t217 + t92;
t21 = qJD(5) * t65 + t209 * t62 + t213 * t73;
t22 = -qJD(5) * t66 - t209 * t73 + t213 * t62;
t220 = t22 * mrSges(6,1) - t21 * mrSges(6,2) + t34;
t44 = -pkin(3) * t233 - t46;
t103 = Ifges(5,5) * t218 - Ifges(5,6) * t217 + Ifges(5,3) * t247;
t206 = Ifges(4,5) * t246;
t202 = -pkin(4) * t214 - pkin(3);
t193 = Ifges(3,5) * t232;
t189 = Ifges(4,1) * t211 + t263;
t187 = Ifges(4,2) * t215 + t264;
t182 = (pkin(4) * t210 + pkin(8)) * t211;
t178 = -mrSges(5,1) * t215 - mrSges(5,3) * t252;
t177 = mrSges(5,2) * t215 - mrSges(5,3) * t253;
t176 = (Ifges(4,1) * t215 - t264) * qJD(3);
t175 = t226 * qJD(4);
t174 = (-Ifges(4,2) * t211 + t263) * qJD(3);
t173 = t225 * qJD(4);
t171 = (mrSges(4,1) * t211 + mrSges(4,2) * t215) * qJD(3);
t170 = t227 * qJD(4);
t165 = (-mrSges(6,1) * t209 - mrSges(6,2) * t213) * qJD(5) * pkin(4);
t163 = t227 * t211;
t150 = -Ifges(5,6) * t215 + t211 * t225;
t149 = -Ifges(5,3) * t215 + (Ifges(5,5) * t214 - t260) * t211;
t140 = pkin(4) * t217 + pkin(8) * t246;
t136 = -mrSges(5,2) * t247 - mrSges(5,3) * t217;
t135 = mrSges(5,1) * t247 - mrSges(5,3) * t218;
t131 = -mrSges(6,1) * t215 + mrSges(6,3) * t156;
t130 = mrSges(6,2) * t215 - mrSges(6,3) * t155;
t127 = -mrSges(4,1) * t254 - t158 * mrSges(4,3);
t126 = mrSges(4,2) * t254 - t157 * mrSges(4,3);
t117 = mrSges(6,1) * t223 + mrSges(6,2) * t169;
t107 = mrSges(5,1) * t217 + mrSges(5,2) * t218;
t106 = mrSges(6,1) * t155 - mrSges(6,2) * t156;
t105 = -t188 * t244 + (Ifges(5,5) * t211 + t215 * t226) * qJD(3);
t104 = -t186 * t244 + (Ifges(5,6) * t211 + t215 * t225) * qJD(3);
t102 = mrSges(4,1) * t233 - mrSges(4,3) * t122;
t101 = -mrSges(4,2) * t233 - mrSges(4,3) * t121;
t100 = Ifges(4,1) * t158 - Ifges(4,4) * t157 - Ifges(4,5) * t254;
t99 = Ifges(4,4) * t158 - Ifges(4,2) * t157 - t238;
t98 = -Ifges(6,1) * t156 - Ifges(6,4) * t155 - Ifges(6,5) * t215;
t97 = -Ifges(6,4) * t156 - Ifges(6,2) * t155 - Ifges(6,6) * t215;
t96 = -Ifges(6,5) * t156 - Ifges(6,6) * t155 - Ifges(6,3) * t215;
t93 = t250 - t294;
t91 = mrSges(5,1) * t157 + mrSges(5,3) * t219;
t90 = -mrSges(5,2) * t157 + mrSges(5,3) * t123;
t77 = -mrSges(5,1) * t123 - mrSges(5,2) * t219;
t74 = mrSges(4,1) * t121 + mrSges(4,2) * t122;
t72 = -Ifges(6,1) * t115 - Ifges(6,4) * t116;
t71 = -Ifges(6,4) * t115 - Ifges(6,2) * t116;
t69 = mrSges(6,1) * t116 - mrSges(6,2) * t115;
t64 = -mrSges(6,2) * t247 + mrSges(6,3) * t79;
t63 = mrSges(6,1) * t247 - mrSges(6,3) * t78;
t61 = Ifges(4,1) * t122 - Ifges(4,4) * t121 + Ifges(4,5) * t233;
t60 = Ifges(4,4) * t122 - Ifges(4,2) * t121 + Ifges(4,6) * t233;
t53 = -Ifges(5,4) * t219 + Ifges(5,2) * t123 + Ifges(5,6) * t157;
t52 = -Ifges(5,5) * t219 + Ifges(5,6) * t123 + Ifges(5,3) * t157;
t49 = -pkin(4) * t123 + t88;
t48 = mrSges(6,1) * t157 - mrSges(6,3) * t68;
t47 = -mrSges(6,2) * t157 + mrSges(6,3) * t67;
t41 = mrSges(5,1) * t121 - mrSges(5,3) * t59;
t40 = -mrSges(5,2) * t121 + mrSges(5,3) * t58;
t37 = -mrSges(6,1) * t79 + mrSges(6,2) * t78;
t36 = Ifges(6,1) * t78 + Ifges(6,4) * t79 + Ifges(6,5) * t247;
t35 = Ifges(6,4) * t78 + Ifges(6,2) * t79 + Ifges(6,6) * t247;
t33 = -mrSges(6,1) * t67 + mrSges(6,2) * t68;
t31 = Ifges(6,1) * t68 + Ifges(6,4) * t67 + Ifges(6,5) * t157;
t30 = Ifges(6,4) * t68 + Ifges(6,2) * t67 + Ifges(6,6) * t157;
t29 = Ifges(6,5) * t68 + Ifges(6,6) * t67 + Ifges(6,3) * t157;
t28 = -mrSges(5,1) * t58 + mrSges(5,2) * t59;
t26 = -pkin(4) * t58 + t44;
t25 = Ifges(5,1) * t59 + Ifges(5,4) * t58 + Ifges(5,5) * t121;
t24 = Ifges(5,4) * t59 + Ifges(5,2) * t58 + Ifges(5,6) * t121;
t15 = -mrSges(6,2) * t121 + mrSges(6,3) * t19;
t14 = mrSges(6,1) * t121 - mrSges(6,3) * t18;
t7 = -mrSges(6,1) * t19 + mrSges(6,2) * t18;
t6 = Ifges(6,1) * t18 + Ifges(6,4) * t19 + Ifges(6,5) * t121;
t5 = Ifges(6,4) * t18 + Ifges(6,2) * t19 + Ifges(6,6) * t121;
t1 = [(t10 * t3 + t11 * t2 + t26 * t49) * t289 + (t12 * t39 + t13 * t38 + t44 * t88) * t290 + (-t216 * t235 + 0.2e1 * (t153 * t216 + t154 * t212) * mrSges(3,3) + ((t161 * t287 + Ifges(3,5) * t208 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t216) * t207) * t216 + (t162 * t287 + Ifges(4,5) * t158 - 0.2e1 * Ifges(3,6) * t208 - Ifges(4,6) * t157 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t212 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t216) * t207) * t212) * qJD(2)) * t207 - t219 * t25 + 0.2e1 * m(4) * (t146 * t154 + t45 * t95 + t46 * t94) + 0.2e1 * m(3) * (t153 * t162 - t154 * t161) + (t29 + t52 - t99) * t121 + 0.2e1 * t10 * t14 + 0.2e1 * t11 * t15 + (t61 + 0.2e1 * t256) * t158 + (t23 + t4 - t60 + 0.2e1 * t257) * t157 + (t193 - 0.2e1 * t258 - 0.2e1 * t259) * t208 + t19 * t30 + t18 * t31 + 0.2e1 * t26 * t33 + 0.2e1 * t39 * t40 + 0.2e1 * t38 * t41 + 0.2e1 * t2 * t47 + 0.2e1 * t3 * t48 + 0.2e1 * t49 * t7 + t58 * t53 + t59 * t54 + t67 * t5 + t68 * t6 + 0.2e1 * t44 * t77 + 0.2e1 * t88 * t28 + 0.2e1 * t12 * t90 + 0.2e1 * t13 * t91 + 0.2e1 * t95 * t101 + 0.2e1 * t94 * t102 + t122 * t100 + t123 * t24 + 0.2e1 * t45 * t126 + 0.2e1 * t46 * t127 + 0.2e1 * t146 * t74; (-m(4) * t154 - t74) * pkin(2) + t6 * t273 + t5 * t274 + t59 * t275 + t105 * t276 + t104 * t277 + t36 * t283 + t35 * t284 + t150 * t285 - t258 - t259 + m(6) * (t10 * t22 + t11 * t21 + t140 * t49 + t182 * t26 + t2 * t66 + t3 * t65) + m(5) * (t12 * t145 + t13 * t144 + t38 * t93 + t39 * t92) + (t34 / 0.2e1 + t103 / 0.2e1 - t174 / 0.2e1) * t157 + (t96 / 0.2e1 + t149 / 0.2e1 - t187 / 0.2e1) * t121 + t193 + (t215 * t101 + (-t102 + t28) * t211 + (-t211 * t126 + (-t127 + t77) * t215) * qJD(3) + m(4) * (-t46 * t211 + t45 * t215 - t246 * t94 - t247 * t95) + m(5) * (t211 * t44 + t246 * t88)) * pkin(8) + (-t216 * t206 / 0.2e1 + (Ifges(4,5) * t211 / 0.2e1 + Ifges(4,6) * t215 / 0.2e1 - Ifges(3,6)) * t248) * t207 + (t45 * mrSges(4,3) - t4 / 0.2e1 - t23 / 0.2e1 + t60 / 0.2e1 - t257) * t215 + ((t100 / 0.2e1 + t237 - t94 * mrSges(4,3) + t53 * t269) * t215 + (t238 / 0.2e1 - t99 / 0.2e1 + t52 / 0.2e1 + t29 / 0.2e1 - t95 * mrSges(4,3)) * t211) * qJD(3) + (t25 * t267 + t24 * t269 - t46 * mrSges(4,3) + t61 / 0.2e1 + t256 + (t268 * t53 + t269 * t54) * qJD(4)) * t211 + t21 * t47 + t22 * t48 + t49 * t37 + t10 * t63 + t11 * t64 + t65 * t14 + t66 * t15 + t78 * t31 / 0.2e1 + t79 * t30 / 0.2e1 + t92 * t90 + t93 * t91 + t19 * t97 / 0.2e1 + t18 * t98 / 0.2e1 + t26 * t106 + t88 * t107 + t2 * t130 + t3 * t131 + t38 * t135 + t39 * t136 + t140 * t33 + t144 * t41 + t145 * t40 + t44 * t163 + t146 * t171 + t158 * t176 / 0.2e1 + t12 * t177 + t13 * t178 + t182 * t7 + t122 * t189 / 0.2e1; -0.2e1 * pkin(2) * t171 + 0.2e1 * t140 * t106 + 0.2e1 * t21 * t130 + 0.2e1 * t22 * t131 + 0.2e1 * t144 * t135 + 0.2e1 * t145 * t136 - t155 * t35 - t156 * t36 + 0.2e1 * t92 * t177 + 0.2e1 * t93 * t178 + 0.2e1 * t182 * t37 + 0.2e1 * t65 * t63 + 0.2e1 * t66 * t64 + t78 * t98 + t79 * t97 + (t144 * t93 + t145 * t92) * t290 + (t140 * t182 + t21 * t66 + t22 * t65) * t289 + (-t103 + t174 - t34 + (-t210 * t150 + t214 * t151 + t163 * t288 + t189) * qJD(3)) * t215 + (t107 * t288 - t210 * t104 + t214 * t105 + t176 + (-t150 * t214 - t151 * t210) * qJD(4) + (pkin(8) ^ 2 * t215 * t290 + t149 - t187 + t96) * qJD(3)) * t211; t292 * t44 + t175 * t276 + t173 * t277 + t18 * t278 + t19 * t279 + t30 * t280 + t31 * t281 + t72 * t283 + t71 * t284 + t186 * t285 + t24 * t267 + t59 * t270 + t6 * t271 + t5 * t272 + ((-t210 * t39 - t214 * t38) * qJD(4) + t224) * mrSges(5,3) + t235 + t25 * t293 + (t10 * t115 - t11 * t116 - t169 * t3 - t2 * t223) * mrSges(6,3) + t229 * t121 + t236 * t157 + (t237 + (-t53 / 0.2e1 + pkin(4) * t33) * t210) * qJD(4) + m(6) * (t10 * t83 + t11 * t82 + t128 * t3 + t129 * t2 + t202 * t26 + t239 * t49) + (m(5) * (-t243 * t38 - t245 * t39 + t224) + t214 * t40 - t210 * t41 - t91 * t243 - t90 * t245) * pkin(9) - pkin(3) * t28 - t45 * mrSges(4,2) + t46 * mrSges(4,1) + t49 * t69 + t82 * t47 + t83 * t48 + t26 * t117 + t128 * t14 + t129 * t15 + t88 * t170 + t202 * t7; t206 + m(6) * (t128 * t22 + t129 * t21 + t140 * t202 + t65 * t83 + t66 * t82) - pkin(3) * t107 + t98 * t281 + t97 * t280 + t79 * t279 + t78 * t278 + t128 * t63 + t129 * t64 + t82 * t130 + t83 * t131 + t140 * t117 + t71 * t274 + t72 * t273 + t35 * t272 + t36 * t271 + t182 * t69 + t202 * t37 + ((-mrSges(4,1) + t292) * qJD(3) * pkin(8) - t236) * t215 + (-t93 * mrSges(5,3) - t186 * t246 / 0.2e1 + t105 / 0.2e1 + (-t150 / 0.2e1 - t145 * mrSges(5,3) + (m(6) * t182 + t106) * pkin(4)) * qJD(4) + (m(5) * (-t93 - t294) - t135 - qJD(4) * t177) * pkin(9)) * t210 + (t115 * t65 - t116 * t66 - t169 * t22 - t21 * t223) * mrSges(6,3) + (qJD(4) * t275 + t246 * t270 + t104 / 0.2e1 + t228 * mrSges(5,3) + (m(5) * t228 - qJD(4) * t178 + t136) * pkin(9)) * t214 + (t175 * t267 + t173 * t269 + pkin(8) * t170 + (t186 * t268 + t188 * t269) * qJD(4) + (pkin(8) * mrSges(4,2) - Ifges(4,6) + t229) * qJD(3)) * t211; (t128 * t83 + t129 * t82 + t202 * t239) * t289 - t115 * t120 + t169 * t72 - t116 * t119 - t223 * t71 + 0.2e1 * t117 * t239 + 0.2e1 * t202 * t69 + t210 * t175 - t186 * t245 - 0.2e1 * pkin(3) * t170 + (qJD(4) * t188 + t173) * t214 + 0.2e1 * (t115 * t128 - t116 * t129 - t169 * t83 - t223 * t82) * mrSges(6,3); t13 * mrSges(5,1) - t12 * mrSges(5,2) + (-t48 * t242 + t213 * t14 + m(6) * (-t10 * t242 + t11 * t241 + t2 * t209 + t213 * t3) + t47 * t241 + t209 * t15) * pkin(4) + t221 + t23; t93 * mrSges(5,1) - t92 * mrSges(5,2) + (t130 * t241 + t209 * t64 + m(6) * (t209 * t21 + t213 * t22 + t241 * t66 - t242 * t65) - t131 * t242 + t213 * t63) * pkin(4) + t103 + t220; t205 + (pkin(9) * t184 - t260) * qJD(4) + (m(6) * (t209 * t82 + t213 * t83 + (-t128 * t209 + t129 * t213) * qJD(5)) + (t213 * t115 - t209 * t116 + (t169 * t209 - t213 * t223) * qJD(5)) * mrSges(6,3)) * pkin(4) + t222; 0.2e1 * t165; t221; t220; t222; t165; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
