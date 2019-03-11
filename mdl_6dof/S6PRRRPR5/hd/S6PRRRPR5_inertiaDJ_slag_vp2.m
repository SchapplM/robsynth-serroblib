% Calculate time derivative of joint inertia matrix for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:23:08
% EndTime: 2019-03-08 23:23:20
% DurationCPUTime: 5.18s
% Computational Cost: add. (7541->585), mult. (21306->879), div. (0->0), fcn. (21255->14), ass. (0->248)
t304 = -Ifges(5,3) - Ifges(6,3);
t200 = cos(qJ(4));
t190 = sin(pkin(13));
t196 = sin(qJ(4));
t253 = t190 * t196;
t259 = cos(pkin(13));
t205 = t259 * t200 - t253;
t151 = t205 * qJD(4);
t218 = t259 * t196;
t157 = t190 * t200 + t218;
t199 = cos(qJ(6));
t195 = sin(qJ(6));
t236 = qJD(6) * t195;
t206 = -t199 * t151 + t157 * t236;
t303 = -m(5) * pkin(3) - mrSges(5,1) * t200 + mrSges(5,2) * t196 - mrSges(4,1);
t281 = t195 / 0.2e1;
t279 = t199 / 0.2e1;
t221 = -t236 / 0.2e1;
t275 = pkin(4) * t190;
t185 = pkin(11) + t275;
t302 = m(7) * t185;
t193 = cos(pkin(7));
t197 = sin(qJ(3));
t191 = sin(pkin(7));
t201 = cos(qJ(3));
t251 = t191 * t201;
t155 = t193 * t197 * pkin(2) + pkin(9) * t251;
t137 = pkin(10) * t193 + t155;
t138 = (-pkin(3) * t201 - pkin(10) * t197 - pkin(2)) * t191;
t94 = t200 * t137 + t196 * t138;
t252 = t191 * t197;
t181 = pkin(9) * t252;
t276 = pkin(2) * t201;
t154 = t193 * t276 - t181;
t150 = t157 * qJD(4);
t237 = qJD(4) * t200;
t238 = qJD(4) * t196;
t301 = Ifges(5,5) * t237 + Ifges(6,5) * t151 - Ifges(5,6) * t238 - Ifges(6,6) * t150;
t202 = cos(qJ(2));
t244 = t201 * t202;
t198 = sin(qJ(2));
t249 = t197 * t198;
t300 = t193 * t244 - t249;
t187 = -pkin(4) * t200 - pkin(3);
t110 = -pkin(5) * t205 - pkin(11) * t157 + t187;
t273 = -qJ(5) - pkin(10);
t171 = t273 * t200;
t124 = -t259 * t171 + t273 * t253;
t67 = t110 * t199 - t124 * t195;
t219 = qJD(4) * t273;
t148 = qJD(5) * t200 + t196 * t219;
t204 = -qJD(5) * t196 + t200 * t219;
t98 = t259 * t148 + t190 * t204;
t233 = pkin(4) * t238;
t99 = pkin(5) * t150 - pkin(11) * t151 + t233;
t20 = t67 * qJD(6) + t195 * t99 + t199 * t98;
t68 = t110 * t195 + t124 * t199;
t21 = -t68 * qJD(6) - t195 * t98 + t199 * t99;
t299 = -t195 * t21 + t199 * t20;
t289 = m(6) * pkin(4);
t298 = t190 * t289 - mrSges(6,2);
t169 = -mrSges(7,1) * t199 + mrSges(7,2) * t195;
t229 = t259 * pkin(4);
t186 = -t229 - pkin(5);
t297 = m(7) * t186 - t259 * t289 - mrSges(6,1) + t169;
t296 = 0.2e1 * m(6);
t295 = 0.2e1 * m(7);
t294 = -2 * mrSges(4,3);
t293 = -2 * mrSges(6,3);
t97 = t148 * t190 - t259 * t204;
t292 = 0.2e1 * t97;
t123 = -t171 * t190 - t273 * t218;
t291 = 0.2e1 * t123;
t239 = qJD(3) * t197;
t227 = t191 * t239;
t152 = t193 * t200 - t196 * t252;
t240 = qJD(3) * t191;
t226 = t201 * t240;
t121 = t152 * qJD(4) + t200 * t226;
t153 = t193 * t196 + t200 * t252;
t122 = -t153 * qJD(4) - t196 * t226;
t75 = t259 * t121 + t190 * t122;
t103 = t190 * t152 + t259 * t153;
t85 = -t103 * t195 - t199 * t251;
t41 = t85 * qJD(6) + t195 * t227 + t199 * t75;
t288 = t41 / 0.2e1;
t209 = -t103 * t199 + t195 * t251;
t42 = t209 * qJD(6) - t195 * t75 + t199 * t227;
t287 = t42 / 0.2e1;
t286 = t85 / 0.2e1;
t235 = qJD(6) * t199;
t188 = Ifges(7,5) * t235;
t285 = Ifges(7,6) * t221 + t188 / 0.2e1;
t268 = Ifges(7,4) * t195;
t213 = Ifges(7,1) * t199 - t268;
t166 = t213 * qJD(6);
t284 = t166 / 0.2e1;
t283 = Ifges(7,5) * t281 + Ifges(7,6) * t279;
t282 = -t195 / 0.2e1;
t280 = t196 / 0.2e1;
t278 = t200 / 0.2e1;
t192 = sin(pkin(6));
t241 = qJD(2) * t192;
t228 = t198 * t241;
t217 = t191 * t228;
t194 = cos(pkin(6));
t79 = t194 * t226 + (t300 * qJD(3) + (-t193 * t249 + t244) * qJD(2)) * t192;
t247 = t198 * t201;
t248 = t197 * t202;
t208 = t193 * t248 + t247;
t115 = t208 * t192 + t194 * t252;
t149 = -t191 * t192 * t202 + t194 * t193;
t84 = t115 * t200 + t149 * t196;
t43 = -t84 * qJD(4) - t196 * t79 + t200 * t217;
t83 = -t115 * t196 + t149 * t200;
t44 = t83 * qJD(4) + t196 * t217 + t200 * t79;
t15 = t190 * t44 - t259 * t43;
t46 = t190 * t84 - t259 * t83;
t274 = t15 * t46;
t14 = -mrSges(7,1) * t42 + mrSges(7,2) * t41;
t62 = mrSges(6,1) * t227 - mrSges(6,3) * t75;
t272 = t14 - t62;
t141 = (pkin(3) * t197 - pkin(10) * t201) * t240;
t142 = t154 * qJD(3);
t55 = -t94 * qJD(4) + t200 * t141 - t142 * t196;
t30 = pkin(4) * t227 - qJ(5) * t121 - qJD(5) * t153 + t55;
t54 = -t137 * t238 + t138 * t237 + t196 * t141 + t200 * t142;
t36 = qJ(5) * t122 + qJD(5) * t152 + t54;
t11 = t190 * t30 + t259 * t36;
t48 = -mrSges(7,1) * t85 - mrSges(7,2) * t209;
t91 = -mrSges(6,1) * t251 - mrSges(6,3) * t103;
t271 = t48 - t91;
t93 = -t137 * t196 + t200 * t138;
t64 = -pkin(4) * t251 - qJ(5) * t153 + t93;
t76 = qJ(5) * t152 + t94;
t35 = t190 * t64 + t259 * t76;
t270 = Ifges(5,4) * t196;
t269 = Ifges(5,4) * t200;
t267 = Ifges(7,4) * t199;
t266 = Ifges(7,6) * t195;
t114 = -t300 * t192 - t194 * t251;
t78 = t194 * t227 + (t208 * qJD(3) + (t193 * t247 + t248) * qJD(2)) * t192;
t58 = t114 * t78;
t265 = t123 * t97;
t264 = t195 * mrSges(7,3);
t262 = t199 * mrSges(7,3);
t143 = t155 * qJD(3);
t258 = t114 * t143;
t257 = t157 * t195;
t256 = t157 * t199;
t255 = t185 * t195;
t254 = t185 * t199;
t173 = Ifges(7,2) * t199 + t268;
t250 = t195 * t173;
t175 = Ifges(7,1) * t195 + t267;
t245 = t199 * t175;
t242 = -mrSges(4,1) * t193 - mrSges(5,1) * t152 + mrSges(5,2) * t153 + mrSges(4,3) * t252;
t74 = t121 * t190 - t259 * t122;
t5 = Ifges(7,5) * t41 + Ifges(7,6) * t42 + Ifges(7,3) * t74;
t232 = mrSges(7,3) * t236;
t231 = mrSges(7,3) * t235;
t225 = t185 * t236;
t224 = t185 * t235;
t222 = t157 * t235;
t39 = t74 * mrSges(6,1) + t75 * mrSges(6,2);
t220 = t235 / 0.2e1;
t105 = t150 * mrSges(6,1) + t151 * mrSges(6,2);
t215 = t123 * t15 + t46 * t97;
t214 = mrSges(7,1) * t195 + mrSges(7,2) * t199;
t212 = -Ifges(7,2) * t195 + t267;
t29 = -pkin(11) * t251 + t35;
t102 = -t259 * t152 + t153 * t190;
t136 = t181 + (-pkin(3) - t276) * t193;
t104 = -pkin(4) * t152 + t136;
t45 = pkin(5) * t102 - pkin(11) * t103 + t104;
t13 = t195 * t45 + t199 * t29;
t12 = -t195 * t29 + t199 * t45;
t211 = -t196 * t55 + t200 * t54;
t47 = t190 * t83 + t259 * t84;
t22 = t114 * t199 - t195 * t47;
t23 = t114 * t195 + t199 * t47;
t210 = -Ifges(5,5) * t121 - Ifges(6,5) * t75 - Ifges(5,6) * t122 + Ifges(6,6) * t74 + t227 * t304;
t10 = -t190 * t36 + t259 * t30;
t34 = -t190 * t76 + t259 * t64;
t207 = t195 * t151 + t222;
t92 = -pkin(4) * t122 + t143;
t49 = -Ifges(7,5) * t206 - t207 * Ifges(7,6) + Ifges(7,3) * t150;
t180 = Ifges(4,5) * t226;
t176 = Ifges(5,1) * t196 + t269;
t174 = Ifges(5,2) * t200 + t270;
t167 = (Ifges(5,1) * t200 - t270) * qJD(4);
t165 = (-Ifges(5,2) * t196 + t269) * qJD(4);
t164 = t212 * qJD(6);
t162 = (mrSges(5,1) * t196 + mrSges(5,2) * t200) * qJD(4);
t161 = t214 * qJD(6);
t160 = -mrSges(4,2) * t193 + mrSges(4,3) * t251;
t140 = (mrSges(4,1) * t197 + mrSges(4,2) * t201) * t240;
t126 = -mrSges(5,1) * t251 - mrSges(5,3) * t153;
t125 = mrSges(5,2) * t251 + mrSges(5,3) * t152;
t118 = Ifges(6,1) * t157 + Ifges(6,4) * t205;
t117 = Ifges(6,4) * t157 + Ifges(6,2) * t205;
t116 = -mrSges(6,1) * t205 + mrSges(6,2) * t157;
t112 = -mrSges(7,1) * t205 - mrSges(7,3) * t256;
t111 = mrSges(7,2) * t205 - mrSges(7,3) * t257;
t108 = t214 * t157;
t107 = Ifges(6,1) * t151 - Ifges(6,4) * t150;
t106 = Ifges(6,4) * t151 - Ifges(6,2) * t150;
t101 = -mrSges(5,2) * t227 + mrSges(5,3) * t122;
t100 = mrSges(5,1) * t227 - mrSges(5,3) * t121;
t96 = Ifges(5,1) * t153 + Ifges(5,4) * t152 - Ifges(5,5) * t251;
t95 = Ifges(5,4) * t153 + Ifges(5,2) * t152 - Ifges(5,6) * t251;
t90 = mrSges(6,2) * t251 - mrSges(6,3) * t102;
t89 = -Ifges(7,5) * t205 + t213 * t157;
t88 = -Ifges(7,6) * t205 + t212 * t157;
t87 = -Ifges(7,3) * t205 + (Ifges(7,5) * t199 - t266) * t157;
t82 = -mrSges(7,2) * t150 - t207 * mrSges(7,3);
t81 = mrSges(7,1) * t150 + t206 * mrSges(7,3);
t77 = -mrSges(5,1) * t122 + mrSges(5,2) * t121;
t66 = Ifges(5,1) * t121 + Ifges(5,4) * t122 + Ifges(5,5) * t227;
t65 = Ifges(5,4) * t121 + Ifges(5,2) * t122 + Ifges(5,6) * t227;
t61 = -mrSges(6,2) * t227 - mrSges(6,3) * t74;
t60 = t207 * mrSges(7,1) - t206 * mrSges(7,2);
t59 = mrSges(6,1) * t102 + mrSges(6,2) * t103;
t57 = Ifges(6,1) * t103 - Ifges(6,4) * t102 - Ifges(6,5) * t251;
t56 = Ifges(6,4) * t103 - Ifges(6,2) * t102 - Ifges(6,6) * t251;
t53 = mrSges(7,1) * t102 + mrSges(7,3) * t209;
t52 = -mrSges(7,2) * t102 + mrSges(7,3) * t85;
t51 = -t206 * Ifges(7,1) - t207 * Ifges(7,4) + Ifges(7,5) * t150;
t50 = -t206 * Ifges(7,4) - t207 * Ifges(7,2) + Ifges(7,6) * t150;
t32 = Ifges(6,1) * t75 - Ifges(6,4) * t74 + Ifges(6,5) * t227;
t31 = Ifges(6,4) * t75 - Ifges(6,2) * t74 + Ifges(6,6) * t227;
t28 = pkin(5) * t251 - t34;
t27 = -Ifges(7,1) * t209 + Ifges(7,4) * t85 + Ifges(7,5) * t102;
t26 = -Ifges(7,4) * t209 + Ifges(7,2) * t85 + Ifges(7,6) * t102;
t25 = -Ifges(7,5) * t209 + Ifges(7,6) * t85 + Ifges(7,3) * t102;
t19 = pkin(5) * t74 - pkin(11) * t75 + t92;
t18 = -mrSges(7,2) * t74 + mrSges(7,3) * t42;
t17 = mrSges(7,1) * t74 - mrSges(7,3) * t41;
t16 = t190 * t43 + t259 * t44;
t9 = pkin(11) * t227 + t11;
t8 = -pkin(5) * t227 - t10;
t7 = Ifges(7,1) * t41 + Ifges(7,4) * t42 + Ifges(7,5) * t74;
t6 = Ifges(7,4) * t41 + Ifges(7,2) * t42 + Ifges(7,6) * t74;
t4 = -t23 * qJD(6) - t16 * t195 + t199 * t78;
t3 = t22 * qJD(6) + t16 * t199 + t195 * t78;
t2 = -t13 * qJD(6) + t19 * t199 - t195 * t9;
t1 = t12 * qJD(6) + t19 * t195 + t199 * t9;
t24 = [0.2e1 * m(7) * (t22 * t4 + t23 * t3 + t274) + 0.2e1 * m(6) * (t16 * t47 + t274 + t58) + 0.2e1 * m(5) * (t43 * t83 + t44 * t84 + t58) + 0.2e1 * m(4) * (t115 * t79 + t149 * t217 + t58); t83 * t100 + t84 * t101 + t44 * t125 + t43 * t126 + t149 * t140 + t16 * t90 + t79 * t160 + t22 * t17 + t23 * t18 + t3 * t52 + t4 * t53 + t47 * t61 + t272 * t46 + t271 * t15 + (t39 + t77) * t114 + (-mrSges(3,1) * t198 - mrSges(3,2) * t202) * t241 + (t59 + t242) * t78 + ((-mrSges(4,1) * t201 + mrSges(4,2) * t197) * t217 + (t114 * t201 - t115 * t197) * qJD(3) * mrSges(4,3)) * t191 + m(4) * (-pkin(2) * t191 ^ 2 * t228 + t115 * t142 - t154 * t78 + t155 * t79 + t258) + m(5) * (t136 * t78 + t43 * t93 + t44 * t94 + t54 * t84 + t55 * t83 + t258) + m(6) * (-t10 * t46 + t104 * t78 + t11 * t47 + t114 * t92 - t15 * t34 + t16 * t35) + m(7) * (t1 * t23 + t12 * t4 + t13 * t3 + t15 * t28 + t2 * t22 + t46 * t8); 0.2e1 * t242 * t143 - t209 * t7 + (t1 * t13 + t12 * t2 + t28 * t8) * t295 + (t10 * t34 + t104 * t92 + t11 * t35) * t296 + (t5 - t31) * t102 + (-0.2e1 * pkin(2) * t140 + t210 * t201 + ((0.2e1 * Ifges(4,4) * t251 + Ifges(4,5) * t193 + t154 * t294) * t201 + (-0.2e1 * Ifges(4,4) * t252 + t155 * t294 + Ifges(5,5) * t153 + Ifges(6,5) * t103 - 0.2e1 * Ifges(4,6) * t193 + Ifges(5,6) * t152 - Ifges(6,6) * t102 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) + t304) * t251) * t197) * qJD(3)) * t191 + 0.2e1 * m(5) * (t136 * t143 + t54 * t94 + t55 * t93) + 0.2e1 * m(4) * (t142 * t155 - t143 * t154) + (t25 - t56) * t74 + t193 * t180 + 0.2e1 * t142 * t160 + t152 * t65 + t153 * t66 + 0.2e1 * t12 * t17 + 0.2e1 * t13 * t18 + 0.2e1 * t28 * t14 + t41 * t27 + t42 * t26 + 0.2e1 * t8 * t48 + 0.2e1 * t1 * t52 + 0.2e1 * t2 * t53 + 0.2e1 * t35 * t61 + 0.2e1 * t34 * t62 + t75 * t57 + t85 * t6 + 0.2e1 * t11 * t90 + 0.2e1 * t10 * t91 + 0.2e1 * t92 * t59 + 0.2e1 * t93 * t100 + 0.2e1 * t94 * t101 + t103 * t32 + 0.2e1 * t104 * t39 + t121 * t96 + t122 * t95 + 0.2e1 * t54 * t125 + 0.2e1 * t55 * t126 + 0.2e1 * t136 * t77; -t79 * mrSges(4,2) + t15 * t108 + t3 * t111 + t4 * t112 + t22 * t81 + t23 * t82 + t46 * t60 + (t105 + t162) * t114 + m(7) * (t20 * t23 + t21 * t22 + t3 * t68 + t4 * t67 + t215) + m(6) * (t114 * t233 + t124 * t16 + t47 * t98 + t215) + (t15 * t157 - t150 * t47 + t151 * t46 + t16 * t205) * mrSges(6,3) + (m(6) * t187 + t116 + t303) * t78 + (m(5) * pkin(10) + mrSges(5,3)) * (-t43 * t196 + t44 * t200 + (-t196 * t84 - t200 * t83) * qJD(4)); t180 + (t96 * t278 + (pkin(4) * t59 - t95 / 0.2e1) * t196) * qJD(4) + (t27 * t279 + t26 * t282 - t34 * mrSges(6,3) + t57 / 0.2e1) * t151 + (t7 * t279 + t6 * t282 - t10 * mrSges(6,3) + t32 / 0.2e1 + (-t199 * t26 / 0.2e1 + t27 * t282) * qJD(6)) * t157 - t209 * t51 / 0.2e1 + ((-Ifges(4,6) + Ifges(5,5) * t280 + Ifges(5,6) * t278 + Ifges(6,5) * t157 / 0.2e1 + Ifges(6,6) * t205 / 0.2e1) * t239 - t301 * t201 / 0.2e1) * t191 - (-t11 * mrSges(6,3) + t5 / 0.2e1 - t31 / 0.2e1) * t205 + (-t126 * t237 - t125 * t238 + m(5) * (-t93 * t237 - t94 * t238 + t211) + t200 * t101 - t196 * t100) * pkin(10) + m(6) * (-t10 * t123 + t104 * t233 + t11 * t124 + t187 * t92 - t34 * t97 + t35 * t98) + ((-t196 * t94 - t200 * t93) * qJD(4) + t211) * mrSges(5,3) + t50 * t286 + t88 * t287 + t89 * t288 + t65 * t278 + t66 * t280 + (-t35 * mrSges(6,3) + t25 / 0.2e1 - t56 / 0.2e1) * t150 + t271 * t97 + t272 * t123 + (t87 / 0.2e1 - t117 / 0.2e1) * t74 + (t49 / 0.2e1 - t106 / 0.2e1) * t102 + m(7) * (t1 * t68 + t12 * t21 + t123 * t8 + t13 * t20 + t2 * t67 + t28 * t97) + t303 * t143 + t187 * t39 + t152 * t165 / 0.2e1 + t153 * t167 / 0.2e1 + t122 * t174 / 0.2e1 + t121 * t176 / 0.2e1 + t136 * t162 - t142 * mrSges(4,2) + t20 * t52 + t21 * t53 + t28 * t60 + t67 * t17 + t68 * t18 - pkin(3) * t77 + t12 * t81 + t13 * t82 + t98 * t90 + t104 * t105 + t103 * t107 / 0.2e1 + t8 * t108 + t1 * t111 + t2 * t112 + t92 * t116 + t75 * t118 / 0.2e1 + t124 * t61; -0.2e1 * pkin(3) * t162 + 0.2e1 * t187 * t105 + t108 * t292 + 0.2e1 * t20 * t111 + 0.2e1 * t21 * t112 + t60 * t291 + t200 * t165 + t196 * t167 + 0.2e1 * t67 * t81 + 0.2e1 * t68 * t82 + (t200 * t176 + (0.2e1 * pkin(4) * t116 - t174) * t196) * qJD(4) + (t124 * t98 + t187 * t233 + t265) * t296 + (t20 * t68 + t21 * t67 + t265) * t295 - (t98 * t293 - t106 + t49) * t205 + (t124 * t293 - t117 + t87) * t150 + (mrSges(6,3) * t291 - t195 * t88 + t199 * t89 + t118) * t151 + (mrSges(6,3) * t292 - t195 * t50 + t199 * t51 + t107 + (-t195 * t89 - t199 * t88) * qJD(6)) * t157; -t4 * t264 - t22 * t231 + (-t4 * t195 + t3 * t199 + (-t195 * t23 - t199 * t22) * qJD(6)) * t302 + t3 * t262 - t23 * t232 + t46 * t161 - t44 * mrSges(5,2) + t43 * mrSges(5,1) + t298 * t16 + t297 * t15; -t210 - t17 * t255 - t209 * t284 - t13 * t232 - t12 * t231 - t53 * t224 - t52 * t225 + t27 * t220 + t26 * t221 + t74 * t283 + t102 * t285 + t164 * t286 + t173 * t287 + t175 * t288 + (t259 * t10 + t11 * t190) * t289 + t61 * t275 + t6 * t279 + t7 * t281 + t1 * t262 + t18 * t254 + t62 * t229 - t2 * t264 + m(7) * (t186 * t8 + (t1 * t199 - t195 * t2 + (-t12 * t199 - t13 * t195) * qJD(6)) * t185) + t186 * t14 + t8 * t169 + t28 * t161 + t10 * mrSges(6,1) - t11 * mrSges(6,2) - t54 * mrSges(5,2) + t55 * mrSges(5,1); t301 + t299 * mrSges(7,3) - t81 * t255 - t164 * t257 / 0.2e1 + t298 * t98 + t297 * t97 - t68 * t232 - t67 * t231 - t173 * t222 / 0.2e1 - t112 * t224 - t111 * t225 + t89 * t220 + t256 * t284 - t205 * t285 + t50 * t279 + t51 * t281 + t82 * t254 + ((-t195 * t68 - t199 * t67) * qJD(6) + t299) * t302 + (t157 * t175 + t88) * t221 + (-t250 / 0.2e1 + t245 / 0.2e1 - mrSges(6,3) * t229) * t151 + (-mrSges(6,3) * t275 + t283) * t150 + (-mrSges(5,1) * t237 + mrSges(5,2) * t238) * pkin(10) + t186 * t60 + t123 * t161; 0.2e1 * t161 * t186 + t164 * t199 + t166 * t195 + (t245 - t250) * qJD(6); m(7) * (t195 * t3 + t199 * t4 + (-t195 * t22 + t199 * t23) * qJD(6)) + m(6) * t78; t199 * t17 + t195 * t18 + (-t195 * t53 + t199 * t52) * qJD(6) + m(7) * (t1 * t195 + t199 * t2 + (-t12 * t195 + t13 * t199) * qJD(6)) + m(6) * t92 + t39; m(7) * (t195 * t20 + t199 * t21 + (-t195 * t67 + t199 * t68) * qJD(6)) + t111 * t235 + t195 * t82 - t112 * t236 + t199 * t81 + m(6) * t233 + t105; 0; 0; mrSges(7,1) * t4 - mrSges(7,2) * t3; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t21 - mrSges(7,2) * t20 + t49; t188 + (t169 * t185 - t266) * qJD(6); -t161; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t24(1) t24(2) t24(4) t24(7) t24(11) t24(16); t24(2) t24(3) t24(5) t24(8) t24(12) t24(17); t24(4) t24(5) t24(6) t24(9) t24(13) t24(18); t24(7) t24(8) t24(9) t24(10) t24(14) t24(19); t24(11) t24(12) t24(13) t24(14) t24(15) t24(20); t24(16) t24(17) t24(18) t24(19) t24(20) t24(21);];
Mq  = res;
