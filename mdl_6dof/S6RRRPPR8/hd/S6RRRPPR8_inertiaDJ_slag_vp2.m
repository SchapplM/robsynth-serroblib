% Calculate time derivative of joint inertia matrix for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:38
% EndTime: 2019-03-09 16:03:50
% DurationCPUTime: 5.20s
% Computational Cost: add. (4276->578), mult. (11114->792), div. (0->0), fcn. (9727->8), ass. (0->231)
t269 = -Ifges(6,5) - Ifges(4,6);
t293 = Ifges(5,6) + t269;
t292 = Ifges(5,4) + Ifges(4,5);
t291 = Ifges(5,2) + Ifges(4,3);
t194 = sin(qJ(6));
t197 = cos(qJ(6));
t198 = cos(qJ(3));
t233 = qJD(6) * t198;
t195 = sin(qJ(3));
t237 = qJD(3) * t195;
t207 = t194 * t233 + t197 * t237;
t290 = -t194 * t237 + t197 * t233;
t274 = -t197 / 0.2e1;
t268 = pkin(9) - qJ(5);
t196 = sin(qJ(2));
t191 = sin(pkin(6));
t199 = cos(qJ(2));
t244 = t191 * t199;
t192 = cos(pkin(6));
t272 = pkin(1) * t192;
t241 = pkin(8) * t244 + t196 * t272;
t117 = t241 * qJD(2);
t245 = t191 * t196;
t124 = t192 * t195 + t198 * t245;
t239 = qJD(2) * t191;
t225 = t199 * t239;
t85 = qJD(3) * t124 + t195 * t225;
t227 = t195 * t245;
t86 = -qJD(3) * t227 + (qJD(3) * t192 + t225) * t198;
t22 = pkin(3) * t85 - qJ(4) * t86 - qJD(4) * t124 + t117;
t289 = -pkin(3) * t198 - qJ(4) * t195;
t288 = 2 * m(6);
t287 = 2 * m(7);
t286 = -2 * mrSges(3,3);
t285 = 2 * mrSges(3,3);
t284 = -2 * mrSges(6,3);
t149 = t268 * t198;
t283 = 0.2e1 * t149;
t193 = qJ(4) + pkin(5);
t282 = 0.2e1 * t193;
t238 = qJD(2) * t196;
t226 = t191 * t238;
t123 = -t192 * t198 + t227;
t88 = t123 * t197 + t194 * t244;
t33 = -qJD(6) * t88 - t194 * t85 - t197 * t226;
t281 = t33 / 0.2e1;
t87 = -t123 * t194 + t197 * t244;
t34 = qJD(6) * t87 - t194 * t226 + t197 * t85;
t280 = t34 / 0.2e1;
t279 = t87 / 0.2e1;
t278 = t88 / 0.2e1;
t200 = -pkin(3) - pkin(4);
t277 = -pkin(4) - pkin(10);
t256 = Ifges(7,6) * t194;
t180 = qJD(6) * t256;
t257 = Ifges(7,5) * t197;
t276 = -qJD(6) * t257 / 0.2e1 + t180 / 0.2e1;
t275 = t194 / 0.2e1;
t273 = t197 / 0.2e1;
t271 = pkin(9) * t195;
t270 = mrSges(6,3) - mrSges(5,2);
t46 = -mrSges(7,1) * t87 + mrSges(7,2) * t88;
t90 = mrSges(6,1) * t244 - mrSges(6,3) * t123;
t267 = t90 - t46;
t266 = mrSges(7,3) * t198;
t265 = Ifges(4,4) * t195;
t264 = Ifges(4,4) * t198;
t263 = Ifges(6,4) * t195;
t262 = Ifges(6,4) * t198;
t261 = Ifges(7,4) * t194;
t260 = Ifges(7,4) * t197;
t259 = Ifges(5,5) * t195;
t258 = Ifges(5,5) * t198;
t125 = -pkin(8) * t245 + t199 * t272;
t116 = t125 * qJD(2);
t255 = t116 * mrSges(3,2);
t254 = t117 * mrSges(3,1);
t253 = t117 * mrSges(4,1);
t252 = t117 * mrSges(4,2);
t65 = mrSges(6,2) * t226 - mrSges(6,3) * t86;
t251 = qJ(5) * t123;
t214 = -Ifges(7,2) * t194 + t260;
t112 = Ifges(7,6) * t195 - t198 * t214;
t250 = t112 * t194;
t216 = Ifges(7,1) * t197 - t261;
t113 = Ifges(7,5) * t195 - t198 * t216;
t249 = t113 * t197;
t120 = qJD(5) * t198 + t237 * t268;
t248 = t120 * t149;
t190 = -pkin(3) + t277;
t247 = t190 * t194;
t246 = t190 * t197;
t213 = Ifges(7,2) * t197 + t261;
t243 = t194 * t213;
t215 = Ifges(7,1) * t194 + t260;
t242 = t197 * t215;
t109 = pkin(9) * t192 + t241;
t110 = (-pkin(2) * t199 - pkin(9) * t196 - pkin(1)) * t191;
t53 = t109 * t198 + t110 * t195;
t236 = qJD(3) * t198;
t240 = qJ(4) * t236 + qJD(4) * t195;
t130 = mrSges(6,1) * t236 + mrSges(6,2) * t237;
t235 = qJD(4) * t149;
t234 = qJD(4) * t199;
t232 = -Ifges(7,5) * t194 / 0.2e1 + Ifges(7,6) * t274 + Ifges(6,6);
t231 = 0.2e1 * t191;
t6 = Ifges(7,5) * t34 + Ifges(7,6) * t33 + Ifges(7,3) * t86;
t144 = -pkin(2) + t289;
t229 = Ifges(6,6) * t244;
t228 = Ifges(6,5) / 0.2e1 + Ifges(4,6) / 0.2e1;
t108 = -pkin(2) * t192 - t125;
t223 = t191 * t234;
t43 = t86 * mrSges(6,1) + mrSges(6,2) * t85;
t52 = -t109 * t195 + t110 * t198;
t128 = pkin(4) * t198 - t144;
t219 = m(5) * pkin(9) - t270;
t47 = pkin(3) * t123 - qJ(4) * t124 + t108;
t18 = pkin(5) * t124 + t123 * t277 - t47;
t49 = pkin(3) * t244 - t52;
t30 = pkin(4) * t244 - qJ(5) * t124 + t49;
t25 = pkin(10) * t244 + t30;
t3 = t18 * t197 - t194 * t25;
t4 = t18 * t194 + t197 * t25;
t218 = t194 * t4 + t197 * t3;
t148 = -t198 * mrSges(5,1) - t195 * mrSges(5,3);
t146 = mrSges(7,1) * t197 - mrSges(7,2) * t194;
t217 = -mrSges(7,1) * t194 - mrSges(7,2) * t197;
t50 = -mrSges(7,2) * t124 + mrSges(7,3) * t87;
t51 = mrSges(7,1) * t124 - mrSges(7,3) * t88;
t212 = -t194 * t50 - t197 * t51;
t100 = pkin(5) * t195 + pkin(10) * t198 + t128;
t147 = t268 * t195;
t69 = t100 * t197 - t147 * t194;
t70 = t100 * t194 + t147 * t197;
t211 = t194 * t70 + t197 * t69;
t142 = -mrSges(7,2) * t195 + t194 * t266;
t143 = mrSges(7,1) * t195 + t197 * t266;
t210 = -t194 * t142 - t197 * t143;
t115 = (pkin(2) * t196 - pkin(9) * t199) * t239;
t24 = -t109 * t236 - t110 * t237 + t115 * t198 - t116 * t195;
t64 = -mrSges(5,1) * t226 + mrSges(5,2) * t86;
t48 = -qJ(4) * t244 + t53;
t27 = Ifges(7,4) * t88 + Ifges(7,2) * t87 + Ifges(7,6) * t124;
t28 = Ifges(7,1) * t88 + Ifges(7,4) * t87 + Ifges(7,5) * t124;
t209 = t27 * t275 + t274 * t28;
t23 = -t109 * t237 + t110 * t236 + t115 * t195 + t116 * t198;
t205 = qJ(4) * t226 + t23;
t204 = t226 * t291 + t292 * t86 + t293 * t85;
t66 = Ifges(7,5) * t207 + Ifges(7,6) * t290 + Ifges(7,3) * t236;
t203 = -qJ(5) * t86 - qJD(5) * t124 - t24;
t201 = qJ(5) * t85 + qJD(5) * t123 + t205;
t183 = Ifges(5,4) * t236;
t182 = Ifges(4,5) * t236;
t181 = Ifges(5,6) * t237;
t164 = Ifges(3,5) * t225;
t159 = Ifges(4,1) * t195 + t264;
t158 = Ifges(5,1) * t195 - t258;
t157 = -Ifges(6,1) * t198 - t263;
t155 = Ifges(4,2) * t198 + t265;
t154 = -Ifges(6,2) * t195 - t262;
t152 = -Ifges(5,3) * t198 + t259;
t150 = mrSges(6,1) * t195 - mrSges(6,2) * t198;
t141 = (Ifges(4,1) * t198 - t265) * qJD(3);
t140 = (Ifges(5,1) * t198 + t259) * qJD(3);
t139 = (Ifges(6,1) * t195 - t262) * qJD(3);
t138 = t216 * qJD(6);
t137 = (-Ifges(4,2) * t195 + t264) * qJD(3);
t136 = (-Ifges(6,2) * t198 + t263) * qJD(3);
t135 = t214 * qJD(6);
t134 = (Ifges(5,3) * t195 + t258) * qJD(3);
t132 = (mrSges(4,1) * t195 + mrSges(4,2) * t198) * qJD(3);
t131 = (mrSges(5,1) * t195 - mrSges(5,3) * t198) * qJD(3);
t129 = t217 * qJD(6);
t127 = t217 * t198;
t122 = -qJD(5) * t195 + t236 * t268;
t121 = pkin(3) * t237 - t240;
t111 = Ifges(7,3) * t195 + (t256 - t257) * t198;
t104 = t200 * t237 + t240;
t99 = -mrSges(7,2) * t236 + mrSges(7,3) * t290;
t98 = mrSges(7,1) * t236 - mrSges(7,3) * t207;
t94 = -mrSges(6,2) * t244 - mrSges(6,3) * t124;
t93 = mrSges(5,1) * t244 + mrSges(5,2) * t124;
t92 = -mrSges(4,1) * t244 - mrSges(4,3) * t124;
t91 = mrSges(4,2) * t244 - mrSges(4,3) * t123;
t89 = -mrSges(5,2) * t123 - mrSges(5,3) * t244;
t74 = (pkin(5) * t198 + t190 * t195) * qJD(3) + t240;
t73 = -mrSges(7,1) * t290 + mrSges(7,2) * t207;
t72 = mrSges(5,1) * t123 - mrSges(5,3) * t124;
t71 = mrSges(6,1) * t124 + mrSges(6,2) * t123;
t68 = t215 * t233 + (Ifges(7,5) * t198 + t195 * t216) * qJD(3);
t67 = t213 * t233 + (Ifges(7,6) * t198 + t195 * t214) * qJD(3);
t63 = mrSges(4,1) * t226 - mrSges(4,3) * t86;
t62 = -mrSges(4,2) * t226 - mrSges(4,3) * t85;
t61 = -mrSges(6,1) * t226 - mrSges(6,3) * t85;
t60 = -mrSges(5,2) * t85 + mrSges(5,3) * t226;
t59 = Ifges(4,1) * t124 - Ifges(4,4) * t123 - Ifges(4,5) * t244;
t58 = Ifges(5,1) * t124 - Ifges(5,4) * t244 + Ifges(5,5) * t123;
t57 = Ifges(6,1) * t123 - Ifges(6,4) * t124 + Ifges(6,5) * t244;
t56 = Ifges(4,4) * t124 - Ifges(4,2) * t123 - Ifges(4,6) * t244;
t55 = Ifges(6,4) * t123 - Ifges(6,2) * t124 + t229;
t54 = Ifges(5,5) * t124 - Ifges(5,6) * t244 + Ifges(5,3) * t123;
t45 = mrSges(4,1) * t85 + mrSges(4,2) * t86;
t44 = mrSges(5,1) * t85 - mrSges(5,3) * t86;
t42 = Ifges(4,1) * t86 - Ifges(4,4) * t85 + Ifges(4,5) * t226;
t41 = Ifges(5,1) * t86 + Ifges(5,4) * t226 + Ifges(5,5) * t85;
t40 = Ifges(6,1) * t85 - Ifges(6,4) * t86 - Ifges(6,5) * t226;
t39 = Ifges(4,4) * t86 - Ifges(4,2) * t85 + Ifges(4,6) * t226;
t38 = Ifges(6,4) * t85 - Ifges(6,2) * t86 - Ifges(6,6) * t226;
t37 = Ifges(5,5) * t86 + Ifges(5,6) * t226 + Ifges(5,3) * t85;
t36 = -t48 - t251;
t35 = -pkin(4) * t123 - t47;
t29 = -t193 * t244 + t251 + t53;
t26 = Ifges(7,5) * t88 + Ifges(7,6) * t87 + Ifges(7,3) * t124;
t21 = -qJD(6) * t70 - t122 * t194 + t197 * t74;
t20 = qJD(6) * t69 + t122 * t197 + t194 * t74;
t19 = -pkin(3) * t226 - t24;
t17 = mrSges(7,1) * t86 - mrSges(7,3) * t34;
t16 = -mrSges(7,2) * t86 + mrSges(7,3) * t33;
t15 = t205 - t223;
t14 = -pkin(4) * t85 - t22;
t13 = -mrSges(7,1) * t33 + mrSges(7,2) * t34;
t12 = -t201 + t223;
t11 = t200 * t226 + t203;
t10 = (pkin(5) * t238 - t234) * t191 + t201;
t9 = t190 * t226 + t203;
t8 = Ifges(7,1) * t34 + Ifges(7,4) * t33 + Ifges(7,5) * t86;
t7 = Ifges(7,4) * t34 + Ifges(7,2) * t33 + Ifges(7,6) * t86;
t5 = pkin(5) * t86 + t277 * t85 - t22;
t2 = -qJD(6) * t4 - t194 * t9 + t197 * t5;
t1 = qJD(6) * t3 + t194 * t5 + t197 * t9;
t31 = [0.2e1 * m(4) * (t108 * t117 + t23 * t53 + t24 * t52) + 0.2e1 * m(5) * (t15 * t48 + t19 * t49 + t22 * t47) + (t1 * t4 + t10 * t29 + t2 * t3) * t287 + (t11 * t30 + t12 * t36 + t14 * t35) * t288 + 0.2e1 * m(3) * (t116 * t241 - t117 * t125) + 0.2e1 * t108 * t45 + t87 * t7 + t88 * t8 + 0.2e1 * t15 * t89 + 0.2e1 * t12 * t90 + 0.2e1 * t23 * t91 + 0.2e1 * t24 * t92 + 0.2e1 * t19 * t93 + 0.2e1 * t11 * t94 + 0.2e1 * t14 * t71 + 0.2e1 * t22 * t72 + 0.2e1 * t48 * t60 + 0.2e1 * t36 * t61 + 0.2e1 * t53 * t62 + 0.2e1 * t52 * t63 + 0.2e1 * t49 * t64 + 0.2e1 * t30 * t65 + 0.2e1 * t1 * t50 + 0.2e1 * t2 * t51 + 0.2e1 * t35 * t43 + 0.2e1 * t10 * t46 + 0.2e1 * t47 * t44 + t33 * t27 + t34 * t28 + 0.2e1 * t29 * t13 + 0.2e1 * t3 * t17 + 0.2e1 * t4 * t16 + (t58 + t59 + t26 - t55) * t86 + (-t56 + t54 + t57) * t85 + (t117 * t196 * t285 + (-Ifges(6,6) * t86 + t116 * t285 - t204) * t199 + ((t125 * t286 + Ifges(3,5) * t192 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t199) * t231) * t199 + (t241 * t286 - 0.2e1 * Ifges(3,6) * t192 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t196) * t231 + (Ifges(6,6) + t292) * t124 + t293 * t123 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(6,3)) - t291) * t244) * t196) * qJD(2)) * t191 + (-t38 + t41 + t42 + t6 + 0.2e1 * t252) * t124 + (t37 - t39 + t40 + 0.2e1 * t253) * t123 + (t164 - 0.2e1 * t254 - 0.2e1 * t255) * t192; (-t61 + t13) * t149 + t164 - t254 - t255 + t68 * t278 + t67 * t279 + t113 * t280 + t112 * t281 + (-t136 / 0.2e1 + t140 / 0.2e1 + t141 / 0.2e1 + t66 / 0.2e1) * t124 + (-t137 / 0.2e1 + t139 / 0.2e1 + t134 / 0.2e1) * t123 + (t152 / 0.2e1 - t155 / 0.2e1 + t157 / 0.2e1) * t85 + m(7) * (t1 * t70 + t10 * t149 - t120 * t29 + t2 * t69 + t20 * t4 + t21 * t3) + m(6) * (t104 * t35 + t11 * t147 - t12 * t149 + t120 * t36 + t122 * t30 + t128 * t14) + (-t154 / 0.2e1 + t158 / 0.2e1 + t159 / 0.2e1 + t111 / 0.2e1) * t86 + t14 * t150 + t1 * t142 + t2 * t143 + t144 * t44 + t147 * t65 + t22 * t148 + t10 * t127 + t128 * t43 + t35 * t130 + t47 * t131 + t108 * t132 + t121 * t72 + t122 * t94 + t3 * t98 + t4 * t99 + t104 * t71 + t70 * t16 + t29 * t73 + t69 * t17 + t20 * t50 + t21 * t51 - pkin(2) * t45 + ((-t182 / 0.2e1 - t183 / 0.2e1 - t181 / 0.2e1) * t199 + (-Ifges(3,6) + (-Ifges(5,6) / 0.2e1 + t228) * t198 + (Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t195) * t238) * t191 + ((-t52 * mrSges(4,3) - t30 * mrSges(6,3) + t49 * mrSges(5,2) - t55 / 0.2e1 + t58 / 0.2e1 + t59 / 0.2e1 + t26 / 0.2e1 - t229 / 0.2e1) * t198 + (-t53 * mrSges(4,3) - t36 * mrSges(6,3) - t48 * mrSges(5,2) + t54 / 0.2e1 - t56 / 0.2e1 + t57 / 0.2e1 + t228 * t244 - t209) * t195 + ((-t92 + t93) * t198 + (-t89 - t91) * t195 + m(5) * (-t195 * t48 + t198 * t49) + m(4) * (-t195 * t53 - t198 * t52)) * pkin(9)) * qJD(3) + (-t38 / 0.2e1 + t41 / 0.2e1 + t42 / 0.2e1 + t6 / 0.2e1 - t24 * mrSges(4,3) + t19 * mrSges(5,2) - t11 * mrSges(6,3) + t252 + (-t63 + t64) * pkin(9)) * t195 + t267 * t120 + m(4) * (-pkin(2) * t117 - t24 * t271) + m(5) * (t121 * t47 + t144 * t22 + t19 * t271) + (-t253 - t37 / 0.2e1 + t39 / 0.2e1 - t40 / 0.2e1 + t8 * t274 + t7 * t275 + t15 * mrSges(5,2) + t23 * mrSges(4,3) + t12 * mrSges(6,3) + (t27 * t273 + t275 * t28) * qJD(6) + (m(4) * t23 + m(5) * t15 + t60 + t62) * pkin(9)) * t198; -0.2e1 * pkin(2) * t132 + 0.2e1 * t104 * t150 - 0.2e1 * t120 * t127 + 0.2e1 * t128 * t130 + 0.2e1 * t144 * t131 + 0.2e1 * t20 * t142 + 0.2e1 * t21 * t143 + t73 * t283 + 0.2e1 * t69 * t98 + 0.2e1 * t70 * t99 + 0.2e1 * (m(5) * t144 + t148) * t121 + (t20 * t70 + t21 * t69 - t248) * t287 + (t104 * t128 + t122 * t147 - t248) * t288 + (t122 * t284 - t136 + t140 + t141 + t66 + (mrSges(6,3) * t283 + t152 - t155 + t157 + t249 - t250) * qJD(3)) * t195 + (0.2e1 * mrSges(6,3) * t120 + t194 * t67 - t197 * t68 - t134 + t137 - t139 + (t112 * t197 + t113 * t194) * qJD(6) + (t147 * t284 + t111 - t154 + t158 + t159) * qJD(3)) * t198; Ifges(6,3) * t226 - t213 * t281 - t215 * t280 - t138 * t278 - t135 * t279 + t124 * t276 + (-t7 / 0.2e1 + t190 * t16 - t1 * mrSges(7,3)) * t197 + (-t8 / 0.2e1 - t190 * t17 + t2 * mrSges(7,3)) * t194 + (t60 - t61) * qJ(4) + m(5) * (-pkin(3) * t19 + qJ(4) * t15 + qJD(4) * t48) + m(6) * (-qJ(4) * t12 - qJD(4) * t36 + t11 * t200) + t200 * t65 + t193 * t13 + t10 * t146 + t29 * t129 - pkin(3) * t64 - t23 * mrSges(4,2) + t24 * mrSges(4,1) - t19 * mrSges(5,1) + t15 * mrSges(5,3) + t11 * mrSges(6,2) - t12 * mrSges(6,1) + (t218 * mrSges(7,3) + (-m(7) * t218 + t212) * t190 + t209) * qJD(6) + t204 + t232 * t86 + m(7) * (qJD(4) * t29 + t1 * t246 + t10 * t193 - t2 * t247) + (t89 - t267) * qJD(4); t181 + t182 + t183 + t195 * t276 + t193 * t73 + t149 * t129 + qJD(4) * t127 + t122 * mrSges(6,2) + (-t146 - mrSges(6,1)) * t120 + (-t67 / 0.2e1 + t190 * t99 - t20 * mrSges(7,3)) * t197 + (-t68 / 0.2e1 - t190 * t98 + t21 * mrSges(7,3)) * t194 + m(6) * (-qJ(4) * t120 + t122 * t200 + t235) + m(7) * (-t120 * t193 + t20 * t246 - t21 * t247 + t235) + (-t249 / 0.2e1 + t250 / 0.2e1 + t211 * mrSges(7,3) + (-m(7) * t211 + t210) * t190) * qJD(6) + (-t138 * t274 - t135 * t275 + (-t213 * t273 - t215 * t275) * qJD(6) + t219 * qJD(4)) * t198 + ((-pkin(3) * mrSges(5,2) - t200 * mrSges(6,3) + t232) * t198 + (-t242 / 0.2e1 + t243 / 0.2e1 + t270 * qJ(4) + t269) * t195 + (m(5) * t289 - t198 * mrSges(4,1) + t195 * mrSges(4,2) + t148) * pkin(9)) * qJD(3); t129 * t282 + t135 * t197 + t138 * t194 + (t242 - t243) * qJD(6) + (m(7) * t282 + 0.2e1 * mrSges(6,1) + 0.2e1 * mrSges(5,3) + 0.2e1 * t146 + 0.2e1 * (m(5) + m(6)) * qJ(4)) * qJD(4); t197 * t16 - t194 * t17 + t212 * qJD(6) + m(7) * (-qJD(6) * t218 + t1 * t197 - t194 * t2) + m(6) * t11 + m(5) * t19 + t64 + t65; -t194 * t98 + t197 * t99 + t210 * qJD(6) + m(7) * (-qJD(6) * t211 - t194 * t21 + t197 * t20) + m(6) * t122 + t219 * t236; 0; 0; t194 * t16 + t197 * t17 + (-t194 * t51 + t197 * t50) * qJD(6) + m(7) * (t1 * t194 + t197 * t2 + (-t194 * t3 + t197 * t4) * qJD(6)) + m(6) * t14 + t43; t194 * t99 + t197 * t98 + (t142 * t197 - t143 * t194) * qJD(6) + m(7) * (t194 * t20 + t197 * t21 + (-t194 * t69 + t197 * t70) * qJD(6)) + m(6) * t104 + t130; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t6; mrSges(7,1) * t21 - mrSges(7,2) * t20 + t66; t180 + (-t146 * t190 - t257) * qJD(6); -t146 * qJD(6); t129; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t31(1) t31(2) t31(4) t31(7) t31(11) t31(16); t31(2) t31(3) t31(5) t31(8) t31(12) t31(17); t31(4) t31(5) t31(6) t31(9) t31(13) t31(18); t31(7) t31(8) t31(9) t31(10) t31(14) t31(19); t31(11) t31(12) t31(13) t31(14) t31(15) t31(20); t31(16) t31(17) t31(18) t31(19) t31(20) t31(21);];
Mq  = res;
