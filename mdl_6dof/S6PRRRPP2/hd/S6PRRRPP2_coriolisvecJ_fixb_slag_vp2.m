% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:21:02
% EndTime: 2018-11-23 15:21:10
% DurationCPUTime: 8.85s
% Computational Cost: add. (3632->551), mult. (9381->684), div. (0->0), fcn. (5944->8), ass. (0->255)
t345 = Ifges(7,4) + Ifges(6,5);
t326 = Ifges(6,4) + Ifges(5,5);
t348 = Ifges(7,5) - t326;
t342 = Ifges(7,2) + Ifges(6,3);
t347 = Ifges(6,6) - Ifges(7,6);
t334 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t171 = cos(qJ(3));
t248 = qJD(2) * t171;
t162 = Ifges(4,4) * t248;
t346 = -t162 / 0.2e1;
t170 = cos(qJ(4));
t344 = t345 * t170;
t167 = sin(qJ(4));
t343 = t345 * t167;
t341 = Ifges(5,6) - Ifges(6,6);
t168 = sin(qJ(3));
t237 = qJD(2) * qJD(3);
t219 = t168 * t237;
t236 = qJD(3) * qJD(4);
t238 = t170 * qJD(3);
t245 = qJD(4) * t167;
t97 = t170 * t236 + (-t168 * t245 + t171 * t238) * qJD(2);
t244 = qJD(4) * t170;
t246 = qJD(3) * t171;
t98 = t167 * t236 + (t167 * t246 + t168 * t244) * qJD(2);
t340 = t347 * t219 + t342 * t98 + t345 * t97;
t290 = pkin(9) - qJ(6);
t148 = t290 * t170;
t259 = t170 * t171;
t308 = pkin(4) + pkin(5);
t169 = sin(qJ(2));
t165 = sin(pkin(6));
t252 = qJD(1) * t165;
t227 = t169 * t252;
t138 = qJD(2) * pkin(8) + t227;
t166 = cos(pkin(6));
t261 = t166 * t171;
t108 = qJD(1) * t261 - t168 * t138;
t213 = pkin(3) * t168 - pkin(9) * t171;
t136 = t213 * qJD(2);
t44 = -t167 * t108 + t136 * t170;
t339 = -(-qJ(6) * t259 - t168 * t308) * qJD(2) + t44 + qJD(4) * t148 - qJD(6) * t167;
t239 = qJD(6) * t170;
t250 = qJD(2) * t168;
t45 = t170 * t108 + t167 * t136;
t37 = qJ(5) * t250 + t45;
t338 = -qJ(6) * t167 * t248 - t245 * t290 - t239 - t37;
t265 = qJ(5) * t170;
t184 = -t167 * t308 + t265;
t251 = qJD(1) * t168;
t225 = t166 * t251;
t242 = qJD(5) * t167;
t337 = t225 - (qJD(2) * t184 - t138) * t171 + qJD(4) * t184 + t242;
t336 = -t342 * t170 + t343;
t172 = cos(qJ(2));
t226 = t172 * t252;
t276 = qJD(2) * pkin(2);
t139 = -t226 - t276;
t221 = Ifges(4,5) * qJD(3) / 0.2e1;
t156 = -qJD(4) + t248;
t150 = t156 * qJ(5);
t133 = t167 * t250 - t238;
t144 = -pkin(3) * t171 - pkin(9) * t168 - pkin(2);
t110 = qJD(2) * t144 - t226;
t109 = t138 * t171 + t225;
t85 = qJD(3) * pkin(9) + t109;
t31 = t167 * t110 + t170 * t85;
t17 = qJ(6) * t133 + t31;
t12 = -t150 + t17;
t134 = qJD(3) * t167 + t170 * t250;
t214 = qJD(3) * pkin(3) + t108;
t180 = qJ(5) * t134 + t214;
t15 = -t133 * t308 + qJD(6) + t180;
t30 = t170 * t110 - t167 * t85;
t186 = t167 * t31 + t170 * t30;
t322 = qJD(5) - t30;
t18 = pkin(4) * t156 + t322;
t19 = -t150 + t31;
t187 = t167 * t19 - t170 * t18;
t282 = Ifges(5,4) * t170;
t198 = -Ifges(5,2) * t167 + t282;
t205 = -mrSges(7,1) * t167 + mrSges(7,2) * t170;
t207 = mrSges(6,1) * t167 - mrSges(6,3) * t170;
t209 = mrSges(5,1) * t167 + mrSges(5,2) * t170;
t298 = t170 / 0.2e1;
t300 = t167 / 0.2e1;
t301 = -t167 / 0.2e1;
t302 = t156 / 0.2e1;
t303 = -t156 / 0.2e1;
t304 = t134 / 0.2e1;
t306 = t133 / 0.2e1;
t307 = -t133 / 0.2e1;
t131 = Ifges(5,4) * t133;
t328 = t345 * t133;
t316 = t334 * t134 + t348 * t156 - t131 + t328;
t32 = pkin(4) * t133 - t180;
t330 = t345 * t134;
t324 = t342 * t133 - t347 * t156 + t330;
t283 = Ifges(5,4) * t167;
t327 = t170 * t334 - t283 + t343;
t329 = t167 * t342 + t344;
t284 = Ifges(5,4) * t134;
t54 = -Ifges(5,2) * t133 - Ifges(5,6) * t156 + t284;
t16 = qJ(6) * t134 + t30;
t323 = qJD(5) - t16;
t9 = t156 * t308 + t323;
t313 = t187 * mrSges(6,2) + t186 * mrSges(5,3) - (t12 * t167 - t170 * t9) * mrSges(7,3) - (Ifges(7,5) * t170 + Ifges(7,6) * t167) * t302 - t32 * t207 - t15 * t205 + t214 * t209 - t54 * t301 - t198 * t307 - t329 * t306 - (-t167 * t341 + t170 * t326) * t303 - t324 * t300 - t327 * t304 - t316 * t298;
t331 = t250 / 0.2e1;
t335 = -t139 * mrSges(4,2) + t108 * mrSges(4,3) - Ifges(4,1) * t331 - t221 + t313 + t346;
t333 = (-Ifges(5,4) + t345) * t98 + t334 * t97 - t348 * t219;
t332 = t334 * t167 + t282 - t344;
t305 = -t134 / 0.2e1;
t220 = -Ifges(4,6) * qJD(3) / 0.2e1;
t247 = qJD(3) * t168;
t321 = qJ(5) * t247 - qJD(5) * t171;
t320 = t167 * t326 + t341 * t170;
t319 = -qJ(5) * t97 - qJD(5) * t134;
t137 = t213 * qJD(3);
t101 = (t137 + t227) * qJD(2);
t262 = t165 * t172;
t223 = qJD(2) * t262;
t182 = qJD(3) * t166 + t223;
t179 = qJD(1) * t182;
t58 = -t138 * t247 + t171 * t179;
t6 = t167 * t101 + t110 * t244 + t170 * t58 - t245 * t85;
t7 = t101 * t170 - t110 * t245 - t167 * t58 - t85 * t244;
t318 = -t167 * t7 + t170 * t6;
t3 = qJ(5) * t219 - t156 * qJD(5) + t6;
t4 = -pkin(4) * t219 - t7;
t317 = t167 * t4 + t170 * t3;
t266 = qJ(5) * t167;
t314 = -t170 * t308 - t266;
t233 = Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t216 = -Ifges(7,6) / 0.2e1 - t233;
t235 = -Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1;
t217 = -Ifges(7,5) / 0.2e1 - t235;
t234 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t312 = t216 * t133 + t217 * t134 - (Ifges(7,3) / 0.2e1 + t234) * t156 + t12 * mrSges(7,2) + t139 * mrSges(4,1) + t19 * mrSges(6,3) + t30 * mrSges(5,1) + t220 - (Ifges(4,4) * t168 + t171 * Ifges(4,2)) * qJD(2) / 0.2e1 + Ifges(7,5) * t305 + Ifges(6,6) * t306 - t109 * mrSges(4,3) - t18 * mrSges(6,1) - t31 * mrSges(5,2) - t9 * mrSges(7,1) + (Ifges(7,6) + Ifges(5,6)) * t307 + t326 * t304 + (Ifges(7,3) + Ifges(5,3) + Ifges(6,2)) * t303;
t311 = t97 / 0.2e1;
t310 = -t98 / 0.2e1;
t309 = t98 / 0.2e1;
t299 = -t170 / 0.2e1;
t297 = pkin(8) * t167;
t292 = t97 * mrSges(7,3);
t61 = -mrSges(6,2) * t98 + mrSges(6,3) * t219;
t66 = -mrSges(5,2) * t219 - mrSges(5,3) * t98;
t289 = t61 + t66;
t63 = mrSges(5,1) * t219 - mrSges(5,3) * t97;
t87 = t97 * mrSges(6,2);
t64 = -mrSges(6,1) * t219 + t87;
t288 = -t63 + t64;
t70 = mrSges(6,1) * t133 - mrSges(6,3) * t134;
t71 = -mrSges(7,1) * t133 + mrSges(7,2) * t134;
t287 = t70 - t71;
t286 = mrSges(5,3) * t133;
t285 = mrSges(5,3) * t134;
t263 = t165 * t169;
t118 = t168 * t263 - t261;
t222 = t138 * t246;
t59 = t168 * t179 + t222;
t274 = t118 * t59;
t270 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t133 + mrSges(5,2) * t134 + mrSges(4,3) * t250;
t267 = qJ(5) * t133;
t264 = qJ(6) * t168;
t260 = t167 * t171;
t258 = t170 * t172;
t102 = -mrSges(7,2) * t156 + mrSges(7,3) * t133;
t107 = -mrSges(6,2) * t133 - mrSges(6,3) * t156;
t257 = t102 + t107;
t103 = mrSges(5,2) * t156 - t286;
t256 = t103 + t107;
t105 = -mrSges(5,1) * t156 - t285;
t106 = mrSges(6,1) * t156 + mrSges(6,2) * t134;
t255 = t105 - t106;
t254 = t167 * t137 + t144 * t244;
t158 = pkin(8) * t259;
t253 = qJD(4) * t158 + t144 * t245;
t113 = t167 * t144 + t158;
t249 = qJD(2) * t169;
t241 = qJD(5) * t170;
t232 = mrSges(4,3) * t248;
t231 = t167 * t262;
t230 = t102 + t256;
t104 = mrSges(7,1) * t156 - mrSges(7,3) * t134;
t229 = t104 - t255;
t228 = -pkin(4) - t297;
t224 = t165 * t249;
t34 = -t98 * mrSges(7,1) + t97 * mrSges(7,2);
t146 = -qJD(3) * mrSges(4,2) + t232;
t218 = -m(4) * t109 - t146;
t157 = pkin(8) * t260;
t112 = t144 * t170 - t157;
t95 = -qJ(5) * t171 + t113;
t212 = -t137 * t170 + t253;
t210 = mrSges(5,1) * t170 - mrSges(5,2) * t167;
t208 = mrSges(6,1) * t170 + mrSges(6,3) * t167;
t206 = mrSges(7,1) * t170 + mrSges(7,2) * t167;
t197 = Ifges(5,2) * t170 + t283;
t190 = Ifges(7,5) * t167 - Ifges(7,6) * t170;
t189 = pkin(4) * t170 + t266;
t188 = pkin(4) * t167 - t265;
t185 = pkin(8) + t188;
t119 = t166 * t168 + t171 * t263;
t77 = t119 * t167 + t165 * t258;
t183 = -m(4) * t108 - m(5) * t214 + t270;
t181 = -pkin(8) + t184;
t41 = (-t168 * t238 - t171 * t245) * pkin(8) + t254;
t1 = qJ(6) * t98 + qJD(6) * t133 + t3;
t2 = -qJ(6) * t97 - qJD(6) * t134 - t219 * t308 - t7;
t178 = -t7 * mrSges(5,1) + t4 * mrSges(6,1) + t2 * mrSges(7,1) + t6 * mrSges(5,2) - t1 * mrSges(7,2) - t3 * mrSges(6,3);
t174 = qJD(2) ^ 2;
t164 = t171 * pkin(4);
t155 = Ifges(6,2) * t219;
t154 = Ifges(5,3) * t219;
t147 = t290 * t167;
t140 = -pkin(3) - t189;
t135 = (-mrSges(4,1) * t171 + mrSges(4,2) * t168) * qJD(2);
t128 = pkin(3) - t314;
t124 = (mrSges(4,1) * t168 + mrSges(4,2) * t171) * t237;
t116 = qJD(4) * t188 - t242;
t114 = t185 * t168;
t100 = (t167 * t169 + t171 * t258) * t252;
t99 = -t170 * t227 + t226 * t260;
t96 = -t112 + t164;
t94 = t181 * t168;
t91 = Ifges(6,4) * t97;
t90 = Ifges(5,5) * t97;
t89 = Ifges(5,6) * t98;
t88 = Ifges(6,6) * t98;
t78 = t119 * t170 - t231;
t76 = qJD(3) * t119 + t168 * t223;
t75 = -qJD(3) * t118 + t171 * t223;
t69 = pkin(4) * t134 + t267;
t67 = t167 * t264 + t95;
t65 = mrSges(7,2) * t219 + mrSges(7,3) * t98;
t62 = -mrSges(7,1) * t219 - t292;
t60 = pkin(5) * t171 + t157 + t164 + (-t144 - t264) * t170;
t48 = t225 + (qJD(2) * t188 + t138) * t171;
t43 = -t134 * t308 - t267;
t42 = t247 * t297 - t212;
t40 = (qJD(4) * t189 - t241) * t168 + t185 * t246;
t39 = -pkin(4) * t250 - t44;
t36 = t228 * t247 + t212;
t35 = mrSges(5,1) * t98 + mrSges(5,2) * t97;
t33 = mrSges(6,1) * t98 - mrSges(6,3) * t97;
t28 = t41 + t321;
t24 = t97 * Ifges(5,4) - t98 * Ifges(5,2) + Ifges(5,6) * t219;
t20 = (qJD(4) * t314 + t241) * t168 + t181 * t246;
t14 = -qJD(4) * t231 + t119 * t244 + t167 * t75 - t170 * t224;
t13 = -qJD(4) * t77 + t167 * t224 + t170 * t75;
t11 = (-pkin(8) * qJD(3) + qJ(6) * qJD(4)) * t170 * t168 + (qJD(6) * t168 + (-pkin(8) * qJD(4) + qJ(6) * qJD(3)) * t171) * t167 + t254 + t321;
t10 = (-qJ(6) * t246 - t137) * t170 + (qJ(6) * t245 - t239 + (-pkin(5) + t228) * qJD(3)) * t168 + t253;
t8 = pkin(4) * t98 + t319 + t59;
t5 = -t182 * t251 - t308 * t98 - t222 - t319;
t21 = [-t119 * mrSges(4,3) * t219 + t75 * t146 + (t65 + t289) * t78 + (t62 + t288) * t77 + t229 * t14 + t230 * t13 + (t270 + t287) * t76 + ((-mrSges(3,2) * t174 - t124) * t172 + (-mrSges(3,1) * t174 + qJD(2) * t135) * t169) * t165 + (qJD(3) * t232 + t33 - t34 + t35) * t118 + m(5) * (t13 * t31 - t14 * t30 - t214 * t76 + t6 * t78 - t7 * t77 + t274) + m(6) * (t118 * t8 + t13 * t19 + t14 * t18 + t3 * t78 + t32 * t76 + t4 * t77) + m(7) * (t1 * t78 - t118 * t5 + t12 * t13 + t14 * t9 - t15 * t76 + t2 * t77) + m(4) * (-t108 * t76 + t109 * t75 + t274 + t119 * t58 + (t139 - t226) * t224); m(6) * (t114 * t8 + t18 * t36 + t19 * t28 + t3 * t95 + t32 * t40 + t4 * t96) + m(7) * (t1 * t67 + t10 * t9 + t11 * t12 + t15 * t20 + t2 * t60 + t5 * t94) + (0.2e1 * (-t276 / 0.2e1 - t139 / 0.2e1) * m(4) - t135) * t227 + (-t154 / 0.2e1 - t155 / 0.2e1 + t89 / 0.2e1 - t90 / 0.2e1 - t91 / 0.2e1 - t88 / 0.2e1 + (m(4) * pkin(8) + mrSges(4,3)) * t58 + (Ifges(7,6) + t233) * t98 + (Ifges(7,5) + t235) * t97 + (-mrSges(4,1) * t249 + t172 * t218) * t252 + (t183 * pkin(8) + t221 + 0.3e1 / 0.2e1 * t162 - t335) * qJD(3) + t178) * t171 + m(5) * (t112 * t7 + t113 * t6 + t30 * t42 + t31 * t41) - m(5) * (t100 * t31 - t30 * t99) - m(6) * (t100 * t19 + t18 * t99) - m(7) * (t100 * t12 + t9 * t99) - t229 * t99 - t230 * t100 + (t198 * t310 + t24 * t301 + t8 * t207 + t5 * t205 + (t1 * t167 - t170 * t2) * mrSges(7,3) + (-t167 * t6 - t170 * t7) * mrSges(5,3) + (-t167 * t3 + t170 * t4) * mrSges(6,2) + (mrSges(4,3) + t209) * t59 + (mrSges(4,2) * t249 + (-m(6) * t32 + m(7) * t15 - t183 - t287) * t172) * t252 + (((-0.3e1 / 0.2e1 * Ifges(4,4) + t217 * t170 + t216 * t167) * t168 + (-0.3e1 / 0.2e1 * Ifges(4,2) - Ifges(7,3) + 0.3e1 / 0.2e1 * Ifges(4,1) - t234) * t171) * qJD(2) + t220 + t312) * qJD(3) + (t54 * t299 + t32 * t208 - t15 * t206 - t214 * t210 + t197 * t306 + t190 * t303 + (t12 * t170 + t167 * t9) * mrSges(7,3) + (t167 * t30 - t170 * t31) * mrSges(5,3) + (-t167 * t18 - t170 * t19) * mrSges(6,2) + t336 * t307 + t320 * t302 + t332 * t305 + t316 * t301) * qJD(4) + t329 * t309 + t340 * t300 + t327 * t311 + (qJD(4) * t324 + t333) * t298 + (t35 + (m(5) + m(4)) * t59 + t218 * qJD(3)) * pkin(8)) * t168 + t60 * t62 + t67 * t65 + t40 * t70 + t20 * t71 + t94 * t34 + t95 * t61 + t96 * t64 + t11 * t102 + t41 * t103 + t10 * t104 + t42 * t105 + t36 * t106 + t28 * t107 + t112 * t63 + t113 * t66 + t114 * t33 - pkin(2) * t124; ((Ifges(4,4) * t331 + t220 - t312 + (-t190 / 0.2e1 + t320 / 0.2e1) * qJD(3)) * t168 + (t346 + t221 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t250 + t335) * t171) * qJD(2) + ((-m(5) * t186 - m(6) * t187 - t167 * t256 - t170 * t255) * pkin(9) - t313) * qJD(4) + t332 * t311 + t333 * t300 + m(6) * (pkin(9) * t317 + t116 * t32 + t140 * t8) + t317 * mrSges(6,2) + m(5) * (-pkin(3) * t59 + pkin(9) * t318) + t318 * mrSges(5,3) + (-t48 + t116) * t70 + t336 * t309 + (-t1 * t170 - t167 * t2) * mrSges(7,3) + t5 * t206 - t8 * t208 + (-mrSges(4,1) - t210) * t59 - m(6) * (t18 * t39 + t19 * t37 + t32 * t48) - m(5) * (-t109 * t214 + t30 * t44 + t31 * t45) + t197 * t310 + t337 * t71 + t338 * t102 + t339 * t104 + (t1 * t148 + t12 * t338 + t128 * t5 + t147 * t2 + t15 * t337 + t339 * t9) * m(7) + t340 * t299 - pkin(3) * t35 - t58 * mrSges(4,2) - t45 * t103 - t44 * t105 - t39 * t106 - t37 * t107 + t128 * t34 + t140 * t33 - t108 * t146 + t147 * t62 + t148 * t65 - t270 * t109 + (t167 * t288 + t170 * t289) * pkin(9) + t24 * t298; (t61 + t65) * qJ(5) + t154 + t155 + (-t133 * t334 - t284 + t324 + t330) * t305 - t178 + Ifges(7,3) * t219 + (-t12 * t134 - t133 * t9) * mrSges(7,3) + (t133 * t18 + t134 * t19) * mrSges(6,2) + t214 * (mrSges(5,1) * t134 - mrSges(5,2) * t133) + (qJ(5) * t1 + t12 * t323 - t15 * t43 - t17 * t9 - t2 * t308) * m(7) - t308 * t62 + (-pkin(4) * t4 + qJ(5) * t3 - t18 * t31 + t19 * t322 - t32 * t69) * m(6) - t89 + t90 + t91 + t88 + t257 * qJD(5) - pkin(4) * t64 - t69 * t70 - t43 * t71 - Ifges(7,5) * t97 - Ifges(7,6) * t98 - t16 * t102 - t17 * t104 - t32 * (mrSges(6,1) * t134 + mrSges(6,3) * t133) - t15 * (-mrSges(7,1) * t134 - mrSges(7,2) * t133) + (-t133 * t326 - t134 * t341) * t302 + (t134 * t342 - t328) * t307 + (t255 + t285) * t31 + (-t256 - t286) * t30 + (-Ifges(5,2) * t134 - t131 + t316) * t306 + (-Ifges(7,5) * t133 + Ifges(7,6) * t134) * t303 + t54 * t304; -t292 + t87 + t257 * t156 + t287 * t134 + (-mrSges(6,1) - mrSges(7,1)) * t219 + (t12 * t156 - t134 * t15 + t2) * m(7) + (t134 * t32 + t156 * t19 + t4) * m(6); -t133 * t102 + t134 * t104 + 0.2e1 * (t5 / 0.2e1 + t12 * t307 + t9 * t304) * m(7) + t34;];
tauc  = t21(:);
