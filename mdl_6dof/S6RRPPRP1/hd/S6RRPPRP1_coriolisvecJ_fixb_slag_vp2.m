% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2018-11-23 16:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:45:10
% EndTime: 2018-11-23 16:45:24
% DurationCPUTime: 14.15s
% Computational Cost: add. (8772->565), mult. (23075->740), div. (0->0), fcn. (17049->8), ass. (0->255)
t361 = -Ifges(6,4) + Ifges(7,5);
t376 = t361 + Ifges(7,5);
t223 = sin(pkin(9));
t226 = sin(qJ(2));
t228 = cos(qJ(2));
t281 = cos(pkin(9));
t204 = t223 * t228 + t226 * t281;
t190 = t204 * qJD(1);
t222 = sin(pkin(10));
t224 = cos(pkin(10));
t161 = qJD(2) * t222 + t190 * t224;
t162 = qJD(2) * t224 - t190 * t222;
t225 = sin(qJ(5));
t227 = cos(qJ(5));
t105 = t225 * t161 - t162 * t227;
t244 = t281 * t228;
t231 = -t223 * t226 + t244;
t191 = t231 * qJD(2);
t180 = qJD(1) * t191;
t203 = t222 * t225 - t227 * t224;
t67 = -qJD(5) * t105 - t180 * t203;
t329 = t67 / 0.2e1;
t205 = t222 * t227 + t224 * t225;
t371 = t227 * t161 + t162 * t225;
t68 = qJD(5) * t371 + t180 * t205;
t327 = t68 / 0.2e1;
t189 = t204 * qJD(2);
t179 = qJD(1) * t189;
t315 = t179 / 0.2e1;
t362 = Ifges(6,1) + Ifges(7,1);
t360 = Ifges(7,4) + Ifges(6,5);
t359 = Ifges(7,2) + Ifges(6,3);
t358 = -Ifges(6,6) + Ifges(7,6);
t261 = qJD(1) * t226;
t188 = -qJD(1) * t244 + t223 * t261;
t128 = t205 * t188;
t193 = t205 * qJD(5);
t350 = t128 + t193;
t129 = t203 * t188;
t192 = t203 * qJD(5);
t349 = t129 + t192;
t300 = -qJ(3) - pkin(7);
t214 = t300 * t228;
t208 = qJD(1) * t214;
t194 = t223 * t208;
t213 = t300 * t226;
t207 = qJD(1) * t213;
t201 = qJD(2) * pkin(2) + t207;
t151 = t201 * t281 + t194;
t138 = -qJD(2) * pkin(3) + qJD(4) - t151;
t99 = -pkin(4) * t162 + t138;
t29 = t105 * pkin(5) - qJ(6) * t371 + t99;
t104 = Ifges(6,4) * t105;
t185 = qJD(5) + t188;
t286 = Ifges(7,5) * t105;
t356 = t185 * t360 + t362 * t371 - t104 + t286;
t375 = -mrSges(6,2) * t99 + mrSges(7,3) * t29 - t356 / 0.2e1;
t312 = -t188 / 0.2e1;
t220 = -pkin(2) * t228 - pkin(1);
t262 = qJD(1) * t220;
t210 = qJD(3) + t262;
t285 = Ifges(5,2) * t222;
t289 = Ifges(5,4) * t224;
t239 = -t285 + t289;
t290 = Ifges(5,4) * t222;
t240 = Ifges(5,1) * t224 - t290;
t241 = mrSges(5,1) * t222 + mrSges(5,2) * t224;
t279 = Ifges(4,5) * qJD(2);
t305 = t224 / 0.2e1;
t306 = -t222 / 0.2e1;
t372 = (Ifges(5,1) * t161 + Ifges(5,4) * t162 + Ifges(5,5) * t188) * t305 + (Ifges(5,4) * t161 + Ifges(5,2) * t162 + Ifges(5,6) * t188) * t306 + t138 * t241 + t210 * mrSges(4,2) + t161 * t240 / 0.2e1 + t162 * t239 / 0.2e1 + t279 / 0.2e1;
t103 = Ifges(7,5) * t371;
t38 = t185 * Ifges(7,6) + t105 * Ifges(7,3) + t103;
t288 = Ifges(6,4) * t371;
t41 = -t105 * Ifges(6,2) + t185 * Ifges(6,6) + t288;
t370 = t99 * mrSges(6,1) + t29 * mrSges(7,1) + t38 / 0.2e1 - t41 / 0.2e1;
t122 = t188 * pkin(3) - t190 * qJ(4) + t210;
t245 = t281 * t208;
t152 = t223 * t201 - t245;
t145 = qJD(2) * qJ(4) + t152;
t76 = t224 * t122 - t145 * t222;
t44 = pkin(4) * t188 - pkin(8) * t161 + t76;
t77 = t222 * t122 + t224 * t145;
t59 = pkin(8) * t162 + t77;
t13 = -t225 * t59 + t227 * t44;
t351 = qJD(6) - t13;
t10 = -pkin(5) * t185 + t351;
t14 = t225 * t44 + t227 * t59;
t11 = qJ(6) * t185 + t14;
t230 = Ifges(5,6) * t162;
t287 = Ifges(5,5) * t161;
t291 = Ifges(4,4) * t190;
t369 = t10 * mrSges(7,1) + t14 * mrSges(6,2) + t77 * mrSges(5,2) + Ifges(4,2) * t312 + Ifges(4,6) * qJD(2) + t291 / 0.2e1 - t11 * mrSges(7,3) - t13 * mrSges(6,1) - t210 * mrSges(4,1) - t76 * mrSges(5,1) - t230 / 0.2e1 - t287 / 0.2e1;
t313 = t185 / 0.2e1;
t318 = t371 / 0.2e1;
t320 = t105 / 0.2e1;
t321 = -t105 / 0.2e1;
t368 = -Ifges(6,2) * t321 + Ifges(7,3) * t320 + t313 * t358 + t318 * t361 + t370;
t367 = -Ifges(6,4) * t321 - Ifges(7,5) * t320 - t360 * t313 - t362 * t318 + t375;
t269 = t180 * t224;
t246 = qJD(2) * t300;
t186 = qJD(3) * t228 + t226 * t246;
t171 = t186 * qJD(1);
t187 = -t226 * qJD(3) + t228 * t246;
t229 = qJD(1) * t187;
t118 = t281 * t171 + t223 * t229;
t114 = qJD(2) * qJD(4) + t118;
t256 = qJD(1) * qJD(2);
t249 = t226 * t256;
t242 = pkin(2) * t249;
t98 = pkin(3) * t179 - qJ(4) * t180 - qJD(4) * t190 + t242;
t50 = -t114 * t222 + t224 * t98;
t31 = pkin(4) * t179 - pkin(8) * t269 + t50;
t270 = t180 * t222;
t51 = t224 * t114 + t222 * t98;
t33 = -pkin(8) * t270 + t51;
t4 = -qJD(5) * t14 - t225 * t33 + t227 * t31;
t2 = -pkin(5) * t179 - t4;
t328 = -t68 / 0.2e1;
t117 = t171 * t223 - t281 * t229;
t88 = pkin(4) * t270 + t117;
t9 = pkin(5) * t68 - qJ(6) * t67 - qJD(6) * t371 + t88;
t366 = mrSges(6,2) * t88 + mrSges(7,2) * t2 - mrSges(6,3) * t4 - mrSges(7,3) * t9 + Ifges(6,4) * t328 + 0.2e1 * t360 * t315 + t376 * t327 + 0.2e1 * t329 * t362;
t365 = -Ifges(6,6) / 0.2e1;
t259 = qJD(2) * t226;
t363 = pkin(2) * t259;
t45 = mrSges(6,1) * t179 - mrSges(6,3) * t67;
t46 = -t179 * mrSges(7,1) + t67 * mrSges(7,2);
t355 = t46 - t45;
t47 = -mrSges(6,2) * t179 - mrSges(6,3) * t68;
t48 = -mrSges(7,2) * t68 + mrSges(7,3) * t179;
t354 = t47 + t48;
t84 = -mrSges(7,2) * t105 + mrSges(7,3) * t185;
t294 = mrSges(6,3) * t105;
t85 = -mrSges(6,2) * t185 - t294;
t298 = t84 + t85;
t293 = mrSges(6,3) * t371;
t86 = mrSges(6,1) * t185 - t293;
t87 = -mrSges(7,1) * t185 + mrSges(7,2) * t371;
t297 = t87 - t86;
t353 = mrSges(4,2) * t190;
t153 = t207 * t223 - t245;
t268 = t188 * t222;
t115 = -pkin(4) * t268 + t153;
t352 = pkin(5) * t350 + qJ(6) * t349 - qJD(6) * t205 - t115;
t295 = mrSges(4,3) * t190;
t263 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t162 + t161 * mrSges(5,2) + t295;
t348 = t179 * t359 + t358 * t68 + t360 * t67;
t123 = -t188 * mrSges(5,2) + mrSges(5,3) * t162;
t124 = mrSges(5,1) * t188 - mrSges(5,3) * t161;
t347 = t123 * t224 - t124 * t222;
t346 = -t222 * t50 + t224 * t51;
t260 = qJD(1) * t228;
t278 = Ifges(3,6) * qJD(2);
t292 = Ifges(3,4) * t226;
t345 = t278 / 0.2e1 + (t228 * Ifges(3,2) + t292) * qJD(1) / 0.2e1 + pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t260);
t264 = t204 * t224;
t150 = -pkin(3) * t231 - t204 * qJ(4) + t220;
t159 = t223 * t213 - t214 * t281;
t92 = t224 * t150 - t159 * t222;
t74 = -pkin(4) * t231 - pkin(8) * t264 + t92;
t265 = t204 * t222;
t93 = t222 * t150 + t224 * t159;
t79 = -pkin(8) * t265 + t93;
t299 = t225 * t74 + t227 * t79;
t110 = pkin(3) * t189 - qJ(4) * t191 - qJD(4) * t204 + t363;
t136 = t186 * t281 + t223 * t187;
t70 = t224 * t110 - t136 * t222;
t36 = -pkin(8) * t191 * t224 + pkin(4) * t189 + t70;
t266 = t191 * t222;
t71 = t222 * t110 + t224 * t136;
t53 = -pkin(8) * t266 + t71;
t8 = -qJD(5) * t299 - t225 * t53 + t227 * t36;
t257 = qJD(5) * t227;
t258 = qJD(5) * t225;
t3 = t225 * t31 + t227 * t33 + t44 * t257 - t258 * t59;
t1 = qJ(6) * t179 + qJD(6) * t185 + t3;
t344 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t337 = mrSges(6,1) * t88 + mrSges(7,1) * t9 + 0.2e1 * Ifges(7,3) * t327 - t67 * Ifges(6,4) / 0.2e1 + t179 * t365 + Ifges(7,6) * t315 + t376 * t329 + (-t328 + t327) * Ifges(6,2);
t335 = -0.2e1 * pkin(1);
t319 = -t371 / 0.2e1;
t314 = -t185 / 0.2e1;
t311 = t188 / 0.2e1;
t310 = -t190 / 0.2e1;
t304 = pkin(2) * t223;
t303 = pkin(7) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t261);
t216 = qJ(4) + t304;
t301 = pkin(8) + t216;
t267 = t188 * t224;
t253 = pkin(2) * t261;
t134 = pkin(3) * t190 + qJ(4) * t188 + t253;
t154 = t207 * t281 + t194;
t82 = t224 * t134 - t154 * t222;
t58 = pkin(4) * t190 + pkin(8) * t267 + t82;
t83 = t222 * t134 + t224 * t154;
t75 = pkin(8) * t268 + t83;
t22 = t225 * t58 + t227 * t75;
t296 = mrSges(4,3) * t188;
t284 = t190 * Ifges(4,1);
t280 = Ifges(3,5) * qJD(2);
t158 = -t281 * t213 - t214 * t223;
t276 = t117 * t158;
t121 = mrSges(5,1) * t270 + mrSges(5,2) * t269;
t250 = t281 * pkin(2);
t24 = t68 * mrSges(6,1) + t67 * mrSges(6,2);
t23 = t68 * mrSges(7,1) - t67 * mrSges(7,3);
t248 = t228 * t256;
t243 = t179 * mrSges(4,1) + t180 * mrSges(4,2);
t135 = t186 * t223 - t281 * t187;
t219 = -t250 - pkin(3);
t100 = pkin(4) * t266 + t135;
t132 = pkin(4) * t265 + t158;
t238 = Ifges(5,5) * t224 - Ifges(5,6) * t222;
t237 = t222 * t51 + t224 * t50;
t236 = t222 * t76 - t224 * t77;
t21 = -t225 * t75 + t227 * t58;
t27 = -t225 * t79 + t227 * t74;
t199 = t301 * t222;
t200 = t301 * t224;
t232 = -t227 * t199 - t200 * t225;
t147 = -t199 * t225 + t200 * t227;
t7 = t225 * t36 + t227 * t53 + t74 * t257 - t258 * t79;
t209 = -t224 * pkin(4) + t219;
t221 = Ifges(3,4) * t260;
t198 = Ifges(3,1) * t261 + t221 + t280;
t184 = Ifges(4,4) * t188;
t169 = -qJD(2) * mrSges(4,2) - t296;
t148 = mrSges(4,1) * t188 + t353;
t142 = -t184 + t279 + t284;
t140 = t203 * t204;
t139 = t205 * t204;
t137 = t203 * pkin(5) - t205 * qJ(6) + t209;
t131 = mrSges(5,1) * t179 - mrSges(5,3) * t269;
t130 = -mrSges(5,2) * t179 - mrSges(5,3) * t270;
t112 = qJD(4) * t205 + qJD(5) * t147;
t111 = -qJD(4) * t203 + qJD(5) * t232;
t97 = t179 * Ifges(5,5) + t180 * t240;
t96 = t179 * Ifges(5,6) + t180 * t239;
t89 = Ifges(5,3) * t188 + t230 + t287;
t81 = t191 * t205 + t257 * t264 - t258 * t265;
t80 = -t191 * t203 - t193 * t204;
t57 = mrSges(6,1) * t105 + mrSges(6,2) * t371;
t56 = mrSges(7,1) * t105 - mrSges(7,3) * t371;
t55 = pkin(5) * t371 + qJ(6) * t105;
t52 = pkin(5) * t139 + qJ(6) * t140 + t132;
t40 = Ifges(7,4) * t371 + t185 * Ifges(7,2) + t105 * Ifges(7,6);
t39 = Ifges(6,5) * t371 - t105 * Ifges(6,6) + t185 * Ifges(6,3);
t26 = pkin(5) * t231 - t27;
t25 = -qJ(6) * t231 + t299;
t16 = -pkin(5) * t190 - t21;
t15 = qJ(6) * t190 + t22;
t12 = pkin(5) * t81 - qJ(6) * t80 + qJD(6) * t140 + t100;
t6 = -pkin(5) * t189 - t8;
t5 = qJ(6) * t189 - qJD(6) * t231 + t7;
t17 = [(-t278 / 0.2e1 + (mrSges(3,1) * t335 - 0.3e1 / 0.2e1 * t292 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t228) * qJD(1) + (t353 + m(4) * (t210 + t262) + t148) * pkin(2) - t345) * t259 + t220 * t243 - t366 * t140 + (t284 / 0.2e1 + t238 * t311 + t142 / 0.2e1 + (-t222 * t77 - t224 * t76) * mrSges(5,3) + t372) * t191 + ((Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t188 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t185 + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t371 + (Ifges(7,6) / 0.2e1 + t365) * t105 + t89 / 0.2e1 + t40 / 0.2e1 + t39 / 0.2e1 - t369) * t189 + (-t151 * t191 - t152 * t189 + t158 * t180 - t159 * t179) * mrSges(4,3) + m(6) * (t100 * t99 + t13 * t8 + t132 * t88 + t14 * t7 + t27 * t4 + t299 * t3) + t299 * t47 + (mrSges(7,2) * t10 - mrSges(6,3) * t13 - t367) * t80 + (-mrSges(7,2) * t11 - mrSges(6,3) * t14 + t368) * t81 + (-t179 * Ifges(4,4) + (Ifges(4,1) + Ifges(5,1) * t224 ^ 2 / 0.2e1 + (-t289 + t285 / 0.2e1) * t222) * t180 - mrSges(5,3) * t237 + t238 * t315 + t305 * t97 + t306 * t96 + (mrSges(4,3) + t241) * t117) * t204 + (-mrSges(7,2) * t1 - mrSges(6,3) * t3 + t358 * t315 + t337) * t139 + (t189 * t310 + t191 * t312) * Ifges(4,4) + m(7) * (t1 * t25 + t10 * t6 + t11 * t5 + t12 * t29 + t2 * t26 + t52 * t9) + (-(Ifges(5,3) + Ifges(4,2)) * t179 - mrSges(4,1) * qJD(1) * t363 - mrSges(5,1) * t50 + mrSges(5,2) * t51 + mrSges(4,3) * t118 - Ifges(6,6) * t328 - Ifges(7,6) * t327 - t360 * t329 - t359 * t315 - (-Ifges(4,4) + t238) * t180 - t344 - t348 / 0.2e1) * t231 + t263 * t135 + (t198 / 0.2e1 - t303 + t280 / 0.2e1 + (mrSges(3,2) * t335 + 0.3e1 / 0.2e1 * Ifges(3,4) * t228) * qJD(1)) * qJD(2) * t228 + m(4) * (t118 * t159 - t135 * t151 + t136 * t152 + t276) + m(5) * (t135 * t138 + t50 * t92 + t51 * t93 + t70 * t76 + t71 * t77 + t276) + t27 * t45 + t26 * t46 + t25 * t48 + t52 * t23 + t12 * t56 + t5 * t84 + t7 * t85 + t8 * t86 + t6 * t87 + t100 * t57 + t71 * t123 + t70 * t124 + t93 * t130 + t92 * t131 + t132 * t24 + t158 * t121 + t136 * t169; (t130 * t224 - t131 * t222) * t216 - t151 * t296 + (Ifges(5,1) * t222 + t289) * t269 / 0.2e1 - (Ifges(5,2) * t224 + t290) * t270 / 0.2e1 - Ifges(3,6) * t249 + (-mrSges(3,1) * t248 + mrSges(3,2) * t249) * pkin(7) + (-t184 + t142) * t311 + t366 * t205 + t352 * t56 - t263 * t153 + (-Ifges(4,2) * t311 + Ifges(6,6) * t320 + Ifges(7,6) * t321 + Ifges(5,3) * t312 + t359 * t314 + t360 * t319 + t369) * t190 + (-t238 * t312 + t372) * t188 + Ifges(3,5) * t248 - (-Ifges(3,2) * t261 + t198 + t221) * t260 / 0.2e1 + (-t226 * (Ifges(3,1) * t228 - t292) / 0.2e1 + pkin(1) * (mrSges(3,1) * t226 + mrSges(3,2) * t228)) * qJD(1) ^ 2 + (-t179 * t304 - t180 * t250) * mrSges(4,3) + (-t1 * t203 - t10 * t349 - t11 * t350) * mrSges(7,2) + (-Ifges(4,1) * t188 - t291 + t39 + t40 + t89) * t310 + (Ifges(6,4) * t320 + Ifges(7,5) * t321 + t314 * t360 + t319 * t362 + t375) * t129 + t367 * t192 + t368 * t193 - (-Ifges(6,2) * t320 + Ifges(7,3) * t321 + t358 * t314 + t361 * t319 - t370) * t128 + t345 * t261 + (-mrSges(5,1) * t224 + mrSges(5,2) * t222 - mrSges(4,1)) * t117 + (t13 * t349 - t14 * t350 - t203 * t3) * mrSges(6,3) + (Ifges(5,5) * t222 + t224 * Ifges(5,6) + t203 * t358) * t315 - t148 * t253 + t297 * t112 + t298 * t111 + t354 * t147 + (-t267 * t76 - t268 * t77 + t346) * mrSges(5,3) + (-t236 * qJD(4) + t117 * t219 - t138 * t153 + t216 * t346 - t76 * t82 - t77 * t83) * m(5) + t347 * qJD(4) - (Ifges(3,5) * t228 - Ifges(3,6) * t226) * t256 / 0.2e1 + ((-t117 * t281 + t118 * t223) * pkin(2) + t151 * t153 - t152 * t154 - t210 * t253) * m(4) + (t1 * t147 + t137 * t9 - t232 * t2 + t352 * t29 + (-t15 + t111) * t11 + (-t16 + t112) * t10) * m(7) + (-t115 * t99 + t232 * t4 + t147 * t3 + t209 * t88 + (-t22 + t111) * t14 + (-t21 - t112) * t13) * m(6) - t355 * t232 + t222 * t97 / 0.2e1 + t219 * t121 + t209 * t24 + t260 * t303 + t96 * t305 + t152 * t295 - t15 * t84 - t22 * t85 - t21 * t86 - t16 * t87 - t115 * t57 - t118 * mrSges(4,2) - t83 * t123 - t82 * t124 + t337 * t203 + t137 * t23 - t154 * t169 - Ifges(4,6) * t179 + Ifges(4,5) * t180; t222 * t130 + t224 * t131 + t354 * t205 + t355 * t203 + (t169 + t347) * t188 + (-t56 - t57 - t263) * t190 + t243 - t349 * t298 + t350 * t297 + (t1 * t205 + t10 * t350 - t11 * t349 - t190 * t29 + t2 * t203) * m(7) + (-t13 * t350 - t14 * t349 - t190 * t99 - t203 * t4 + t205 * t3) * m(6) + (-t138 * t190 - t188 * t236 + t237) * m(5) + (t151 * t190 + t152 * t188 + t242) * m(4); -t297 * t371 + t298 * t105 - t162 * t123 + t161 * t124 + t121 + t23 + t24 + (-t10 * t371 + t105 * t11 + t9) * m(7) + (t105 * t14 + t13 * t371 + t88) * m(6) + (t161 * t76 - t162 * t77 + t117) * m(5); t344 + (-Ifges(6,2) * t371 - t104 + t356) * t320 + (-t105 * t360 + t358 * t371) * t314 + (-t105 * t362 + t103 - t288 + t38) * t319 + (-pkin(5) * t2 + qJ(6) * t1 - t10 * t14 + t11 * t351 - t29 * t55) * m(7) + (t10 * t105 + t11 * t371) * mrSges(7,2) + (Ifges(7,3) * t371 - t286) * t321 + t41 * t318 + (t293 - t297) * t14 + (-t294 - t298) * t13 - pkin(5) * t46 + qJ(6) * t48 - t55 * t56 + qJD(6) * t84 - t29 * (mrSges(7,1) * t371 + mrSges(7,3) * t105) - t99 * (mrSges(6,1) * t371 - mrSges(6,2) * t105) + t348; t371 * t56 - t185 * t84 + 0.2e1 * (t2 / 0.2e1 + t29 * t318 + t11 * t314) * m(7) + t46;];
tauc  = t17(:);
