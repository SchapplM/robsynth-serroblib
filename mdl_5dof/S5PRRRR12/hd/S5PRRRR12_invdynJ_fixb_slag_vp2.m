% Calculate vector of inverse dynamics joint torques for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR12_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:21
% EndTime: 2024-09-28 18:07:28
% DurationCPUTime: 3.91s
% Computational Cost: add. (4911->438), mult. (10287->672), div. (0->0), fcn. (8162->26), ass. (0->243)
t212 = sin(qJ(3));
t213 = sin(qJ(2));
t216 = cos(qJ(3));
t217 = cos(qJ(2));
t232 = t212 * t217 + t213 * t216;
t206 = sin(pkin(5));
t271 = qJD(1) * t206;
t115 = t232 * t271;
t276 = t212 * t213;
t231 = t216 * t217 - t276;
t132 = t231 * t206;
t116 = qJD(1) * t132;
t211 = sin(qJ(4));
t215 = cos(qJ(4));
t52 = -t115 * t215 - t116 * t211;
t324 = pkin(2) * t216;
t189 = pkin(3) + t324;
t267 = qJD(4) * t215;
t268 = qJD(4) * t211;
t275 = t212 * t215;
t82 = -t189 * t268 + (-t212 * t267 + (-t211 * t216 - t275) * qJD(3)) * pkin(2);
t354 = t52 - t82;
t277 = t211 * t212;
t341 = -t115 * t211 + t116 * t215 - t189 * t267 - (-t212 * t268 + (t215 * t216 - t277) * qJD(3)) * pkin(2);
t202 = qJD(2) + qJD(3);
t193 = qJD(4) + t202;
t205 = sin(pkin(6));
t210 = sin(qJ(5));
t214 = cos(qJ(5));
t225 = (-mrSges(6,1) * t214 + mrSges(6,2) * t210) * t205;
t288 = t205 * t193 * t225;
t331 = -t193 * mrSges(5,1) + t288;
t269 = qJD(1) * t217;
t163 = qJD(2) * pkin(2) + t206 * t269;
t248 = t213 * t271;
t105 = t163 * t212 + t216 * t248;
t300 = t105 * t211;
t104 = t216 * t163 - t212 * t248;
t91 = pkin(3) * t202 + t104;
t42 = t215 * t91 - t300;
t353 = -m(5) * t42 + t331;
t204 = sin(pkin(11));
t209 = cos(pkin(5));
t289 = t204 * t209;
t207 = cos(pkin(11));
t199 = t205 * pkin(10);
t261 = t204 * t199;
t147 = pkin(4) * t207 + t209 * t261;
t260 = t207 * t199;
t149 = pkin(4) * t289 - t260;
t335 = t147 * t211 + t149 * t215;
t64 = -pkin(3) * t289 - t335;
t76 = -t147 * t215 + t149 * t211;
t68 = pkin(3) * t207 - t76;
t352 = -t212 * t68 + t64 * t216;
t148 = -t204 * pkin(4) + t209 * t260;
t284 = t207 * t209;
t150 = pkin(4) * t284 + t261;
t242 = t148 * t211 + t150 * t215;
t65 = pkin(3) * t284 + t242;
t77 = t148 * t215 - t150 * t211;
t71 = pkin(3) * t204 - t77;
t351 = -t212 * t71 + t65 * t216;
t350 = -t212 * t65 - t216 * t71;
t349 = -t212 * t64 - t216 * t68;
t347 = t193 * t205 ^ 2;
t208 = cos(pkin(6));
t283 = t208 * t210;
t144 = pkin(2) * t275 + t211 * t189;
t117 = t199 + t144;
t143 = -pkin(2) * t277 + t215 * t189;
t140 = pkin(4) + t143;
t282 = t208 * t214;
t56 = -t117 * t210 + t140 * t282;
t346 = qJD(5) * t56 - t341 * t214 - t354 * t283;
t57 = t117 * t214 + t140 * t283;
t345 = -qJD(5) * t57 + t341 * t210 - t354 * t282;
t323 = pkin(3) * t211;
t167 = t199 + t323;
t322 = pkin(3) * t215;
t187 = pkin(4) + t322;
t107 = t167 * t214 + t187 * t283;
t311 = pkin(3) * qJD(4);
t299 = t105 * t215;
t44 = -t104 * t211 - t299;
t45 = t104 * t215 - t300;
t344 = t210 * t45 - t282 * t44 - t107 * qJD(5) + (-t210 * t215 - t211 * t282) * t311;
t106 = -t167 * t210 + t187 * t282;
t343 = -t214 * t45 - t283 * t44 + t106 * qJD(5) + (-t211 * t283 + t214 * t215) * t311;
t326 = m(4) * pkin(2);
t340 = -mrSges(3,1) - t326;
t286 = t205 * t214;
t153 = pkin(4) * t283 + pkin(10) * t286;
t43 = t211 * t91 + t299;
t339 = -t153 * qJD(5) + t210 * t42 + t282 * t43;
t287 = t205 * t210;
t152 = pkin(4) * t282 - pkin(10) * t287;
t338 = t152 * qJD(5) - t214 * t42 + t283 * t43;
t270 = qJD(1) * t209;
t41 = pkin(4) * t193 + t42;
t229 = t205 * t270 + t208 * t41;
t39 = t193 * t199 + t43;
t14 = -t210 * t39 + t214 * t229;
t337 = qJD(5) * t14;
t313 = Ifges(6,4) * t210;
t336 = t210 * (Ifges(6,1) * t214 - t313) * t347;
t133 = t232 * t206;
t203 = qJ(2) + qJ(3);
t195 = pkin(5) - t203;
t184 = -qJ(4) + t195;
t171 = sin(t184) / 0.2e1;
t194 = pkin(5) + t203;
t183 = qJ(4) + t194;
t172 = cos(t183) / 0.2e1;
t177 = sin(t183);
t178 = cos(t184);
t200 = qJ(4) + t203;
t185 = sin(t200);
t186 = cos(t200);
t314 = mrSges(6,3) * t205;
t334 = -((-t185 * t283 + t186 * t214) * mrSges(6,1) + (-t185 * t282 - t186 * t210) * mrSges(6,2) + t185 * t314) * t206 - (t177 / 0.2e1 + t171) * mrSges(5,1) - (t172 - t178 / 0.2e1) * mrSges(5,2);
t281 = t209 * t210;
t251 = t208 * t281;
t122 = t204 * t214 + t207 * t251;
t279 = t209 * t214;
t250 = t208 * t279;
t123 = -t204 * t210 + t207 * t250;
t124 = t204 * t283 - t207 * t279;
t125 = -t204 * t282 - t207 * t281;
t141 = t171 - t177 / 0.2e1;
t142 = t178 / 0.2e1 + t172;
t293 = t186 * t204;
t333 = -(-t122 * t185 - t124 * t186) * mrSges(6,1) - (-t123 * t185 + t125 * t186) * mrSges(6,2) - (t185 * t284 + t293) * t314 - (t142 * t207 - t185 * t204) * mrSges(5,1) - (t141 * t207 - t293) * mrSges(5,2);
t120 = t204 * t251 - t207 * t214;
t121 = -t204 * t250 - t207 * t210;
t126 = t204 * t279 + t207 * t283;
t127 = t204 * t281 - t207 * t282;
t292 = t186 * t207;
t332 = -(t120 * t185 - t126 * t186) * mrSges(6,1) - (-t121 * t185 + t127 * t186) * mrSges(6,2) - (-t185 * t289 + t292) * t314 - (-t142 * t204 - t185 * t207) * mrSges(5,1) - (-t141 * t204 - t292) * mrSges(5,2);
t174 = sin(t195) / 0.2e1;
t175 = cos(t194) / 0.2e1;
t181 = sin(t194);
t182 = cos(t195);
t330 = -(t181 / 0.2e1 + t174) * mrSges(4,1) - (t175 - t182 / 0.2e1) * mrSges(4,2) + t334;
t155 = t174 - t181 / 0.2e1;
t156 = t182 / 0.2e1 + t175;
t197 = cos(t203);
t196 = sin(t203);
t291 = t196 * t204;
t329 = -(t156 * t207 - t291) * mrSges(4,1) - (t155 * t207 - t197 * t204) * mrSges(4,2) + t333;
t290 = t196 * t207;
t328 = -(-t156 * t204 - t290) * mrSges(4,1) - (-t155 * t204 - t197 * t207) * mrSges(4,2) + t332;
t15 = t210 * t229 + t214 * t39;
t66 = t132 * t215 - t133 * t211;
t234 = t205 * t209 + t208 * t66;
t67 = t132 * t211 + t133 * t215;
t27 = t210 * t234 + t214 * t67;
t325 = m(5) * pkin(3);
t317 = mrSges(4,1) * t202;
t316 = mrSges(4,2) * t202;
t315 = mrSges(5,2) * t193;
t312 = Ifges(6,4) * t214;
t201 = qJDD(2) + qJDD(3);
t191 = qJDD(4) + t201;
t310 = t191 * mrSges(5,1);
t309 = t191 * mrSges(5,2);
t29 = -t205 * t41 + t208 * t270;
t307 = t205 * t29;
t266 = qJD(5) * t210;
t111 = (t191 * t214 - t193 * t266) * t205;
t265 = qJD(5) * t214;
t112 = (t191 * t210 + t193 * t265) * t205;
t47 = -mrSges(6,1) * t111 + mrSges(6,2) * t112;
t306 = t205 * t47;
t305 = t210 * t67;
t162 = -pkin(4) * t211 + t199 * t215;
t294 = t162 * t212;
t285 = t206 * t217;
t280 = t209 * t213;
t278 = t209 * t217;
t264 = qJDD(1) * t209;
t258 = mrSges(6,3) * t286;
t257 = pkin(2) * qJD(3) * t202;
t256 = pkin(3) * t268;
t255 = t193 * t287;
t254 = t193 * t286;
t253 = t206 * t287;
t252 = t206 * t286;
t164 = t191 * t208 + qJDD(5);
t249 = Ifges(6,5) * t112 + Ifges(6,6) * t111 + Ifges(6,3) * t164;
t246 = t287 / 0.2e1;
t245 = qJD(2) * t269;
t241 = mrSges(6,3) * t255;
t240 = mrSges(6,3) * t254;
t237 = qJD(2) * t248;
t161 = -pkin(4) * t215 - t199 * t211;
t159 = pkin(3) - t161;
t233 = -t159 * t216 - t294;
t230 = -pkin(3) * t276 + (pkin(3) * t216 + pkin(2)) * t217;
t138 = qJDD(1) * t285 - t237;
t118 = qJDD(2) * pkin(2) + t138;
t139 = (qJDD(1) * t213 + t245) * t206;
t38 = -t105 * qJD(3) + t216 * t118 - t139 * t212;
t32 = pkin(3) * t201 + t38;
t37 = t104 * qJD(3) + t118 * t212 + t139 * t216;
t13 = -qJD(4) * t43 - t211 * t37 + t215 * t32;
t10 = pkin(4) * t191 + t13;
t226 = t10 * t208 + t205 * t264;
t12 = -t105 * t268 + t211 * t32 + t215 * t37 + t91 * t267;
t224 = (Ifges(6,2) * t214 + t313) * t205;
t223 = (mrSges(6,1) * t210 + mrSges(6,2) * t214) * t307;
t166 = t193 * t208 + qJD(5);
t222 = t166 * t205 * (Ifges(6,5) * t214 - Ifges(6,6) * t210);
t9 = t191 * t199 + t12;
t3 = t210 * t226 + t214 * t9 + t337;
t4 = -qJD(5) * t15 - t210 * t9 + t226 * t214;
t79 = Ifges(6,6) * t166 + t193 * t224;
t8 = -t10 * t205 + t208 * t264;
t160 = Ifges(6,4) * t254;
t80 = Ifges(6,1) * t255 + Ifges(6,5) * t166 + t160;
t220 = -t258 * t337 - t12 * mrSges(5,2) + t3 * (-mrSges(6,2) * t208 + t258) + t13 * mrSges(5,1) + Ifges(5,3) * t191 + t4 * (mrSges(6,1) * t208 - mrSges(6,3) * t287) + qJD(5) * t223 + (Ifges(6,1) * t112 + Ifges(6,4) * t111 + Ifges(6,5) * t164) * t246 + (Ifges(6,4) * t112 + Ifges(6,2) * t111 + Ifges(6,6) * t164) * t286 / 0.2e1 + t208 * t249 / 0.2e1 + t111 * (Ifges(6,6) * t208 + t224) / 0.2e1 + t112 * (Ifges(6,5) * t208 + (Ifges(6,1) * t210 + t312) * t205) / 0.2e1 + t8 * t225 + t164 * (Ifges(6,3) * t208 + (Ifges(6,5) * t210 + Ifges(6,6) * t214) * t205) / 0.2e1 + (t222 + t336) * qJD(5) / 0.2e1 + (-mrSges(6,3) * t15 - t79 / 0.2e1) * t205 * t266 + (t205 * t80 + (-Ifges(6,2) * t210 + t312) * t347) * t265 / 0.2e1;
t219 = t38 * mrSges(4,1) - t37 * mrSges(4,2) + Ifges(4,3) * t201 + t220;
t218 = qJD(2) ^ 2;
t192 = t209 ^ 2 * qJDD(1);
t168 = t205 * t279;
t165 = -pkin(2) * t213 - pkin(3) * t196;
t151 = t162 * t216;
t119 = t231 * t209 * pkin(3);
t114 = t230 * t209;
t109 = -mrSges(6,2) * t166 + t240;
t108 = mrSges(6,1) * t166 - t241;
t78 = (-t159 * t212 + t151) * t206 * t213;
t75 = mrSges(6,1) * t164 - mrSges(6,3) * t112;
t74 = -mrSges(6,2) * t164 + mrSges(6,3) * t111;
t73 = t202 * t133;
t72 = t202 * t132;
t46 = -t205 * t66 + t208 * t209;
t26 = t282 * t66 + t168 - t305;
t19 = -qJD(4) * t67 - t211 * t72 - t215 * t73;
t18 = qJD(4) * t66 - t211 * t73 + t215 * t72;
t6 = -qJD(5) * t27 - t18 * t210 + t19 * t282;
t5 = t19 * t283 + t18 * t214 + (t214 * t234 - t305) * qJD(5);
t1 = [-t19 * t288 + m(2) * qJDD(1) + t6 * t108 + t5 * t109 + t26 * t75 + t27 * t74 + t46 * t47 + (-t18 * t193 - t191 * t67) * mrSges(5,2) + (-t133 * t201 - t202 * t72) * mrSges(4,2) + (t19 * t193 + t191 * t66) * mrSges(5,1) + (t132 * t201 - t202 * t73) * mrSges(4,1) + ((-qJDD(2) * t213 - t217 * t218) * mrSges(3,2) + (qJDD(2) * t217 - t213 * t218) * mrSges(3,1)) * t206 + (-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3) + m(5) * (t12 * t67 + t13 * t66 + t18 * t43 + t19 * t42 + t192) + m(4) * (-t104 * t73 + t105 * t72 + t132 * t38 + t133 * t37 + t192) + m(6) * (t14 * t6 + t15 * t5 - t19 * t307 + t26 * t4 + t27 * t3 + t46 * t8) + m(3) * (t192 + (t138 * t217 + t139 * t213) * t206); t57 * t74 + t56 * t75 + t143 * t310 + t116 * t316 + (t212 * t37 + t216 * t38 + (-t104 * t212 + t105 * t216) * qJD(3)) * t326 - t140 * t306 - t144 * t309 + t115 * t317 - m(4) * (-t104 * t115 + t105 * t116) + t219 + t341 * t315 + (t12 * t144 + t13 * t143 - t341 * t43) * m(5) + t345 * t108 + (t307 * t52 + t3 * t57 + t4 * t56 + (-t140 * t8 - t29 * t82) * t205 + t346 * t15 + t345 * t14) * m(6) + t346 * t109 + (-(-t204 * t217 - t207 * t280) * mrSges(3,2) - m(6) * (-(pkin(2) * t204 - t350) * t213 + (pkin(2) * t284 + t351) * t217) - m(5) * (t114 * t207 + t165 * t204) + t340 * (-t204 * t213 + t207 * t278) + t329) * g(2) + (-(t204 * t280 - t207 * t217) * mrSges(3,2) - m(6) * (-(pkin(2) * t207 - t349) * t213 + (-pkin(2) * t289 + t352) * t217) - m(5) * (-t114 * t204 + t165 * t207) + t340 * (-t204 * t278 - t207 * t213) + t328) * g(1) + (-m(6) * ((pkin(2) - t233) * t285 + t78) - t285 * t326 + (-m(5) * t230 - mrSges(3,1) * t217 + mrSges(3,2) * t213) * t206 + t330) * g(3) + (-t212 * pkin(2) * t201 - t216 * t257) * mrSges(4,2) + (t201 * t324 - t212 * t257) * mrSges(4,1) + Ifges(3,3) * qJDD(2) + (t206 * t245 - t139) * mrSges(3,2) + (t237 + t138) * mrSges(3,1) + t353 * t354; t106 * t75 + t107 * t74 + t104 * t316 + t105 * t317 + t310 * t322 + (t12 * t211 + t13 * t215 + (-t211 * t42 + t215 * t43) * qJD(4)) * t325 - t187 * t306 + t219 - m(5) * t43 * t45 - t309 * t323 + t353 * t44 + (-pkin(3) * t267 + t45) * t315 + t331 * t256 + t343 * t109 + t344 * t108 + (t44 * t307 + t106 * t4 + t107 * t3 + (-t187 * t8 + t256 * t29) * t205 + t343 * t15 + t344 * t14) * m(6) + (-m(6) * (-t233 * t285 + t78) - t132 * t325 + t330) * g(3) + (-m(6) * (t350 * t213 + t351 * t217) - m(5) * (-pkin(3) * t291 + t119 * t207) + t329) * g(2) + (-m(6) * (t349 * t213 + t352 * t217) - m(5) * (-pkin(3) * t290 - t119 * t204) + t328) * g(1); -pkin(4) * t306 + t152 * t75 + t153 * t74 + t42 * t315 + t220 - t331 * t43 + t338 * t109 + t339 * t108 + t334 * g(3) + t333 * g(2) + t332 * g(1) + (-g(1) * ((t76 * t212 - t216 * t335) * t217 + (t212 * t335 + t76 * t216) * t213) - t43 * t307 - pkin(4) * t205 * t8 + t152 * t4 + t153 * t3 - g(3) * ((t161 * t212 + t151) * t213 - (t161 * t216 - t294) * t217) * t206 - g(2) * ((t77 * t212 + t216 * t242) * t217 + (-t212 * t242 + t77 * t216) * t213) + t338 * t15 + t339 * t14) * m(6); -t3 * mrSges(6,2) + t4 * mrSges(6,1) - g(1) * ((t121 * t186 + t127 * t185 + t204 * t252) * mrSges(6,1) + (t120 * t186 + t126 * t185 - t204 * t253) * mrSges(6,2)) - g(2) * ((t123 * t186 + t125 * t185 - t207 * t252) * mrSges(6,1) + (-t122 * t186 + t124 * t185 + t207 * t253) * mrSges(6,2)) - g(3) * (-t205 * mrSges(6,2) * t281 + t168 * mrSges(6,1) + ((-t185 * t210 + t186 * t282) * mrSges(6,1) + (-t185 * t214 - t186 * t283) * mrSges(6,2)) * t206) + t249 - (-Ifges(6,2) * t255 + t160 + t80) * t254 / 0.2e1 + (t241 + t108) * t15 + (-t109 + t240) * t14 + (-t223 + t79 * t246 - t222 / 0.2e1 - t336 / 0.2e1) * t193;];
tau = t1;
