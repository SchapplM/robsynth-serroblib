% Calculate vector of inverse dynamics joint torques for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:26
% EndTime: 2019-12-31 17:02:31
% DurationCPUTime: 2.95s
% Computational Cost: add. (6708->347), mult. (5357->438), div. (0->0), fcn. (4060->8), ass. (0->212)
t183 = sin(qJ(1));
t274 = pkin(1) * qJD(1);
t237 = t183 * t274;
t179 = qJ(1) + qJ(2);
t173 = sin(t179);
t174 = cos(t179);
t124 = rSges(3,1) * t173 + rSges(3,2) * t174;
t178 = qJD(1) + qJD(2);
t261 = t124 * t178;
t90 = -t237 - t261;
t244 = qJD(4) * t178;
t105 = -qJDD(4) * t174 + t173 * t244;
t177 = pkin(7) + qJ(4);
t171 = sin(t177);
t172 = cos(t177);
t122 = rSges(5,1) * t172 - rSges(5,2) * t171;
t112 = t122 * qJD(4);
t275 = rSges(5,2) * t172;
t277 = rSges(5,1) * t171;
t121 = t275 + t277;
t182 = -pkin(6) - qJ(3);
t255 = t173 * t182;
t134 = t178 * t255;
t181 = cos(pkin(7));
t167 = pkin(3) * t181 + pkin(2);
t139 = t174 * t167;
t176 = qJDD(1) + qJDD(2);
t184 = cos(qJ(1));
t185 = qJD(1) ^ 2;
t202 = (-qJDD(1) * t183 - t184 * t185) * pkin(1);
t247 = qJD(3) * t178;
t199 = qJDD(3) * t173 + t174 * t247 + t202;
t158 = t174 * qJ(3);
t284 = pkin(2) * t173;
t123 = -t158 + t284;
t298 = -t173 * t167 - t174 * t182;
t75 = t123 + t298;
t260 = t171 * t173;
t140 = rSges(5,2) * t260;
t258 = t172 * t173;
t240 = rSges(5,1) * t258;
t85 = -t174 * rSges(5,3) - t140 + t240;
t243 = -t123 + t75 - t85;
t245 = qJD(4) * t174;
t296 = t174 * pkin(2) + t173 * qJ(3);
t257 = t172 * t174;
t207 = rSges(5,1) * t257 + t173 * rSges(5,3);
t235 = qJD(4) * t275;
t259 = t171 * t174;
t238 = rSges(5,2) * t259;
t246 = qJD(4) * t173;
t234 = -t173 * t235 - t178 * t238 - t246 * t277;
t53 = t178 * t207 + t234;
t156 = qJD(3) * t174;
t87 = t178 * t296 - t156;
t9 = -t112 * t245 + t105 * t121 + t243 * t176 + (-t87 + t134 - t53 + (t296 - t139) * t178) * t178 + t199;
t305 = -g(1) + t9;
t180 = sin(pkin(7));
t276 = rSges(4,2) * t180;
t239 = t174 * t276;
t129 = t178 * t239;
t278 = rSges(4,1) * t181;
t241 = t173 * t278;
t150 = t173 * t276;
t248 = t174 * rSges(4,3) + t150;
t88 = t241 - t248;
t268 = -t123 - t88;
t297 = -t173 * rSges(4,3) - t174 * t278;
t22 = t268 * t176 + (t297 * t178 + t129 - t87) * t178 + t199;
t304 = -g(1) + t22;
t254 = t174 * t178;
t256 = t173 * t178;
t102 = rSges(3,1) * t254 - rSges(3,2) * t256;
t303 = -t102 * t178 - t124 * t176 - g(1) + t202;
t104 = qJDD(4) * t173 + t174 * t244;
t144 = qJ(3) * t254;
t175 = t184 * pkin(1);
t285 = pkin(1) * t183;
t224 = qJDD(1) * t175 - t185 * t285;
t155 = qJD(3) * t173;
t249 = t144 + t155;
t194 = -qJDD(3) * t174 + t176 * t296 + t173 * t247 + t224 + t178 * (-pkin(2) * t256 + t249);
t86 = -t238 + t207;
t61 = t139 - t255 + t86;
t281 = -t296 + t61;
t206 = rSges(5,3) * t254 + t178 * t140 - t174 * t235;
t233 = t171 * t245;
t52 = (-t172 * t256 - t233) * rSges(5,1) + t206;
t10 = -t112 * t246 - t104 * t121 + t281 * t176 + (-t144 + t52 + (t298 + t284) * t178) * t178 + t194;
t302 = -g(2) + t10;
t251 = rSges(4,3) * t254 + t178 * t150;
t89 = -t239 - t297;
t23 = t176 * t89 + t178 * (-t178 * t241 + t251) + t194;
t301 = -g(2) + t23;
t126 = t174 * rSges(3,1) - rSges(3,2) * t173;
t300 = t126 * t176 - t178 * t261 - g(2) + t224;
t69 = t296 + t89;
t299 = t178 * t69;
t154 = Icges(5,4) * t172;
t211 = -Icges(5,2) * t171 + t154;
t119 = Icges(5,1) * t171 + t154;
t108 = t178 * t123;
t295 = t178 * t88 + t108 + t251;
t116 = Icges(5,5) * t172 - Icges(5,6) * t171;
t115 = Icges(5,5) * t171 + Icges(5,6) * t172;
t195 = Icges(5,3) * t178 - qJD(4) * t115;
t204 = t211 * t174;
t81 = Icges(5,6) * t173 + t204;
t272 = t171 * t81;
t266 = Icges(5,4) * t171;
t120 = Icges(5,1) * t172 - t266;
t205 = t120 * t174;
t83 = Icges(5,5) * t173 + t205;
t218 = -t172 * t83 + t272;
t294 = -t116 * t256 + t174 * t195 + t178 * t218;
t203 = t116 * t174;
t80 = Icges(5,4) * t258 - Icges(5,2) * t260 - Icges(5,6) * t174;
t273 = t171 * t80;
t137 = Icges(5,4) * t260;
t82 = Icges(5,1) * t258 - Icges(5,5) * t174 - t137;
t219 = -t172 * t82 + t273;
t293 = t173 * t195 + (t203 + t219) * t178;
t117 = Icges(5,2) * t172 + t266;
t209 = t171 * t117 - t172 * t119;
t292 = t116 * qJD(4) + t178 * t209;
t78 = Icges(5,5) * t258 - Icges(5,6) * t260 - Icges(5,3) * t174;
t29 = -t173 * t219 - t174 * t78;
t74 = t178 * t85;
t291 = -rSges(5,1) * t233 - t178 * t75 + t108 + t155 + t206 + t74;
t290 = t173 * (-t117 * t174 + t83) - t174 * (-Icges(5,2) * t258 - t137 + t82);
t289 = t104 / 0.2e1;
t288 = t105 / 0.2e1;
t287 = t173 / 0.2e1;
t286 = -t174 / 0.2e1;
t283 = -t173 * t78 - t82 * t257;
t79 = Icges(5,3) * t173 + t203;
t282 = t173 * t79 + t83 * t257;
t221 = -t121 * t245 + t155;
t200 = t221 - t237;
t27 = t178 * t243 + t200;
t270 = t178 * t27;
t263 = t115 * t174;
t42 = -t173 * t209 - t263;
t269 = t42 * t178;
t264 = t115 * t173;
t262 = t116 * t178;
t253 = -t117 + t120;
t252 = t119 + t211;
t242 = t296 + t281;
t236 = t184 * t274;
t232 = -pkin(2) - t278;
t231 = -t246 / 0.2e1;
t230 = t246 / 0.2e1;
t229 = -t245 / 0.2e1;
t228 = t245 / 0.2e1;
t227 = -t78 + t272;
t65 = t83 * t258;
t226 = t174 * t79 - t65;
t222 = -t156 + t236;
t149 = rSges(2,1) * t184 - rSges(2,2) * t183;
t148 = rSges(2,1) * t183 + rSges(2,2) * t184;
t40 = t171 * t82 + t172 * t80;
t41 = t171 * t83 + t172 * t81;
t100 = t121 * t246;
t28 = t178 * t242 - t100 + t222;
t217 = -t28 * t173 - t27 * t174;
t30 = -t260 * t81 - t226;
t216 = t173 * t30 - t174 * t29;
t31 = -t259 * t80 - t283;
t32 = -t259 * t81 + t282;
t215 = t173 * t32 - t174 * t31;
t214 = t173 * t85 + t174 * t86;
t213 = t134 + t156 - t234;
t210 = t117 * t172 + t119 * t171;
t60 = -t85 + t298;
t201 = -t81 * t173 + t80 * t174;
t68 = t173 * t232 + t158 + t248;
t198 = (-t171 * t252 + t172 * t253) * t178;
t197 = Icges(5,5) * t178 - qJD(4) * t119;
t196 = Icges(5,6) * t178 - qJD(4) * t117;
t11 = qJD(4) * t216 + t269;
t110 = t211 * qJD(4);
t111 = t120 * qJD(4);
t43 = -t174 * t209 + t264;
t38 = t43 * t178;
t12 = qJD(4) * t215 + t38;
t49 = t173 * t196 + t178 * t204;
t51 = t173 * t197 + t178 * t205;
t15 = -qJD(4) * t219 + t171 * t51 + t172 * t49;
t48 = t174 * t196 - t211 * t256;
t50 = -t120 * t256 + t174 * t197;
t16 = -qJD(4) * t218 + t171 * t50 + t172 * t48;
t189 = -qJD(4) * t210 - t110 * t171 + t111 * t172 + t115 * t178;
t19 = t292 * t173 + t189 * t174;
t20 = t189 * t173 - t292 * t174;
t192 = (t38 + ((t30 - t65 + (t79 + t273) * t174 + t283) * t174 + t282 * t173) * qJD(4)) * t228 + (-qJD(4) * t209 + t110 * t172 + t111 * t171) * t178 + (t41 + t43) * t289 + (t40 + t42) * t288 + (-t269 + ((t174 * t227 - t282 + t32) * t174 + (t173 * t227 + t226 + t31) * t173) * qJD(4) + t11) * t231 + (t16 + t19) * t230 + (t15 + t20 + t12) * t229 + (Icges(3,3) + t210 + Icges(4,2) * t181 ^ 2 + (Icges(4,1) * t180 + 0.2e1 * Icges(4,4) * t181) * t180) * t176;
t191 = -qJD(4) * t41 - t171 * t48 + t172 * t50 + t178 * t79;
t190 = -qJD(4) * t40 - t171 * t49 + t172 * t51 + t178 * t78;
t54 = t178 * t268 + t155 - t237;
t55 = t222 + t299;
t188 = (t54 * t232 * t174 + (t54 * (-rSges(4,3) - qJ(3)) + t55 * t232) * t173) * t178;
t187 = -t290 * t171 + t201 * t172;
t186 = (t27 * (-t207 - t139) + t28 * (t298 - t240)) * t178;
t99 = t121 * t174;
t98 = t121 * t173;
t91 = t126 * t178 + t236;
t39 = t214 * qJD(4);
t6 = t191 * t173 - t294 * t174;
t5 = t190 * t173 - t293 * t174;
t4 = t294 * t173 + t191 * t174;
t3 = t293 * t173 + t190 * t174;
t1 = [Icges(2,3) * qJDD(1) + t192 + (t300 * (t126 + t175) + t303 * (-t124 - t285) + (-t102 - t236 + t91) * t90) * m(3) + (g(1) * t148 - g(2) * t149 + (t148 ^ 2 + t149 ^ 2) * qJDD(1)) * m(2) + (t27 * (t213 - t236) + t186 + t302 * (t175 + t61) + t305 * (t60 - t285) + (-t200 + t27 - t237 + t291) * t28) * m(5) + (t54 * (t129 - t222) + t188 + t301 * (t175 + t69) + t304 * (t68 - t285) + (t54 + t144 + t295) * t55) * m(4); t192 + (t242 * t270 + t186 + t302 * t61 + t305 * t60 + (-t221 + t291) * t28 + (-t100 - t156 + t213) * t27) * m(5) + (t188 + t301 * t69 + t304 * t68 + (-t155 + t249 + t295) * t55 + (t129 + t299) * t54) * m(4) + (-t261 * t91 - t102 * t90 + (t178 * t90 + t300) * t126 + (t178 * t91 - t303) * t124) * m(3); (-m(4) - m(5)) * (g(1) * t173 - g(2) * t174) + 0.2e1 * (t10 * t286 + t287 * t9) * m(5) + 0.2e1 * (t22 * t287 + t23 * t286) * m(4); t12 * t254 / 0.2e1 + (t104 * t32 + t105 * t31 + t176 * t43 + t178 * t19 + (t173 * t4 - t174 * t3) * qJD(4)) * t287 + t215 * t289 + ((t178 * t32 - t3) * t174 + (t178 * t31 + t4) * t173) * t230 + t11 * t256 / 0.2e1 + (t104 * t30 + t105 * t29 + t176 * t42 + t178 * t20 + (t173 * t6 - t174 * t5) * qJD(4)) * t286 + t216 * t288 + ((t178 * t30 - t5) * t174 + (t178 * t29 + t6) * t173) * t229 + t176 * (t173 * t41 - t174 * t40) / 0.2e1 + t178 * ((t178 * t41 - t15) * t174 + (t178 * t40 + t16) * t173) / 0.2e1 + ((-t246 * t263 + t262) * t173 + (t198 + (t173 * t264 + t187) * qJD(4)) * t174) * t231 + ((-t245 * t264 - t262) * t174 + (t198 + (t174 * t263 + t187) * qJD(4)) * t173) * t228 - t178 * ((t171 * t253 + t172 * t252) * t178 + (t201 * t171 + t290 * t172) * qJD(4)) / 0.2e1 + ((t104 * t85 - t105 * t86 + (t173 * t53 + t174 * t52) * qJD(4)) * t214 + t39 * ((t52 + t74) * t174 + (-t178 * t86 + t53) * t173) + t217 * t112 + ((-t178 * t28 - t9) * t174 + (-t10 + t270) * t173) * t121 - (t27 * t98 - t28 * t99) * t178 - (t39 * (-t173 * t98 - t174 * t99) + t217 * t122) * qJD(4) + g(1) * t99 + g(2) * t98 - g(3) * t122) * m(5);];
tau = t1;
