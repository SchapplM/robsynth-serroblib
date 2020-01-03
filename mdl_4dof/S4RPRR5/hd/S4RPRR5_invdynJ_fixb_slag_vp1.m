% Calculate vector of inverse dynamics joint torques for
% S4RPRR5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:32
% EndTime: 2019-12-31 16:51:39
% DurationCPUTime: 5.14s
% Computational Cost: add. (4751->350), mult. (10372->457), div. (0->0), fcn. (10714->6), ass. (0->182)
t299 = sin(qJ(3));
t300 = sin(qJ(1));
t301 = cos(qJ(3));
t302 = cos(qJ(1));
t128 = t302 * t299 - t300 * t301;
t171 = qJDD(1) - qJDD(3);
t172 = qJD(1) - qJD(3);
t169 = t302 * pkin(2);
t175 = qJD(1) ^ 2;
t162 = qJD(2) * t300;
t231 = qJD(1) * t300;
t255 = t302 * pkin(1) + t300 * qJ(2);
t232 = qJD(1) * t302;
t258 = qJ(2) * t232 + t162;
t194 = qJDD(1) * t255 - qJDD(2) * t302 + (-pkin(1) * t231 + t162 + t258) * qJD(1);
t245 = t300 * pkin(2);
t179 = qJDD(1) * t169 - t175 * t245 + t194;
t173 = sin(qJ(4));
t174 = cos(qJ(4));
t221 = rSges(5,1) * t173 + rSges(5,2) * t174;
t144 = -rSges(5,1) * t174 + rSges(5,2) * t173;
t129 = t144 * qJD(4);
t252 = qJD(4) * t129;
t127 = -t300 * t299 - t302 * t301;
t268 = t127 * t174;
t269 = t127 * t173;
t75 = -rSges(5,1) * t268 + rSges(5,2) * t269 + t128 * rSges(5,3);
t279 = -t127 * pkin(3) + pkin(6) * t128 + t75;
t250 = qJD(4) * t173;
t96 = t172 * t128;
t202 = t127 * t250 + t174 * t96;
t249 = qJD(4) * t174;
t203 = t127 * t249 - t173 * t96;
t95 = t172 * t127;
t27 = rSges(5,1) * t202 + t203 * rSges(5,2) + t95 * rSges(5,3);
t310 = t96 * pkin(3) + t95 * pkin(6) + t27;
t76 = qJD(4) * t95 + qJDD(4) * t128;
t8 = -t128 * t252 + t279 * t171 + t172 * t310 + t76 * t221 + t179;
t341 = g(2) - t8;
t165 = t302 * qJ(2);
t246 = t300 * pkin(1);
t141 = t246 - t165;
t163 = qJD(2) * t302;
t191 = -qJDD(1) * t141 + qJDD(2) * t300 + (-qJD(1) * t255 + 0.2e1 * t163) * qJD(1);
t177 = (-qJDD(1) * t300 - t175 * t302) * pkin(2) + t191;
t200 = t128 * t250 - t174 * t95;
t201 = t128 * t249 + t173 * t95;
t26 = rSges(5,1) * t200 + t201 * rSges(5,2) + t96 * rSges(5,3);
t320 = t95 * pkin(3) - t96 * pkin(6) - t26;
t265 = t128 * t174;
t266 = t128 * t173;
t74 = -rSges(5,1) * t265 + rSges(5,2) * t266 - t127 * rSges(5,3);
t40 = -t128 * pkin(3) - t127 * pkin(6) + t74;
t77 = qJD(4) * t96 - qJDD(4) * t127;
t7 = -t127 * t252 - t40 * t171 + t172 * t320 - t77 * t221 + t177;
t340 = -g(1) + t7;
t261 = t128 * rSges(4,1) - t127 * rSges(4,2);
t49 = -t95 * rSges(4,1) - t96 * rSges(4,2);
t18 = t171 * t261 - t172 * t49 + t177;
t339 = -g(1) + t18;
t338 = t172 * t40;
t277 = Icges(5,4) * t174;
t209 = -Icges(5,2) * t173 + t277;
t64 = -Icges(5,6) * t128 + t127 * t209;
t282 = t173 * t64;
t278 = Icges(5,4) * t173;
t211 = Icges(5,1) * t174 - t278;
t68 = -Icges(5,5) * t128 + t127 * t211;
t212 = t174 * t68 - t282;
t65 = Icges(5,6) * t127 + t128 * t209;
t283 = t173 * t65;
t69 = Icges(5,5) * t127 + t128 * t211;
t336 = -t174 * t69 + t283;
t207 = Icges(5,5) * t174 - Icges(5,6) * t173;
t61 = Icges(5,3) * t127 + t128 * t207;
t293 = -t127 * t61 - t265 * t69;
t60 = -Icges(5,3) * t128 + t127 * t207;
t335 = -(t60 - t283) * t128 + t293;
t251 = qJD(4) * t221;
t334 = t128 * t251 + t172 * t279;
t213 = -t173 * t68 - t174 * t64;
t215 = -t173 * t69 - t174 * t65;
t253 = qJD(4) * t128;
t254 = qJD(4) * t127;
t303 = -t172 / 0.2e1;
t333 = ((-t213 * t127 + t215 * t128) * qJD(4) + t213 * t254 - t215 * t253) * t303;
t50 = t96 * rSges(4,1) - t95 * rSges(4,2);
t97 = rSges(4,1) * t127 + rSges(4,2) * t128;
t19 = -t171 * t97 + t172 * t50 + t179;
t332 = -g(2) + t19;
t331 = t128 * t60;
t330 = t128 * t61;
t326 = t61 + t282;
t325 = t172 * t261;
t321 = t172 * t97;
t206 = Icges(5,5) * t173 + Icges(5,6) * t174;
t80 = t206 * t127;
t79 = t206 * t128;
t319 = qJD(1) * t141 - t162 + t258;
t239 = t169 + t255;
t188 = qJD(1) * t239 - t163;
t208 = Icges(5,2) * t174 + t278;
t210 = Icges(5,1) * t173 + t277;
t205 = -t173 * t208 + t174 * t210;
t41 = t128 * t205 + t80;
t38 = t41 * t172;
t42 = t127 * t205 - t79;
t39 = t42 * t172;
t256 = t302 * rSges(3,1) + t300 * rSges(3,3);
t313 = t255 + t256;
t309 = t127 * t251 - t338;
t196 = -t246 - t245;
t307 = pkin(2) * t231 + t196 * qJD(1) + t319;
t287 = t208 * t128 - t69;
t289 = -t210 * t128 - t65;
t306 = t173 * t287 + t174 * t289;
t292 = t127 * t60 + t265 * t68;
t291 = t268 * t69 - t330;
t290 = t268 * t68 - t331;
t288 = -t210 * t127 - t64;
t286 = t208 * t127 - t68;
t264 = t207 * t172;
t260 = -t208 + t211;
t259 = -t209 - t210;
t248 = -t302 / 0.2e1;
t247 = t300 / 0.2e1;
t240 = t300 * rSges(3,1);
t230 = t254 / 0.2e1;
t229 = -t253 / 0.2e1;
t228 = t253 / 0.2e1;
t23 = Icges(5,4) * t202 + Icges(5,2) * t203 + Icges(5,6) * t95;
t25 = Icges(5,1) * t202 + Icges(5,4) * t203 + Icges(5,5) * t95;
t180 = qJD(4) * t213 + t173 * t23 - t174 * t25;
t22 = Icges(5,4) * t200 + Icges(5,2) * t201 + Icges(5,6) * t96;
t24 = Icges(5,1) * t200 + Icges(5,4) * t201 + Icges(5,5) * t96;
t181 = qJD(4) * t215 + t173 * t22 - t174 * t24;
t20 = Icges(5,5) * t200 + Icges(5,6) * t201 + Icges(5,3) * t96;
t21 = Icges(5,5) * t202 + Icges(5,6) * t203 + Icges(5,3) * t95;
t223 = -(-t127 * t20 + t128 * t181 - t336 * t95 - t61 * t96) * t127 + t128 * (-t127 * t21 + t128 * t180 + t212 * t95 - t60 * t96);
t222 = -t127 * (t127 * t181 + t128 * t20 + t336 * t96 - t61 * t95) + t128 * (t127 * t180 + t128 * t21 - t212 * t96 - t60 * t95);
t14 = -t266 * t65 - t293;
t15 = -t266 * t64 + t292;
t220 = -t127 * t14 + t128 * t15;
t16 = -t269 * t65 + t291;
t17 = -t269 * t64 + t290;
t219 = -t127 * t16 + t128 * t17;
t218 = t127 * t27 + t128 * t26;
t190 = t162 + (-t245 - t141) * qJD(1);
t30 = t190 + t309;
t31 = t188 + t334;
t217 = -t127 * t30 - t128 * t31;
t216 = t127 * t75 + t128 * t74;
t104 = -t173 * t210 - t174 * t208;
t195 = -t240 - t246;
t147 = rSges(2,1) * t302 - rSges(2,2) * t300;
t143 = rSges(2,1) * t300 + rSges(2,2) * t302;
t189 = t173 * t286 + t174 * t288;
t187 = t165 + t196;
t184 = (t173 * t259 + t174 * t260) * t172;
t125 = t209 * qJD(4);
t126 = t211 * qJD(4);
t178 = qJD(4) * t104 - t125 * t173 + t126 * t174;
t10 = qJD(4) * t212 - t173 * t25 - t174 * t23;
t124 = t207 * qJD(4);
t12 = t124 * t127 + t128 * t178 + t205 * t95 - t206 * t96;
t13 = -t124 * t128 + t127 * t178 - t205 * t96 - t206 * t95;
t5 = qJD(4) * t220 + t38;
t6 = qJD(4) * t219 + t39;
t9 = -t336 * qJD(4) - t173 * t24 - t174 * t22;
t176 = t172 * (-t126 * t173 + t208 * t250 + (-qJD(4) * t210 - t125) * t174) + t228 * t5 - (-t213 + t42) * t76 / 0.2e1 - (-t215 + t41) * t77 / 0.2e1 + (t10 + t13) * t229 + (-Icges(4,3) + t104) * t171 + (t12 + t6 + t9) * t230;
t167 = t302 * rSges(3,3);
t160 = rSges(3,3) * t232;
t142 = t240 - t167;
t86 = t221 * t127;
t85 = t221 * t128;
t59 = t188 - t321;
t58 = t190 + t325;
t44 = qJDD(1) * t256 + qJD(1) * (-rSges(3,1) * t231 + t160) + t194;
t43 = -qJDD(1) * t142 - t175 * t256 + t191;
t32 = t216 * qJD(4);
t11 = (-t127 * t65 - t128 * t64) * t173 + t291 + t292;
t1 = [-t176 + (-t38 + ((t16 + (-t212 + t61) * t128) * t128 + (t326 * t127 + t17 - t290 - t331) * t127) * qJD(4)) * t229 + (t39 + ((t336 * t128 + t14 + t290) * t128 + (-t326 * t128 - t11 + t15) * t127) * qJD(4)) * t230 - m(2) * (-g(1) * t143 + g(2) * t147) + t333 + (m(2) * (t143 ^ 2 + t147 ^ 2) + Icges(2,3) + Icges(3,2)) * qJDD(1) + (-t341 * (t239 + t279) + (-t188 + t320) * t30 + (-t221 * t254 + t30 + t307 + t310 + t338) * t31 + t340 * (t187 - t40)) * m(5) + (t332 * (-t97 + t239) + (-t188 - t49) * t58 + (t307 + t50 + t58 - t325) * t59 + t339 * (t187 + t261)) * m(4) + ((-g(2) + t44) * t313 + (-g(1) + t43) * (t165 + t167 + t195) + (t160 + (t142 + t195) * qJD(1) + t319) * (qJD(1) * t313 - t163)) * m(3); (-m(3) - m(4) - m(5)) * (g(1) * t300 - g(2) * t302) + 0.2e1 * (t247 * t7 + t248 * t8) * m(5) + 0.2e1 * (t18 * t247 + t19 * t248) * m(4) + 0.2e1 * (t247 * t43 + t248 * t44) * m(3); t176 + (t38 + ((t11 - t16) * t128 + (t212 * t127 - t17 + t335) * t127) * qJD(4)) * t229 + (-t39 + ((-t14 - t335) * t128 + (-t15 + (-t336 + t60) * t127 - t330) * t127) * qJD(4)) * t230 + t333 + (t340 * t40 + (t309 - t310) * t31 + (-t320 - t334) * t30 + t341 * t279) * m(5) + (-t50 * t59 + (t49 + t321) * t58 - (-t172 * t59 + t339) * t261 + t332 * t97) * m(4); t95 * t6 / 0.2e1 + t128 * (t222 * qJD(4) + t13 * t172 + t16 * t77 + t17 * t76 + t171 * t42) / 0.2e1 + t76 * t219 / 0.2e1 + (t16 * t96 + t17 * t95 + t222) * t228 + t96 * t5 / 0.2e1 - t127 * (t223 * qJD(4) + t12 * t172 + t14 * t77 + t15 * t76 + t171 * t41) / 0.2e1 + t77 * t220 / 0.2e1 - (t14 * t96 + t15 * t95 + t223) * t254 / 0.2e1 + t171 * (t127 * t215 - t128 * t213) / 0.2e1 + t172 * (t10 * t128 - t127 * t9 - t213 * t95 - t215 * t96) / 0.2e1 + ((t80 * t253 - t264) * t128 + (t184 + (-t306 * t127 + (-t79 + t189) * t128) * qJD(4)) * t127) * t229 + ((t79 * t254 + t264) * t127 + (t184 + (t189 * t128 + (-t306 - t80) * t127) * qJD(4)) * t128) * t230 + ((t173 * t260 - t174 * t259) * t172 + ((t127 * t287 - t128 * t286) * t174 + (-t127 * t289 + t128 * t288) * t173) * qJD(4)) * t303 + ((qJD(4) * t218 + t74 * t76 - t75 * t77) * t216 + t32 * (t74 * t95 - t75 * t96 + t218) + t217 * t129 - (-t127 * t7 - t128 * t8 + t30 * t96 - t31 * t95) * t221 - (-t30 * t85 + t31 * t86) * t172 - (t32 * (t127 * t86 + t128 * t85) + t217 * t144) * qJD(4) - g(1) * t86 - g(2) * t85 - g(3) * t144) * m(5);];
tau = t1;
