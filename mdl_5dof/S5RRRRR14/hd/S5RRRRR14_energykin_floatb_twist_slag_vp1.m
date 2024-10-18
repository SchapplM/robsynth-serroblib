% Calculate kinetic energy for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR14_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:17
% EndTime: 2024-09-27 18:42:19
% DurationCPUTime: 0.69s
% Computational Cost: add. (2795->340), mult. (1978->478), div. (0->0), fcn. (1752->22), ass. (0->166)
t301 = -pkin(6) - pkin(7);
t272 = pkin(8) + pkin(9);
t269 = sin(qJ(1));
t300 = pkin(1) * t269;
t271 = cos(qJ(1));
t299 = pkin(1) * t271;
t266 = cos(pkin(5));
t298 = t266 * pkin(8);
t270 = cos(qJ(3));
t297 = t270 * pkin(3);
t296 = Icges(2,4) * t269;
t264 = qJ(1) + qJ(2);
t254 = sin(t264);
t295 = Icges(3,4) * t254;
t265 = sin(pkin(5));
t294 = t254 * t265;
t256 = cos(t264);
t293 = t256 * t265;
t268 = sin(qJ(3));
t292 = t266 * t268;
t291 = t266 * t270;
t263 = qJ(3) + qJ(4);
t203 = pkin(3) * t292 - t265 * t272;
t246 = cos(qJ(4)) * pkin(4) + pkin(3);
t262 = pkin(10) + t272;
t283 = pkin(4) * sin(qJ(4)) * t270;
t290 = -t262 * t265 + (t246 * t268 + t283) * t266 - t203;
t255 = cos(t263);
t289 = pkin(4) * t255;
t288 = qJD(3) * t265;
t287 = -qJD(3) - qJD(4);
t252 = V_base(6) + qJD(1);
t286 = t252 * t299 + V_base(2);
t285 = V_base(4) * t300 + V_base(3);
t284 = V_base(5) * pkin(6) + V_base(1);
t213 = t254 * t288 + V_base(4);
t251 = pkin(5) - t263;
t250 = pkin(5) + t263;
t280 = pkin(8) * t265 + t203;
t180 = qJD(4) * t294 + t213;
t241 = qJD(2) + t252;
t218 = qJD(3) * t266 + t241;
t196 = qJD(4) * t266 + t218;
t279 = V_base(5) * pkin(7) - t252 * t300 + t284;
t197 = pkin(2) * t254 - pkin(8) * t293;
t198 = pkin(2) * t256 + pkin(8) * t294;
t278 = V_base(4) * t197 + (-t198 - t299) * V_base(5) + t285;
t277 = -t197 * t241 + V_base(5) * t298 + t279;
t276 = t241 * t198 + (-t298 + t301) * V_base(4) + t286;
t150 = t254 * t297 + t256 * t280;
t151 = -t254 * t280 + t256 * t297;
t212 = -t256 * t288 + V_base(5);
t275 = t213 * t150 - t151 * t212 + t278;
t191 = t265 * t268 * pkin(3) + (-pkin(8) + t272) * t266;
t274 = -t150 * t218 + t212 * t191 + t277;
t273 = t218 * t151 - t191 * t213 + t276;
t258 = qJ(5) + t263;
t257 = Icges(2,4) * t271;
t253 = sin(t263);
t245 = cos(t258);
t244 = sin(t258);
t243 = -qJ(5) + t251;
t242 = qJ(5) + t250;
t240 = cos(t250);
t239 = sin(t251);
t238 = Icges(3,4) * t256;
t235 = cos(t242);
t234 = sin(t243);
t233 = cos(t251) / 0.2e1;
t232 = sin(t250) / 0.2e1;
t231 = cos(t243) / 0.2e1;
t230 = sin(t242) / 0.2e1;
t228 = rSges(2,1) * t271 - rSges(2,2) * t269;
t227 = rSges(2,1) * t269 + rSges(2,2) * t271;
t224 = Icges(2,1) * t271 - t296;
t223 = Icges(2,1) * t269 + t257;
t222 = -Icges(2,2) * t269 + t257;
t221 = Icges(2,2) * t271 + t296;
t217 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t216 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t215 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t211 = rSges(3,1) * t256 - rSges(3,2) * t254;
t210 = rSges(3,1) * t254 + rSges(3,2) * t256;
t209 = Icges(3,1) * t256 - t295;
t208 = Icges(3,1) * t254 + t238;
t207 = -Icges(3,2) * t254 + t238;
t206 = Icges(3,2) * t256 + t295;
t202 = t233 - t240 / 0.2e1;
t201 = t233 + t240 / 0.2e1;
t200 = t232 - t239 / 0.2e1;
t199 = t232 + t239 / 0.2e1;
t195 = t231 - t235 / 0.2e1;
t194 = t231 + t235 / 0.2e1;
t193 = t230 - t234 / 0.2e1;
t192 = t230 + t234 / 0.2e1;
t189 = -t254 * t292 + t256 * t270;
t188 = -t254 * t291 - t256 * t268;
t187 = t254 * t270 + t256 * t292;
t186 = -t254 * t268 + t256 * t291;
t185 = rSges(4,3) * t266 + (rSges(4,1) * t268 + rSges(4,2) * t270) * t265;
t184 = Icges(4,5) * t266 + (Icges(4,1) * t268 + Icges(4,4) * t270) * t265;
t183 = Icges(4,6) * t266 + (Icges(4,4) * t268 + Icges(4,2) * t270) * t265;
t182 = Icges(4,3) * t266 + (Icges(4,5) * t268 + Icges(4,6) * t270) * t265;
t181 = qJD(5) * t266 + t196;
t179 = t287 * t293 + V_base(5);
t177 = V_base(5) * rSges(2,3) - t227 * t252 + t284;
t176 = t228 * t252 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t175 = t227 * V_base(4) - t228 * V_base(5) + V_base(3);
t173 = -t200 * t254 + t255 * t256;
t172 = -t201 * t254 - t253 * t256;
t171 = t200 * t256 + t254 * t255;
t170 = t201 * t256 - t253 * t254;
t168 = qJD(5) * t294 + t180;
t167 = V_base(5) + (-qJD(5) + t287) * t293;
t166 = -t193 * t254 + t245 * t256;
t165 = -t194 * t254 - t244 * t256;
t164 = t193 * t256 + t245 * t254;
t163 = t194 * t256 - t244 * t254;
t162 = rSges(5,1) * t202 + rSges(5,2) * t199 + rSges(5,3) * t266;
t161 = Icges(5,1) * t202 + Icges(5,4) * t199 + Icges(5,5) * t266;
t160 = Icges(5,4) * t202 + Icges(5,2) * t199 + Icges(5,6) * t266;
t159 = Icges(5,5) * t202 + Icges(5,6) * t199 + Icges(5,3) * t266;
t158 = rSges(6,1) * t195 + rSges(6,2) * t192 + rSges(6,3) * t266;
t157 = (t262 - t272) * t266 + (t283 + (-pkin(3) + t246) * t268) * t265;
t156 = Icges(6,1) * t195 + Icges(6,4) * t192 + Icges(6,5) * t266;
t155 = Icges(6,4) * t195 + Icges(6,2) * t192 + Icges(6,6) * t266;
t154 = Icges(6,5) * t195 + Icges(6,6) * t192 + Icges(6,3) * t266;
t153 = V_base(5) * rSges(3,3) - t210 * t241 + t279;
t152 = t211 * t241 + (-rSges(3,3) + t301) * V_base(4) + t286;
t149 = t210 * V_base(4) + (-t211 - t299) * V_base(5) + t285;
t148 = rSges(4,1) * t189 + rSges(4,2) * t188 + rSges(4,3) * t294;
t147 = rSges(4,1) * t187 + rSges(4,2) * t186 - rSges(4,3) * t293;
t146 = Icges(4,1) * t189 + Icges(4,4) * t188 + Icges(4,5) * t294;
t145 = Icges(4,1) * t187 + Icges(4,4) * t186 - Icges(4,5) * t293;
t144 = Icges(4,4) * t189 + Icges(4,2) * t188 + Icges(4,6) * t294;
t143 = Icges(4,4) * t187 + Icges(4,2) * t186 - Icges(4,6) * t293;
t142 = Icges(4,5) * t189 + Icges(4,6) * t188 + Icges(4,3) * t294;
t141 = Icges(4,5) * t187 + Icges(4,6) * t186 - Icges(4,3) * t293;
t138 = rSges(5,1) * t173 + rSges(5,2) * t172 + rSges(5,3) * t294;
t137 = rSges(5,1) * t171 + rSges(5,2) * t170 - rSges(5,3) * t293;
t136 = Icges(5,1) * t173 + Icges(5,4) * t172 + Icges(5,5) * t294;
t135 = Icges(5,1) * t171 + Icges(5,4) * t170 - Icges(5,5) * t293;
t134 = Icges(5,4) * t173 + Icges(5,2) * t172 + Icges(5,6) * t294;
t133 = Icges(5,4) * t171 + Icges(5,2) * t170 - Icges(5,6) * t293;
t132 = Icges(5,5) * t173 + Icges(5,6) * t172 + Icges(5,3) * t294;
t131 = Icges(5,5) * t171 + Icges(5,6) * t170 - Icges(5,3) * t293;
t130 = rSges(6,1) * t166 + rSges(6,2) * t165 + rSges(6,3) * t294;
t129 = rSges(6,1) * t164 + rSges(6,2) * t163 - rSges(6,3) * t293;
t128 = Icges(6,1) * t166 + Icges(6,4) * t165 + Icges(6,5) * t294;
t127 = Icges(6,1) * t164 + Icges(6,4) * t163 - Icges(6,5) * t293;
t126 = Icges(6,4) * t166 + Icges(6,2) * t165 + Icges(6,6) * t294;
t125 = Icges(6,4) * t164 + Icges(6,2) * t163 - Icges(6,6) * t293;
t124 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t294;
t123 = Icges(6,5) * t164 + Icges(6,6) * t163 - Icges(6,3) * t293;
t122 = -t254 * t290 + t256 * t289;
t121 = t254 * t289 + t256 * t290;
t120 = -t147 * t218 + t185 * t212 + t277;
t119 = t148 * t218 - t185 * t213 + t276;
t118 = t147 * t213 - t148 * t212 + t278;
t117 = -t137 * t196 + t162 * t179 + t274;
t116 = t138 * t196 - t162 * t180 + t273;
t115 = t137 * t180 - t138 * t179 + t275;
t114 = -t121 * t196 - t129 * t181 + t157 * t179 + t158 * t167 + t274;
t113 = t122 * t196 + t130 * t181 - t157 * t180 - t158 * t168 + t273;
t112 = t121 * t180 - t122 * t179 + t129 * t168 - t130 * t167 + t275;
t1 = m(1) * (t215 ^ 2 + t216 ^ 2 + t217 ^ 2) / 0.2e1 + m(2) * (t175 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + m(3) * (t149 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(4) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + t213 * ((t142 * t294 + t144 * t188 + t146 * t189) * t213 + (t141 * t294 + t143 * t188 + t145 * t189) * t212 + (t182 * t294 + t183 * t188 + t184 * t189) * t218) / 0.2e1 + t212 * ((-t142 * t293 + t144 * t186 + t146 * t187) * t213 + (-t141 * t293 + t143 * t186 + t145 * t187) * t212 + (-t182 * t293 + t183 * t186 + t184 * t187) * t218) / 0.2e1 + t218 * ((t141 * t212 + t142 * t213 + t182 * t218) * t266 + ((t144 * t270 + t146 * t268) * t213 + (t143 * t270 + t145 * t268) * t212 + (t183 * t270 + t184 * t268) * t218) * t265) / 0.2e1 + m(5) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + t180 * ((t132 * t294 + t172 * t134 + t173 * t136) * t180 + (t131 * t294 + t133 * t172 + t135 * t173) * t179 + (t159 * t294 + t160 * t172 + t161 * t173) * t196) / 0.2e1 + t179 * ((-t132 * t293 + t134 * t170 + t136 * t171) * t180 + (-t131 * t293 + t170 * t133 + t171 * t135) * t179 + (-t159 * t293 + t160 * t170 + t161 * t171) * t196) / 0.2e1 + t196 * ((t132 * t266 + t134 * t199 + t136 * t202) * t180 + (t131 * t266 + t133 * t199 + t135 * t202) * t179 + (t159 * t266 + t160 * t199 + t161 * t202) * t196) / 0.2e1 + m(6) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + t168 * ((t124 * t294 + t165 * t126 + t166 * t128) * t168 + (t123 * t294 + t125 * t165 + t127 * t166) * t167 + (t154 * t294 + t155 * t165 + t156 * t166) * t181) / 0.2e1 + t167 * ((-t124 * t293 + t126 * t163 + t128 * t164) * t168 + (-t123 * t293 + t163 * t125 + t164 * t127) * t167 + (-t154 * t293 + t155 * t163 + t156 * t164) * t181) / 0.2e1 + t181 * ((t124 * t266 + t126 * t192 + t128 * t195) * t168 + (t123 * t266 + t125 * t192 + t127 * t195) * t167 + (t266 * t154 + t192 * t155 + t195 * t156) * t181) / 0.2e1 + ((-t206 * t254 + t208 * t256 - t221 * t269 + t223 * t271 + Icges(1,4)) * V_base(5) + (-t207 * t254 + t209 * t256 - t222 * t269 + t224 * t271 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t206 * t256 + t208 * t254 + t221 * t271 + t223 * t269 + Icges(1,2)) * V_base(5) + (t207 * t256 + t209 * t254 + t222 * t271 + t224 * t269 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t269 + Icges(2,6) * t271) * V_base(5) + (Icges(2,5) * t271 - Icges(2,6) * t269) * V_base(4) + Icges(2,3) * t252 / 0.2e1) * t252 + ((Icges(3,5) * t254 + Icges(3,6) * t256) * V_base(5) + (Icges(3,5) * t256 - Icges(3,6) * t254) * V_base(4) + Icges(3,3) * t241 / 0.2e1) * t241;
T = t1;
