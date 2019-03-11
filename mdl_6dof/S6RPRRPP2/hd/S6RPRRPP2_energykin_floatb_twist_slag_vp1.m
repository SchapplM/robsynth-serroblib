% Calculate kinetic energy for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:24
% EndTime: 2019-03-09 04:31:26
% DurationCPUTime: 2.61s
% Computational Cost: add. (1812->268), mult. (2001->365), div. (0->0), fcn. (1917->8), ass. (0->131)
t290 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t289 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t288 = Icges(7,5) - Icges(6,4) - Icges(5,5);
t287 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t286 = -Icges(6,6) + Icges(7,6) + Icges(5,6);
t285 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t284 = rSges(7,1) + pkin(5);
t283 = rSges(7,3) + qJ(6);
t223 = qJ(1) + pkin(9);
t217 = sin(t223);
t218 = cos(t223);
t227 = cos(qJ(4));
t224 = sin(qJ(4));
t228 = cos(qJ(3));
t258 = t224 * t228;
t160 = t217 * t258 + t218 * t227;
t257 = t227 * t228;
t161 = t217 * t257 - t218 * t224;
t225 = sin(qJ(3));
t260 = t217 * t225;
t280 = t286 * t160 + t288 * t161 - t285 * t260;
t162 = -t217 * t227 + t218 * t258;
t163 = t217 * t224 + t218 * t257;
t259 = t218 * t225;
t279 = t286 * t162 + t288 * t163 - t285 * t259;
t278 = t287 * t160 + t289 * t161 - t286 * t260;
t277 = t287 * t162 + t289 * t163 - t286 * t259;
t276 = t289 * t160 + t290 * t161 - t288 * t260;
t275 = t289 * t162 + t290 * t163 - t288 * t259;
t274 = t285 * t228 + (t286 * t224 + t288 * t227) * t225;
t273 = t286 * t228 + (t287 * t224 + t289 * t227) * t225;
t272 = t288 * t228 + (t289 * t224 + t290 * t227) * t225;
t226 = sin(qJ(1));
t267 = pkin(1) * t226;
t229 = cos(qJ(1));
t266 = pkin(1) * t229;
t265 = -pkin(6) - qJ(2);
t264 = Icges(2,4) * t226;
t263 = Icges(3,4) * t217;
t262 = Icges(4,4) * t225;
t261 = Icges(4,4) * t228;
t256 = rSges(7,2) * t160 + t284 * t161 - t283 * t260;
t255 = rSges(7,2) * t162 + t284 * t163 - t283 * t259;
t254 = t283 * t228 + (rSges(7,2) * t224 + t284 * t227) * t225;
t253 = qJD(4) * t225;
t252 = qJD(6) * t225;
t219 = V_base(6) + qJD(1);
t251 = t219 * t266 + V_base(2);
t250 = V_base(5) * pkin(6) + V_base(1);
t196 = qJD(3) * t217 + V_base(4);
t189 = pkin(2) * t217 - pkin(7) * t218;
t247 = -t189 - t267;
t246 = V_base(5) * qJ(2) + t250;
t245 = V_base(4) * t267 + qJD(2) + V_base(3);
t244 = pkin(3) * t228 + pkin(8) * t225;
t195 = -qJD(3) * t218 + V_base(5);
t243 = rSges(4,1) * t228 - rSges(4,2) * t225;
t242 = Icges(4,1) * t228 - t262;
t241 = -Icges(4,2) * t225 + t261;
t240 = Icges(4,5) * t228 - Icges(4,6) * t225;
t239 = (-Icges(4,3) * t218 + t217 * t240) * t195 + (Icges(4,3) * t217 + t218 * t240) * t196 + (Icges(4,5) * t225 + Icges(4,6) * t228) * t219;
t190 = pkin(2) * t218 + pkin(7) * t217;
t238 = t219 * t190 + t265 * V_base(4) + t251;
t177 = t244 * t217;
t210 = pkin(3) * t225 - pkin(8) * t228;
t237 = t195 * t210 + (-t177 + t247) * t219 + t246;
t236 = V_base(4) * t189 + (-t190 - t266) * V_base(5) + t245;
t178 = t244 * t218;
t235 = t219 * t178 - t196 * t210 + t238;
t158 = t217 * t253 + t195;
t180 = (pkin(4) * t227 + qJ(5) * t224) * t225;
t234 = qJD(5) * t162 + t158 * t180 + t237;
t138 = pkin(4) * t163 + qJ(5) * t162;
t206 = -qJD(4) * t228 + t219;
t233 = qJD(5) * t160 + t206 * t138 + t235;
t232 = t196 * t177 - t195 * t178 + t236;
t137 = pkin(4) * t161 + qJ(5) * t160;
t159 = t218 * t253 + t196;
t231 = qJD(5) * t225 * t224 + t159 * t137 + t232;
t146 = -Icges(4,6) * t218 + t217 * t241;
t147 = Icges(4,6) * t217 + t218 * t241;
t148 = -Icges(4,5) * t218 + t217 * t242;
t149 = Icges(4,5) * t217 + t218 * t242;
t200 = Icges(4,2) * t228 + t262;
t203 = Icges(4,1) * t225 + t261;
t230 = (-t147 * t225 + t149 * t228) * t196 + (-t146 * t225 + t148 * t228) * t195 + (-t200 * t225 + t203 * t228) * t219;
t221 = Icges(2,4) * t229;
t216 = Icges(3,4) * t218;
t209 = rSges(2,1) * t229 - t226 * rSges(2,2);
t208 = t226 * rSges(2,1) + rSges(2,2) * t229;
t207 = rSges(4,1) * t225 + rSges(4,2) * t228;
t205 = Icges(2,1) * t229 - t264;
t204 = Icges(2,1) * t226 + t221;
t202 = -Icges(2,2) * t226 + t221;
t201 = Icges(2,2) * t229 + t264;
t194 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t193 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t192 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t188 = rSges(3,1) * t218 - rSges(3,2) * t217;
t187 = rSges(3,1) * t217 + rSges(3,2) * t218;
t186 = Icges(3,1) * t218 - t263;
t185 = Icges(3,1) * t217 + t216;
t184 = -Icges(3,2) * t217 + t216;
t183 = Icges(3,2) * t218 + t263;
t175 = -rSges(5,3) * t228 + (rSges(5,1) * t227 - rSges(5,2) * t224) * t225;
t174 = -rSges(6,2) * t228 + (rSges(6,1) * t227 + rSges(6,3) * t224) * t225;
t153 = rSges(4,3) * t217 + t218 * t243;
t152 = -rSges(4,3) * t218 + t217 * t243;
t151 = V_base(5) * rSges(2,3) - t208 * t219 + t250;
t150 = t209 * t219 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t143 = t208 * V_base(4) - t209 * V_base(5) + V_base(3);
t136 = V_base(5) * rSges(3,3) + (-t187 - t267) * t219 + t246;
t135 = t188 * t219 + (-rSges(3,3) + t265) * V_base(4) + t251;
t134 = V_base(4) * t187 + (-t188 - t266) * V_base(5) + t245;
t133 = rSges(5,1) * t163 - rSges(5,2) * t162 + rSges(5,3) * t259;
t132 = rSges(6,1) * t163 + rSges(6,2) * t259 + rSges(6,3) * t162;
t130 = rSges(5,1) * t161 - rSges(5,2) * t160 + rSges(5,3) * t260;
t129 = rSges(6,1) * t161 + rSges(6,2) * t260 + rSges(6,3) * t160;
t107 = t195 * t207 + (-t152 + t247) * t219 + t246;
t106 = t153 * t219 - t196 * t207 + t238;
t105 = t196 * t152 - t195 * t153 + t236;
t104 = -t130 * t206 + t158 * t175 + t237;
t103 = t133 * t206 - t159 * t175 + t235;
t102 = t159 * t130 - t158 * t133 + t232;
t101 = t158 * t174 + (-t129 - t137) * t206 + t234;
t100 = t132 * t206 + (-t174 - t180) * t159 + t233;
t99 = -t218 * t252 + t254 * t158 + (-t137 - t256) * t206 + t234;
t98 = -t217 * t252 + t255 * t206 + (-t180 - t254) * t159 + t233;
t97 = t159 * t129 + (-t132 - t138) * t158 + t231;
t96 = qJD(6) * t228 + t256 * t159 + (-t138 - t255) * t158 + t231;
t1 = t196 * (t239 * t217 + t230 * t218) / 0.2e1 + t195 * (t230 * t217 - t239 * t218) / 0.2e1 + m(7) * (t96 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t101 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(3) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(2) * (t143 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(1) * (t192 ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + ((t160 * t273 + t161 * t272 - t260 * t274) * t206 + (t160 * t277 + t161 * t275 - t260 * t279) * t159 + (t278 * t160 + t276 * t161 - t280 * t260) * t158) * t158 / 0.2e1 + ((t162 * t273 + t163 * t272 - t259 * t274) * t206 + (t277 * t162 + t275 * t163 - t279 * t259) * t159 + (t162 * t278 + t163 * t276 - t259 * t280) * t158) * t159 / 0.2e1 + ((t158 * t280 + t159 * t279 + t274 * t206) * t228 + ((t224 * t273 + t227 * t272) * t206 + (t224 * t277 + t227 * t275) * t159 + (t224 * t278 + t227 * t276) * t158) * t225) * t206 / 0.2e1 + ((t147 * t228 + t149 * t225) * t196 + (t146 * t228 + t148 * t225) * t195 + (t228 * t200 + t225 * t203 + Icges(2,3) + Icges(3,3)) * t219) * t219 / 0.2e1 + ((-t183 * t217 + t185 * t218 - t226 * t201 + t204 * t229 + Icges(1,4)) * V_base(5) + (-t217 * t184 + t218 * t186 - t226 * t202 + t229 * t205 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t218 * t183 + t217 * t185 + t229 * t201 + t226 * t204 + Icges(1,2)) * V_base(5) + (t184 * t218 + t186 * t217 + t202 * t229 + t226 * t205 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t219 * (Icges(2,5) * t226 + Icges(3,5) * t217 + Icges(2,6) * t229 + Icges(3,6) * t218) + V_base(4) * t219 * (Icges(2,5) * t229 + Icges(3,5) * t218 - Icges(2,6) * t226 - Icges(3,6) * t217) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
