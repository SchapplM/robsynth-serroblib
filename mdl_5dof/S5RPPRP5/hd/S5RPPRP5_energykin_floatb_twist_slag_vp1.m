% Calculate kinetic energy for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:17
% EndTime: 2019-12-31 17:53:19
% DurationCPUTime: 2.45s
% Computational Cost: add. (766->221), mult. (1534->283), div. (0->0), fcn. (1487->6), ass. (0->117)
t277 = Icges(3,4) - Icges(4,5);
t276 = Icges(3,1) + Icges(4,1);
t275 = Icges(3,2) + Icges(4,3);
t192 = cos(pkin(7));
t274 = t277 * t192;
t191 = sin(pkin(7));
t273 = t277 * t191;
t272 = Icges(4,4) + Icges(3,5);
t271 = Icges(3,6) - Icges(4,6);
t270 = t275 * t191 - t274;
t269 = t276 * t192 - t273;
t268 = Icges(5,1) + Icges(6,1);
t267 = -Icges(5,4) + Icges(6,5);
t266 = Icges(6,4) + Icges(5,5);
t265 = Icges(5,2) + Icges(6,3);
t264 = Icges(6,2) + Icges(5,3);
t263 = Icges(5,6) - Icges(6,6);
t194 = sin(qJ(1));
t195 = cos(qJ(1));
t262 = t270 * t194 + t271 * t195;
t261 = -t271 * t194 + t270 * t195;
t260 = t269 * t194 - t272 * t195;
t259 = t272 * t194 + t269 * t195;
t258 = -t275 * t192 - t273;
t257 = t276 * t191 + t274;
t256 = t271 * t191 - t272 * t192;
t255 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t254 = rSges(6,1) + pkin(4);
t253 = rSges(6,3) + qJ(5);
t193 = sin(qJ(4));
t238 = cos(qJ(4));
t217 = t191 * t238;
t230 = t192 * t194;
t143 = t193 * t230 - t194 * t217;
t153 = t191 * t193 + t192 * t238;
t144 = t153 * t194;
t252 = -t263 * t143 + t266 * t144 + t264 * t195;
t229 = t192 * t195;
t145 = t193 * t229 - t195 * t217;
t146 = t153 * t195;
t251 = -t263 * t145 + t266 * t146 - t264 * t194;
t250 = t265 * t143 + t267 * t144 - t263 * t195;
t249 = t265 * t145 + t267 * t146 + t263 * t194;
t248 = t267 * t143 + t268 * t144 + t266 * t195;
t247 = t267 * t145 + t268 * t146 - t266 * t194;
t154 = -t192 * t193 + t217;
t246 = t265 * t153 + t267 * t154;
t245 = -t263 * t153 + t266 * t154;
t244 = t267 * t153 + t268 * t154;
t187 = V_base(6) + qJD(1);
t189 = Icges(2,4) * t195;
t235 = Icges(2,4) * t194;
t243 = (-t272 * t191 - t271 * t192) * t187 + (t256 * t194 + t255 * t195 + t235) * V_base(5) + (-t255 * t194 + t256 * t195 + t189) * V_base(4);
t242 = (t258 * t191 + t257 * t192) * t187 + (Icges(2,1) * t194 + t262 * t191 + t260 * t192 + t189) * V_base(5) + (Icges(2,1) * t195 + t261 * t191 + t259 * t192 - t235) * V_base(4);
t237 = pkin(3) * t191;
t167 = pkin(2) * t191 - qJ(3) * t192;
t236 = -pkin(5) - t167;
t228 = rSges(6,2) * t195 + t253 * t143 + t254 * t144;
t227 = -rSges(6,2) * t194 + t253 * t145 + t254 * t146;
t226 = t253 * t153 + t254 * t154;
t213 = pkin(2) * t192 + qJ(3) * t191;
t148 = t213 * t194;
t176 = t194 * pkin(1) - qJ(2) * t195;
t225 = -t148 - t176;
t149 = t213 * t195;
t178 = pkin(1) * t195 + t194 * qJ(2);
t224 = -t149 - t178;
t223 = qJD(3) * t191;
t222 = V_base(4) * t176 + V_base(3);
t221 = V_base(5) * pkin(5) + V_base(1);
t156 = pkin(3) * t230 + pkin(6) * t195;
t218 = -t156 + t225;
t216 = qJD(2) * t194 + t221;
t215 = rSges(3,1) * t192 - rSges(3,2) * t191;
t214 = rSges(4,1) * t192 + rSges(4,3) * t191;
t206 = -qJD(2) * t195 + t187 * t178 + V_base(2);
t205 = V_base(5) * t167 + t195 * t223 + t216;
t204 = -qJD(3) * t192 + V_base(4) * t148 + t222;
t203 = V_base(5) * t237 + t205;
t202 = t187 * t149 + t194 * t223 + t206;
t157 = pkin(3) * t229 - t194 * pkin(6);
t199 = V_base(4) * t156 + (-t157 + t224) * V_base(5) + t204;
t198 = t187 * t157 + (t236 - t237) * V_base(4) + t202;
t183 = -qJD(4) * t194 + V_base(4);
t182 = qJD(4) * t195 + V_base(5);
t179 = rSges(2,1) * t195 - t194 * rSges(2,2);
t177 = t194 * rSges(2,1) + rSges(2,2) * t195;
t171 = Icges(2,5) * t195 - Icges(2,6) * t194;
t170 = Icges(2,5) * t194 + Icges(2,6) * t195;
t169 = rSges(3,1) * t191 + rSges(3,2) * t192;
t168 = rSges(4,1) * t191 - rSges(4,3) * t192;
t160 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t159 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t158 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t141 = t194 * rSges(3,3) + t195 * t215;
t140 = t194 * rSges(4,2) + t195 * t214;
t139 = -rSges(3,3) * t195 + t194 * t215;
t138 = -rSges(4,2) * t195 + t194 * t214;
t124 = V_base(5) * rSges(2,3) - t177 * t187 + t221;
t123 = t179 * t187 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t122 = t177 * V_base(4) - t179 * V_base(5) + V_base(3);
t121 = rSges(5,1) * t154 - rSges(5,2) * t153;
t110 = rSges(5,1) * t146 - rSges(5,2) * t145 - rSges(5,3) * t194;
t108 = t144 * rSges(5,1) - t143 * rSges(5,2) + rSges(5,3) * t195;
t94 = t169 * V_base(5) + (-t139 - t176) * t187 + t216;
t93 = t187 * t141 + (-pkin(5) - t169) * V_base(4) + t206;
t92 = t139 * V_base(4) + (-t141 - t178) * V_base(5) + t222;
t91 = t168 * V_base(5) + (-t138 + t225) * t187 + t205;
t90 = t187 * t140 + (-t168 + t236) * V_base(4) + t202;
t89 = t138 * V_base(4) + (-t140 + t224) * V_base(5) + t204;
t88 = t121 * t182 + (-t108 + t218) * t187 + t203;
t87 = t187 * t110 - t183 * t121 + t198;
t86 = t108 * t183 - t110 * t182 + t199;
t85 = qJD(5) * t145 + t226 * t182 + (t218 - t228) * t187 + t203;
t84 = qJD(5) * t143 - t183 * t226 + t187 * t227 + t198;
t83 = qJD(5) * t153 - t182 * t227 + t183 * t228 + t199;
t1 = m(1) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(2) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + ((t246 * t143 + t244 * t144 + t245 * t195) * t187 + (t249 * t143 + t247 * t144 + t251 * t195) * t183 + (t250 * t143 + t248 * t144 + t252 * t195) * t182) * t182 / 0.2e1 + ((t246 * t145 + t244 * t146 - t245 * t194) * t187 + (t249 * t145 + t247 * t146 - t251 * t194) * t183 + (t250 * t145 + t248 * t146 - t252 * t194) * t182) * t183 / 0.2e1 + (Icges(1,1) * V_base(4) + t171 * t187 - t243 * t194 + t242 * t195) * V_base(4) / 0.2e1 + (Icges(1,2) * V_base(5) + t170 * t187 + t242 * t194 + t243 * t195) * V_base(5) / 0.2e1 + ((t249 * t153 + t247 * t154) * t183 + (t250 * t153 + t248 * t154) * t182 + (t260 * t191 - t262 * t192 + t170) * V_base(5) + (t259 * t191 - t261 * t192 + t171) * V_base(4) + (t246 * t153 + t244 * t154 + t257 * t191 - t258 * t192 + Icges(2,3)) * t187) * t187 / 0.2e1 + V_base(5) * V_base(4) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
