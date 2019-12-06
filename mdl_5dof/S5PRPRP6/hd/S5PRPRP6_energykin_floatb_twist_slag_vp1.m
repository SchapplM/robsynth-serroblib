% Calculate kinetic energy for
% S5PRPRP6
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:16
% EndTime: 2019-12-05 15:40:19
% DurationCPUTime: 2.78s
% Computational Cost: add. (788->234), mult. (1597->318), div. (0->0), fcn. (1501->6), ass. (0->119)
t279 = Icges(3,4) + Icges(4,6);
t278 = Icges(3,1) + Icges(4,2);
t277 = -Icges(3,2) - Icges(4,3);
t193 = cos(qJ(2));
t276 = t279 * t193;
t191 = sin(qJ(2));
t275 = t279 * t191;
t274 = Icges(4,4) - Icges(3,5);
t273 = Icges(4,5) - Icges(3,6);
t272 = t277 * t191 + t276;
t271 = t278 * t193 - t275;
t270 = Icges(4,1) + Icges(3,3);
t269 = Icges(5,1) + Icges(6,1);
t268 = Icges(5,4) - Icges(6,5);
t267 = Icges(6,4) + Icges(5,5);
t266 = Icges(5,2) + Icges(6,3);
t265 = Icges(6,6) - Icges(5,6);
t264 = Icges(5,3) + Icges(6,2);
t188 = sin(pkin(7));
t189 = cos(pkin(7));
t263 = t272 * t188 + t273 * t189;
t262 = -t273 * t188 + t272 * t189;
t261 = t271 * t188 + t274 * t189;
t260 = -t274 * t188 + t271 * t189;
t259 = t277 * t193 - t275;
t258 = t278 * t191 + t276;
t257 = t273 * t191 - t274 * t193;
t256 = rSges(6,1) + pkin(4);
t255 = rSges(6,3) + qJ(5);
t190 = sin(qJ(4));
t192 = cos(qJ(4));
t227 = t191 * t192;
t144 = t188 * t190 - t189 * t227;
t228 = t190 * t191;
t145 = t188 * t192 + t189 * t228;
t229 = t189 * t193;
t254 = t266 * t144 - t268 * t145 + t265 * t229;
t253 = t265 * t144 + t267 * t145 + t264 * t229;
t146 = t188 * t227 + t189 * t190;
t147 = t188 * t228 - t189 * t192;
t230 = t188 * t193;
t252 = -t265 * t146 + t267 * t147 + t264 * t230;
t251 = t266 * t146 + t268 * t147 - t265 * t230;
t248 = -t268 * t144 + t269 * t145 + t267 * t229;
t247 = t268 * t146 + t269 * t147 + t267 * t230;
t246 = (t268 * t190 + t266 * t192) * t193 + t265 * t191;
t245 = (-t267 * t190 + t265 * t192) * t193 + t264 * t191;
t244 = (-t269 * t190 - t268 * t192) * t193 + t267 * t191;
t180 = -qJD(2) * t189 + V_base(5);
t181 = qJD(2) * t188 + V_base(4);
t243 = (t259 * t191 + t258 * t193) * V_base(6) + (-t262 * t191 + t260 * t193) * t181 + (-t263 * t191 + t261 * t193) * t180;
t242 = (-t274 * t191 - t273 * t193) * V_base(6) + (t270 * t188 + t257 * t189) * t181 + (t257 * t188 - t270 * t189) * t180;
t236 = pkin(6) * t191;
t235 = Icges(2,4) * t188;
t226 = rSges(6,2) * t229 + t255 * t144 + t256 * t145;
t225 = rSges(6,2) * t230 - t255 * t146 + t256 * t147;
t224 = t191 * rSges(6,2) + (-t256 * t190 + t255 * t192) * t193;
t211 = pkin(2) * t193 + qJ(3) * t191;
t149 = t211 * t188;
t168 = pkin(1) * t188 - pkin(5) * t189;
t223 = -t149 - t168;
t222 = qJD(3) * t191;
t221 = qJD(3) * t193;
t220 = qJD(4) * t193;
t219 = V_base(5) * qJ(1) + V_base(1);
t215 = qJD(1) + V_base(3);
t176 = t191 * pkin(2) - qJ(3) * t193;
t214 = t180 * t176 + t189 * t222 + t219;
t213 = rSges(3,1) * t193 - rSges(3,2) * t191;
t212 = -rSges(4,2) * t193 + rSges(4,3) * t191;
t169 = pkin(1) * t189 + pkin(5) * t188;
t204 = -V_base(4) * qJ(1) + V_base(6) * t169 + V_base(2);
t203 = V_base(4) * t168 - V_base(5) * t169 + t215;
t150 = t211 * t189;
t202 = V_base(6) * t150 + t188 * t222 + t204;
t201 = t181 * t149 + t203;
t156 = -t189 * pkin(3) + pkin(6) * t230;
t198 = t180 * t236 + (-t156 + t223) * V_base(6) + t214;
t155 = t188 * pkin(3) + pkin(6) * t229;
t197 = t181 * t156 + (-t150 - t155) * t180 + t201;
t196 = V_base(6) * t155 + (-t176 - t236) * t181 + t202;
t186 = Icges(2,4) * t189;
t182 = qJD(4) * t191 + V_base(6);
t178 = t191 * rSges(3,1) + rSges(3,2) * t193;
t177 = -t191 * rSges(4,2) - rSges(4,3) * t193;
t167 = rSges(2,1) * t189 - rSges(2,2) * t188;
t166 = rSges(2,1) * t188 + rSges(2,2) * t189;
t165 = Icges(2,1) * t189 - t235;
t164 = Icges(2,1) * t188 + t186;
t163 = -Icges(2,2) * t188 + t186;
t162 = Icges(2,2) * t189 + t235;
t159 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t158 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t157 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t143 = t189 * t220 + t181;
t142 = t188 * t220 + t180;
t139 = t191 * rSges(5,3) + (-rSges(5,1) * t190 - rSges(5,2) * t192) * t193;
t131 = -t189 * rSges(4,1) + t188 * t212;
t130 = t188 * rSges(4,1) + t189 * t212;
t129 = t188 * rSges(3,3) + t189 * t213;
t128 = -t189 * rSges(3,3) + t188 * t213;
t113 = V_base(5) * rSges(2,3) - t166 * V_base(6) + t219;
t112 = t167 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t111 = t166 * V_base(4) - t167 * V_base(5) + t215;
t108 = t147 * rSges(5,1) + t146 * rSges(5,2) + rSges(5,3) * t230;
t106 = t145 * rSges(5,1) - t144 * rSges(5,2) + rSges(5,3) * t229;
t92 = t178 * t180 + (-t128 - t168) * V_base(6) + t219;
t91 = t129 * V_base(6) - t178 * t181 + t204;
t90 = t128 * t181 - t129 * t180 + t203;
t89 = t177 * t180 + (-t131 + t223) * V_base(6) + t214;
t88 = t130 * V_base(6) + (-t176 - t177) * t181 + t202;
t87 = -t221 + t181 * t131 + (-t130 - t150) * t180 + t201;
t86 = -t108 * t182 + t139 * t142 + t198;
t85 = t106 * t182 - t139 * t143 + t196;
t84 = -t142 * t106 + t143 * t108 + t197 - t221;
t83 = qJD(5) * t144 + t142 * t224 - t182 * t225 + t198;
t82 = -qJD(5) * t146 - t143 * t224 + t182 * t226 + t196;
t81 = (qJD(5) * t192 - qJD(3)) * t193 + t225 * t143 - t226 * t142 + t197;
t1 = m(1) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(2) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(3) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + ((-t146 * t246 + t147 * t244 + t230 * t245) * t182 + (-t146 * t254 + t248 * t147 + t253 * t230) * t143 + (t251 * t146 + t247 * t147 + t252 * t230) * t142) * t142 / 0.2e1 + ((t144 * t246 + t145 * t244 + t229 * t245) * t182 + (t254 * t144 + t248 * t145 + t253 * t229) * t143 + (-t251 * t144 + t145 * t247 + t229 * t252) * t142) * t143 / 0.2e1 + (t243 * t188 - t242 * t189) * t180 / 0.2e1 + (t242 * t188 + t243 * t189) * t181 / 0.2e1 + (((-t190 * t244 + t192 * t246) * t182 + (-t248 * t190 + t192 * t254) * t143 + (-t190 * t247 - t192 * t251) * t142) * t193 + (t252 * t142 + t143 * t253 + t245 * t182) * t191) * t182 / 0.2e1 + ((-t162 * t188 + t164 * t189 + Icges(1,4)) * V_base(5) + (-t188 * t163 + t189 * t165 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t189 * t162 + t188 * t164 + Icges(1,2)) * V_base(5) + (t163 * t189 + t165 * t188 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t260 * t191 + t262 * t193) * t181 + (t261 * t191 + t263 * t193) * t180 + (t258 * t191 - t259 * t193 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t189 - Icges(2,6) * t188 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t188 + Icges(2,6) * t189 + Icges(1,6));
T = t1;
