% Calculate kinetic energy for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:00
% EndTime: 2019-12-05 15:35:02
% DurationCPUTime: 2.31s
% Computational Cost: add. (1242->251), mult. (1606->345), div. (0->0), fcn. (1510->8), ass. (0->128)
t274 = Icges(5,1) + Icges(6,1);
t273 = -Icges(5,4) + Icges(6,5);
t272 = Icges(6,4) + Icges(5,5);
t271 = Icges(5,2) + Icges(6,3);
t270 = -Icges(6,6) + Icges(5,6);
t269 = Icges(3,3) + Icges(4,3);
t268 = -Icges(5,3) - Icges(6,2);
t194 = qJ(2) + pkin(8);
t190 = sin(t194);
t191 = cos(t194);
t199 = sin(qJ(2));
t201 = cos(qJ(2));
t267 = Icges(3,5) * t201 + Icges(4,5) * t191 - Icges(3,6) * t199 - Icges(4,6) * t190;
t266 = rSges(6,1) + pkin(4);
t265 = rSges(6,3) + qJ(5);
t196 = cos(pkin(7));
t200 = cos(qJ(4));
t233 = t196 * t200;
t195 = sin(pkin(7));
t198 = sin(qJ(4));
t236 = t195 * t198;
t154 = t191 * t236 + t233;
t234 = t196 * t198;
t235 = t195 * t200;
t155 = t191 * t235 - t234;
t238 = t190 * t195;
t264 = -t154 * t270 + t155 * t272 - t238 * t268;
t263 = t154 * t271 + t155 * t273 - t238 * t270;
t156 = t191 * t234 - t235;
t157 = t191 * t233 + t236;
t237 = t190 * t196;
t262 = t156 * t271 + t157 * t273 - t237 * t270;
t259 = -t156 * t270 + t157 * t272 - t237 * t268;
t258 = t273 * t154 + t155 * t274 + t272 * t238;
t257 = t273 * t156 + t157 * t274 + t272 * t237;
t256 = t270 * t191 + (t198 * t271 + t200 * t273) * t190;
t255 = t268 * t191 + (-t198 * t270 + t200 * t272) * t190;
t254 = -t272 * t191 + (t273 * t198 + t200 * t274) * t190;
t239 = Icges(4,4) * t191;
t215 = -Icges(4,2) * t190 + t239;
t131 = -Icges(4,6) * t196 + t195 * t215;
t132 = Icges(4,6) * t195 + t196 * t215;
t240 = Icges(4,4) * t190;
t217 = Icges(4,1) * t191 - t240;
t133 = -Icges(4,5) * t196 + t195 * t217;
t134 = Icges(4,5) * t195 + t196 * t217;
t241 = Icges(3,4) * t201;
t216 = -Icges(3,2) * t199 + t241;
t143 = -Icges(3,6) * t196 + t195 * t216;
t144 = Icges(3,6) * t195 + t196 * t216;
t242 = Icges(3,4) * t199;
t218 = Icges(3,1) * t201 - t242;
t145 = -Icges(3,5) * t196 + t195 * t218;
t146 = Icges(3,5) * t195 + t196 * t218;
t159 = Icges(4,2) * t191 + t240;
t160 = Icges(4,1) * t190 + t239;
t181 = Icges(3,2) * t201 + t242;
t182 = Icges(3,1) * t199 + t241;
t184 = -qJD(2) * t196 + V_base(5);
t185 = qJD(2) * t195 + V_base(4);
t253 = (-t159 * t190 + t160 * t191 - t181 * t199 + t182 * t201) * V_base(6) + (-t132 * t190 + t134 * t191 - t144 * t199 + t146 * t201) * t185 + (-t131 * t190 + t133 * t191 - t143 * t199 + t145 * t201) * t184;
t252 = (Icges(3,5) * t199 + Icges(4,5) * t190 + Icges(3,6) * t201 + Icges(4,6) * t191) * V_base(6) + (t195 * t269 + t196 * t267) * t185 + (t195 * t267 - t196 * t269) * t184;
t246 = pkin(2) * t199;
t245 = pkin(2) * t201;
t243 = Icges(2,4) * t195;
t232 = rSges(6,2) * t238 + t265 * t154 + t266 * t155;
t231 = rSges(6,2) * t237 + t265 * t156 + t266 * t157;
t230 = -rSges(6,2) * t191 + (t265 * t198 + t266 * t200) * t190;
t127 = -qJ(3) * t196 + t195 * t245;
t178 = pkin(1) * t195 - pkin(5) * t196;
t229 = -t127 - t178;
t228 = qJD(4) * t190;
t227 = V_base(5) * qJ(1) + V_base(1);
t223 = qJD(1) + V_base(3);
t222 = qJD(3) * t195 + t184 * t246 + t227;
t221 = pkin(3) * t191 + pkin(6) * t190;
t220 = rSges(3,1) * t201 - rSges(3,2) * t199;
t219 = rSges(4,1) * t191 - rSges(4,2) * t190;
t179 = pkin(1) * t196 + pkin(5) * t195;
t212 = -V_base(4) * qJ(1) + V_base(6) * t179 + V_base(2);
t211 = V_base(4) * t178 - t179 * V_base(5) + t223;
t210 = t185 * t127 + t211;
t128 = qJ(3) * t195 + t196 * t245;
t207 = -qJD(3) * t196 + V_base(6) * t128 + t212;
t149 = t221 * t195;
t162 = pkin(3) * t190 - pkin(6) * t191;
t206 = t184 * t162 + (-t149 + t229) * V_base(6) + t222;
t150 = t221 * t196;
t205 = t185 * t149 + (-t128 - t150) * t184 + t210;
t204 = V_base(6) * t150 + (-t162 - t246) * t185 + t207;
t192 = Icges(2,4) * t196;
t183 = t199 * rSges(3,1) + rSges(3,2) * t201;
t177 = rSges(2,1) * t196 - rSges(2,2) * t195;
t176 = rSges(2,1) * t195 + rSges(2,2) * t196;
t175 = -qJD(4) * t191 + V_base(6);
t174 = Icges(2,1) * t196 - t243;
t173 = Icges(2,1) * t195 + t192;
t172 = -Icges(2,2) * t195 + t192;
t171 = Icges(2,2) * t196 + t243;
t168 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t167 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t166 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t161 = rSges(4,1) * t190 + rSges(4,2) * t191;
t153 = t196 * t228 + t185;
t152 = t195 * t228 + t184;
t148 = t195 * rSges(3,3) + t196 * t220;
t147 = -t196 * rSges(3,3) + t195 * t220;
t138 = rSges(4,3) * t195 + t196 * t219;
t137 = -rSges(4,3) * t196 + t195 * t219;
t136 = V_base(5) * rSges(2,3) - t176 * V_base(6) + t227;
t135 = t177 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t126 = -rSges(5,3) * t191 + (rSges(5,1) * t200 - rSges(5,2) * t198) * t190;
t116 = t176 * V_base(4) - t177 * V_base(5) + t223;
t112 = rSges(5,1) * t157 - rSges(5,2) * t156 + rSges(5,3) * t237;
t110 = rSges(5,1) * t155 - rSges(5,2) * t154 + rSges(5,3) * t238;
t96 = t183 * t184 + (-t147 - t178) * V_base(6) + t227;
t95 = t148 * V_base(6) - t183 * t185 + t212;
t94 = t147 * t185 - t148 * t184 + t211;
t93 = t161 * t184 + (-t137 + t229) * V_base(6) + t222;
t92 = t138 * V_base(6) + (-t161 - t246) * t185 + t207;
t91 = t137 * t185 + (-t128 - t138) * t184 + t210;
t90 = -t110 * t175 + t126 * t152 + t206;
t89 = t112 * t175 - t126 * t153 + t204;
t88 = t110 * t153 - t112 * t152 + t205;
t87 = qJD(5) * t156 + t152 * t230 - t175 * t232 + t206;
t86 = qJD(5) * t154 - t153 * t230 + t175 * t231 + t204;
t85 = qJD(5) * t190 * t198 - t152 * t231 + t153 * t232 + t205;
t1 = m(1) * (t166 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + m(2) * (t116 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(3) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + ((t154 * t256 + t155 * t254 + t238 * t255) * t175 + (t154 * t262 + t155 * t257 + t238 * t259) * t153 + (t263 * t154 + t258 * t155 + t264 * t238) * t152) * t152 / 0.2e1 + ((t156 * t256 + t157 * t254 + t237 * t255) * t175 + (t262 * t156 + t257 * t157 + t259 * t237) * t153 + (t263 * t156 + t258 * t157 + t237 * t264) * t152) * t153 / 0.2e1 + ((-t152 * t264 - t259 * t153 - t255 * t175) * t191 + ((t198 * t256 + t200 * t254) * t175 + (t198 * t262 + t200 * t257) * t153 + (t198 * t263 + t200 * t258) * t152) * t190) * t175 / 0.2e1 + (t253 * t195 - t252 * t196) * t184 / 0.2e1 + (t252 * t195 + t253 * t196) * t185 / 0.2e1 + ((-t171 * t195 + t173 * t196 + Icges(1,4)) * V_base(5) + (-t195 * t172 + t196 * t174 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t196 * t171 + t195 * t173 + Icges(1,2)) * V_base(5) + (t172 * t196 + t174 * t195 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t132 * t191 + t134 * t190 + t144 * t201 + t199 * t146) * t185 + (t131 * t191 + t133 * t190 + t143 * t201 + t199 * t145) * t184 + (t191 * t159 + t190 * t160 + t201 * t181 + t199 * t182 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t196 - Icges(2,6) * t195 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t195 + Icges(2,6) * t196 + Icges(1,6));
T = t1;
