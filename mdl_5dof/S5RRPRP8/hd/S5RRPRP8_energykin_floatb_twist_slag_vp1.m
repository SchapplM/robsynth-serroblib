% Calculate kinetic energy for
% S5RRPRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:16
% EndTime: 2019-12-31 20:03:19
% DurationCPUTime: 2.46s
% Computational Cost: add. (838->226), mult. (1585->304), div. (0->0), fcn. (1523->6), ass. (0->119)
t269 = Icges(3,4) - Icges(4,5);
t268 = Icges(3,1) + Icges(4,1);
t267 = Icges(3,2) + Icges(4,3);
t185 = cos(qJ(2));
t266 = t269 * t185;
t182 = sin(qJ(2));
t265 = t269 * t182;
t264 = Icges(4,4) + Icges(3,5);
t263 = Icges(3,6) - Icges(4,6);
t262 = t267 * t182 - t266;
t261 = t268 * t185 - t265;
t260 = Icges(5,1) + Icges(6,1);
t259 = Icges(5,4) + Icges(6,4);
t258 = Icges(5,5) + Icges(6,5);
t257 = Icges(4,2) + Icges(3,3);
t256 = Icges(5,2) + Icges(6,2);
t255 = Icges(5,6) + Icges(6,6);
t254 = Icges(5,3) + Icges(6,3);
t183 = sin(qJ(1));
t186 = cos(qJ(1));
t253 = t262 * t183 + t263 * t186;
t252 = -t263 * t183 + t262 * t186;
t251 = t261 * t183 - t264 * t186;
t250 = t264 * t183 + t261 * t186;
t249 = -t267 * t185 - t265;
t248 = t268 * t182 + t266;
t247 = -t263 * t182 + t264 * t185;
t246 = rSges(6,3) + qJ(5);
t184 = cos(qJ(4));
t181 = sin(qJ(4));
t217 = t181 * t185;
t144 = t182 * t184 - t217;
t133 = t144 * t183;
t218 = t181 * t182;
t199 = t184 * t185 + t218;
t134 = t199 * t183;
t245 = t255 * t133 + t258 * t134 + t254 * t186;
t135 = t144 * t186;
t136 = t199 * t186;
t244 = t255 * t135 + t258 * t136 - t254 * t183;
t243 = t256 * t133 + t259 * t134 + t255 * t186;
t242 = t256 * t135 + t259 * t136 - t255 * t183;
t241 = t259 * t133 + t260 * t134 + t258 * t186;
t240 = t259 * t135 + t260 * t136 - t258 * t183;
t239 = t258 * t144 - t255 * t199;
t238 = t259 * t144 - t256 * t199;
t237 = t260 * t144 - t259 * t199;
t171 = -qJD(2) * t186 + V_base(5);
t172 = qJD(2) * t183 + V_base(4);
t176 = V_base(6) + qJD(1);
t236 = (t249 * t182 + t248 * t185) * t176 + (t252 * t182 + t250 * t185) * t172 + (t253 * t182 + t251 * t185) * t171;
t235 = (t264 * t182 + t263 * t185) * t176 + (t257 * t183 + t247 * t186) * t172 + (t247 * t183 - t257 * t186) * t171;
t229 = pkin(3) * t182;
t228 = pkin(3) * t185;
t227 = pkin(4) * t184;
t195 = pkin(4) * t218 + t185 * t227;
t225 = t134 * rSges(6,1) + t133 * rSges(6,2) + t183 * t195 + t186 * t246;
t224 = rSges(6,1) * t136 + rSges(6,2) * t135 - t183 * t246 + t186 * t195;
t223 = Icges(2,4) * t183;
t216 = rSges(6,1) * t144 - rSges(6,2) * t199 - pkin(4) * t217 + t182 * t227;
t206 = pkin(2) * t185 + qJ(3) * t182;
t138 = t206 * t183;
t169 = t183 * pkin(1) - pkin(6) * t186;
t215 = -t138 - t169;
t214 = qJD(3) * t182;
t213 = V_base(5) * pkin(5) + V_base(1);
t147 = pkin(7) * t186 + t183 * t228;
t210 = -t147 + t215;
t164 = pkin(2) * t182 - qJ(3) * t185;
t209 = t171 * t164 + t186 * t214 + t213;
t208 = rSges(3,1) * t185 - rSges(3,2) * t182;
t207 = rSges(4,1) * t185 + rSges(4,3) * t182;
t198 = t171 * t229 + t209;
t170 = pkin(1) * t186 + t183 * pkin(6);
t197 = -V_base(4) * pkin(5) + t176 * t170 + V_base(2);
t196 = V_base(4) * t169 - t170 * V_base(5) + V_base(3);
t139 = t206 * t186;
t192 = t176 * t139 + t183 * t214 + t197;
t191 = -qJD(3) * t185 + t172 * t138 + t196;
t148 = -t183 * pkin(7) + t186 * t228;
t190 = t176 * t148 + (-t164 - t229) * t172 + t192;
t189 = t172 * t147 + (-t139 - t148) * t171 + t191;
t178 = Icges(2,4) * t186;
t168 = rSges(2,1) * t186 - t183 * rSges(2,2);
t167 = t183 * rSges(2,1) + rSges(2,2) * t186;
t166 = rSges(3,1) * t182 + rSges(3,2) * t185;
t165 = rSges(4,1) * t182 - rSges(4,3) * t185;
t163 = Icges(2,1) * t186 - t223;
t162 = Icges(2,1) * t183 + t178;
t159 = -Icges(2,2) * t183 + t178;
t158 = Icges(2,2) * t186 + t223;
t151 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t150 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t149 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t142 = -qJD(4) * t183 + t172;
t141 = V_base(5) + (-qJD(2) + qJD(4)) * t186;
t131 = t183 * rSges(3,3) + t186 * t208;
t130 = t183 * rSges(4,2) + t186 * t207;
t129 = -rSges(3,3) * t186 + t183 * t208;
t128 = -rSges(4,2) * t186 + t183 * t207;
t111 = V_base(5) * rSges(2,3) - t167 * t176 + t213;
t110 = t168 * t176 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t109 = t167 * V_base(4) - t168 * V_base(5) + V_base(3);
t108 = rSges(5,1) * t144 - rSges(5,2) * t199;
t98 = rSges(5,1) * t136 + rSges(5,2) * t135 - rSges(5,3) * t183;
t96 = t134 * rSges(5,1) + t133 * rSges(5,2) + rSges(5,3) * t186;
t82 = t166 * t171 + (-t129 - t169) * t176 + t213;
t81 = t131 * t176 - t166 * t172 + t197;
t80 = t129 * t172 - t131 * t171 + t196;
t79 = t165 * t171 + (-t128 + t215) * t176 + t209;
t78 = t130 * t176 + (-t164 - t165) * t172 + t192;
t77 = t128 * t172 + (-t130 - t139) * t171 + t191;
t76 = t108 * t141 + (-t96 + t210) * t176 + t198;
t75 = -t108 * t142 + t176 * t98 + t190;
t74 = -t141 * t98 + t142 * t96 + t189;
t73 = -qJD(5) * t183 + t216 * t141 + (t210 - t225) * t176 + t198;
t72 = qJD(5) * t186 - t142 * t216 + t176 * t224 + t190;
t71 = -t141 * t224 + t142 * t225 + t189;
t1 = m(1) * (t149 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(2) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(3) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(4) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(5) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + ((t133 * t238 + t134 * t237 + t186 * t239) * t176 + (t133 * t242 + t134 * t240 + t186 * t244) * t142 + (t243 * t133 + t241 * t134 + t245 * t186) * t141) * t141 / 0.2e1 + ((t135 * t238 + t136 * t237 - t183 * t239) * t176 + (t242 * t135 + t240 * t136 - t244 * t183) * t142 + (t243 * t135 + t241 * t136 - t245 * t183) * t141) * t142 / 0.2e1 + (t236 * t183 - t235 * t186) * t171 / 0.2e1 + (t235 * t183 + t236 * t186) * t172 / 0.2e1 + ((-t183 * t158 + t162 * t186 + Icges(1,4)) * V_base(5) + (-t183 * t159 + t163 * t186 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t158 * t186 + t183 * t162 + Icges(1,2)) * V_base(5) + (t159 * t186 + t183 * t163 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t250 * t182 - t252 * t185) * t172 + (t251 * t182 - t253 * t185) * t171 + (t144 * t240 - t199 * t242) * t142 + (t144 * t241 - t199 * t243) * t141 + (t237 * t144 + t248 * t182 - t249 * t185 - t238 * t199 + Icges(2,3)) * t176) * t176 / 0.2e1 + V_base(4) * t176 * (Icges(2,5) * t186 - Icges(2,6) * t183) + V_base(5) * t176 * (Icges(2,5) * t183 + Icges(2,6) * t186) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
