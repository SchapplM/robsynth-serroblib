% Calculate kinetic energy for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:33
% EndTime: 2019-03-09 01:50:35
% DurationCPUTime: 2.71s
% Computational Cost: add. (739->260), mult. (1284->342), div. (0->0), fcn. (1066->6), ass. (0->130)
t282 = -Icges(5,4) - Icges(6,6);
t281 = Icges(5,1) + Icges(6,2);
t280 = Icges(5,2) + Icges(6,3);
t195 = cos(qJ(4));
t279 = t282 * t195;
t192 = sin(qJ(4));
t278 = t282 * t192;
t277 = Icges(6,4) - Icges(5,5);
t276 = Icges(6,5) - Icges(5,6);
t275 = t280 * t195 - t278;
t274 = t281 * t192 - t279;
t273 = Icges(2,4) + Icges(3,6) - Icges(4,6);
t272 = Icges(6,1) + Icges(5,3);
t193 = sin(qJ(1));
t196 = cos(qJ(1));
t271 = t275 * t193 - t276 * t196;
t270 = t276 * t193 + t275 * t196;
t269 = t274 * t193 - t277 * t196;
t268 = t277 * t193 + t274 * t196;
t267 = t280 * t192 + t279;
t266 = t281 * t195 + t278;
t265 = -t277 * t192 - t276 * t195;
t264 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t263 = -Icges(3,4) + Icges(2,5) + Icges(4,5);
t262 = Icges(4,4) + Icges(3,5) - Icges(2,6);
t261 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t260 = t273 * t193;
t259 = t273 * t196;
t175 = qJD(4) * t196 + V_base(5);
t176 = -qJD(4) * t193 + V_base(4);
t180 = V_base(6) + qJD(1);
t258 = (t266 * t192 - t267 * t195) * t180 + (t268 * t192 + t270 * t195) * t176 + (t269 * t192 + t271 * t195) * t175;
t257 = (t276 * t192 - t277 * t195) * t180 + (-t272 * t193 + t265 * t196) * t176 + (t265 * t193 + t272 * t196) * t175;
t256 = t264 * t196 - t260;
t255 = t264 * t193 + t259;
t254 = -t261 * t196 - t260;
t253 = t261 * t193 - t259;
t250 = -pkin(2) - pkin(6);
t246 = pkin(7) * t193;
t245 = pkin(7) * t196;
t236 = qJ(3) * t193;
t235 = qJ(3) * t196;
t234 = t192 * t193;
t233 = t192 * t196;
t232 = t193 * t195;
t231 = t195 * t196;
t230 = qJD(5) * t196;
t229 = qJD(6) * t192;
t163 = t193 * pkin(1) - qJ(2) * t196;
t228 = V_base(4) * t163 + V_base(3);
t227 = V_base(5) * pkin(6) + V_base(1);
t224 = V_base(4) * t236 + t228;
t223 = qJD(2) * t193 + t227;
t222 = -t163 - t236;
t170 = pkin(1) * t196 + t193 * qJ(2);
t221 = -t170 - t235;
t220 = rSges(5,1) * t192 + rSges(5,2) * t195;
t219 = -rSges(6,2) * t192 - rSges(6,3) * t195;
t218 = pkin(4) * t192 - qJ(5) * t195;
t211 = -qJD(2) * t196 + t180 * t170 + V_base(2);
t210 = V_base(5) * pkin(2) + qJD(3) * t196 + t223;
t209 = t222 - t245;
t208 = V_base(5) * pkin(3) + t210;
t130 = t218 * t193;
t207 = -t130 + t209;
t206 = qJD(3) * t193 + t180 * t235 + t211;
t167 = pkin(4) * t195 + qJ(5) * t192;
t205 = t175 * t167 + t208;
t202 = V_base(4) * t245 + t224 + (t221 + t246) * V_base(5);
t201 = (-pkin(3) + t250) * V_base(4) + t206;
t200 = qJD(5) * t192 + t176 * t130 + t202;
t131 = t218 * t196;
t199 = -qJD(5) * t232 + t180 * t131 + t201;
t194 = cos(qJ(6));
t191 = sin(qJ(6));
t173 = rSges(2,1) * t196 - t193 * rSges(2,2);
t172 = -rSges(3,2) * t196 + t193 * rSges(3,3);
t171 = -rSges(4,2) * t196 + t193 * rSges(4,3);
t169 = rSges(5,1) * t195 - rSges(5,2) * t192;
t168 = -rSges(6,2) * t195 + rSges(6,3) * t192;
t166 = t193 * rSges(2,1) + rSges(2,2) * t196;
t165 = -t193 * rSges(3,2) - rSges(3,3) * t196;
t164 = t193 * rSges(4,2) + rSges(4,3) * t196;
t162 = qJD(6) * t195 + t180;
t137 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t136 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t135 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t134 = pkin(5) * t196 + pkin(8) * t234;
t133 = -t193 * pkin(5) + pkin(8) * t233;
t128 = -t191 * t232 + t194 * t196;
t127 = -t191 * t196 - t194 * t232;
t126 = -t191 * t231 - t193 * t194;
t125 = t193 * t191 - t194 * t231;
t124 = t196 * t229 + t176;
t123 = t193 * t229 + t175;
t121 = rSges(6,1) * t196 + t193 * t219;
t120 = -t193 * rSges(6,1) + t196 * t219;
t119 = -t193 * rSges(5,3) + t196 * t220;
t118 = rSges(5,3) * t196 + t193 * t220;
t117 = rSges(7,3) * t195 + (rSges(7,1) * t191 + rSges(7,2) * t194) * t192;
t108 = Icges(7,5) * t195 + (Icges(7,1) * t191 + Icges(7,4) * t194) * t192;
t105 = Icges(7,6) * t195 + (Icges(7,4) * t191 + Icges(7,2) * t194) * t192;
t102 = Icges(7,3) * t195 + (Icges(7,5) * t191 + Icges(7,6) * t194) * t192;
t99 = V_base(5) * rSges(2,3) - t166 * t180 + t227;
t98 = t173 * t180 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t97 = t166 * V_base(4) - t173 * V_base(5) + V_base(3);
t96 = rSges(7,1) * t128 + rSges(7,2) * t127 + rSges(7,3) * t234;
t95 = t126 * rSges(7,1) + t125 * rSges(7,2) + rSges(7,3) * t233;
t94 = Icges(7,1) * t128 + Icges(7,4) * t127 + Icges(7,5) * t234;
t93 = Icges(7,1) * t126 + Icges(7,4) * t125 + Icges(7,5) * t233;
t92 = Icges(7,4) * t128 + Icges(7,2) * t127 + Icges(7,6) * t234;
t91 = Icges(7,4) * t126 + Icges(7,2) * t125 + Icges(7,6) * t233;
t90 = Icges(7,5) * t128 + Icges(7,6) * t127 + Icges(7,3) * t234;
t89 = Icges(7,5) * t126 + Icges(7,6) * t125 + Icges(7,3) * t233;
t88 = V_base(5) * rSges(3,1) + (-t163 - t165) * t180 + t223;
t87 = t180 * t172 + (-rSges(3,1) - pkin(6)) * V_base(4) + t211;
t86 = t165 * V_base(4) + (-t170 - t172) * V_base(5) + t228;
t85 = V_base(5) * rSges(4,1) + (-t171 + t222) * t180 + t210;
t84 = t180 * t164 + (-rSges(4,1) + t250) * V_base(4) + t206;
t83 = V_base(4) * t171 + (-t164 + t221) * V_base(5) + t224;
t82 = t175 * t169 + (-t118 + t209) * t180 + t208;
t81 = -t176 * t169 + (t119 - t246) * t180 + t201;
t80 = t176 * t118 - t175 * t119 + t202;
t79 = -t195 * t230 + t175 * t168 + (-t121 + t207) * t180 + t205;
t78 = (t120 - t246) * t180 + (-t167 - t168) * t176 + t199;
t77 = t176 * t121 + (-t120 - t131) * t175 + t200;
t76 = t123 * t117 - t162 * t96 + (pkin(8) * t175 - t230) * t195 + (-t134 + t207) * t180 + t205;
t75 = -t124 * t117 + t162 * t95 + (t133 - t246) * t180 + (-pkin(8) * t195 - t167) * t176 + t199;
t74 = -t123 * t95 + t124 * t96 + t176 * t134 + (-t131 - t133) * t175 + t200;
t1 = t124 * ((t125 * t91 + t126 * t93 + t89 * t233) * t124 + (t125 * t92 + t126 * t94 + t233 * t90) * t123 + (t102 * t233 + t125 * t105 + t126 * t108) * t162) / 0.2e1 + t123 * ((t127 * t91 + t128 * t93 + t234 * t89) * t124 + (t127 * t92 + t128 * t94 + t90 * t234) * t123 + (t102 * t234 + t105 * t127 + t108 * t128) * t162) / 0.2e1 + t162 * ((t102 * t162 + t90 * t123 + t89 * t124) * t195 + ((t191 * t93 + t194 * t91) * t124 + (t191 * t94 + t194 * t92) * t123 + (t105 * t194 + t108 * t191) * t162) * t192) / 0.2e1 + m(1) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(2) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(3) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(7) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(5) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + (t258 * t193 + t257 * t196) * t175 / 0.2e1 + (-t257 * t193 + t258 * t196) * t176 / 0.2e1 + ((t193 * t254 + t196 * t255 + Icges(1,4)) * V_base(5) + (t253 * t193 + t256 * t196 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t255 * t193 - t254 * t196 + Icges(1,2)) * V_base(5) + (t193 * t256 - t196 * t253 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t270 * t192 + t268 * t195) * t176 + (-t271 * t192 + t269 * t195) * t175 + (t267 * t192 + t266 * t195 + Icges(3,1) + Icges(4,1) + Icges(2,3)) * t180) * t180 / 0.2e1 + t180 * V_base(4) * (t262 * t193 + t263 * t196) + t180 * V_base(5) * (t263 * t193 - t262 * t196) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
