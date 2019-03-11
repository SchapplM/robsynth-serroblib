% Calculate kinetic energy for
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:10
% EndTime: 2019-03-09 01:37:12
% DurationCPUTime: 2.27s
% Computational Cost: add. (987->265), mult. (1760->355), div. (0->0), fcn. (1868->8), ass. (0->131)
t256 = Icges(2,4) - Icges(4,5) + Icges(3,6);
t255 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t254 = Icges(4,1) + Icges(2,2) + Icges(3,3);
t253 = -Icges(3,4) + Icges(2,5) - Icges(4,6);
t252 = Icges(4,4) - Icges(3,5) + Icges(2,6);
t192 = sin(qJ(1));
t251 = t256 * t192;
t195 = cos(qJ(1));
t250 = t256 * t195;
t249 = -t254 * t195 - t251;
t248 = t254 * t192 - t250;
t247 = t255 * t192 + t250;
t246 = t255 * t195 - t251;
t189 = cos(pkin(9));
t238 = sin(pkin(9));
t138 = t192 * t189 + t195 * t238;
t131 = -qJD(5) * t138 + V_base(5);
t137 = t189 * t195 - t192 * t238;
t132 = -qJD(5) * t137 + V_base(4);
t194 = cos(qJ(5));
t191 = sin(qJ(5));
t235 = Icges(6,4) * t191;
t152 = Icges(6,2) * t194 + t235;
t234 = Icges(6,4) * t194;
t159 = Icges(6,1) * t191 + t234;
t180 = V_base(6) + qJD(1);
t207 = -Icges(6,2) * t191 + t234;
t96 = -Icges(6,6) * t137 + t138 * t207;
t97 = -Icges(6,6) * t138 - t137 * t207;
t208 = Icges(6,1) * t194 - t235;
t98 = -Icges(6,5) * t137 + t138 * t208;
t99 = -Icges(6,5) * t138 - t137 * t208;
t243 = (t191 * t96 - t194 * t98) * t132 + (t191 * t97 - t194 * t99) * t131 + (t152 * t191 - t159 * t194) * t180;
t242 = -pkin(2) - pkin(6);
t240 = pkin(3) * t192;
t239 = pkin(3) * t195;
t236 = Icges(5,4) * t137;
t230 = qJ(3) * t192;
t229 = t137 * t191;
t228 = t138 * t191;
t227 = t180 * t195;
t190 = sin(qJ(6));
t226 = t190 * t194;
t193 = cos(qJ(6));
t225 = t193 * t194;
t224 = qJD(6) * t191;
t223 = -qJ(4) + t242;
t166 = t192 * pkin(1) - qJ(2) * t195;
t222 = V_base(4) * t166 + V_base(3);
t221 = V_base(5) * pkin(6) + V_base(1);
t218 = V_base(4) * t230 + t222;
t217 = qJD(2) * t192 + t221;
t216 = -t166 - t230;
t170 = pkin(1) * t195 + t192 * qJ(2);
t215 = -qJ(3) * t195 - t170;
t214 = pkin(5) * t194 + pkin(8) * t191;
t123 = -pkin(4) * t137 - pkin(7) * t138;
t213 = -t123 + t216;
t212 = qJD(4) + t218;
t211 = rSges(6,1) * t194 - rSges(6,2) * t191;
t206 = Icges(6,5) * t194 - Icges(6,6) * t191;
t204 = -qJD(2) * t195 + t180 * t170 + V_base(2);
t203 = V_base(5) * pkin(2) + qJD(3) * t195 + t217;
t202 = t215 - t240;
t201 = qJ(3) * t227 + qJD(3) * t192 + t204;
t200 = -t131 * (-Icges(6,3) * t138 - t137 * t206) - t132 * (-Icges(6,3) * t137 + t138 * t206) - (Icges(6,5) * t191 + Icges(6,6) * t194) * t180;
t199 = pkin(3) * t227 + V_base(5) * qJ(4) + t203;
t198 = t180 * t240 + t201;
t122 = pkin(4) * t138 - pkin(7) * t137;
t197 = t180 * t122 + t223 * V_base(4) + t198;
t196 = (-t122 + t202) * V_base(5) + t212 + (t123 - t239) * V_base(4);
t175 = pkin(5) * t191 - pkin(8) * t194;
t173 = rSges(2,1) * t195 - t192 * rSges(2,2);
t172 = -rSges(4,1) * t195 + t192 * rSges(4,3);
t171 = -rSges(3,2) * t195 + t192 * rSges(3,3);
t169 = t192 * rSges(2,1) + rSges(2,2) * t195;
t168 = t192 * rSges(4,1) + rSges(4,3) * t195;
t167 = -t192 * rSges(3,2) - rSges(3,3) * t195;
t165 = rSges(6,1) * t191 + rSges(6,2) * t194;
t164 = -qJD(6) * t194 + t180;
t142 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t141 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t140 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t135 = Icges(5,4) * t138;
t130 = -rSges(7,3) * t194 + (rSges(7,1) * t193 - rSges(7,2) * t190) * t191;
t129 = -Icges(7,5) * t194 + (Icges(7,1) * t193 - Icges(7,4) * t190) * t191;
t128 = -Icges(7,6) * t194 + (Icges(7,4) * t193 - Icges(7,2) * t190) * t191;
t127 = -Icges(7,3) * t194 + (Icges(7,5) * t193 - Icges(7,6) * t190) * t191;
t126 = V_base(5) * rSges(2,3) - t169 * t180 + t221;
t125 = t173 * t180 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t124 = t169 * V_base(4) - t173 * V_base(5) + V_base(3);
t121 = -rSges(5,1) * t137 + rSges(5,2) * t138;
t120 = rSges(5,1) * t138 + rSges(5,2) * t137;
t119 = -Icges(5,1) * t137 + t135;
t118 = Icges(5,1) * t138 + t236;
t117 = Icges(5,2) * t138 - t236;
t116 = Icges(5,2) * t137 + t135;
t112 = t214 * t137;
t111 = t214 * t138;
t109 = -t137 * t225 - t138 * t190;
t108 = t137 * t226 - t138 * t193;
t107 = -t137 * t190 + t138 * t225;
t106 = -t137 * t193 - t138 * t226;
t105 = t138 * t224 + t132;
t104 = -t137 * t224 + t131;
t103 = V_base(5) * rSges(3,1) + (-t166 - t167) * t180 + t217;
t102 = t180 * t171 + (-rSges(3,1) - pkin(6)) * V_base(4) + t204;
t101 = -rSges(6,3) * t138 - t137 * t211;
t100 = -rSges(6,3) * t137 + t138 * t211;
t93 = t167 * V_base(4) + (-t170 - t171) * V_base(5) + t222;
t92 = -V_base(5) * rSges(4,2) + (-t172 + t216) * t180 + t203;
t91 = t180 * t168 + (rSges(4,2) + t242) * V_base(4) + t201;
t90 = V_base(4) * t172 + (-t168 + t215) * V_base(5) + t218;
t89 = rSges(7,1) * t109 + rSges(7,2) * t108 - rSges(7,3) * t229;
t88 = rSges(7,1) * t107 + rSges(7,2) * t106 + rSges(7,3) * t228;
t87 = Icges(7,1) * t109 + Icges(7,4) * t108 - Icges(7,5) * t229;
t86 = Icges(7,1) * t107 + Icges(7,4) * t106 + Icges(7,5) * t228;
t85 = Icges(7,4) * t109 + Icges(7,2) * t108 - Icges(7,6) * t229;
t84 = Icges(7,4) * t107 + Icges(7,2) * t106 + Icges(7,6) * t228;
t83 = Icges(7,5) * t109 + Icges(7,6) * t108 - Icges(7,3) * t229;
t82 = Icges(7,5) * t107 + Icges(7,6) * t106 + Icges(7,3) * t228;
t81 = V_base(5) * rSges(5,3) + (-t121 + t216) * t180 + t199;
t80 = t180 * t120 + (-rSges(5,3) + t223) * V_base(4) + t198;
t79 = (t121 - t239) * V_base(4) + (-t120 + t202) * V_base(5) + t212;
t78 = t131 * t165 + (-t101 + t213) * t180 + t199;
t77 = t180 * t100 - t132 * t165 + t197;
t76 = -t131 * t100 + t132 * t101 + t196;
t75 = t104 * t130 + t131 * t175 - t164 * t89 + (t112 + t213) * t180 + t199;
t74 = -t105 * t130 + t180 * t111 - t132 * t175 + t164 * t88 + t197;
t73 = -t104 * t88 + t105 * t89 - t131 * t111 - t132 * t112 + t196;
t1 = m(1) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(2) * (t124 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t102 ^ 2 + t103 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(5) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(7) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + t164 * ((-t104 * t83 - t105 * t82 - t127 * t164) * t194 + ((-t190 * t84 + t193 * t86) * t105 + (-t190 * t85 + t193 * t87) * t104 + (-t128 * t190 + t129 * t193) * t164) * t191) / 0.2e1 + t104 * ((t108 * t84 + t109 * t86 - t229 * t82) * t105 + (t108 * t85 + t109 * t87 - t83 * t229) * t104 + (t108 * t128 + t109 * t129 - t127 * t229) * t164) / 0.2e1 + t105 * ((t106 * t84 + t107 * t86 + t82 * t228) * t105 + (t106 * t85 + t107 * t87 + t228 * t83) * t104 + (t106 * t128 + t107 * t129 + t127 * t228) * t164) / 0.2e1 + t131 * (t243 * t137 + t200 * t138) / 0.2e1 + t132 * (t200 * t137 - t243 * t138) / 0.2e1 + ((t191 * t98 + t194 * t96) * t132 + (t191 * t99 + t194 * t97) * t131 + (t152 * t194 + t159 * t191 + Icges(3,1) + Icges(4,2) + Icges(2,3) + Icges(5,3)) * t180) * t180 / 0.2e1 + ((t117 * t137 + t119 * t138 + t192 * t249 + t247 * t195 + Icges(1,4)) * V_base(5) + (t116 * t137 + t118 * t138 + t192 * t248 + t195 * t246 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t117 * t138 - t119 * t137 + t247 * t192 - t195 * t249 + Icges(1,2)) * V_base(5) + (t116 * t138 - t118 * t137 + t192 * t246 - t195 * t248 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t180 * (-Icges(5,5) * t137 + Icges(5,6) * t138 + t253 * t192 + t252 * t195) + V_base(4) * t180 * (Icges(5,5) * t138 + Icges(5,6) * t137 - t252 * t192 + t253 * t195) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
