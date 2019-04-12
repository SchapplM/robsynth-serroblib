% Calculate kinetic energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR14V3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:18
% EndTime: 2019-04-12 15:03:23
% DurationCPUTime: 5.09s
% Computational Cost: add. (1419->322), mult. (3126->489), div. (0->0), fcn. (3548->10), ass. (0->145)
t265 = Icges(3,4) - Icges(4,5);
t264 = Icges(3,1) + Icges(4,1);
t263 = Icges(3,2) + Icges(4,3);
t204 = sin(qJ(2));
t262 = t265 * t204;
t208 = cos(qJ(2));
t261 = t265 * t208;
t260 = Icges(4,4) + Icges(3,5);
t259 = Icges(3,6) - Icges(4,6);
t258 = t263 * t204 - t261;
t257 = t264 * t208 - t262;
t256 = Icges(4,2) + Icges(3,3);
t205 = sin(qJ(1));
t209 = cos(qJ(1));
t255 = t258 * t205 + t259 * t209;
t254 = -t259 * t205 + t258 * t209;
t253 = t257 * t205 - t260 * t209;
t252 = t260 * t205 + t257 * t209;
t251 = -t263 * t208 - t262;
t250 = t264 * t204 + t261;
t249 = -t259 * t204 + t260 * t208;
t190 = -qJD(2) * t209 + V_base(5);
t191 = qJD(2) * t205 + V_base(4);
t198 = V_base(6) + qJD(1);
t248 = (t204 * t251 + t208 * t250) * t198 + (t204 * t254 + t208 * t252) * t191 + (t204 * t255 + t208 * t253) * t190;
t247 = (t260 * t204 + t259 * t208) * t198 + (t205 * t256 + t249 * t209) * t191 + (t249 * t205 - t209 * t256) * t190;
t243 = cos(qJ(5));
t242 = Icges(2,4) * t205;
t237 = qJ(3) * t191;
t203 = sin(qJ(4));
t236 = t203 * t204;
t235 = t204 * t205;
t234 = t204 * t209;
t233 = t205 * t208;
t232 = t208 * t209;
t231 = qJD(3) * t204;
t230 = qJD(4) * t204;
t227 = qJ(3) * t234;
t226 = t204 * t243;
t158 = t209 * t230 + t191;
t225 = t198 * t227 + t205 * t231 + t208 * t237 + V_base(2);
t224 = rSges(3,1) * t208 - rSges(3,2) * t204;
t223 = rSges(4,1) * t208 + rSges(4,3) * t204;
t207 = cos(qJ(4));
t163 = t203 * t232 - t205 * t207;
t124 = qJD(5) * t163 + t158;
t216 = -qJD(3) * t208 + t235 * t237 + V_base(3);
t157 = t205 * t230 + t190;
t183 = -qJD(4) * t208 + t198;
t161 = t203 * t233 + t207 * t209;
t123 = qJD(5) * t161 + t157;
t156 = qJD(5) * t236 + t183;
t213 = -t190 * t227 + t216;
t212 = t209 * t231 + V_base(1) + (-t190 * t208 - t198 * t235) * qJ(3);
t206 = cos(qJ(6));
t202 = sin(qJ(5));
t201 = sin(qJ(6));
t200 = Icges(2,4) * t209;
t187 = rSges(2,1) * t209 - t205 * rSges(2,2);
t186 = t205 * rSges(2,1) + rSges(2,2) * t209;
t185 = rSges(3,1) * t204 + rSges(3,2) * t208;
t184 = rSges(4,1) * t204 - rSges(4,3) * t208;
t182 = Icges(2,1) * t209 - t242;
t181 = Icges(2,1) * t205 + t200;
t178 = -Icges(2,2) * t205 + t200;
t177 = Icges(2,2) * t209 + t242;
t170 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t169 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t168 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t164 = t205 * t203 + t207 * t232;
t162 = -t203 * t209 + t207 * t233;
t160 = -t208 * t202 + t207 * t226;
t159 = t204 * t207 * t202 + t208 * t243;
t153 = t205 * rSges(3,3) + t209 * t224;
t152 = t205 * rSges(4,2) + t209 * t223;
t151 = -rSges(3,3) * t209 + t205 * t224;
t150 = -rSges(4,2) * t209 + t205 * t223;
t149 = -rSges(5,3) * t208 + (rSges(5,1) * t207 - rSges(5,2) * t203) * t204;
t144 = -Icges(5,5) * t208 + (Icges(5,1) * t207 - Icges(5,4) * t203) * t204;
t139 = -Icges(5,6) * t208 + (Icges(5,4) * t207 - Icges(5,2) * t203) * t204;
t134 = -Icges(5,3) * t208 + (Icges(5,5) * t207 - Icges(5,6) * t203) * t204;
t133 = V_base(5) * rSges(2,3) - t186 * t198 + V_base(1);
t132 = -V_base(4) * rSges(2,3) + t187 * t198 + V_base(2);
t131 = t164 * t243 + t202 * t234;
t130 = t164 * t202 - t209 * t226;
t129 = t162 * t243 + t202 * t235;
t128 = t162 * t202 - t205 * t226;
t127 = t160 * t206 + t201 * t236;
t126 = -t160 * t201 + t206 * t236;
t125 = t186 * V_base(4) - t187 * V_base(5) + V_base(3);
t122 = qJD(6) * t159 + t156;
t121 = t164 * rSges(5,1) - t163 * rSges(5,2) + rSges(5,3) * t234;
t120 = rSges(5,1) * t162 - rSges(5,2) * t161 + rSges(5,3) * t235;
t119 = rSges(6,1) * t160 - rSges(6,2) * t159 + rSges(6,3) * t236;
t118 = Icges(5,1) * t164 - Icges(5,4) * t163 + Icges(5,5) * t234;
t117 = Icges(5,1) * t162 - Icges(5,4) * t161 + Icges(5,5) * t235;
t116 = Icges(6,1) * t160 - Icges(6,4) * t159 + Icges(6,5) * t236;
t115 = Icges(5,4) * t164 - Icges(5,2) * t163 + Icges(5,6) * t234;
t114 = Icges(5,4) * t162 - Icges(5,2) * t161 + Icges(5,6) * t235;
t113 = Icges(6,4) * t160 - Icges(6,2) * t159 + Icges(6,6) * t236;
t112 = Icges(5,5) * t164 - Icges(5,6) * t163 + Icges(5,3) * t234;
t111 = Icges(5,5) * t162 - Icges(5,6) * t161 + Icges(5,3) * t235;
t110 = Icges(6,5) * t160 - Icges(6,6) * t159 + Icges(6,3) * t236;
t109 = -t151 * t198 + t185 * t190 + V_base(1);
t108 = t153 * t198 - t185 * t191 + V_base(2);
t107 = t131 * t206 + t163 * t201;
t106 = -t131 * t201 + t163 * t206;
t105 = t129 * t206 + t161 * t201;
t104 = -t129 * t201 + t161 * t206;
t103 = t151 * t191 - t153 * t190 + V_base(3);
t102 = qJD(6) * t130 + t124;
t101 = qJD(6) * t128 + t123;
t100 = rSges(6,1) * t131 - rSges(6,2) * t130 + rSges(6,3) * t163;
t99 = rSges(6,1) * t129 - rSges(6,2) * t128 + rSges(6,3) * t161;
t98 = rSges(7,1) * t127 + rSges(7,2) * t126 + rSges(7,3) * t159;
t97 = -t150 * t198 + t184 * t190 + t212;
t96 = t152 * t198 - t184 * t191 + t225;
t95 = Icges(6,1) * t131 - Icges(6,4) * t130 + Icges(6,5) * t163;
t94 = Icges(6,1) * t129 - Icges(6,4) * t128 + Icges(6,5) * t161;
t93 = Icges(7,1) * t127 + Icges(7,4) * t126 + Icges(7,5) * t159;
t92 = Icges(6,4) * t131 - Icges(6,2) * t130 + Icges(6,6) * t163;
t91 = Icges(6,4) * t129 - Icges(6,2) * t128 + Icges(6,6) * t161;
t90 = Icges(7,4) * t127 + Icges(7,2) * t126 + Icges(7,6) * t159;
t89 = Icges(6,5) * t131 - Icges(6,6) * t130 + Icges(6,3) * t163;
t88 = Icges(6,5) * t129 - Icges(6,6) * t128 + Icges(6,3) * t161;
t87 = Icges(7,5) * t127 + Icges(7,6) * t126 + Icges(7,3) * t159;
t86 = t191 * t150 + (-t152 - t227) * t190 + t216;
t85 = -t120 * t183 + t149 * t157 + t212;
t84 = t121 * t183 - t149 * t158 + t225;
t83 = rSges(7,1) * t107 + rSges(7,2) * t106 + rSges(7,3) * t130;
t82 = rSges(7,1) * t105 + rSges(7,2) * t104 + rSges(7,3) * t128;
t81 = Icges(7,1) * t107 + Icges(7,4) * t106 + Icges(7,5) * t130;
t80 = Icges(7,1) * t105 + Icges(7,4) * t104 + Icges(7,5) * t128;
t79 = Icges(7,4) * t107 + Icges(7,2) * t106 + Icges(7,6) * t130;
t78 = Icges(7,4) * t105 + Icges(7,2) * t104 + Icges(7,6) * t128;
t77 = Icges(7,5) * t107 + Icges(7,6) * t106 + Icges(7,3) * t130;
t76 = Icges(7,5) * t105 + Icges(7,6) * t104 + Icges(7,3) * t128;
t75 = t158 * t120 - t157 * t121 + t213;
t74 = t119 * t123 - t156 * t99 + t212;
t73 = t100 * t156 - t119 * t124 + t225;
t72 = -t123 * t100 + t124 * t99 + t213;
t71 = t101 * t98 - t122 * t82 + t212;
t70 = -t102 * t98 + t122 * t83 + t225;
t69 = -t101 * t83 + t102 * t82 + t213;
t1 = m(5) * (t75 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + t101 * ((t104 * t79 + t105 * t81 + t128 * t77) * t102 + (t104 * t78 + t105 * t80 + t128 * t76) * t101 + (t104 * t90 + t105 * t93 + t128 * t87) * t122) / 0.2e1 + t102 * ((t106 * t79 + t107 * t81 + t130 * t77) * t102 + (t106 * t78 + t107 * t80 + t130 * t76) * t101 + (t106 * t90 + t107 * t93 + t130 * t87) * t122) / 0.2e1 + m(3) * (t103 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(4) * (t86 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + t158 * ((t112 * t234 - t163 * t115 + t164 * t118) * t158 + (t111 * t234 - t163 * t114 + t164 * t117) * t157 + (t134 * t234 - t163 * t139 + t164 * t144) * t183) / 0.2e1 + t157 * ((t112 * t235 - t115 * t161 + t118 * t162) * t158 + (t111 * t235 - t114 * t161 + t117 * t162) * t157 + (t134 * t235 - t139 * t161 + t144 * t162) * t183) / 0.2e1 + t156 * ((-t159 * t92 + t160 * t95 + t236 * t89) * t124 + (-t159 * t91 + t160 * t94 + t236 * t88) * t123 + (t110 * t236 - t113 * t159 + t116 * t160) * t156) / 0.2e1 + t183 * ((-t111 * t157 - t112 * t158 - t134 * t183) * t208 + ((-t115 * t203 + t118 * t207) * t158 + (-t114 * t203 + t117 * t207) * t157 + (-t139 * t203 + t144 * t207) * t183) * t204) / 0.2e1 + t124 * ((-t130 * t92 + t131 * t95 + t163 * t89) * t124 + (-t130 * t91 + t131 * t94 + t163 * t88) * t123 + (t110 * t163 - t113 * t130 + t116 * t131) * t156) / 0.2e1 + m(1) * (t168 ^ 2 + t169 ^ 2 + t170 ^ 2) / 0.2e1 + t123 * ((-t128 * t92 + t129 * t95 + t161 * t89) * t124 + (-t128 * t91 + t129 * t94 + t161 * t88) * t123 + (t110 * t161 - t113 * t128 + t116 * t129) * t156) / 0.2e1 + t122 * ((t126 * t79 + t127 * t81 + t159 * t77) * t102 + (t126 * t78 + t127 * t80 + t159 * t76) * t101 + (t126 * t90 + t127 * t93 + t159 * t87) * t122) / 0.2e1 + m(2) * (t125 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + (t248 * t205 - t247 * t209) * t190 / 0.2e1 + (t247 * t205 + t248 * t209) * t191 / 0.2e1 + ((-t205 * t177 + t181 * t209 + Icges(1,4)) * V_base(5) + (-t205 * t178 + t182 * t209 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t177 * t209 + t205 * t181 + Icges(1,2)) * V_base(5) + (t178 * t209 + t205 * t182 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t204 * t252 - t208 * t254) * t191 + (t204 * t253 - t208 * t255) * t190 + (t204 * t250 - t208 * t251 + Icges(2,3)) * t198) * t198 / 0.2e1 + t198 * V_base(4) * (Icges(2,5) * t209 - Icges(2,6) * t205) + t198 * V_base(5) * (Icges(2,5) * t205 + Icges(2,6) * t209) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
