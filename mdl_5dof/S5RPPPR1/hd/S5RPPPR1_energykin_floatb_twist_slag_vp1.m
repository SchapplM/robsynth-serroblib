% Calculate kinetic energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:42
% EndTime: 2019-12-05 17:28:44
% DurationCPUTime: 2.09s
% Computational Cost: add. (1318->299), mult. (1336->395), div. (0->0), fcn. (1228->10), ass. (0->139)
t188 = qJ(1) + pkin(7);
t181 = sin(t188);
t183 = cos(t188);
t194 = sin(qJ(1));
t195 = cos(qJ(1));
t246 = -Icges(2,5) * t194 - Icges(3,5) * t181 - Icges(2,6) * t195 - Icges(3,6) * t183;
t245 = Icges(2,5) * t195 + Icges(3,5) * t183 - Icges(2,6) * t194 - Icges(3,6) * t181;
t190 = sin(pkin(8));
t192 = cos(pkin(8));
t231 = Icges(4,4) * t192;
t205 = -Icges(4,2) * t190 + t231;
t113 = Icges(4,6) * t183 - t181 * t205;
t232 = Icges(4,4) * t190;
t206 = Icges(4,1) * t192 - t232;
t115 = Icges(4,5) * t183 - t181 * t206;
t233 = Icges(3,4) * t183;
t244 = -Icges(3,1) * t181 - t113 * t190 + t115 * t192 - t233;
t114 = Icges(4,6) * t181 + t183 * t205;
t116 = Icges(4,5) * t181 + t183 * t206;
t234 = Icges(3,4) * t181;
t243 = Icges(3,1) * t183 - t114 * t190 + t116 * t192 - t234;
t191 = cos(pkin(9));
t238 = pkin(4) * t191;
t242 = pkin(6) * t190 + t192 * t238;
t240 = pkin(1) * t194;
t239 = pkin(1) * t195;
t237 = -pkin(5) - qJ(2);
t236 = Icges(2,4) * t194;
t235 = Icges(2,4) * t195;
t189 = sin(pkin(9));
t230 = t181 * t189;
t229 = t181 * t190;
t228 = t181 * t192;
t227 = t183 * t189;
t226 = t183 * t190;
t225 = t183 * t192;
t224 = t189 * t192;
t223 = t191 * t192;
t207 = pkin(3) * t192 + qJ(4) * t190;
t140 = t207 * t181;
t150 = -pkin(2) * t181 + qJ(3) * t183;
t221 = t140 - t150;
t220 = qJD(4) * t190;
t219 = qJD(5) * t190;
t218 = V_base(6) * pkin(5) + V_base(2);
t163 = pkin(3) * t190 - qJ(4) * t192;
t215 = -t163 + t237;
t184 = V_base(4) + qJD(1);
t152 = pkin(2) * t183 + qJ(3) * t181;
t214 = -t152 - t239;
t213 = qJD(3) * t181 + t150 * t184 + V_base(3);
t212 = V_base(6) * qJ(2) + t218;
t141 = t207 * t183;
t211 = -t141 + t214;
t210 = qJD(3) * t183 + t212;
t209 = V_base(5) * t239 + V_base(6) * t240 + qJD(2) + V_base(1);
t208 = rSges(4,1) * t192 - rSges(4,2) * t190;
t204 = Icges(4,5) * t192 - Icges(4,6) * t190;
t161 = Icges(4,2) * t192 + t232;
t162 = Icges(4,1) * t190 + t231;
t201 = t161 * t190 - t162 * t192;
t200 = -t140 * t184 + t183 * t220 + t213;
t199 = V_base(5) * t152 + t209;
t198 = t163 * V_base(6) - t181 * t220 + t210;
t197 = -qJD(4) * t192 + V_base(5) * t141 + t199;
t196 = (Icges(4,3) * t183 - t181 * t204) * V_base(5) + (Icges(4,3) * t181 + t183 * t204) * V_base(6) + (Icges(4,5) * t190 + Icges(4,6) * t192) * t184;
t187 = pkin(9) + qJ(5);
t182 = cos(t187);
t180 = sin(t187);
t174 = rSges(2,1) * t195 - rSges(2,2) * t194;
t173 = -rSges(2,1) * t194 - rSges(2,2) * t195;
t172 = Icges(2,1) * t195 - t236;
t171 = -Icges(2,1) * t194 - t235;
t170 = -Icges(2,2) * t194 + t235;
t169 = -Icges(2,2) * t195 - t236;
t166 = -qJD(5) * t192 + t184;
t164 = rSges(4,1) * t190 + rSges(4,2) * t192;
t159 = -V_base(5) * rSges(1,1) + rSges(1,2) * V_base(4) + V_base(3);
t158 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t157 = -rSges(1,2) * V_base(6) + V_base(5) * rSges(1,3) + V_base(1);
t156 = -t181 * t219 + V_base(5);
t155 = t183 * t219 + V_base(6);
t153 = rSges(3,1) * t183 - rSges(3,2) * t181;
t151 = -rSges(3,1) * t181 - rSges(3,2) * t183;
t147 = -Icges(3,2) * t181 + t233;
t146 = -Icges(3,2) * t183 - t234;
t139 = t183 * t223 + t230;
t138 = t181 * t191 - t183 * t224;
t137 = -t181 * t223 + t227;
t136 = t181 * t224 + t183 * t191;
t135 = -rSges(5,3) * t192 + (rSges(5,1) * t191 - rSges(5,2) * t189) * t190;
t134 = -Icges(5,5) * t192 + (Icges(5,1) * t191 - Icges(5,4) * t189) * t190;
t133 = -Icges(5,6) * t192 + (Icges(5,4) * t191 - Icges(5,2) * t189) * t190;
t132 = -Icges(5,3) * t192 + (Icges(5,5) * t191 - Icges(5,6) * t189) * t190;
t130 = t180 * t181 + t182 * t225;
t129 = -t180 * t225 + t181 * t182;
t128 = t180 * t183 - t182 * t228;
t127 = t180 * t228 + t182 * t183;
t126 = -rSges(6,3) * t192 + (rSges(6,1) * t182 - rSges(6,2) * t180) * t190;
t125 = -Icges(6,5) * t192 + (Icges(6,1) * t182 - Icges(6,4) * t180) * t190;
t124 = -Icges(6,6) * t192 + (Icges(6,4) * t182 - Icges(6,2) * t180) * t190;
t123 = -Icges(6,3) * t192 + (Icges(6,5) * t182 - Icges(6,6) * t180) * t190;
t121 = V_base(6) * rSges(2,3) - t174 * t184 + t218;
t120 = t173 * t184 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t119 = -pkin(6) * t192 + t190 * t238;
t118 = rSges(4,3) * t181 + t183 * t208;
t117 = rSges(4,3) * t183 - t181 * t208;
t110 = -t173 * V_base(6) + t174 * V_base(5) + V_base(1);
t109 = V_base(6) * rSges(3,3) + (-t153 - t239) * t184 + t212;
t108 = V_base(3) + (t151 - t240) * t184 + (-rSges(3,3) + t237) * V_base(5);
t107 = -t151 * V_base(6) + t153 * V_base(5) + t209;
t106 = rSges(5,1) * t139 + rSges(5,2) * t138 + rSges(5,3) * t226;
t105 = rSges(5,1) * t137 + rSges(5,2) * t136 - rSges(5,3) * t229;
t104 = pkin(4) * t230 + t183 * t242;
t103 = pkin(4) * t227 - t181 * t242;
t102 = Icges(5,1) * t139 + Icges(5,4) * t138 + Icges(5,5) * t226;
t101 = Icges(5,1) * t137 + Icges(5,4) * t136 - Icges(5,5) * t229;
t100 = Icges(5,4) * t139 + Icges(5,2) * t138 + Icges(5,6) * t226;
t99 = Icges(5,4) * t137 + Icges(5,2) * t136 - Icges(5,6) * t229;
t98 = Icges(5,5) * t139 + Icges(5,6) * t138 + Icges(5,3) * t226;
t97 = Icges(5,5) * t137 + Icges(5,6) * t136 - Icges(5,3) * t229;
t96 = rSges(6,1) * t130 + rSges(6,2) * t129 + rSges(6,3) * t226;
t95 = rSges(6,1) * t128 + rSges(6,2) * t127 - rSges(6,3) * t229;
t94 = Icges(6,1) * t130 + Icges(6,4) * t129 + Icges(6,5) * t226;
t93 = Icges(6,1) * t128 + Icges(6,4) * t127 - Icges(6,5) * t229;
t92 = Icges(6,4) * t130 + Icges(6,2) * t129 + Icges(6,6) * t226;
t91 = Icges(6,4) * t128 + Icges(6,2) * t127 - Icges(6,6) * t229;
t90 = Icges(6,5) * t130 + Icges(6,6) * t129 + Icges(6,3) * t226;
t89 = Icges(6,5) * t128 + Icges(6,6) * t127 - Icges(6,3) * t229;
t88 = V_base(6) * t164 + (-t118 + t214) * t184 + t210;
t87 = (t117 - t240) * t184 + (-t164 + t237) * V_base(5) + t213;
t86 = t118 * V_base(5) + (-t117 - t150) * V_base(6) + t199;
t85 = V_base(6) * t135 + (-t106 + t211) * t184 + t198;
t84 = (t105 - t240) * t184 + (-t135 + t215) * V_base(5) + t200;
t83 = t106 * V_base(5) + (-t105 + t221) * V_base(6) + t197;
t82 = V_base(6) * t119 + t155 * t126 - t166 * t96 + (-t104 + t211) * t184 + t198;
t81 = -t126 * t156 + t166 * t95 + (t103 - t240) * t184 + (-t119 + t215) * V_base(5) + t200;
t80 = t104 * V_base(5) - t155 * t95 + t156 * t96 + (-t103 + t221) * V_base(6) + t197;
t1 = m(1) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(2) * (t110 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(3) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(4) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + t166 * ((-t123 * t166 - t155 * t90 - t156 * t89) * t192 + ((-t124 * t180 + t125 * t182) * t166 + (-t180 * t91 + t182 * t93) * t156 + (-t180 * t92 + t182 * t94) * t155) * t190) / 0.2e1 + t156 * ((-t123 * t229 + t124 * t127 + t125 * t128) * t166 + (t127 * t91 + t128 * t93 - t89 * t229) * t156 + (t127 * t92 + t128 * t94 - t229 * t90) * t155) / 0.2e1 + t155 * ((t123 * t226 + t124 * t129 + t125 * t130) * t166 + (t129 * t91 + t130 * t93 + t226 * t89) * t156 + (t129 * t92 + t130 * t94 + t90 * t226) * t155) / 0.2e1 + (((t114 - t98) * t192 + (-t100 * t189 + t102 * t191 + t116) * t190 + t245) * V_base(6) + ((t113 - t97) * t192 + (t101 * t191 - t189 * t99 + t115) * t190 + t246) * V_base(5) + (Icges(2,3) + Icges(3,3) + (t161 - t132) * t192 + (-t133 * t189 + t134 * t191 + t162) * t190) * t184) * t184 / 0.2e1 + (t196 * t183 + (-t132 * t229 + t133 * t136 + t134 * t137 + t181 * t201 + t246) * t184 + (t100 * t136 + t102 * t137 - t147 * t183 - t170 * t195 - t194 * t172 - t181 * t243 - t229 * t98 + Icges(1,6)) * V_base(6) + (t101 * t137 + t136 * t99 - t146 * t183 - t169 * t195 - t194 * t171 - t181 * t244 - t229 * t97 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t196 * t181 + (t132 * t226 + t133 * t138 + t134 * t139 - t183 * t201 + t245) * t184 + (t100 * t138 + t102 * t139 - t147 * t181 - t194 * t170 + t172 * t195 + t183 * t243 + t226 * t98 + Icges(1,3)) * V_base(6) + (t101 * t139 + t138 * t99 - t146 * t181 - t194 * t169 + t171 * t195 + t183 * t244 + t226 * t97 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
