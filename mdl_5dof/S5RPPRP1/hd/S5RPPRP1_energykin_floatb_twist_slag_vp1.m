% Calculate kinetic energy for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:04
% EndTime: 2020-01-03 11:25:07
% DurationCPUTime: 2.18s
% Computational Cost: add. (1213->248), mult. (1371->330), div. (0->0), fcn. (1263->8), ass. (0->123)
t241 = Icges(5,1) + Icges(6,1);
t240 = Icges(5,4) + Icges(6,4);
t239 = -Icges(6,5) - Icges(5,5);
t238 = Icges(5,2) + Icges(6,2);
t237 = -Icges(6,6) - Icges(5,6);
t236 = -Icges(6,3) - Icges(5,3);
t167 = qJ(1) + pkin(7);
t161 = sin(t167);
t162 = cos(t167);
t173 = cos(qJ(4));
t169 = cos(pkin(8));
t171 = sin(qJ(4));
t202 = t169 * t171;
t119 = -t161 * t202 - t162 * t173;
t201 = t169 * t173;
t203 = t162 * t171;
t120 = t161 * t201 - t203;
t168 = sin(pkin(8));
t206 = t161 * t168;
t235 = -t237 * t119 - t239 * t120 - t236 * t206;
t121 = -t161 * t173 + t162 * t202;
t205 = t161 * t171;
t122 = -t162 * t201 - t205;
t204 = t162 * t168;
t234 = -t237 * t121 - t239 * t122 + t236 * t204;
t233 = t238 * t119 + t240 * t120 - t237 * t206;
t232 = t238 * t121 + t240 * t122 + t237 * t204;
t231 = t240 * t119 + t241 * t120 - t239 * t206;
t230 = t240 * t121 + t241 * t122 + t239 * t204;
t229 = t236 * t169 + (t237 * t171 - t239 * t173) * t168;
t228 = t237 * t169 + (-t238 * t171 + t240 * t173) * t168;
t227 = t239 * t169 + (-t240 * t171 + t241 * t173) * t168;
t172 = sin(qJ(1));
t174 = cos(qJ(1));
t226 = Icges(2,5) * t172 + Icges(3,5) * t161 + Icges(2,6) * t174 + Icges(3,6) * t162;
t225 = -Icges(2,5) * t174 - Icges(3,5) * t162 + Icges(2,6) * t172 + Icges(3,6) * t161;
t207 = Icges(4,4) * t169;
t185 = -Icges(4,2) * t168 + t207;
t100 = -Icges(4,6) * t162 + t161 * t185;
t208 = Icges(4,4) * t168;
t186 = Icges(4,1) * t169 - t208;
t102 = -Icges(4,5) * t162 + t161 * t186;
t209 = Icges(3,4) * t162;
t224 = Icges(3,1) * t161 - t100 * t168 + t102 * t169 + t209;
t101 = -Icges(4,6) * t161 - t162 * t185;
t103 = -Icges(4,5) * t161 - t162 * t186;
t159 = Icges(3,4) * t161;
t223 = -Icges(3,1) * t162 - t101 * t168 + t103 * t169 + t159;
t215 = pkin(4) * t173;
t222 = qJ(5) * t168 + t169 * t215;
t217 = pkin(1) * t172;
t216 = pkin(1) * t174;
t214 = -pkin(5) - qJ(2);
t212 = rSges(6,1) * t120 + rSges(6,2) * t119 + rSges(6,3) * t206 - pkin(4) * t203 + t161 * t222;
t211 = rSges(6,1) * t122 + rSges(6,2) * t121 - rSges(6,3) * t204 - pkin(4) * t205 - t162 * t222;
t210 = Icges(2,4) * t174;
t200 = (-qJ(5) - rSges(6,3)) * t169 + (rSges(6,1) * t173 - rSges(6,2) * t171 + t215) * t168;
t199 = qJD(4) * t168;
t198 = qJD(5) * t168;
t163 = V_base(4) + qJD(1);
t197 = t163 * t217 + V_base(3);
t196 = V_base(6) * pkin(5) + V_base(2);
t193 = qJD(2) + V_base(1);
t192 = t174 * V_base(5);
t133 = pkin(2) * t161 - qJ(3) * t162;
t191 = -t133 - t217;
t135 = -pkin(2) * t162 - qJ(3) * t161;
t190 = V_base(5) * t135 + t193;
t189 = V_base(6) * qJ(2) + t163 * t216 + t196;
t188 = pkin(3) * t169 + pkin(6) * t168;
t187 = rSges(4,1) * t169 - rSges(4,2) * t168;
t184 = Icges(4,5) * t169 - Icges(4,6) * t168;
t144 = Icges(4,2) * t169 + t208;
t145 = Icges(4,1) * t168 + t207;
t181 = t144 * t168 - t145 * t169;
t180 = -qJD(3) * t161 + t163 * t133 + t197;
t179 = -qJD(3) * t162 + t189;
t178 = -(Icges(4,5) * t168 + Icges(4,6) * t169) * t163 - (-Icges(4,3) * t162 + t161 * t184) * V_base(5) - (-Icges(4,3) * t161 - t162 * t184) * V_base(6);
t124 = t188 * t162;
t148 = pkin(3) * t168 - pkin(6) * t169;
t177 = V_base(6) * t148 + (t124 - t135) * t163 + t179;
t123 = t188 * t161;
t176 = t163 * t123 + (-t148 + t214) * V_base(5) + t180;
t175 = -V_base(5) * t124 + (-t123 + t191) * V_base(6) - pkin(1) * t192 + t190;
t165 = Icges(2,4) * t172;
t156 = -rSges(2,1) * t174 + t172 * rSges(2,2);
t155 = t172 * rSges(2,1) + rSges(2,2) * t174;
t154 = -Icges(2,1) * t174 + t165;
t153 = Icges(2,1) * t172 + t210;
t152 = Icges(2,2) * t172 - t210;
t151 = Icges(2,2) * t174 + t165;
t147 = -qJD(4) * t169 + t163;
t146 = rSges(4,1) * t168 + rSges(4,2) * t169;
t142 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t141 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t140 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t138 = t161 * t199 + V_base(5);
t137 = -t162 * t199 + V_base(6);
t136 = -rSges(3,1) * t162 + rSges(3,2) * t161;
t134 = rSges(3,1) * t161 + rSges(3,2) * t162;
t130 = Icges(3,2) * t161 - t209;
t129 = Icges(3,2) * t162 + t159;
t118 = -rSges(5,3) * t169 + (rSges(5,1) * t173 - rSges(5,2) * t171) * t168;
t107 = V_base(6) * rSges(2,3) - t156 * t163 + t196;
t106 = t155 * t163 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t105 = -rSges(4,3) * t161 - t162 * t187;
t104 = -rSges(4,3) * t162 + t161 * t187;
t97 = -t155 * V_base(6) + t156 * V_base(5) + V_base(1);
t96 = V_base(6) * rSges(3,3) - t136 * t163 + t189;
t95 = t134 * t163 + (-rSges(3,3) + t214) * V_base(5) + t197;
t94 = -V_base(6) * t134 + V_base(5) * t136 + (-t172 * V_base(6) - t192) * pkin(1) + t193;
t93 = rSges(5,1) * t122 + rSges(5,2) * t121 - rSges(5,3) * t204;
t91 = rSges(5,1) * t120 + rSges(5,2) * t119 + rSges(5,3) * t206;
t75 = t146 * V_base(6) + (-t105 - t135) * t163 + t179;
t74 = t104 * t163 + (-t146 + t214) * V_base(5) + t180;
t73 = (t105 - t216) * V_base(5) + (-t104 + t191) * V_base(6) + t190;
t72 = t118 * t137 - t147 * t93 + t177;
t71 = -t118 * t138 + t147 * t91 + t176;
t70 = -t137 * t91 + t138 * t93 + t175;
t69 = t137 * t200 - t147 * t211 + t161 * t198 + t177;
t68 = -t138 * t200 + t147 * t212 - t162 * t198 + t176;
t67 = -qJD(5) * t169 - t137 * t212 + t138 * t211 + t175;
t1 = m(1) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(2) * (t106 ^ 2 + t107 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(4) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(5) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(6) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + ((t121 * t228 + t122 * t227 - t204 * t229) * t147 + (t233 * t121 + t231 * t122 - t204 * t235) * t138 + (t232 * t121 + t230 * t122 - t234 * t204) * t137) * t137 / 0.2e1 + ((t119 * t228 + t120 * t227 + t206 * t229) * t147 + (t233 * t119 + t231 * t120 + t235 * t206) * t138 + (t119 * t232 + t120 * t230 + t206 * t234) * t137) * t138 / 0.2e1 + ((-t234 * t137 - t138 * t235 - t229 * t147) * t169 + ((-t171 * t228 + t173 * t227) * t147 + (-t171 * t233 + t173 * t231) * t138 + (-t171 * t232 + t173 * t230) * t137) * t168) * t147 / 0.2e1 + ((t101 * t169 + t103 * t168 + t225) * V_base(6) + (t100 * t169 + t102 * t168 + t226) * V_base(5) + (t169 * t144 + t168 * t145 + Icges(2,3) + Icges(3,3)) * t163) * t163 / 0.2e1 + (t178 * t162 + (-t181 * t161 + t226) * t163 + (t130 * t162 + t152 * t174 + t172 * t154 + t161 * t223 + Icges(1,6)) * V_base(6) + (t129 * t162 + t151 * t174 + t172 * t153 + t161 * t224 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t178 * t161 + (t181 * t162 + t225) * t163 + (t130 * t161 + t172 * t152 - t154 * t174 - t162 * t223 + Icges(1,3)) * V_base(6) + (t129 * t161 + t172 * t151 - t153 * t174 - t162 * t224 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
