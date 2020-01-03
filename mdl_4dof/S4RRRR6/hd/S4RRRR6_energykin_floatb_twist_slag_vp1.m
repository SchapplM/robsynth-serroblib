% Calculate kinetic energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:19
% EndTime: 2019-12-31 17:29:21
% DurationCPUTime: 2.03s
% Computational Cost: add. (1271->278), mult. (2890->419), div. (0->0), fcn. (3384->10), ass. (0->121)
t232 = cos(qJ(3));
t204 = cos(pkin(4));
t231 = pkin(6) * t204;
t208 = sin(qJ(1));
t230 = Icges(2,4) * t208;
t203 = sin(pkin(4));
t229 = t203 * t208;
t210 = cos(qJ(2));
t228 = t203 * t210;
t211 = cos(qJ(1));
t227 = t203 * t211;
t207 = sin(qJ(2));
t226 = t207 * t208;
t225 = t207 * t211;
t224 = t208 * t210;
t223 = t210 * t211;
t222 = qJD(2) * t203;
t221 = V_base(5) * pkin(5) + V_base(1);
t218 = t203 * t232;
t184 = t208 * t222 + V_base(4);
t200 = V_base(6) + qJD(1);
t174 = t204 * t224 + t225;
t153 = qJD(3) * t174 + t184;
t185 = qJD(2) * t204 + t200;
t183 = -t211 * t222 + V_base(5);
t178 = t208 * pkin(1) - pkin(6) * t227;
t217 = -t178 * t200 + V_base(5) * t231 + t221;
t179 = pkin(1) * t211 + pkin(6) * t229;
t216 = V_base(4) * t178 - t179 * V_base(5) + V_base(3);
t172 = -t204 * t223 + t226;
t152 = qJD(3) * t172 + t183;
t168 = -qJD(3) * t228 + t185;
t215 = t200 * t179 + V_base(2) + (-pkin(5) - t231) * V_base(4);
t173 = t204 * t225 + t224;
t148 = pkin(2) * t173 + pkin(7) * t172;
t177 = (pkin(2) * t207 - pkin(7) * t210) * t203;
t214 = -t148 * t185 + t183 * t177 + t217;
t175 = -t204 * t226 + t223;
t149 = pkin(2) * t175 + pkin(7) * t174;
t213 = t184 * t148 - t149 * t183 + t216;
t212 = t185 * t149 - t177 * t184 + t215;
t209 = cos(qJ(4));
t206 = sin(qJ(3));
t205 = sin(qJ(4));
t201 = Icges(2,4) * t211;
t193 = rSges(2,1) * t211 - t208 * rSges(2,2);
t192 = t208 * rSges(2,1) + rSges(2,2) * t211;
t191 = Icges(2,1) * t211 - t230;
t190 = Icges(2,1) * t208 + t201;
t189 = -Icges(2,2) * t208 + t201;
t188 = Icges(2,2) * t211 + t230;
t182 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t181 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t180 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t171 = t204 * t206 + t207 * t218;
t170 = t203 * t206 * t207 - t204 * t232;
t165 = rSges(3,3) * t204 + (rSges(3,1) * t207 + rSges(3,2) * t210) * t203;
t164 = Icges(3,5) * t204 + (Icges(3,1) * t207 + Icges(3,4) * t210) * t203;
t163 = Icges(3,6) * t204 + (Icges(3,4) * t207 + Icges(3,2) * t210) * t203;
t162 = Icges(3,3) * t204 + (Icges(3,5) * t207 + Icges(3,6) * t210) * t203;
t161 = V_base(5) * rSges(2,3) - t192 * t200 + t221;
t160 = t193 * t200 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t158 = t192 * V_base(4) - t193 * V_base(5) + V_base(3);
t157 = t175 * t232 + t206 * t229;
t156 = t175 * t206 - t208 * t218;
t155 = t173 * t232 - t206 * t227;
t154 = t173 * t206 + t211 * t218;
t151 = t171 * t209 - t205 * t228;
t150 = -t171 * t205 - t209 * t228;
t147 = qJD(4) * t170 + t168;
t146 = pkin(3) * t171 + pkin(8) * t170;
t145 = rSges(3,1) * t175 - rSges(3,2) * t174 + rSges(3,3) * t229;
t144 = t173 * rSges(3,1) - t172 * rSges(3,2) - rSges(3,3) * t227;
t143 = Icges(3,1) * t175 - Icges(3,4) * t174 + Icges(3,5) * t229;
t142 = Icges(3,1) * t173 - Icges(3,4) * t172 - Icges(3,5) * t227;
t141 = Icges(3,4) * t175 - Icges(3,2) * t174 + Icges(3,6) * t229;
t140 = Icges(3,4) * t173 - Icges(3,2) * t172 - Icges(3,6) * t227;
t139 = Icges(3,5) * t175 - Icges(3,6) * t174 + Icges(3,3) * t229;
t138 = Icges(3,5) * t173 - Icges(3,6) * t172 - Icges(3,3) * t227;
t137 = rSges(4,1) * t171 - rSges(4,2) * t170 - rSges(4,3) * t228;
t136 = Icges(4,1) * t171 - Icges(4,4) * t170 - Icges(4,5) * t228;
t135 = Icges(4,4) * t171 - Icges(4,2) * t170 - Icges(4,6) * t228;
t134 = Icges(4,5) * t171 - Icges(4,6) * t170 - Icges(4,3) * t228;
t131 = t157 * t209 + t174 * t205;
t130 = -t157 * t205 + t174 * t209;
t129 = t155 * t209 + t172 * t205;
t128 = -t155 * t205 + t172 * t209;
t127 = qJD(4) * t156 + t153;
t126 = qJD(4) * t154 + t152;
t125 = pkin(3) * t157 + pkin(8) * t156;
t124 = pkin(3) * t155 + pkin(8) * t154;
t123 = rSges(4,1) * t157 - rSges(4,2) * t156 + rSges(4,3) * t174;
t122 = rSges(4,1) * t155 - rSges(4,2) * t154 + rSges(4,3) * t172;
t121 = Icges(4,1) * t157 - Icges(4,4) * t156 + Icges(4,5) * t174;
t120 = Icges(4,1) * t155 - Icges(4,4) * t154 + Icges(4,5) * t172;
t119 = Icges(4,4) * t157 - Icges(4,2) * t156 + Icges(4,6) * t174;
t118 = Icges(4,4) * t155 - Icges(4,2) * t154 + Icges(4,6) * t172;
t117 = Icges(4,5) * t157 - Icges(4,6) * t156 + Icges(4,3) * t174;
t116 = Icges(4,5) * t155 - Icges(4,6) * t154 + Icges(4,3) * t172;
t115 = rSges(5,1) * t151 + rSges(5,2) * t150 + rSges(5,3) * t170;
t114 = Icges(5,1) * t151 + Icges(5,4) * t150 + Icges(5,5) * t170;
t113 = Icges(5,4) * t151 + Icges(5,2) * t150 + Icges(5,6) * t170;
t112 = Icges(5,5) * t151 + Icges(5,6) * t150 + Icges(5,3) * t170;
t111 = -t144 * t185 + t165 * t183 + t217;
t110 = t145 * t185 - t165 * t184 + t215;
t109 = rSges(5,1) * t131 + rSges(5,2) * t130 + rSges(5,3) * t156;
t108 = rSges(5,1) * t129 + rSges(5,2) * t128 + rSges(5,3) * t154;
t107 = Icges(5,1) * t131 + Icges(5,4) * t130 + Icges(5,5) * t156;
t106 = Icges(5,1) * t129 + Icges(5,4) * t128 + Icges(5,5) * t154;
t105 = Icges(5,4) * t131 + Icges(5,2) * t130 + Icges(5,6) * t156;
t104 = Icges(5,4) * t129 + Icges(5,2) * t128 + Icges(5,6) * t154;
t103 = Icges(5,5) * t131 + Icges(5,6) * t130 + Icges(5,3) * t156;
t102 = Icges(5,5) * t129 + Icges(5,6) * t128 + Icges(5,3) * t154;
t101 = t144 * t184 - t145 * t183 + t216;
t100 = -t122 * t168 + t137 * t152 + t214;
t99 = t123 * t168 - t137 * t153 + t212;
t98 = t122 * t153 - t123 * t152 + t213;
t97 = -t108 * t147 + t115 * t126 - t124 * t168 + t146 * t152 + t214;
t96 = t109 * t147 - t115 * t127 + t125 * t168 - t146 * t153 + t212;
t95 = t108 * t127 - t109 * t126 + t124 * t153 - t125 * t152 + t213;
t1 = m(1) * (t180 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + m(2) * (t158 ^ 2 + t160 ^ 2 + t161 ^ 2) / 0.2e1 + m(3) * (t101 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + t184 * ((t139 * t229 - t141 * t174 + t143 * t175) * t184 + (t138 * t229 - t140 * t174 + t142 * t175) * t183 + (t162 * t229 - t163 * t174 + t164 * t175) * t185) / 0.2e1 + t183 * ((-t139 * t227 - t172 * t141 + t173 * t143) * t184 + (-t138 * t227 - t172 * t140 + t173 * t142) * t183 + (-t162 * t227 - t172 * t163 + t173 * t164) * t185) / 0.2e1 + t185 * ((t138 * t183 + t139 * t184 + t162 * t185) * t204 + ((t141 * t210 + t143 * t207) * t184 + (t140 * t210 + t142 * t207) * t183 + (t163 * t210 + t164 * t207) * t185) * t203) / 0.2e1 + m(4) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + t153 * ((t174 * t117 - t156 * t119 + t157 * t121) * t153 + (t116 * t174 - t118 * t156 + t120 * t157) * t152 + (t134 * t174 - t135 * t156 + t136 * t157) * t168) / 0.2e1 + t152 * ((t117 * t172 - t119 * t154 + t121 * t155) * t153 + (t172 * t116 - t154 * t118 + t155 * t120) * t152 + (t134 * t172 - t135 * t154 + t136 * t155) * t168) / 0.2e1 + t168 * ((-t117 * t228 - t119 * t170 + t121 * t171) * t153 + (-t116 * t228 - t118 * t170 + t120 * t171) * t152 + (-t134 * t228 - t170 * t135 + t171 * t136) * t168) / 0.2e1 + m(5) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + t127 * ((t156 * t103 + t130 * t105 + t131 * t107) * t127 + (t102 * t156 + t104 * t130 + t106 * t131) * t126 + (t112 * t156 + t113 * t130 + t114 * t131) * t147) / 0.2e1 + t126 * ((t103 * t154 + t105 * t128 + t107 * t129) * t127 + (t154 * t102 + t128 * t104 + t129 * t106) * t126 + (t112 * t154 + t113 * t128 + t114 * t129) * t147) / 0.2e1 + t147 * ((t103 * t170 + t105 * t150 + t107 * t151) * t127 + (t102 * t170 + t104 * t150 + t106 * t151) * t126 + (t170 * t112 + t150 * t113 + t151 * t114) * t147) / 0.2e1 + ((-t208 * t188 + t190 * t211 + Icges(1,4)) * V_base(5) + (-t208 * t189 + t191 * t211 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t188 * t211 + t208 * t190 + Icges(1,2)) * V_base(5) + (t189 * t211 + t208 * t191 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t208 + Icges(2,6) * t211) * V_base(5) + (Icges(2,5) * t211 - Icges(2,6) * t208) * V_base(4) + Icges(2,3) * t200 / 0.2e1) * t200;
T = t1;
