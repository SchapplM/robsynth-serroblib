% Calculate kinetic energy for
% S6RPPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:48
% EndTime: 2019-03-09 01:29:50
% DurationCPUTime: 2.14s
% Computational Cost: add. (1109->264), mult. (1106->348), div. (0->0), fcn. (886->8), ass. (0->132)
t246 = Icges(3,4) + Icges(4,6) - Icges(5,6);
t245 = Icges(3,1) + Icges(4,2) + Icges(5,3);
t244 = -Icges(4,4) + Icges(3,5) + Icges(5,5);
t243 = Icges(5,4) + Icges(4,5) - Icges(3,6);
t242 = Icges(3,2) + Icges(5,2) + Icges(4,3);
t179 = qJ(1) + pkin(9);
t171 = sin(t179);
t241 = t246 * t171;
t172 = cos(t179);
t240 = t246 * t172;
t239 = t245 * t172 - t241;
t238 = t245 * t171 + t240;
t237 = -t242 * t172 - t241;
t236 = t242 * t171 - t240;
t182 = sin(qJ(1));
t232 = pkin(1) * t182;
t185 = cos(qJ(1));
t231 = pkin(1) * t185;
t230 = pkin(7) * t171;
t229 = pkin(7) * t172;
t228 = -pkin(6) - qJ(2);
t227 = Icges(2,4) * t182;
t181 = sin(qJ(5));
t225 = Icges(6,4) * t181;
t184 = cos(qJ(5));
t224 = Icges(6,4) * t184;
t220 = qJ(4) * t171;
t219 = qJ(4) * t172;
t218 = t171 * t184;
t217 = t172 * t184;
t180 = sin(qJ(6));
t216 = t180 * t181;
t183 = cos(qJ(6));
t215 = t181 * t183;
t214 = qJD(6) * t184;
t213 = -pkin(3) + t228;
t173 = V_base(6) + qJD(1);
t212 = t173 * t231 + V_base(2);
t211 = V_base(5) * pkin(6) + V_base(1);
t144 = qJD(5) * t172 + V_base(5);
t132 = pkin(2) * t171 - qJ(3) * t172;
t208 = -t132 - t232;
t136 = pkin(2) * t172 + qJ(3) * t171;
t207 = -t136 - t231;
t206 = V_base(5) * qJ(2) + t211;
t205 = V_base(4) * t232 + qJD(2) + V_base(3);
t204 = qJD(3) * t171 + t206;
t203 = pkin(5) * t181 - pkin(8) * t184;
t202 = V_base(4) * t132 + t205;
t145 = -qJD(5) * t171 + V_base(4);
t201 = rSges(6,1) * t181 + rSges(6,2) * t184;
t200 = Icges(6,1) * t181 + t224;
t199 = Icges(6,2) * t184 + t225;
t198 = Icges(6,5) * t181 + Icges(6,6) * t184;
t197 = V_base(4) * t220 + t202;
t196 = t208 - t220;
t195 = t207 - t219;
t194 = -qJD(3) * t172 + t173 * t136 + t212;
t193 = V_base(5) * pkin(3) + qJD(4) * t172 + t204;
t192 = (Icges(6,3) * t172 + t171 * t198) * t144 + (-Icges(6,3) * t171 + t172 * t198) * t145 + (Icges(6,5) * t184 - Icges(6,6) * t181) * t173;
t191 = V_base(5) * pkin(4) + t193;
t190 = t196 - t229;
t189 = qJD(4) * t171 + t173 * t219 + t194;
t188 = (-pkin(4) + t213) * V_base(4) + t189;
t187 = V_base(4) * t229 + t197 + (t195 + t230) * V_base(5);
t149 = -Icges(6,2) * t181 + t224;
t152 = Icges(6,1) * t184 - t225;
t92 = Icges(6,6) * t172 + t171 * t199;
t93 = -Icges(6,6) * t171 + t172 * t199;
t94 = Icges(6,5) * t172 + t171 * t200;
t95 = -Icges(6,5) * t171 + t172 * t200;
t186 = (t181 * t95 + t184 * t93) * t145 + (t181 * t94 + t184 * t92) * t144 + (t149 * t184 + t152 * t181) * t173;
t175 = Icges(2,4) * t185;
t160 = pkin(5) * t184 + pkin(8) * t181;
t158 = rSges(2,1) * t185 - t182 * rSges(2,2);
t157 = rSges(6,1) * t184 - rSges(6,2) * t181;
t156 = t182 * rSges(2,1) + rSges(2,2) * t185;
t155 = qJD(6) * t181 + t173;
t154 = Icges(2,1) * t185 - t227;
t153 = Icges(2,1) * t182 + t175;
t151 = -Icges(2,2) * t182 + t175;
t150 = Icges(2,2) * t185 + t227;
t143 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t142 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t141 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t139 = rSges(3,1) * t172 - rSges(3,2) * t171;
t138 = -rSges(4,2) * t172 + rSges(4,3) * t171;
t137 = -rSges(5,2) * t172 + rSges(5,3) * t171;
t135 = rSges(3,1) * t171 + rSges(3,2) * t172;
t134 = -rSges(4,2) * t171 - rSges(4,3) * t172;
t133 = rSges(5,2) * t171 + rSges(5,3) * t172;
t112 = t203 * t172;
t111 = t203 * t171;
t110 = rSges(7,3) * t181 + (rSges(7,1) * t183 - rSges(7,2) * t180) * t184;
t108 = Icges(7,5) * t181 + (Icges(7,1) * t183 - Icges(7,4) * t180) * t184;
t107 = Icges(7,6) * t181 + (Icges(7,4) * t183 - Icges(7,2) * t180) * t184;
t106 = Icges(7,3) * t181 + (Icges(7,5) * t183 - Icges(7,6) * t180) * t184;
t105 = -t171 * t180 + t172 * t215;
t104 = -t171 * t183 - t172 * t216;
t103 = t171 * t215 + t172 * t180;
t102 = -t171 * t216 + t172 * t183;
t101 = -t172 * t214 + t145;
t100 = -t171 * t214 + t144;
t99 = -rSges(6,3) * t171 + t172 * t201;
t98 = rSges(6,3) * t172 + t171 * t201;
t97 = V_base(5) * rSges(2,3) - t156 * t173 + t211;
t96 = t158 * t173 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t89 = t156 * V_base(4) - t158 * V_base(5) + V_base(3);
t88 = V_base(5) * rSges(3,3) + (-t135 - t232) * t173 + t206;
t87 = t139 * t173 + (-rSges(3,3) + t228) * V_base(4) + t212;
t86 = V_base(4) * t135 + (-t139 - t231) * V_base(5) + t205;
t85 = rSges(7,1) * t105 + rSges(7,2) * t104 - rSges(7,3) * t217;
t84 = rSges(7,1) * t103 + rSges(7,2) * t102 - rSges(7,3) * t218;
t83 = Icges(7,1) * t105 + Icges(7,4) * t104 - Icges(7,5) * t217;
t82 = Icges(7,1) * t103 + Icges(7,4) * t102 - Icges(7,5) * t218;
t81 = Icges(7,4) * t105 + Icges(7,2) * t104 - Icges(7,6) * t217;
t80 = Icges(7,4) * t103 + Icges(7,2) * t102 - Icges(7,6) * t218;
t79 = Icges(7,5) * t105 + Icges(7,6) * t104 - Icges(7,3) * t217;
t78 = Icges(7,5) * t103 + Icges(7,6) * t102 - Icges(7,3) * t218;
t77 = V_base(5) * rSges(4,1) + (-t134 + t208) * t173 + t204;
t76 = t138 * t173 + (-rSges(4,1) + t228) * V_base(4) + t194;
t75 = V_base(4) * t134 + (-t138 + t207) * V_base(5) + t202;
t74 = V_base(5) * rSges(5,1) + (-t137 + t196) * t173 + t193;
t73 = t133 * t173 + (-rSges(5,1) + t213) * V_base(4) + t189;
t72 = V_base(4) * t137 + (-t133 + t195) * V_base(5) + t197;
t71 = t144 * t157 + (t190 - t98) * t173 + t191;
t70 = -t145 * t157 + (t99 - t230) * t173 + t188;
t69 = -t144 * t99 + t145 * t98 + t187;
t68 = t100 * t110 + t144 * t160 - t155 * t84 + (-t111 + t190) * t173 + t191;
t67 = -t101 * t110 - t145 * t160 + t155 * t85 + (t112 - t230) * t173 + t188;
t66 = -t100 * t85 + t101 * t84 + t145 * t111 - t144 * t112 + t187;
t1 = t101 * ((t104 * t81 + t105 * t83 - t79 * t217) * t101 + (t104 * t80 + t105 * t82 - t217 * t78) * t100 + (t104 * t107 + t105 * t108 - t106 * t217) * t155) / 0.2e1 + t100 * ((t102 * t81 + t103 * t83 - t218 * t79) * t101 + (t102 * t80 + t103 * t82 - t78 * t218) * t100 + (t102 * t107 + t103 * t108 - t106 * t218) * t155) / 0.2e1 + t144 * (t186 * t171 + t192 * t172) / 0.2e1 + t145 * (-t192 * t171 + t186 * t172) / 0.2e1 + m(7) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + m(4) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(3) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(2) * (t89 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(1) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + t155 * ((t78 * t100 + t79 * t101 + t106 * t155) * t181 + ((-t180 * t81 + t183 * t83) * t101 + (-t180 * t80 + t183 * t82) * t100 + (-t107 * t180 + t108 * t183) * t155) * t184) / 0.2e1 + ((-t181 * t93 + t184 * t95) * t145 + (-t181 * t92 + t184 * t94) * t144 + (-t181 * t149 + t184 * t152 + Icges(4,1) + Icges(5,1) + Icges(2,3) + Icges(3,3)) * t173) * t173 / 0.2e1 + ((-t182 * t150 + t153 * t185 + t171 * t237 + t172 * t238 + Icges(1,4)) * V_base(5) + (-t182 * t151 + t185 * t154 + t236 * t171 + t172 * t239 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t185 * t150 + t182 * t153 + t171 * t238 - t172 * t237 + Icges(1,2)) * V_base(5) + (t151 * t185 + t182 * t154 + t171 * t239 - t236 * t172 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t173 * (Icges(2,5) * t182 + Icges(2,6) * t185 + t244 * t171 - t243 * t172) + V_base(4) * t173 * (Icges(2,5) * t185 - Icges(2,6) * t182 + t243 * t171 + t244 * t172) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
