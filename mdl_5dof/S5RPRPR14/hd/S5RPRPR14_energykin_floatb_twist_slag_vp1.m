% Calculate kinetic energy for
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR14_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR14_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR14_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR14_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:27
% EndTime: 2019-12-31 18:34:29
% DurationCPUTime: 2.16s
% Computational Cost: add. (887->253), mult. (1150->347), div. (0->0), fcn. (988->8), ass. (0->128)
t249 = Icges(2,4) + Icges(3,6);
t248 = Icges(2,1) + Icges(3,2);
t247 = -Icges(3,4) + Icges(2,5);
t246 = Icges(3,5) - Icges(2,6);
t245 = Icges(2,2) + Icges(3,3);
t244 = Icges(4,3) + Icges(5,3);
t169 = qJ(3) + pkin(8);
t158 = sin(t169);
t159 = cos(t169);
t172 = sin(qJ(3));
t175 = cos(qJ(3));
t243 = Icges(4,5) * t172 + Icges(5,5) * t158 + Icges(4,6) * t175 + Icges(5,6) * t159;
t176 = cos(qJ(1));
t242 = t249 * t176;
t173 = sin(qJ(1));
t241 = t249 * t173;
t220 = Icges(4,4) * t172;
t192 = Icges(4,2) * t175 + t220;
t104 = Icges(4,6) * t176 + t173 * t192;
t105 = Icges(4,6) * t173 - t176 * t192;
t219 = Icges(4,4) * t175;
t194 = Icges(4,1) * t172 + t219;
t106 = Icges(4,5) * t176 + t173 * t194;
t107 = Icges(4,5) * t173 - t176 * t194;
t217 = Icges(5,4) * t159;
t122 = -Icges(5,2) * t158 + t217;
t218 = Icges(5,4) * t158;
t123 = Icges(5,1) * t159 - t218;
t139 = -Icges(4,2) * t172 + t219;
t144 = Icges(4,1) * t175 - t220;
t155 = qJD(3) * t173 + V_base(5);
t156 = qJD(3) * t176 + V_base(4);
t160 = V_base(6) + qJD(1);
t191 = Icges(5,2) * t159 + t218;
t95 = Icges(5,6) * t176 + t173 * t191;
t96 = Icges(5,6) * t173 - t176 * t191;
t193 = Icges(5,1) * t158 + t217;
t97 = Icges(5,5) * t176 + t173 * t193;
t98 = Icges(5,5) * t173 - t176 * t193;
t240 = t155 * (t105 * t175 + t107 * t172 + t158 * t98 + t159 * t96) + t156 * (t104 * t175 + t106 * t172 + t158 * t97 + t159 * t95) + t160 * (t122 * t159 + t123 * t158 + t139 * t175 + t144 * t172);
t239 = -t245 * t176 - t241;
t238 = t245 * t173 - t242;
t237 = t248 * t173 + t242;
t236 = t248 * t176 - t241;
t233 = (Icges(4,5) * t175 + Icges(5,5) * t159 - Icges(4,6) * t172 - Icges(5,6) * t158) * t160 + (t243 * t173 + t244 * t176) * t156 + (t244 * t173 - t243 * t176) * t155;
t226 = pkin(3) * t172;
t225 = pkin(3) * t175;
t224 = pkin(6) * t176;
t223 = t173 * pkin(6);
t214 = t159 * t173;
t213 = t159 * t176;
t171 = sin(qJ(5));
t212 = t171 * t176;
t211 = t173 * t171;
t174 = cos(qJ(5));
t210 = t173 * t174;
t209 = t174 * t176;
t208 = qJD(5) * t159;
t147 = t173 * pkin(1) - qJ(2) * t176;
t207 = V_base(4) * t147 + V_base(3);
t206 = V_base(5) * pkin(5) + V_base(1);
t203 = -t147 - t223;
t202 = qJD(2) * t173 + t206;
t112 = qJ(4) * t173 - t176 * t226;
t201 = -t112 + t203;
t200 = V_base(5) * pkin(2) + t202;
t199 = pkin(4) * t158 - pkin(7) * t159;
t198 = rSges(4,1) * t172 + rSges(4,2) * t175;
t197 = rSges(5,1) * t158 + rSges(5,2) * t159;
t151 = pkin(1) * t176 + t173 * qJ(2);
t184 = -qJD(2) * t176 + t160 * t151 + V_base(2);
t183 = qJD(4) * t176 + t155 * t225 + t200;
t180 = V_base(4) * t223 + (-t151 - t224) * V_base(5) + t207;
t179 = t160 * t224 + (-pkin(2) - pkin(5)) * V_base(4) + t184;
t178 = t156 * t112 + t180;
t113 = qJ(4) * t176 + t173 * t226;
t177 = qJD(4) * t173 + t160 * t113 + t179;
t153 = rSges(2,1) * t176 - t173 * rSges(2,2);
t152 = -rSges(3,2) * t176 + t173 * rSges(3,3);
t150 = rSges(4,1) * t175 - rSges(4,2) * t172;
t149 = t173 * rSges(2,1) + rSges(2,2) * t176;
t148 = -t173 * rSges(3,2) - rSges(3,3) * t176;
t131 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t130 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t129 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t128 = qJD(5) * t158 + t160;
t125 = pkin(4) * t159 + pkin(7) * t158;
t124 = rSges(5,1) * t159 - rSges(5,2) * t158;
t119 = -t158 * t209 + t211;
t118 = t158 * t212 + t210;
t117 = t158 * t210 + t212;
t116 = -t158 * t211 + t209;
t115 = -t173 * t208 + t156;
t114 = t176 * t208 + t155;
t111 = t199 * t176;
t110 = t199 * t173;
t109 = t173 * rSges(4,3) - t176 * t198;
t108 = rSges(4,3) * t176 + t173 * t198;
t100 = t173 * rSges(5,3) - t176 * t197;
t99 = rSges(5,3) * t176 + t173 * t197;
t92 = rSges(6,3) * t158 + (rSges(6,1) * t174 - rSges(6,2) * t171) * t159;
t91 = V_base(5) * rSges(2,3) - t149 * t160 + t206;
t90 = t153 * t160 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t89 = Icges(6,5) * t158 + (Icges(6,1) * t174 - Icges(6,4) * t171) * t159;
t88 = Icges(6,6) * t158 + (Icges(6,4) * t174 - Icges(6,2) * t171) * t159;
t87 = Icges(6,3) * t158 + (Icges(6,5) * t174 - Icges(6,6) * t171) * t159;
t85 = t149 * V_base(4) - t153 * V_base(5) + V_base(3);
t84 = V_base(5) * rSges(3,1) + (-t147 - t148) * t160 + t202;
t83 = t160 * t152 + (-rSges(3,1) - pkin(5)) * V_base(4) + t184;
t82 = t119 * rSges(6,1) + t118 * rSges(6,2) + rSges(6,3) * t213;
t81 = rSges(6,1) * t117 + rSges(6,2) * t116 - rSges(6,3) * t214;
t80 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t213;
t79 = Icges(6,1) * t117 + Icges(6,4) * t116 - Icges(6,5) * t214;
t78 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t213;
t77 = Icges(6,4) * t117 + Icges(6,2) * t116 - Icges(6,6) * t214;
t76 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t213;
t75 = Icges(6,5) * t117 + Icges(6,6) * t116 - Icges(6,3) * t214;
t74 = t148 * V_base(4) + (-t151 - t152) * V_base(5) + t207;
t73 = t150 * t155 + (-t109 + t203) * t160 + t200;
t72 = t160 * t108 - t156 * t150 + t179;
t71 = -t155 * t108 + t156 * t109 + t180;
t70 = t124 * t155 + (-t100 + t201) * t160 + t183;
t69 = t160 * t99 + (-t124 - t225) * t156 + t177;
t68 = t156 * t100 + (-t113 - t99) * t155 + t178;
t67 = t114 * t92 + t125 * t155 - t128 * t82 + (t111 + t201) * t160 + t183;
t66 = t160 * t110 - t115 * t92 + t128 * t81 + (-t125 - t225) * t156 + t177;
t65 = -t156 * t111 - t114 * t81 + t115 * t82 + (-t110 - t113) * t155 + t178;
t1 = m(1) * (t129 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(2) * (t85 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(3) * (t74 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(4) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + t115 * ((t116 * t77 + t117 * t79 - t75 * t214) * t115 + (t116 * t78 + t117 * t80 - t214 * t76) * t114 + (t116 * t88 + t117 * t89 - t214 * t87) * t128) / 0.2e1 + t114 * ((t118 * t77 + t119 * t79 + t213 * t75) * t115 + (t118 * t78 + t119 * t80 + t76 * t213) * t114 + (t118 * t88 + t119 * t89 + t213 * t87) * t128) / 0.2e1 + t128 * ((t76 * t114 + t75 * t115 + t87 * t128) * t158 + ((-t171 * t77 + t174 * t79) * t115 + (-t171 * t78 + t174 * t80) * t114 + (-t171 * t88 + t174 * t89) * t128) * t159) / 0.2e1 + (t233 * t173 - t240 * t176) * t155 / 0.2e1 + (t240 * t173 + t233 * t176) * t156 / 0.2e1 + ((t173 * t239 + t237 * t176 + Icges(1,4)) * V_base(5) + (t173 * t238 + t176 * t236 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t237 * t173 - t176 * t239 + Icges(1,2)) * V_base(5) + (t173 * t236 - t176 * t238 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t104 * t172 + t106 * t175 - t158 * t95 + t159 * t97) * t156 + (-t105 * t172 + t107 * t175 - t158 * t96 + t159 * t98) * t155 + (-t122 * t158 + t123 * t159 - t139 * t172 + t144 * t175 + Icges(3,1) + Icges(2,3)) * t160) * t160 / 0.2e1 + t160 * V_base(5) * (t247 * t173 - t246 * t176) + t160 * V_base(4) * (t246 * t173 + t247 * t176) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
