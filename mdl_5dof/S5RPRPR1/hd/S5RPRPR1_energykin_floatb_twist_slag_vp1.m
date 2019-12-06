% Calculate kinetic energy for
% S5RPRPR1
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:46:59
% EndTime: 2019-12-05 17:47:01
% DurationCPUTime: 1.88s
% Computational Cost: add. (834->231), mult. (948->307), div. (0->0), fcn. (730->8), ass. (0->122)
t246 = Icges(2,4) + Icges(3,6);
t245 = Icges(2,1) + Icges(3,2);
t244 = -Icges(3,4) + Icges(2,5);
t243 = Icges(3,5) - Icges(2,6);
t242 = Icges(2,2) + Icges(3,3);
t241 = Icges(4,3) + Icges(5,3);
t162 = qJ(3) + pkin(8);
t149 = sin(t162);
t150 = cos(t162);
t164 = sin(qJ(3));
t166 = cos(qJ(3));
t240 = Icges(4,5) * t164 + Icges(5,5) * t149 + Icges(4,6) * t166 + Icges(5,6) * t150;
t167 = cos(qJ(1));
t239 = t246 * t167;
t165 = sin(qJ(1));
t238 = t246 * t165;
t212 = Icges(5,4) * t150;
t110 = -Icges(5,2) * t149 + t212;
t213 = Icges(5,4) * t149;
t111 = Icges(5,1) * t150 - t213;
t214 = Icges(4,4) * t166;
t128 = -Icges(4,2) * t164 + t214;
t215 = Icges(4,4) * t164;
t133 = Icges(4,1) * t166 - t215;
t144 = qJD(3) * t165 + V_base(5);
t145 = qJD(3) * t167 + V_base(4);
t152 = V_base(6) + qJD(1);
t184 = Icges(5,2) * t150 + t213;
t87 = Icges(5,6) * t167 + t165 * t184;
t88 = Icges(5,6) * t165 - t167 * t184;
t187 = Icges(5,1) * t149 + t212;
t89 = Icges(5,5) * t167 + t165 * t187;
t90 = Icges(5,5) * t165 - t167 * t187;
t185 = Icges(4,2) * t166 + t215;
t96 = Icges(4,6) * t167 + t165 * t185;
t97 = Icges(4,6) * t165 - t167 * t185;
t188 = Icges(4,1) * t164 + t214;
t98 = Icges(4,5) * t167 + t165 * t188;
t99 = Icges(4,5) * t165 - t167 * t188;
t237 = t144 * (t149 * t90 + t150 * t88 + t164 * t99 + t166 * t97) + t145 * (t149 * t89 + t150 * t87 + t164 * t98 + t166 * t96) + t152 * (t110 * t150 + t111 * t149 + t128 * t166 + t133 * t164);
t236 = -t242 * t167 - t238;
t235 = t242 * t165 - t239;
t234 = t245 * t165 + t239;
t233 = t245 * t167 - t238;
t230 = (Icges(4,5) * t166 + Icges(5,5) * t150 - Icges(4,6) * t164 - Icges(5,6) * t149) * t152 + (t240 * t165 + t241 * t167) * t145 + (t241 * t165 - t240 * t167) * t144;
t151 = qJ(5) + t162;
t147 = sin(t151);
t148 = cos(t151);
t210 = Icges(6,4) * t148;
t105 = -Icges(6,2) * t147 + t210;
t211 = Icges(6,4) * t147;
t106 = Icges(6,1) * t148 - t211;
t113 = qJD(5) * t165 + t144;
t114 = qJD(5) * t167 + t145;
t183 = Icges(6,2) * t148 + t211;
t76 = Icges(6,6) * t167 + t165 * t183;
t77 = Icges(6,6) * t165 - t167 * t183;
t186 = Icges(6,1) * t147 + t210;
t78 = Icges(6,5) * t167 + t165 * t186;
t79 = Icges(6,5) * t165 - t167 * t186;
t226 = (t147 * t78 + t148 * t76) * t114 + (t147 * t79 + t148 * t77) * t113 + (t105 * t148 + t106 * t147) * t152;
t222 = pkin(3) * t164;
t221 = pkin(3) * t166;
t220 = pkin(4) * t150;
t219 = pkin(6) * t167;
t218 = t165 * pkin(6);
t136 = t165 * pkin(1) - qJ(2) * t167;
t206 = V_base(4) * t136 + V_base(3);
t205 = V_base(5) * pkin(5) + V_base(1);
t202 = pkin(4) * t149;
t201 = -t136 - t218;
t200 = qJD(2) * t165 + t205;
t102 = qJ(4) * t165 - t167 * t222;
t199 = -t102 + t201;
t198 = V_base(5) * pkin(2) + t200;
t197 = rSges(4,1) * t164 + rSges(4,2) * t166;
t196 = rSges(5,1) * t149 + rSges(5,2) * t150;
t195 = rSges(6,1) * t147 + rSges(6,2) * t148;
t180 = Icges(6,5) * t147 + Icges(6,6) * t148;
t140 = pkin(1) * t167 + t165 * qJ(2);
t176 = -qJD(2) * t167 + t152 * t140 + V_base(2);
t175 = qJD(4) * t167 + t144 * t221 + t198;
t174 = (Icges(6,5) * t148 - Icges(6,6) * t147) * t152 + (Icges(6,3) * t165 - t167 * t180) * t113 + (Icges(6,3) * t167 + t165 * t180) * t114;
t171 = V_base(4) * t218 + (-t140 - t219) * V_base(5) + t206;
t170 = t152 * t219 + (-pkin(2) - pkin(5)) * V_base(4) + t176;
t169 = t145 * t102 + t171;
t103 = qJ(4) * t167 + t165 * t222;
t168 = qJD(4) * t165 + t152 * t103 + t170;
t142 = rSges(2,1) * t167 - t165 * rSges(2,2);
t141 = -rSges(3,2) * t167 + t165 * rSges(3,3);
t139 = rSges(4,1) * t166 - rSges(4,2) * t164;
t138 = t165 * rSges(2,1) + rSges(2,2) * t167;
t137 = -t165 * rSges(3,2) - rSges(3,3) * t167;
t119 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t118 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t117 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t112 = rSges(5,1) * t150 - rSges(5,2) * t149;
t107 = rSges(6,1) * t148 - rSges(6,2) * t147;
t101 = t165 * rSges(4,3) - t167 * t197;
t100 = rSges(4,3) * t167 + t165 * t197;
t92 = t165 * rSges(5,3) - t167 * t196;
t91 = rSges(5,3) * t167 + t165 * t196;
t84 = t165 * rSges(6,3) - t167 * t195;
t83 = rSges(6,3) * t167 + t165 * t195;
t82 = V_base(5) * rSges(2,3) - t138 * t152 + t205;
t81 = t142 * t152 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t73 = t138 * V_base(4) - t142 * V_base(5) + V_base(3);
t72 = pkin(7) * t167 + t165 * t202;
t71 = pkin(7) * t165 - t167 * t202;
t70 = V_base(5) * rSges(3,1) + (-t136 - t137) * t152 + t200;
t69 = t152 * t141 + (-rSges(3,1) - pkin(5)) * V_base(4) + t176;
t68 = t137 * V_base(4) + (-t140 - t141) * V_base(5) + t206;
t67 = t139 * t144 + (-t101 + t201) * t152 + t198;
t66 = t152 * t100 - t145 * t139 + t170;
t65 = -t144 * t100 + t145 * t101 + t171;
t64 = t112 * t144 + (t199 - t92) * t152 + t175;
t63 = t152 * t91 + (-t112 - t221) * t145 + t168;
t62 = t145 * t92 + (-t103 - t91) * t144 + t169;
t61 = t144 * t220 + t107 * t113 + (t199 - t71 - t84) * t152 + t175;
t60 = -t114 * t107 + (t72 + t83) * t152 + (-t220 - t221) * t145 + t168;
t59 = -t113 * t83 + t114 * t84 + t145 * t71 + (-t103 - t72) * t144 + t169;
t1 = m(1) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(2) * (t73 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(3) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + t114 * (t226 * t165 + t174 * t167) / 0.2e1 + t113 * (t174 * t165 - t226 * t167) / 0.2e1 + (t230 * t165 - t167 * t237) * t144 / 0.2e1 + (t165 * t237 + t230 * t167) * t145 / 0.2e1 + ((t165 * t236 + t234 * t167 + Icges(1,4)) * V_base(5) + (t165 * t235 + t167 * t233 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t234 * t165 - t167 * t236 + Icges(1,2)) * V_base(5) + (t165 * t233 - t167 * t235 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t147 * t76 + t148 * t78) * t114 + (-t147 * t77 + t148 * t79) * t113 + (-t149 * t87 + t150 * t89 - t164 * t96 + t166 * t98) * t145 + (-t149 * t88 + t150 * t90 - t164 * t97 + t166 * t99) * t144 + (-t147 * t105 + t106 * t148 - t110 * t149 + t111 * t150 - t128 * t164 + t133 * t166 + Icges(3,1) + Icges(2,3)) * t152) * t152 / 0.2e1 + t152 * V_base(5) * (t244 * t165 - t243 * t167) + t152 * V_base(4) * (t243 * t165 + t244 * t167) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
