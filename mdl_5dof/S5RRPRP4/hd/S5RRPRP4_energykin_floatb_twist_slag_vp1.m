% Calculate kinetic energy for
% S5RRPRP4
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:34
% EndTime: 2019-12-31 19:52:36
% DurationCPUTime: 1.60s
% Computational Cost: add. (814->190), mult. (759->233), div. (0->0), fcn. (539->6), ass. (0->103)
t239 = -Icges(5,4) + Icges(6,5);
t238 = Icges(5,1) + Icges(6,1);
t237 = Icges(5,2) + Icges(6,3);
t148 = cos(qJ(4));
t236 = t239 * t148;
t146 = sin(qJ(4));
t235 = t239 * t146;
t234 = Icges(6,4) + Icges(5,5);
t233 = Icges(5,6) - Icges(6,6);
t232 = -t148 * t237 + t235;
t231 = t146 * t238 - t236;
t230 = rSges(6,1) + pkin(4);
t229 = rSges(6,3) + qJ(5);
t228 = Icges(3,4) + Icges(4,6);
t145 = qJ(1) + qJ(2);
t139 = sin(t145);
t140 = cos(t145);
t227 = t139 * t232 - t140 * t233;
t226 = -t139 * t233 - t140 * t232;
t225 = t139 * t231 + t140 * t234;
t224 = t139 * t234 - t140 * t231;
t223 = Icges(3,1) + Icges(4,2);
t222 = -Icges(4,4) + Icges(3,5);
t221 = Icges(4,5) - Icges(3,6);
t220 = -Icges(3,2) - Icges(4,3);
t219 = Icges(6,2) + Icges(5,3);
t218 = t146 * t237 + t236;
t217 = t148 * t238 + t235;
t216 = t146 * t234 + t148 * t233;
t215 = t228 * t140;
t214 = t146 * t230 - t148 * t229;
t213 = t228 * t139;
t113 = qJD(4) * t139 + V_base(5);
t114 = qJD(4) * t140 + V_base(4);
t138 = V_base(6) + qJD(1);
t137 = qJD(2) + t138;
t212 = (t146 * t224 - t148 * t226) * t113 + (t146 * t225 - t148 * t227) * t114 + (t146 * t217 - t148 * t218) * t137;
t211 = t140 * t220 - t213;
t208 = t139 * t220 + t215;
t207 = t139 * t223 + t215;
t206 = t140 * t223 - t213;
t205 = (-t146 * t233 + t148 * t234) * t137 + (t139 * t216 + t140 * t219) * t114 + (t139 * t219 - t140 * t216) * t113;
t201 = -pkin(5) - pkin(6);
t147 = sin(qJ(1));
t197 = pkin(1) * t147;
t149 = cos(qJ(1));
t196 = pkin(1) * t149;
t195 = pkin(7) * t139;
t194 = pkin(7) * t140;
t193 = t140 * rSges(6,2) + t139 * t214;
t192 = t139 * rSges(6,2) - t140 * t214;
t191 = Icges(2,4) * t147;
t183 = t146 * t229 + t148 * t230;
t182 = qJD(3) * t140;
t181 = qJD(5) * t148;
t180 = t138 * t196 + V_base(2);
t179 = t197 * V_base(4) + V_base(3);
t178 = V_base(5) * pkin(5) + V_base(1);
t106 = pkin(2) * t140 + qJ(3) * t139;
t175 = -t106 - t196;
t103 = pkin(2) * t139 - qJ(3) * t140;
t174 = -t103 - t195;
t173 = t106 * t137 + t180;
t172 = t103 * V_base(4) + t179;
t171 = rSges(5,1) * t146 + rSges(5,2) * t148;
t156 = V_base(5) * pkin(6) - t138 * t197 + t178;
t153 = qJD(3) * t139 + t156;
t152 = V_base(5) * pkin(3) + t153;
t151 = t137 * t194 + (-pkin(3) + t201) * V_base(4) + t173;
t150 = V_base(4) * t195 + (t175 - t194) * V_base(5) + t172;
t141 = Icges(2,4) * t149;
t131 = rSges(2,1) * t149 - rSges(2,2) * t147;
t130 = rSges(5,1) * t148 - rSges(5,2) * t146;
t127 = rSges(2,1) * t147 + rSges(2,2) * t149;
t126 = Icges(2,1) * t149 - t191;
t125 = Icges(2,1) * t147 + t141;
t122 = -Icges(2,2) * t147 + t141;
t121 = Icges(2,2) * t149 + t191;
t112 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t111 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t110 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t108 = rSges(3,1) * t140 - rSges(3,2) * t139;
t107 = -rSges(4,2) * t140 + rSges(4,3) * t139;
t105 = rSges(3,1) * t139 + rSges(3,2) * t140;
t104 = -rSges(4,2) * t139 - rSges(4,3) * t140;
t86 = t139 * rSges(5,3) - t140 * t171;
t84 = rSges(5,3) * t140 + t139 * t171;
t70 = V_base(5) * rSges(2,3) - t127 * t138 + t178;
t69 = t138 * t131 + V_base(2) + (-pkin(5) - rSges(2,3)) * V_base(4);
t68 = t127 * V_base(4) - t131 * V_base(5) + V_base(3);
t67 = V_base(5) * rSges(3,3) - t105 * t137 + t156;
t66 = t137 * t108 + (-rSges(3,3) + t201) * V_base(4) + t180;
t65 = V_base(4) * t105 + (-t108 - t196) * V_base(5) + t179;
t64 = V_base(5) * rSges(4,1) + (-t103 - t104) * t137 + t153;
t63 = -t182 + t137 * t107 + (-rSges(4,1) + t201) * V_base(4) + t173;
t62 = V_base(4) * t104 + (-t107 + t175) * V_base(5) + t172;
t61 = t113 * t130 + (t174 - t86) * t137 + t152;
t60 = -t114 * t130 + t137 * t84 + t151 - t182;
t59 = -t113 * t84 + t114 * t86 + t150;
t58 = -t139 * t181 + t183 * t113 + (t174 - t192) * t137 + t152;
t57 = (-qJD(3) + t181) * t140 + t193 * t137 - t183 * t114 + t151;
t56 = qJD(5) * t146 - t113 * t193 + t114 * t192 + t150;
t1 = m(1) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(2) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(3) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(4) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(5) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + (t205 * t139 - t212 * t140) * t113 / 0.2e1 + (t212 * t139 + t205 * t140) * t114 / 0.2e1 + ((t146 * t227 + t148 * t225) * t114 + (t146 * t226 + t148 * t224) * t113 + (t218 * t146 + t217 * t148 + Icges(4,1) + Icges(3,3)) * t137) * t137 / 0.2e1 + ((-t147 * t121 + t149 * t125 + t139 * t211 + t140 * t207 + Icges(1,4)) * V_base(5) + (-t147 * t122 + t149 * t126 - t208 * t139 + t206 * t140 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t149 * t121 + t147 * t125 + t207 * t139 - t211 * t140 + Icges(1,2)) * V_base(5) + (t149 * t122 + t147 * t126 + t139 * t206 + t140 * t208 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t137 * (t139 * t222 - t140 * t221) + V_base(4) * t137 * (t139 * t221 + t140 * t222) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t147 + Icges(2,6) * t149) * V_base(5) + (Icges(2,5) * t149 - Icges(2,6) * t147) * V_base(4) + Icges(2,3) * t138 / 0.2e1) * t138;
T = t1;
