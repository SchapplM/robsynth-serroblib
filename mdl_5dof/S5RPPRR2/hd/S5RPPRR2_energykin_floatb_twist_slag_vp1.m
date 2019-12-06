% Calculate kinetic energy for
% S5RPPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:30
% EndTime: 2019-12-05 17:39:32
% DurationCPUTime: 2.20s
% Computational Cost: add. (816->237), mult. (930->318), div. (0->0), fcn. (712->8), ass. (0->128)
t240 = Icges(2,4) + Icges(3,6);
t239 = Icges(2,1) + Icges(3,2);
t238 = -Icges(3,4) + Icges(2,5);
t237 = Icges(3,5) - Icges(2,6);
t236 = Icges(2,2) + Icges(3,3);
t167 = cos(qJ(1));
t235 = t240 * t167;
t166 = sin(qJ(1));
t234 = t240 * t166;
t233 = t239 * t166 + t235;
t232 = t239 * t167 - t234;
t231 = t238 * t166 - t237 * t167;
t230 = t237 * t166 + t238 * t167;
t163 = sin(pkin(8));
t164 = cos(pkin(8));
t219 = Icges(4,4) * t163;
t185 = Icges(4,2) * t164 + t219;
t97 = Icges(4,6) * t166 - t167 * t185;
t218 = Icges(4,4) * t164;
t188 = Icges(4,1) * t163 + t218;
t99 = Icges(4,5) * t166 - t167 * t188;
t229 = t163 * t99 + t164 * t97 - t236 * t167 - t234;
t96 = Icges(4,6) * t167 + t166 * t185;
t98 = Icges(4,5) * t167 + t166 * t188;
t228 = t163 * t98 + t164 * t96 + t236 * t166 - t235;
t162 = pkin(8) + qJ(4);
t151 = qJ(5) + t162;
t147 = sin(t151);
t148 = cos(t151);
t214 = Icges(6,4) * t148;
t105 = -Icges(6,2) * t147 + t214;
t215 = Icges(6,4) * t147;
t106 = Icges(6,1) * t148 - t215;
t143 = qJD(4) * t166 + V_base(5);
t113 = qJD(5) * t166 + t143;
t144 = qJD(4) * t167 + V_base(4);
t114 = qJD(5) * t167 + t144;
t152 = V_base(6) + qJD(1);
t183 = Icges(6,2) * t148 + t215;
t76 = Icges(6,6) * t167 + t166 * t183;
t77 = Icges(6,6) * t166 - t167 * t183;
t186 = Icges(6,1) * t147 + t214;
t78 = Icges(6,5) * t167 + t166 * t186;
t79 = Icges(6,5) * t166 - t167 * t186;
t227 = (t147 * t78 + t148 * t76) * t114 + (t147 * t79 + t148 * t77) * t113 + (t105 * t148 + t106 * t147) * t152;
t149 = sin(t162);
t150 = cos(t162);
t216 = Icges(5,4) * t150;
t110 = -Icges(5,2) * t149 + t216;
t217 = Icges(5,4) * t149;
t111 = Icges(5,1) * t150 - t217;
t184 = Icges(5,2) * t150 + t217;
t86 = Icges(5,6) * t167 + t166 * t184;
t87 = Icges(5,6) * t166 - t167 * t184;
t187 = Icges(5,1) * t149 + t216;
t88 = Icges(5,5) * t167 + t166 * t187;
t89 = Icges(5,5) * t166 - t167 * t187;
t226 = (t149 * t88 + t150 * t86) * t144 + (t149 * t89 + t150 * t87) * t143 + (t110 * t150 + t111 * t149) * t152;
t225 = -pkin(2) - pkin(5);
t223 = pkin(3) * t163;
t222 = pkin(3) * t164;
t221 = pkin(4) * t150;
t211 = qJ(3) * t167;
t210 = t166 * qJ(3);
t136 = t166 * pkin(1) - qJ(2) * t167;
t207 = V_base(4) * t136 + V_base(3);
t206 = V_base(5) * pkin(5) + V_base(1);
t203 = pkin(4) * t149;
t202 = V_base(4) * t210 + t207;
t201 = qJD(2) * t166 + t206;
t200 = -t136 - t210;
t139 = pkin(1) * t167 + t166 * qJ(2);
t199 = -t139 - t211;
t102 = pkin(6) * t166 - t167 * t223;
t198 = -t102 + t200;
t197 = rSges(4,1) * t163 + rSges(4,2) * t164;
t196 = rSges(5,1) * t149 + rSges(5,2) * t150;
t195 = rSges(6,1) * t147 + rSges(6,2) * t148;
t182 = Icges(4,5) * t163 + Icges(4,6) * t164;
t181 = Icges(5,5) * t149 + Icges(5,6) * t150;
t180 = Icges(6,5) * t147 + Icges(6,6) * t148;
t121 = -Icges(4,2) * t163 + t218;
t122 = Icges(4,1) * t164 - t219;
t177 = t121 * t164 + t122 * t163;
t176 = -qJD(2) * t167 + t152 * t139 + V_base(2);
t175 = V_base(5) * pkin(2) + qJD(3) * t167 + t201;
t174 = V_base(5) * t222 + t175;
t173 = qJD(3) * t166 + t152 * t211 + t176;
t172 = (Icges(6,5) * t148 - Icges(6,6) * t147) * t152 + t113 * (Icges(6,3) * t166 - t167 * t180) + t114 * (Icges(6,3) * t167 + t166 * t180);
t171 = (Icges(5,5) * t150 - Icges(5,6) * t149) * t152 + t143 * (Icges(5,3) * t166 - t167 * t181) + t144 * (Icges(5,3) * t167 + t166 * t181);
t170 = (Icges(4,5) * t164 - Icges(4,6) * t163) * t152 + (Icges(4,3) * t167 + t166 * t182) * V_base(4) + (Icges(4,3) * t166 - t167 * t182) * V_base(5);
t103 = pkin(6) * t167 + t166 * t223;
t169 = V_base(4) * t102 + (-t103 + t199) * V_base(5) + t202;
t168 = t152 * t103 + (-t222 + t225) * V_base(4) + t173;
t141 = rSges(2,1) * t167 - t166 * rSges(2,2);
t140 = -rSges(3,2) * t167 + t166 * rSges(3,3);
t138 = t166 * rSges(2,1) + rSges(2,2) * t167;
t137 = -t166 * rSges(3,2) - rSges(3,3) * t167;
t123 = rSges(4,1) * t164 - rSges(4,2) * t163;
t119 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t118 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t117 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t112 = rSges(5,1) * t150 - rSges(5,2) * t149;
t107 = rSges(6,1) * t148 - rSges(6,2) * t147;
t101 = t166 * rSges(4,3) - t167 * t197;
t100 = rSges(4,3) * t167 + t166 * t197;
t92 = t166 * rSges(5,3) - t167 * t196;
t91 = rSges(5,3) * t167 + t166 * t196;
t83 = t166 * rSges(6,3) - t167 * t195;
t82 = rSges(6,3) * t167 + t166 * t195;
t81 = V_base(5) * rSges(2,3) - t138 * t152 + t206;
t80 = t141 * t152 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t73 = t138 * V_base(4) - t141 * V_base(5) + V_base(3);
t72 = pkin(7) * t167 + t166 * t203;
t71 = pkin(7) * t166 - t167 * t203;
t70 = V_base(5) * rSges(3,1) + (-t136 - t137) * t152 + t201;
t69 = t152 * t140 + (-rSges(3,1) - pkin(5)) * V_base(4) + t176;
t68 = t137 * V_base(4) + (-t139 - t140) * V_base(5) + t207;
t67 = t123 * V_base(5) + (-t101 + t200) * t152 + t175;
t66 = t152 * t100 + (-t123 + t225) * V_base(4) + t173;
t65 = V_base(4) * t101 + (-t100 + t199) * V_base(5) + t202;
t64 = t112 * t143 + (t198 - t92) * t152 + t174;
t63 = -t144 * t112 + t152 * t91 + t168;
t62 = -t143 * t91 + t144 * t92 + t169;
t61 = t143 * t221 + t107 * t113 + (t198 - t71 - t83) * t152 + t174;
t60 = -t144 * t221 - t114 * t107 + (t72 + t82) * t152 + t168;
t59 = -t113 * t82 + t114 * t83 - t143 * t72 + t144 * t71 + t169;
t1 = m(1) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(2) * (t73 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(3) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + t144 * (t226 * t166 + t171 * t167) / 0.2e1 + t143 * (t171 * t166 - t226 * t167) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + t114 * (t227 * t166 + t172 * t167) / 0.2e1 + t113 * (t172 * t166 - t227 * t167) / 0.2e1 + (t170 * t167 + (t177 * t166 + t230) * t152 + (t229 * t166 + t233 * t167 + Icges(1,4)) * V_base(5) + (t228 * t166 + t232 * t167 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t170 * t166 + (-t177 * t167 + t231) * t152 + (t233 * t166 - t229 * t167 + Icges(1,2)) * V_base(5) + (t232 * t166 - t228 * t167 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t149 * t86 + t150 * t88) * t144 + (-t149 * t87 + t150 * t89) * t143 + (-t147 * t76 + t148 * t78) * t114 + (-t147 * t77 + t148 * t79) * t113 + (-t163 * t97 + t164 * t99 + t231) * V_base(5) + (-t163 * t96 + t164 * t98 + t230) * V_base(4) + (-t147 * t105 + t106 * t148 - t110 * t149 + t111 * t150 - t121 * t163 + t122 * t164 + Icges(3,1) + Icges(2,3)) * t152) * t152 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
