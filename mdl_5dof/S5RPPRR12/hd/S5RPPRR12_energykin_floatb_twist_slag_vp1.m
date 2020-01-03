% Calculate kinetic energy for
% S5RPPRR12
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR12_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR12_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:54
% EndTime: 2019-12-31 18:06:56
% DurationCPUTime: 2.05s
% Computational Cost: add. (869->259), mult. (1132->358), div. (0->0), fcn. (970->8), ass. (0->134)
t243 = Icges(2,4) + Icges(3,6);
t242 = Icges(2,1) + Icges(3,2);
t241 = -Icges(3,4) + Icges(2,5);
t240 = Icges(3,5) - Icges(2,6);
t239 = Icges(2,2) + Icges(3,3);
t176 = cos(qJ(1));
t238 = t243 * t176;
t174 = sin(qJ(1));
t237 = t243 * t174;
t236 = t242 * t174 + t238;
t235 = t242 * t176 - t237;
t234 = t241 * t174 - t240 * t176;
t233 = t240 * t174 + t241 * t176;
t171 = cos(pkin(8));
t170 = sin(pkin(8));
t224 = Icges(4,4) * t170;
t192 = Icges(4,2) * t171 + t224;
t105 = Icges(4,6) * t174 - t176 * t192;
t223 = Icges(4,4) * t171;
t194 = Icges(4,1) * t170 + t223;
t107 = Icges(4,5) * t174 - t176 * t194;
t232 = t105 * t171 + t107 * t170 - t239 * t176 - t237;
t104 = Icges(4,6) * t176 + t174 * t192;
t106 = Icges(4,5) * t176 + t174 * t194;
t231 = t104 * t171 + t106 * t170 + t239 * t174 - t238;
t169 = pkin(8) + qJ(4);
t158 = sin(t169);
t159 = cos(t169);
t221 = Icges(5,4) * t159;
t122 = -Icges(5,2) * t158 + t221;
t222 = Icges(5,4) * t158;
t123 = Icges(5,1) * t159 - t222;
t154 = qJD(4) * t174 + V_base(5);
t155 = qJD(4) * t176 + V_base(4);
t160 = V_base(6) + qJD(1);
t191 = Icges(5,2) * t159 + t222;
t94 = Icges(5,6) * t176 + t174 * t191;
t95 = Icges(5,6) * t174 - t176 * t191;
t193 = Icges(5,1) * t158 + t221;
t96 = Icges(5,5) * t176 + t174 * t193;
t97 = Icges(5,5) * t174 - t176 * t193;
t230 = (t158 * t96 + t159 * t94) * t155 + (t158 * t97 + t159 * t95) * t154 + (t122 * t159 + t123 * t158) * t160;
t229 = -pkin(2) - pkin(5);
t227 = pkin(3) * t170;
t226 = pkin(3) * t171;
t218 = qJ(3) * t176;
t217 = t159 * t174;
t216 = t159 * t176;
t173 = sin(qJ(5));
t215 = t173 * t176;
t214 = t174 * qJ(3);
t213 = t174 * t173;
t175 = cos(qJ(5));
t212 = t174 * t175;
t211 = t175 * t176;
t209 = qJD(5) * t159;
t147 = t174 * pkin(1) - qJ(2) * t176;
t208 = V_base(4) * t147 + V_base(3);
t207 = V_base(5) * pkin(5) + V_base(1);
t204 = V_base(4) * t214 + t208;
t203 = qJD(2) * t174 + t207;
t202 = -t147 - t214;
t150 = pkin(1) * t176 + t174 * qJ(2);
t201 = -t150 - t218;
t200 = pkin(4) * t158 - pkin(7) * t159;
t112 = pkin(6) * t174 - t176 * t227;
t199 = -t112 + t202;
t198 = rSges(4,1) * t170 + rSges(4,2) * t171;
t197 = rSges(5,1) * t158 + rSges(5,2) * t159;
t190 = Icges(4,5) * t170 + Icges(4,6) * t171;
t189 = Icges(5,5) * t158 + Icges(5,6) * t159;
t132 = -Icges(4,2) * t170 + t223;
t133 = Icges(4,1) * t171 - t224;
t185 = t132 * t171 + t133 * t170;
t184 = -qJD(2) * t176 + t160 * t150 + V_base(2);
t183 = V_base(5) * pkin(2) + qJD(3) * t176 + t203;
t182 = V_base(5) * t226 + t183;
t181 = qJD(3) * t174 + t160 * t218 + t184;
t180 = (Icges(5,5) * t159 - Icges(5,6) * t158) * t160 + t154 * (Icges(5,3) * t174 - t176 * t189) + t155 * (Icges(5,3) * t176 + t174 * t189);
t179 = (Icges(4,3) * t176 + t174 * t190) * V_base(4) + (Icges(4,3) * t174 - t176 * t190) * V_base(5) + (Icges(4,5) * t171 - Icges(4,6) * t170) * t160;
t113 = pkin(6) * t176 + t174 * t227;
t178 = V_base(4) * t112 + (-t113 + t201) * V_base(5) + t204;
t177 = t160 * t113 + (-t226 + t229) * V_base(4) + t181;
t152 = rSges(2,1) * t176 - t174 * rSges(2,2);
t151 = -rSges(3,2) * t176 + t174 * rSges(3,3);
t149 = t174 * rSges(2,1) + rSges(2,2) * t176;
t148 = -t174 * rSges(3,2) - rSges(3,3) * t176;
t134 = rSges(4,1) * t171 - rSges(4,2) * t170;
t130 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t129 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t128 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t127 = qJD(5) * t158 + t160;
t125 = pkin(4) * t159 + pkin(7) * t158;
t124 = rSges(5,1) * t159 - rSges(5,2) * t158;
t119 = -t158 * t211 + t213;
t118 = t158 * t215 + t212;
t117 = t158 * t212 + t215;
t116 = -t158 * t213 + t211;
t115 = -t174 * t209 + t155;
t114 = t176 * t209 + t154;
t111 = t200 * t176;
t110 = t200 * t174;
t109 = t174 * rSges(4,3) - t176 * t198;
t108 = rSges(4,3) * t176 + t174 * t198;
t100 = t174 * rSges(5,3) - t176 * t197;
t99 = rSges(5,3) * t176 + t174 * t197;
t91 = rSges(6,3) * t158 + (rSges(6,1) * t175 - rSges(6,2) * t173) * t159;
t90 = V_base(5) * rSges(2,3) - t149 * t160 + t207;
t89 = t152 * t160 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t88 = Icges(6,5) * t158 + (Icges(6,1) * t175 - Icges(6,4) * t173) * t159;
t87 = Icges(6,6) * t158 + (Icges(6,4) * t175 - Icges(6,2) * t173) * t159;
t86 = Icges(6,3) * t158 + (Icges(6,5) * t175 - Icges(6,6) * t173) * t159;
t85 = t149 * V_base(4) - t152 * V_base(5) + V_base(3);
t84 = V_base(5) * rSges(3,1) + (-t147 - t148) * t160 + t203;
t83 = t160 * t151 + (-rSges(3,1) - pkin(5)) * V_base(4) + t184;
t82 = t119 * rSges(6,1) + t118 * rSges(6,2) + rSges(6,3) * t216;
t81 = rSges(6,1) * t117 + rSges(6,2) * t116 - rSges(6,3) * t217;
t80 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t216;
t79 = Icges(6,1) * t117 + Icges(6,4) * t116 - Icges(6,5) * t217;
t78 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t216;
t77 = Icges(6,4) * t117 + Icges(6,2) * t116 - Icges(6,6) * t217;
t76 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t216;
t75 = Icges(6,5) * t117 + Icges(6,6) * t116 - Icges(6,3) * t217;
t74 = t148 * V_base(4) + (-t150 - t151) * V_base(5) + t208;
t73 = t134 * V_base(5) + (-t109 + t202) * t160 + t183;
t72 = t160 * t108 + (-t134 + t229) * V_base(4) + t181;
t71 = V_base(4) * t109 + (-t108 + t201) * V_base(5) + t204;
t70 = t124 * t154 + (-t100 + t199) * t160 + t182;
t69 = -t155 * t124 + t160 * t99 + t177;
t68 = t155 * t100 - t154 * t99 + t178;
t67 = t114 * t91 + t125 * t154 - t127 * t82 + (t111 + t199) * t160 + t182;
t66 = t160 * t110 - t115 * t91 - t155 * t125 + t127 * t81 + t177;
t65 = -t154 * t110 - t155 * t111 - t114 * t81 + t115 * t82 + t178;
t1 = m(1) * (t128 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(2) * (t85 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(3) * (t74 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(4) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + t155 * (t230 * t174 + t180 * t176) / 0.2e1 + t154 * (t180 * t174 - t230 * t176) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + t115 * ((t116 * t77 + t117 * t79 - t217 * t75) * t115 + (t116 * t78 + t117 * t80 - t217 * t76) * t114 + (t116 * t87 + t117 * t88 - t217 * t86) * t127) / 0.2e1 + t114 * ((t118 * t77 + t119 * t79 + t216 * t75) * t115 + (t118 * t78 + t119 * t80 + t216 * t76) * t114 + (t118 * t87 + t119 * t88 + t216 * t86) * t127) / 0.2e1 + t127 * ((t114 * t76 + t115 * t75 + t127 * t86) * t158 + ((-t173 * t77 + t175 * t79) * t115 + (-t173 * t78 + t175 * t80) * t114 + (-t173 * t87 + t175 * t88) * t127) * t159) / 0.2e1 + ((-t158 * t94 + t159 * t96) * t155 + (-t158 * t95 + t159 * t97) * t154 + (-t105 * t170 + t107 * t171 + t234) * V_base(5) + (-t104 * t170 + t106 * t171 + t233) * V_base(4) + (-t122 * t158 + t123 * t159 - t132 * t170 + t133 * t171 + Icges(3,1) + Icges(2,3)) * t160) * t160 / 0.2e1 + (t179 * t176 + (t185 * t174 + t233) * t160 + (t232 * t174 + t176 * t236 + Icges(1,4)) * V_base(5) + (t174 * t231 + t176 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t179 * t174 + (-t185 * t176 + t234) * t160 + (t174 * t236 - t232 * t176 + Icges(1,2)) * V_base(5) + (t174 * t235 - t231 * t176 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
