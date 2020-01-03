% Calculate kinetic energy for
% S5RPRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:37
% EndTime: 2019-12-31 18:10:39
% DurationCPUTime: 1.64s
% Computational Cost: add. (920->187), mult. (922->235), div. (0->0), fcn. (702->6), ass. (0->99)
t242 = Icges(4,4) - Icges(6,4) - Icges(5,5);
t241 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t240 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t166 = cos(qJ(3));
t239 = t242 * t166;
t164 = sin(qJ(3));
t238 = t242 * t164;
t237 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t236 = Icges(4,6) - Icges(5,6) + Icges(6,6);
t235 = t240 * t164 - t239;
t234 = t241 * t166 - t238;
t233 = rSges(6,1) + pkin(4);
t163 = qJ(1) + pkin(7);
t157 = sin(t163);
t158 = cos(t163);
t232 = t235 * t157 + t236 * t158;
t231 = -t236 * t157 + t235 * t158;
t230 = t234 * t157 - t237 * t158;
t229 = t237 * t157 + t234 * t158;
t228 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t227 = -t240 * t166 - t238;
t226 = t241 * t164 + t239;
t225 = -t236 * t164 + t237 * t166;
t224 = rSges(6,3) + qJ(5);
t223 = rSges(6,2) * t164 + t233 * t166;
t129 = -qJD(3) * t158 + V_base(5);
t130 = qJD(3) * t157 + V_base(4);
t159 = V_base(6) + qJD(1);
t220 = (t227 * t164 + t226 * t166) * t159 + (t231 * t164 + t229 * t166) * t130 + (t232 * t164 + t230 * t166) * t129;
t219 = (t237 * t164 + t236 * t166) * t159 + (t228 * t157 + t225 * t158) * t130 + (t225 * t157 - t228 * t158) * t129;
t165 = sin(qJ(1));
t215 = pkin(1) * t165;
t167 = cos(qJ(1));
t214 = pkin(1) * t167;
t212 = -pkin(5) - qJ(2);
t211 = Icges(2,4) * t165;
t210 = Icges(3,4) * t157;
t203 = t223 * t157 + t224 * t158;
t202 = -t224 * t157 + t223 * t158;
t201 = qJD(4) * t164;
t200 = t159 * t214 + V_base(2);
t199 = V_base(5) * pkin(5) + V_base(1);
t124 = pkin(2) * t157 - pkin(6) * t158;
t196 = -t124 - t215;
t195 = -rSges(6,2) * t166 + t233 * t164;
t194 = V_base(5) * qJ(2) + t199;
t193 = V_base(4) * t215 + qJD(2) + V_base(3);
t188 = pkin(3) * t166 + qJ(4) * t164;
t110 = t188 * t157;
t192 = -t110 + t196;
t191 = rSges(4,1) * t166 - rSges(4,2) * t164;
t190 = rSges(5,1) * t166 + rSges(5,3) * t164;
t148 = pkin(3) * t164 - qJ(4) * t166;
t178 = t129 * t148 + t158 * t201 + t194;
t125 = pkin(2) * t158 + pkin(6) * t157;
t174 = t159 * t125 + t212 * V_base(4) + t200;
t111 = t188 * t158;
t173 = t159 * t111 + t157 * t201 + t174;
t172 = V_base(4) * t124 + (-t125 - t214) * V_base(5) + t193;
t171 = -qJD(4) * t166 + t130 * t110 + t172;
t161 = Icges(2,4) * t167;
t156 = Icges(3,4) * t158;
t153 = rSges(2,1) * t167 - t165 * rSges(2,2);
t152 = t165 * rSges(2,1) + rSges(2,2) * t167;
t151 = rSges(4,1) * t164 + rSges(4,2) * t166;
t150 = rSges(5,1) * t164 - rSges(5,3) * t166;
t145 = Icges(2,1) * t167 - t211;
t144 = Icges(2,1) * t165 + t161;
t140 = -Icges(2,2) * t165 + t161;
t139 = Icges(2,2) * t167 + t211;
t128 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t127 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t126 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t123 = rSges(3,1) * t158 - rSges(3,2) * t157;
t122 = rSges(3,1) * t157 + rSges(3,2) * t158;
t121 = Icges(3,1) * t158 - t210;
t120 = Icges(3,1) * t157 + t156;
t119 = -Icges(3,2) * t157 + t156;
t118 = Icges(3,2) * t158 + t210;
t107 = rSges(4,3) * t157 + t158 * t191;
t106 = rSges(5,2) * t157 + t158 * t190;
t104 = -rSges(4,3) * t158 + t157 * t191;
t103 = -rSges(5,2) * t158 + t157 * t190;
t101 = V_base(5) * rSges(2,3) - t152 * t159 + t199;
t100 = t153 * t159 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t81 = t152 * V_base(4) - t153 * V_base(5) + V_base(3);
t79 = V_base(5) * rSges(3,3) + (-t122 - t215) * t159 + t194;
t78 = t123 * t159 + (-rSges(3,3) + t212) * V_base(4) + t200;
t77 = V_base(4) * t122 + (-t123 - t214) * V_base(5) + t193;
t76 = t129 * t151 + (-t104 + t196) * t159 + t194;
t75 = t107 * t159 - t130 * t151 + t174;
t74 = t130 * t104 - t129 * t107 + t172;
t73 = t129 * t150 + (-t103 + t192) * t159 + t178;
t72 = t106 * t159 + (-t148 - t150) * t130 + t173;
t71 = -qJD(5) * t157 + t195 * t129 + (t192 - t203) * t159 + t178;
t70 = qJD(5) * t158 + t202 * t159 + (-t148 - t195) * t130 + t173;
t69 = t130 * t103 + (-t106 - t111) * t129 + t171;
t68 = t203 * t130 + (-t111 - t202) * t129 + t171;
t1 = m(1) * (t126 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t101 ^ 2 + t81 ^ 2) / 0.2e1 + m(3) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(4) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(5) * (t69 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t68 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + (t220 * t157 - t219 * t158) * t129 / 0.2e1 + (t219 * t157 + t220 * t158) * t130 / 0.2e1 + ((-t118 * t157 + t120 * t158 - t165 * t139 + t144 * t167 + Icges(1,4)) * V_base(5) + (-t119 * t157 + t121 * t158 - t165 * t140 + t145 * t167 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t118 * t158 + t120 * t157 + t139 * t167 + t165 * t144 + Icges(1,2)) * V_base(5) + (t119 * t158 + t121 * t157 + t140 * t167 + t165 * t145 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t229 * t164 - t231 * t166) * t130 + (t230 * t164 - t232 * t166) * t129 + (t226 * t164 - t227 * t166 + Icges(2,3) + Icges(3,3)) * t159) * t159 / 0.2e1 + t159 * V_base(5) * (Icges(2,5) * t165 + Icges(3,5) * t157 + Icges(2,6) * t167 + Icges(3,6) * t158) + t159 * V_base(4) * (Icges(2,5) * t167 + Icges(3,5) * t158 - Icges(2,6) * t165 - Icges(3,6) * t157) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
