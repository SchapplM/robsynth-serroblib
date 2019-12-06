% Calculate kinetic energy for
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:08
% EndTime: 2019-12-05 16:06:10
% DurationCPUTime: 1.94s
% Computational Cost: add. (1123->211), mult. (929->268), div. (0->0), fcn. (709->8), ass. (0->112)
t245 = Icges(5,4) - Icges(6,5);
t244 = Icges(5,1) + Icges(6,1);
t243 = Icges(5,2) + Icges(6,3);
t162 = qJ(3) + pkin(8);
t156 = cos(t162);
t242 = t245 * t156;
t154 = sin(t162);
t241 = t245 * t154;
t240 = Icges(6,4) + Icges(5,5);
t239 = Icges(5,6) - Icges(6,6);
t238 = t243 * t154 - t242;
t237 = t244 * t156 - t241;
t236 = rSges(6,1) + pkin(4);
t235 = rSges(6,3) + qJ(5);
t161 = pkin(7) + qJ(2);
t153 = sin(t161);
t155 = cos(t161);
t234 = t238 * t153 + t239 * t155;
t233 = -t239 * t153 + t238 * t155;
t232 = t237 * t153 - t240 * t155;
t231 = t240 * t153 + t237 * t155;
t230 = -t243 * t156 - t241;
t229 = t244 * t154 + t242;
t228 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t166 = sin(qJ(3));
t167 = cos(qJ(3));
t227 = Icges(4,5) * t167 - Icges(4,6) * t166 - t239 * t154 + t240 * t156;
t226 = t235 * t154 + t236 * t156;
t207 = Icges(4,4) * t167;
t185 = -Icges(4,2) * t166 + t207;
t100 = Icges(4,6) * t153 + t185 * t155;
t208 = Icges(4,4) * t166;
t188 = Icges(4,1) * t167 - t208;
t101 = -Icges(4,5) * t155 + t188 * t153;
t102 = Icges(4,5) * t153 + t188 * t155;
t140 = -qJD(3) * t155 + V_base(5);
t141 = qJD(3) * t153 + V_base(4);
t145 = Icges(4,2) * t167 + t208;
t146 = Icges(4,1) * t166 + t207;
t158 = V_base(6) + qJD(2);
t99 = -Icges(4,6) * t155 + t185 * t153;
t223 = (-t145 * t166 + t146 * t167 + t230 * t154 + t229 * t156) * t158 + (-t100 * t166 + t102 * t167 + t233 * t154 + t231 * t156) * t141 + (t101 * t167 + t234 * t154 + t232 * t156 - t166 * t99) * t140;
t222 = (Icges(4,5) * t166 + Icges(4,6) * t167 + t240 * t154 + t239 * t156) * t158 + (t228 * t153 + t227 * t155) * t141 + (t227 * t153 - t228 * t155) * t140;
t164 = cos(pkin(7));
t218 = pkin(1) * t164;
t217 = pkin(3) * t166;
t216 = pkin(3) * t167;
t215 = -pkin(5) - qJ(1);
t213 = -rSges(6,2) * t155 + t226 * t153;
t212 = rSges(6,2) * t153 + t226 * t155;
t129 = pkin(2) * t153 - pkin(6) * t155;
t78 = -qJ(4) * t155 + t216 * t153;
t211 = -t129 - t78;
t163 = sin(pkin(7));
t210 = Icges(2,4) * t163;
t209 = Icges(3,4) * t153;
t202 = t236 * t154 - t235 * t156;
t201 = qJD(5) * t154;
t194 = pkin(1) * V_base(6);
t200 = t164 * t194 + V_base(2);
t199 = V_base(5) * qJ(1) + V_base(1);
t195 = qJD(1) + V_base(3);
t193 = V_base(4) * t163 * pkin(1) + t195;
t192 = rSges(4,1) * t167 - rSges(4,2) * t166;
t191 = rSges(5,1) * t156 - rSges(5,2) * t154;
t176 = V_base(5) * pkin(5) - t163 * t194 + t199;
t130 = pkin(2) * t155 + pkin(6) * t153;
t175 = t158 * t130 + t215 * V_base(4) + t200;
t174 = qJD(4) * t153 + t140 * t217 + t176;
t173 = V_base(4) * t129 + (-t130 - t218) * V_base(5) + t193;
t172 = t141 * t78 + t173;
t79 = qJ(4) * t153 + t216 * t155;
t171 = -qJD(4) * t155 + t158 * t79 + t175;
t157 = Icges(2,4) * t164;
t151 = Icges(3,4) * t155;
t147 = t166 * rSges(4,1) + rSges(4,2) * t167;
t143 = rSges(2,1) * t164 - rSges(2,2) * t163;
t142 = rSges(2,1) * t163 + rSges(2,2) * t164;
t139 = Icges(2,1) * t164 - t210;
t138 = Icges(2,1) * t163 + t157;
t137 = -Icges(2,2) * t163 + t157;
t136 = Icges(2,2) * t164 + t210;
t133 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t132 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t131 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t128 = rSges(3,1) * t155 - rSges(3,2) * t153;
t127 = rSges(5,1) * t154 + rSges(5,2) * t156;
t124 = rSges(3,1) * t153 + rSges(3,2) * t155;
t123 = Icges(3,1) * t155 - t209;
t122 = Icges(3,1) * t153 + t151;
t119 = -Icges(3,2) * t153 + t151;
t118 = Icges(3,2) * t155 + t209;
t106 = V_base(5) * rSges(2,3) - t142 * V_base(6) + t199;
t105 = t143 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t104 = t153 * rSges(4,3) + t192 * t155;
t103 = -t155 * rSges(4,3) + t192 * t153;
t96 = rSges(5,3) * t153 + t191 * t155;
t94 = -rSges(5,3) * t155 + t191 * t153;
t80 = t142 * V_base(4) - t143 * V_base(5) + t195;
t76 = V_base(5) * rSges(3,3) - t124 * t158 + t176;
t75 = t128 * t158 + (-rSges(3,3) + t215) * V_base(4) + t200;
t73 = t124 * V_base(4) + (-t128 - t218) * V_base(5) + t193;
t72 = t140 * t147 + (-t103 - t129) * t158 + t176;
t71 = t104 * t158 - t141 * t147 + t175;
t70 = t103 * t141 - t104 * t140 + t173;
t69 = t127 * t140 + (-t94 + t211) * t158 + t174;
t68 = t158 * t96 + (-t127 - t217) * t141 + t171;
t67 = t141 * t94 + (-t79 - t96) * t140 + t172;
t66 = t155 * t201 + t202 * t140 + (t211 - t213) * t158 + t174;
t65 = t153 * t201 + t212 * t158 + (-t202 - t217) * t141 + t171;
t64 = -qJD(5) * t156 + t213 * t141 + (-t79 - t212) * t140 + t172;
t1 = m(1) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(2) * (t105 ^ 2 + t106 ^ 2 + t80 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (t223 * t153 - t222 * t155) * t140 / 0.2e1 + (t222 * t153 + t223 * t155) * t141 / 0.2e1 + ((-t118 * t153 + t122 * t155 - t136 * t163 + t138 * t164 + Icges(1,4)) * V_base(5) + (-t153 * t119 + t155 * t123 - t163 * t137 + t164 * t139 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t155 * t118 + t153 * t122 + t164 * t136 + t163 * t138 + Icges(1,2)) * V_base(5) + (t119 * t155 + t123 * t153 + t137 * t164 + t139 * t163 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t100 * t167 + t166 * t102 + t231 * t154 - t233 * t156) * t141 + (t166 * t101 + t232 * t154 - t234 * t156 + t167 * t99) * t140 + (t167 * t145 + t166 * t146 + t229 * t154 - t230 * t156 + Icges(3,3)) * t158) * t158 / 0.2e1 + V_base(4) * t158 * (Icges(3,5) * t155 - Icges(3,6) * t153) + V_base(5) * t158 * (Icges(3,5) * t153 + Icges(3,6) * t155) + ((Icges(2,5) * t163 + Icges(2,6) * t164 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t164 - Icges(2,6) * t163 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
