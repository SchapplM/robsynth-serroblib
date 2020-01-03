% Calculate kinetic energy for
% S5RPPRP2
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:02
% EndTime: 2019-12-31 17:49:04
% DurationCPUTime: 2.17s
% Computational Cost: add. (1100->219), mult. (911->279), div. (0->0), fcn. (691->8), ass. (0->119)
t246 = Icges(5,4) - Icges(6,5);
t245 = Icges(5,1) + Icges(6,1);
t244 = Icges(5,2) + Icges(6,3);
t161 = pkin(8) + qJ(4);
t155 = cos(t161);
t243 = t246 * t155;
t153 = sin(t161);
t242 = t246 * t153;
t241 = Icges(6,4) + Icges(5,5);
t240 = Icges(5,6) - Icges(6,6);
t239 = t244 * t153 - t243;
t238 = t245 * t155 - t242;
t237 = rSges(6,1) + pkin(4);
t236 = rSges(6,3) + qJ(5);
t162 = qJ(1) + pkin(7);
t154 = sin(t162);
t156 = cos(t162);
t235 = t239 * t154 + t240 * t156;
t234 = -t240 * t154 + t239 * t156;
t233 = t238 * t154 - t241 * t156;
t232 = t241 * t154 + t238 * t156;
t231 = Icges(6,2) + Icges(5,3);
t230 = -t244 * t155 - t242;
t229 = t245 * t153 + t243;
t228 = -t240 * t153 + t241 * t155;
t227 = t236 * t153 + t237 * t155;
t166 = sin(qJ(1));
t167 = cos(qJ(1));
t226 = Icges(2,5) * t166 + Icges(3,5) * t154 + Icges(2,6) * t167 + Icges(3,6) * t156;
t225 = Icges(2,5) * t167 + Icges(3,5) * t156 - Icges(2,6) * t166 - Icges(3,6) * t154;
t136 = -qJD(4) * t156 + V_base(5);
t137 = qJD(4) * t154 + V_base(4);
t157 = V_base(6) + qJD(1);
t224 = (t230 * t153 + t229 * t155) * t157 + (t234 * t153 + t232 * t155) * t137 + (t235 * t153 + t233 * t155) * t136;
t223 = (t241 * t153 + t240 * t155) * t157 + (t231 * t154 + t228 * t156) * t137 + (t228 * t154 - t231 * t156) * t136;
t219 = pkin(1) * t166;
t218 = pkin(1) * t167;
t163 = sin(pkin(8));
t217 = pkin(3) * t163;
t164 = cos(pkin(8));
t216 = pkin(3) * t164;
t215 = -pkin(5) - qJ(2);
t214 = -rSges(6,2) * t156 + t227 * t154;
t213 = rSges(6,2) * t154 + t227 * t156;
t212 = Icges(2,4) * t166;
t211 = Icges(3,4) * t154;
t210 = Icges(4,4) * t163;
t209 = Icges(4,4) * t164;
t203 = t237 * t153 - t236 * t155;
t202 = qJD(5) * t153;
t201 = t157 * t218 + V_base(2);
t200 = V_base(5) * pkin(5) + V_base(1);
t126 = pkin(2) * t154 - qJ(3) * t156;
t197 = -t126 - t219;
t128 = pkin(2) * t156 + qJ(3) * t154;
t196 = -t128 - t218;
t195 = V_base(5) * qJ(2) + t200;
t194 = V_base(4) * t219 + qJD(2) + V_base(3);
t78 = -pkin(6) * t156 + t154 * t216;
t193 = t197 - t78;
t192 = qJD(3) * t154 + t195;
t191 = V_base(4) * t126 + t194;
t190 = rSges(4,1) * t164 - rSges(4,2) * t163;
t189 = rSges(5,1) * t155 - rSges(5,2) * t153;
t186 = Icges(4,1) * t164 - t210;
t183 = -Icges(4,2) * t163 + t209;
t180 = Icges(4,5) * t164 - Icges(4,6) * t163;
t177 = V_base(5) * t217 + t192;
t176 = -qJD(3) * t156 + t157 * t128 + t201;
t173 = (Icges(4,5) * t163 + Icges(4,6) * t164) * t157 + (-Icges(4,3) * t156 + t154 * t180) * V_base(5) + (Icges(4,3) * t154 + t156 * t180) * V_base(4);
t79 = pkin(6) * t154 + t156 * t216;
t172 = V_base(4) * t78 + (t196 - t79) * V_base(5) + t191;
t171 = t157 * t79 + (t215 - t217) * V_base(4) + t176;
t100 = Icges(4,6) * t154 + t156 * t183;
t101 = -Icges(4,5) * t156 + t154 * t186;
t102 = Icges(4,5) * t154 + t156 * t186;
t134 = Icges(4,2) * t164 + t210;
t135 = Icges(4,1) * t163 + t209;
t99 = -Icges(4,6) * t156 + t154 * t183;
t168 = (-t100 * t163 + t102 * t164) * V_base(4) + (t101 * t164 - t163 * t99) * V_base(5) + (-t134 * t163 + t135 * t164) * t157;
t159 = Icges(2,4) * t167;
t151 = Icges(3,4) * t156;
t146 = rSges(2,1) * t167 - t166 * rSges(2,2);
t145 = t166 * rSges(2,1) + rSges(2,2) * t167;
t144 = Icges(2,1) * t167 - t212;
t143 = Icges(2,1) * t166 + t159;
t142 = -Icges(2,2) * t166 + t159;
t141 = Icges(2,2) * t167 + t212;
t138 = rSges(4,1) * t163 + rSges(4,2) * t164;
t132 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t131 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t130 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = rSges(3,1) * t156 - rSges(3,2) * t154;
t127 = rSges(3,1) * t154 + rSges(3,2) * t156;
t125 = rSges(5,1) * t153 + rSges(5,2) * t155;
t122 = Icges(3,1) * t156 - t211;
t121 = Icges(3,1) * t154 + t151;
t118 = -Icges(3,2) * t154 + t151;
t117 = Icges(3,2) * t156 + t211;
t106 = V_base(5) * rSges(2,3) - t145 * t157 + t200;
t105 = t146 * t157 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t104 = rSges(4,3) * t154 + t156 * t190;
t103 = -rSges(4,3) * t156 + t154 * t190;
t96 = t145 * V_base(4) - t146 * V_base(5) + V_base(3);
t95 = rSges(5,3) * t154 + t156 * t189;
t93 = -rSges(5,3) * t156 + t154 * t189;
t75 = V_base(5) * rSges(3,3) + (-t127 - t219) * t157 + t195;
t74 = t129 * t157 + (-rSges(3,3) + t215) * V_base(4) + t201;
t73 = V_base(4) * t127 + (-t129 - t218) * V_base(5) + t194;
t72 = t138 * V_base(5) + (-t103 + t197) * t157 + t192;
t71 = t104 * t157 + (-t138 + t215) * V_base(4) + t176;
t70 = V_base(4) * t103 + (-t104 + t196) * V_base(5) + t191;
t69 = t125 * t136 + (t193 - t93) * t157 + t177;
t68 = -t125 * t137 + t157 * t95 + t171;
t67 = -t136 * t95 + t137 * t93 + t172;
t66 = t156 * t202 + t203 * t136 + (t193 - t214) * t157 + t177;
t65 = -t137 * t203 + t154 * t202 + t157 * t213 + t171;
t64 = -qJD(5) * t155 - t136 * t213 + t137 * t214 + t172;
t1 = m(1) * (t130 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(2) * (t105 ^ 2 + t106 ^ 2 + t96 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (t224 * t154 - t223 * t156) * t136 / 0.2e1 + (t223 * t154 + t224 * t156) * t137 / 0.2e1 + (t173 * t154 + t168 * t156 + t225 * t157 + (-t117 * t154 + t121 * t156 - t166 * t141 + t143 * t167 + Icges(1,4)) * V_base(5) + (-t154 * t118 + t156 * t122 - t166 * t142 + t167 * t144 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t168 * t154 - t173 * t156 + t226 * t157 + (t156 * t117 + t154 * t121 + t167 * t141 + t166 * t143 + Icges(1,2)) * V_base(5) + (t118 * t156 + t122 * t154 + t142 * t167 + t166 * t144 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t232 * t153 - t234 * t155) * t137 + (t233 * t153 - t235 * t155) * t136 + (t101 * t163 + t164 * t99 + t226) * V_base(5) + (t100 * t164 + t102 * t163 + t225) * V_base(4) + (t164 * t134 + t163 * t135 + t229 * t153 - t230 * t155 + Icges(2,3) + Icges(3,3)) * t157) * t157 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
