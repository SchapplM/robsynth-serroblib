% Calculate kinetic energy for
% S5RPRRP5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:38
% EndTime: 2019-12-31 18:40:40
% DurationCPUTime: 1.92s
% Computational Cost: add. (1038->199), mult. (751->250), div. (0->0), fcn. (529->8), ass. (0->105)
t228 = Icges(5,4) - Icges(6,5);
t227 = Icges(5,1) + Icges(6,1);
t226 = Icges(5,2) + Icges(6,3);
t157 = cos(qJ(4));
t225 = t228 * t157;
t155 = sin(qJ(4));
t224 = t228 * t155;
t223 = Icges(6,4) + Icges(5,5);
t222 = Icges(5,6) - Icges(6,6);
t221 = t226 * t155 - t225;
t220 = t227 * t157 - t224;
t219 = rSges(6,1) + pkin(4);
t218 = rSges(6,3) + qJ(5);
t154 = qJ(1) + pkin(8);
t148 = qJ(3) + t154;
t143 = sin(t148);
t144 = cos(t148);
t217 = t221 * t143 + t222 * t144;
t216 = -t222 * t143 + t221 * t144;
t215 = t220 * t143 - t223 * t144;
t214 = t223 * t143 + t220 * t144;
t213 = Icges(6,2) + Icges(5,3);
t212 = -t226 * t157 - t224;
t211 = t227 * t155 + t225;
t210 = -t222 * t155 + t223 * t157;
t209 = t218 * t155 + t219 * t157;
t115 = -qJD(4) * t144 + V_base(5);
t116 = qJD(4) * t143 + V_base(4);
t149 = V_base(6) + qJD(1);
t145 = qJD(3) + t149;
t206 = (t212 * t155 + t211 * t157) * t145 + (t216 * t155 + t214 * t157) * t116 + (t217 * t155 + t215 * t157) * t115;
t205 = (t223 * t155 + t222 * t157) * t145 + (t213 * t143 + t210 * t144) * t116 + (t210 * t143 - t213 * t144) * t115;
t156 = sin(qJ(1));
t200 = pkin(1) * t156;
t158 = cos(qJ(1));
t199 = pkin(1) * t158;
t146 = sin(t154);
t198 = pkin(2) * t146;
t147 = cos(t154);
t197 = pkin(2) * t147;
t196 = -pkin(5) - qJ(2);
t195 = -rSges(6,2) * t144 + t209 * t143;
t194 = rSges(6,2) * t143 + t209 * t144;
t193 = Icges(2,4) * t156;
t192 = Icges(3,4) * t146;
t191 = Icges(4,4) * t143;
t186 = t219 * t155 - t218 * t157;
t185 = qJD(5) * t155;
t184 = -pkin(6) + t196;
t183 = t149 * t199 + V_base(2);
t182 = V_base(5) * pkin(5) + V_base(1);
t179 = t149 * t197 + t183;
t178 = V_base(5) * qJ(2) + t182;
t177 = V_base(4) * t200 + qJD(2) + V_base(3);
t176 = -t197 - t199;
t175 = V_base(4) * t198 + t177;
t174 = rSges(5,1) * t157 - rSges(5,2) * t155;
t106 = pkin(3) * t144 + pkin(7) * t143;
t163 = t145 * t106 + t184 * V_base(4) + t179;
t162 = V_base(5) * pkin(6) + (-t198 - t200) * t149 + t178;
t105 = pkin(3) * t143 - pkin(7) * t144;
t161 = V_base(4) * t105 + (-t106 + t176) * V_base(5) + t175;
t151 = Icges(2,4) * t158;
t142 = Icges(3,4) * t147;
t140 = Icges(4,4) * t144;
t137 = rSges(2,1) * t158 - t156 * rSges(2,2);
t136 = t156 * rSges(2,1) + rSges(2,2) * t158;
t135 = rSges(5,1) * t155 + rSges(5,2) * t157;
t132 = Icges(2,1) * t158 - t193;
t131 = Icges(2,1) * t156 + t151;
t128 = -Icges(2,2) * t156 + t151;
t127 = Icges(2,2) * t158 + t193;
t119 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t118 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t117 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t114 = rSges(3,1) * t147 - rSges(3,2) * t146;
t113 = rSges(3,1) * t146 + rSges(3,2) * t147;
t112 = Icges(3,1) * t147 - t192;
t111 = Icges(3,1) * t146 + t142;
t110 = -Icges(3,2) * t146 + t142;
t109 = Icges(3,2) * t147 + t192;
t104 = rSges(4,1) * t144 - rSges(4,2) * t143;
t103 = rSges(4,1) * t143 + rSges(4,2) * t144;
t102 = Icges(4,1) * t144 - t191;
t101 = Icges(4,1) * t143 + t140;
t100 = -Icges(4,2) * t143 + t140;
t99 = Icges(4,2) * t144 + t191;
t92 = V_base(5) * rSges(2,3) - t136 * t149 + t182;
t91 = t137 * t149 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t90 = t136 * V_base(4) - t137 * V_base(5) + V_base(3);
t89 = rSges(5,3) * t143 + t174 * t144;
t87 = -rSges(5,3) * t144 + t174 * t143;
t73 = V_base(5) * rSges(3,3) + (-t113 - t200) * t149 + t178;
t72 = t114 * t149 + (-rSges(3,3) + t196) * V_base(4) + t183;
t71 = V_base(4) * t113 + (-t114 - t199) * V_base(5) + t177;
t70 = V_base(5) * rSges(4,3) - t103 * t145 + t162;
t69 = t104 * t145 + (-rSges(4,3) + t184) * V_base(4) + t179;
t68 = V_base(4) * t103 + (-t104 + t176) * V_base(5) + t175;
t67 = t115 * t135 + (-t105 - t87) * t145 + t162;
t66 = -t116 * t135 + t145 * t89 + t163;
t65 = -t115 * t89 + t116 * t87 + t161;
t64 = t144 * t185 + t186 * t115 + (-t105 - t195) * t145 + t162;
t63 = -t186 * t116 + t143 * t185 + t194 * t145 + t163;
t62 = -qJD(5) * t157 - t194 * t115 + t195 * t116 + t161;
t1 = m(1) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(2) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(3) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(4) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + (t206 * t143 - t205 * t144) * t115 / 0.2e1 + (t205 * t143 + t206 * t144) * t116 / 0.2e1 + ((t214 * t155 - t216 * t157) * t116 + (t215 * t155 - t217 * t157) * t115 + (t211 * t155 - t212 * t157 + Icges(4,3)) * t145) * t145 / 0.2e1 + ((t101 * t144 - t109 * t146 + t111 * t147 - t156 * t127 + t131 * t158 - t143 * t99 + Icges(1,4)) * V_base(5) + (-t143 * t100 + t144 * t102 - t146 * t110 + t147 * t112 - t156 * t128 + t158 * t132 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t143 * t101 + t147 * t109 + t146 * t111 + t158 * t127 + t156 * t131 + t144 * t99 + Icges(1,2)) * V_base(5) + (t100 * t144 + t102 * t143 + t110 * t147 + t112 * t146 + t128 * t158 + t156 * t132 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t145 * (Icges(4,5) * t144 - Icges(4,6) * t143) + V_base(5) * t145 * (Icges(4,5) * t143 + Icges(4,6) * t144) + ((Icges(2,5) * t158 + Icges(3,5) * t147 - Icges(2,6) * t156 - Icges(3,6) * t146) * V_base(4) + (Icges(2,5) * t156 + Icges(3,5) * t146 + Icges(2,6) * t158 + Icges(3,6) * t147) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * t149) * t149 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
