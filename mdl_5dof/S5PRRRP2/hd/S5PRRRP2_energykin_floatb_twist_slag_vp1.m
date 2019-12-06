% Calculate kinetic energy for
% S5PRRRP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:34
% EndTime: 2019-12-05 16:41:36
% DurationCPUTime: 2.03s
% Computational Cost: add. (1025->198), mult. (751->253), div. (0->0), fcn. (529->8), ass. (0->105)
t228 = Icges(5,4) - Icges(6,5);
t227 = Icges(5,1) + Icges(6,1);
t226 = Icges(5,2) + Icges(6,3);
t158 = cos(qJ(4));
t225 = t228 * t158;
t157 = sin(qJ(4));
t224 = t228 * t157;
t223 = Icges(6,4) + Icges(5,5);
t222 = Icges(5,6) - Icges(6,6);
t221 = t226 * t157 - t225;
t220 = t227 * t158 - t224;
t219 = rSges(6,1) + pkin(4);
t218 = rSges(6,3) + qJ(5);
t154 = pkin(8) + qJ(2);
t149 = qJ(3) + t154;
t143 = sin(t149);
t144 = cos(t149);
t217 = t221 * t143 + t222 * t144;
t216 = -t222 * t143 + t221 * t144;
t215 = t220 * t143 - t223 * t144;
t214 = t223 * t143 + t220 * t144;
t213 = Icges(6,2) + Icges(5,3);
t212 = -t226 * t158 - t224;
t211 = t227 * t157 + t225;
t210 = -t222 * t157 + t223 * t158;
t209 = t218 * t157 + t219 * t158;
t115 = -qJD(4) * t144 + V_base(5);
t116 = qJD(4) * t143 + V_base(4);
t150 = V_base(6) + qJD(2);
t145 = qJD(3) + t150;
t206 = (t212 * t157 + t211 * t158) * t145 + (t216 * t157 + t214 * t158) * t116 + (t217 * t157 + t215 * t158) * t115;
t205 = (t223 * t157 + t222 * t158) * t145 + (t213 * t143 + t210 * t144) * t116 + (t210 * t143 - t213 * t144) * t115;
t156 = cos(pkin(8));
t201 = pkin(1) * t156;
t200 = pkin(2) * t150;
t199 = -pkin(5) - qJ(1);
t198 = -t144 * rSges(6,2) + t143 * t209;
t197 = t143 * rSges(6,2) + t144 * t209;
t155 = sin(pkin(8));
t196 = Icges(2,4) * t155;
t146 = sin(t154);
t195 = Icges(3,4) * t146;
t194 = Icges(4,4) * t143;
t189 = t219 * t157 - t218 * t158;
t188 = qJD(5) * t157;
t187 = -pkin(6) + t199;
t180 = pkin(1) * V_base(6);
t186 = t156 * t180 + V_base(2);
t185 = V_base(5) * qJ(1) + V_base(1);
t181 = qJD(1) + V_base(3);
t147 = cos(t154);
t179 = t147 * t200 + t186;
t178 = V_base(4) * t155 * pkin(1) + t181;
t177 = -pkin(2) * t147 - t201;
t176 = V_base(4) * pkin(2) * t146 + t178;
t175 = rSges(5,1) * t158 - rSges(5,2) * t157;
t164 = V_base(5) * pkin(5) - t155 * t180 + t185;
t106 = pkin(3) * t144 + pkin(7) * t143;
t163 = t145 * t106 + t187 * V_base(4) + t179;
t162 = V_base(5) * pkin(6) - t146 * t200 + t164;
t105 = pkin(3) * t143 - pkin(7) * t144;
t161 = V_base(4) * t105 + (-t106 + t177) * V_base(5) + t176;
t148 = Icges(2,4) * t156;
t142 = Icges(3,4) * t147;
t139 = Icges(4,4) * t144;
t137 = t157 * rSges(5,1) + rSges(5,2) * t158;
t128 = rSges(2,1) * t156 - rSges(2,2) * t155;
t127 = rSges(2,1) * t155 + rSges(2,2) * t156;
t126 = Icges(2,1) * t156 - t196;
t125 = Icges(2,1) * t155 + t148;
t124 = -Icges(2,2) * t155 + t148;
t123 = Icges(2,2) * t156 + t196;
t119 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t118 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t117 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t114 = rSges(3,1) * t147 - rSges(3,2) * t146;
t113 = rSges(3,1) * t146 + rSges(3,2) * t147;
t112 = Icges(3,1) * t147 - t195;
t111 = Icges(3,1) * t146 + t142;
t110 = -Icges(3,2) * t146 + t142;
t109 = Icges(3,2) * t147 + t195;
t104 = rSges(4,1) * t144 - rSges(4,2) * t143;
t103 = rSges(4,1) * t143 + rSges(4,2) * t144;
t102 = Icges(4,1) * t144 - t194;
t101 = Icges(4,1) * t143 + t139;
t100 = -Icges(4,2) * t143 + t139;
t99 = Icges(4,2) * t144 + t194;
t92 = V_base(5) * rSges(2,3) - t127 * V_base(6) + t185;
t91 = t128 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t90 = t143 * rSges(5,3) + t144 * t175;
t88 = -t144 * rSges(5,3) + t143 * t175;
t74 = t127 * V_base(4) - t128 * V_base(5) + t181;
t73 = V_base(5) * rSges(3,3) - t113 * t150 + t164;
t72 = t114 * t150 + (-rSges(3,3) + t199) * V_base(4) + t186;
t71 = t113 * V_base(4) + (-t114 - t201) * V_base(5) + t178;
t70 = V_base(5) * rSges(4,3) - t103 * t145 + t162;
t69 = t104 * t145 + (-rSges(4,3) + t187) * V_base(4) + t179;
t68 = t103 * V_base(4) + (-t104 + t177) * V_base(5) + t176;
t67 = t115 * t137 + (-t105 - t88) * t145 + t162;
t66 = -t116 * t137 + t145 * t90 + t163;
t65 = -t115 * t90 + t116 * t88 + t161;
t64 = t144 * t188 + t189 * t115 + (-t105 - t198) * t145 + t162;
t63 = -t116 * t189 + t143 * t188 + t145 * t197 + t163;
t62 = -qJD(5) * t158 - t115 * t197 + t116 * t198 + t161;
t1 = m(1) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(2) * (t74 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(3) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(4) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + (t206 * t143 - t205 * t144) * t115 / 0.2e1 + (t205 * t143 + t206 * t144) * t116 / 0.2e1 + ((t214 * t157 - t216 * t158) * t116 + (t215 * t157 - t217 * t158) * t115 + (t211 * t157 - t212 * t158 + Icges(4,3)) * t145) * t145 / 0.2e1 + ((t101 * t144 - t109 * t146 + t111 * t147 - t123 * t155 + t125 * t156 - t143 * t99 + Icges(1,4)) * V_base(5) + (-t143 * t100 + t144 * t102 - t146 * t110 + t147 * t112 - t155 * t124 + t156 * t126 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t143 * t101 + t147 * t109 + t146 * t111 + t156 * t123 + t155 * t125 + t144 * t99 + Icges(1,2)) * V_base(5) + (t100 * t144 + t102 * t143 + t110 * t147 + t112 * t146 + t124 * t156 + t126 * t155 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t145 * (Icges(4,5) * t144 - Icges(4,6) * t143) + V_base(5) * t145 * (Icges(4,5) * t143 + Icges(4,6) * t144) + ((Icges(2,5) * t155 + Icges(2,6) * t156 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t156 - Icges(2,6) * t155 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + ((Icges(3,5) * t146 + Icges(3,6) * t147) * V_base(5) + (Icges(3,5) * t147 - Icges(3,6) * t146) * V_base(4) + Icges(3,3) * t150 / 0.2e1) * t150;
T = t1;
