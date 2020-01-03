% Calculate kinetic energy for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:00
% EndTime: 2019-12-31 16:47:02
% DurationCPUTime: 1.23s
% Computational Cost: add. (388->152), mult. (651->189), div. (0->0), fcn. (483->4), ass. (0->81)
t197 = -Icges(4,4) + Icges(5,5);
t196 = Icges(4,1) + Icges(5,1);
t195 = Icges(4,2) + Icges(5,3);
t114 = cos(qJ(3));
t194 = t197 * t114;
t112 = sin(qJ(3));
t193 = t197 * t112;
t192 = Icges(5,4) + Icges(4,5);
t191 = Icges(4,6) - Icges(5,6);
t190 = -t195 * t114 + t193;
t189 = t196 * t112 - t194;
t188 = rSges(5,1) + pkin(3);
t187 = rSges(5,3) + qJ(4);
t186 = Icges(2,4) + Icges(3,6);
t113 = sin(qJ(1));
t115 = cos(qJ(1));
t185 = t190 * t113 - t191 * t115;
t184 = -t191 * t113 - t190 * t115;
t183 = t189 * t113 + t192 * t115;
t182 = t192 * t113 - t189 * t115;
t181 = t195 * t112 + t194;
t180 = t196 * t114 + t193;
t179 = Icges(2,1) + Icges(3,2);
t178 = -Icges(3,4) + Icges(2,5);
t177 = Icges(3,5) - Icges(2,6);
t176 = Icges(2,2) + Icges(3,3);
t175 = Icges(5,2) + Icges(4,3);
t174 = t192 * t112 + t191 * t114;
t173 = t186 * t115;
t172 = t188 * t112 - t187 * t114;
t171 = t186 * t113;
t104 = qJD(3) * t113 + V_base(5);
t105 = qJD(3) * t115 + V_base(4);
t107 = V_base(6) + qJD(1);
t170 = (t182 * t112 - t184 * t114) * t104 + (t183 * t112 - t185 * t114) * t105 + (t180 * t112 - t181 * t114) * t107;
t169 = -t176 * t115 - t171;
t168 = t176 * t113 - t173;
t167 = t179 * t113 + t173;
t166 = t179 * t115 - t171;
t163 = (-t191 * t112 + t192 * t114) * t107 + (t174 * t113 + t175 * t115) * t105 + (t175 * t113 - t174 * t115) * t104;
t156 = pkin(5) * t113;
t155 = pkin(5) * t115;
t154 = rSges(5,2) * t115 + t172 * t113;
t153 = t113 * rSges(5,2) - t172 * t115;
t152 = t187 * t112 + t188 * t114;
t144 = qJD(2) * t115;
t143 = qJD(4) * t114;
t100 = pkin(1) * t115 + t113 * qJ(2);
t142 = t107 * t100 + V_base(2);
t94 = t113 * pkin(1) - qJ(2) * t115;
t141 = V_base(4) * t94 + V_base(3);
t140 = V_base(5) * pkin(4) + V_base(1);
t137 = -t94 - t156;
t136 = qJD(2) * t113 + t140;
t135 = V_base(5) * pkin(2) + t136;
t134 = rSges(4,1) * t112 + rSges(4,2) * t114;
t117 = t107 * t155 + (-pkin(2) - pkin(4)) * V_base(4) + t142;
t116 = V_base(4) * t156 + (-t100 - t155) * V_base(5) + t141;
t102 = rSges(2,1) * t115 - t113 * rSges(2,2);
t101 = -rSges(3,2) * t115 + t113 * rSges(3,3);
t99 = rSges(4,1) * t114 - rSges(4,2) * t112;
t96 = t113 * rSges(2,1) + rSges(2,2) * t115;
t95 = -t113 * rSges(3,2) - rSges(3,3) * t115;
t75 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t74 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t73 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t68 = t113 * rSges(4,3) - t115 * t134;
t66 = rSges(4,3) * t115 + t113 * t134;
t52 = V_base(5) * rSges(2,3) - t107 * t96 + t140;
t51 = t102 * t107 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t50 = -t102 * V_base(5) + t96 * V_base(4) + V_base(3);
t49 = V_base(5) * rSges(3,1) + (-t94 - t95) * t107 + t136;
t48 = -t144 + t107 * t101 + (-rSges(3,1) - pkin(4)) * V_base(4) + t142;
t47 = t95 * V_base(4) + (-t100 - t101) * V_base(5) + t141;
t46 = t104 * t99 + (t137 - t68) * t107 + t135;
t45 = -t105 * t99 + t107 * t66 + t117 - t144;
t44 = -t104 * t66 + t105 * t68 + t116;
t43 = -t113 * t143 + t152 * t104 + (t137 - t153) * t107 + t135;
t42 = (-qJD(2) + t143) * t115 + t154 * t107 - t152 * t105 + t117;
t41 = qJD(4) * t112 - t154 * t104 + t153 * t105 + t116;
t1 = m(1) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(2) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(3) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + m(4) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + m(5) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + (t163 * t113 - t170 * t115) * t104 / 0.2e1 + (t170 * t113 + t163 * t115) * t105 / 0.2e1 + ((t169 * t113 + t167 * t115 + Icges(1,4)) * V_base(5) + (t168 * t113 + t166 * t115 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t167 * t113 - t169 * t115 + Icges(1,2)) * V_base(5) + (t166 * t113 - t168 * t115 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t185 * t112 + t183 * t114) * t105 + (t184 * t112 + t182 * t114) * t104 + (t181 * t112 + t180 * t114 + Icges(3,1) + Icges(2,3)) * t107) * t107 / 0.2e1 + t107 * V_base(5) * (t178 * t113 - t177 * t115) + t107 * V_base(4) * (t177 * t113 + t178 * t115) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
