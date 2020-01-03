% Calculate kinetic energy for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:05
% EndTime: 2019-12-31 18:16:06
% DurationCPUTime: 1.69s
% Computational Cost: add. (550->176), mult. (932->219), div. (0->0), fcn. (714->4), ass. (0->95)
t247 = -Icges(4,4) + Icges(6,4) + Icges(5,5);
t246 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t245 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t150 = cos(qJ(3));
t244 = t247 * t150;
t148 = sin(qJ(3));
t243 = t247 * t148;
t242 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t241 = -Icges(4,6) + Icges(5,6) - Icges(6,6);
t240 = -t245 * t150 + t243;
t239 = t246 * t148 - t244;
t238 = rSges(6,1) + pkin(4);
t237 = Icges(2,4) + Icges(3,6);
t236 = Icges(2,1) + Icges(3,2);
t235 = -Icges(3,4) + Icges(2,5);
t234 = Icges(3,5) - Icges(2,6);
t233 = Icges(2,2) + Icges(3,3);
t149 = sin(qJ(1));
t151 = cos(qJ(1));
t232 = t240 * t149 + t241 * t151;
t231 = t241 * t149 - t240 * t151;
t230 = t239 * t149 + t242 * t151;
t229 = t242 * t149 - t239 * t151;
t228 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t227 = t245 * t148 + t244;
t226 = t246 * t150 + t243;
t225 = -t242 * t148 + t241 * t150;
t224 = -rSges(6,3) - qJ(5);
t223 = t237 * t151;
t222 = -rSges(6,2) * t150 + t238 * t148;
t221 = t237 * t149;
t138 = qJD(3) * t149 + V_base(5);
t139 = qJD(3) * t151 + V_base(4);
t142 = V_base(6) + qJD(1);
t220 = (t229 * t148 - t231 * t150) * t138 + (t230 * t148 - t232 * t150) * t139 + (t226 * t148 - t227 * t150) * t142;
t219 = -t233 * t151 - t221;
t218 = t233 * t149 - t223;
t217 = t236 * t149 + t223;
t216 = t236 * t151 - t221;
t213 = (t241 * t148 + t242 * t150) * t142 + (-t225 * t149 + t228 * t151) * t139 + (t228 * t149 + t225 * t151) * t138;
t205 = pkin(6) * t149;
t204 = pkin(6) * t151;
t203 = t222 * t149 + t224 * t151;
t202 = t224 * t149 - t222 * t151;
t192 = qJD(4) * t150;
t127 = t149 * pkin(1) - qJ(2) * t151;
t191 = V_base(4) * t127 + V_base(3);
t190 = V_base(5) * pkin(5) + V_base(1);
t187 = rSges(6,2) * t148 + t238 * t150;
t186 = -t127 - t205;
t185 = qJD(2) * t149 + t190;
t179 = pkin(3) * t148 - qJ(4) * t150;
t99 = t179 * t151;
t184 = t186 + t99;
t183 = V_base(5) * pkin(2) + t185;
t182 = rSges(4,1) * t148 + rSges(4,2) * t150;
t181 = rSges(5,1) * t148 - rSges(5,3) * t150;
t134 = pkin(1) * t151 + t149 * qJ(2);
t160 = -qJD(2) * t151 + t142 * t134 + V_base(2);
t130 = pkin(3) * t150 + qJ(4) * t148;
t156 = t138 * t130 - t149 * t192 + t183;
t155 = V_base(4) * t205 + (-t134 - t204) * V_base(5) + t191;
t154 = t142 * t204 + (-pkin(2) - pkin(5)) * V_base(4) + t160;
t153 = qJD(4) * t148 - t139 * t99 + t155;
t98 = t179 * t149;
t152 = t142 * t98 + t151 * t192 + t154;
t136 = rSges(2,1) * t151 - t149 * rSges(2,2);
t135 = -rSges(3,2) * t151 + t149 * rSges(3,3);
t133 = rSges(4,1) * t150 - rSges(4,2) * t148;
t132 = rSges(5,1) * t150 + rSges(5,3) * t148;
t129 = t149 * rSges(2,1) + rSges(2,2) * t151;
t128 = -t149 * rSges(3,2) - rSges(3,3) * t151;
t105 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t104 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t103 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t95 = t149 * rSges(4,3) - t151 * t182;
t94 = t149 * rSges(5,2) - t151 * t181;
t92 = rSges(4,3) * t151 + t149 * t182;
t91 = rSges(5,2) * t151 + t149 * t181;
t69 = V_base(5) * rSges(2,3) - t129 * t142 + t190;
t68 = t136 * t142 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t67 = t129 * V_base(4) - t136 * V_base(5) + V_base(3);
t66 = V_base(5) * rSges(3,1) + (-t127 - t128) * t142 + t185;
t65 = t142 * t135 + (-rSges(3,1) - pkin(5)) * V_base(4) + t160;
t64 = t128 * V_base(4) + (-t134 - t135) * V_base(5) + t191;
t63 = t133 * t138 + (t186 - t95) * t142 + t183;
t62 = -t139 * t133 + t142 * t92 + t154;
t61 = -t138 * t92 + t139 * t95 + t155;
t60 = t132 * t138 + (t184 - t94) * t142 + t156;
t59 = t142 * t91 + (-t130 - t132) * t139 + t152;
t58 = t139 * t94 + (-t91 - t98) * t138 + t153;
t57 = -qJD(5) * t151 + t187 * t138 + (t184 - t202) * t142 + t156;
t56 = -qJD(5) * t149 + t203 * t142 + (-t130 - t187) * t139 + t152;
t55 = t202 * t139 + (-t98 - t203) * t138 + t153;
t1 = m(1) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(2) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(3) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(4) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(5) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + (t213 * t149 - t220 * t151) * t138 / 0.2e1 + (t220 * t149 + t213 * t151) * t139 / 0.2e1 + ((t219 * t149 + t217 * t151 + Icges(1,4)) * V_base(5) + (t218 * t149 + t216 * t151 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t217 * t149 - t219 * t151 + Icges(1,2)) * V_base(5) + (t216 * t149 - t218 * t151 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t232 * t148 + t230 * t150) * t139 + (t231 * t148 + t229 * t150) * t138 + (t227 * t148 + t226 * t150 + Icges(3,1) + Icges(2,3)) * t142) * t142 / 0.2e1 + t142 * V_base(5) * (t235 * t149 - t234 * t151) + t142 * V_base(4) * (t234 * t149 + t235 * t151) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
