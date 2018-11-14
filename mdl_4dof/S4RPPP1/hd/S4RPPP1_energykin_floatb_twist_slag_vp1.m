% Calculate kinetic energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RPPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:23
% EndTime: 2018-11-14 13:45:24
% DurationCPUTime: 1.51s
% Computational Cost: add. (1556->186), mult. (1752->231), div. (0->0), fcn. (1481->10), ass. (0->93)
t237 = Icges(3,1) + Icges(4,2) + Icges(5,3);
t236 = Icges(4,1) + Icges(5,1) + Icges(3,3);
t235 = Icges(3,4) + Icges(4,6) - Icges(5,6);
t234 = Icges(4,4) - Icges(3,5) - Icges(5,5);
t233 = Icges(5,4) + Icges(4,5) - Icges(3,6);
t232 = Icges(3,2) + Icges(5,2) + Icges(4,3);
t231 = rSges(5,1) + pkin(3);
t230 = rSges(5,3) + qJ(4);
t181 = V_base(6) + qJD(1);
t219 = t181 / 0.2e1;
t229 = V_base(4) / 0.2e1;
t184 = sin(pkin(6));
t188 = sin(qJ(1));
t189 = cos(qJ(1));
t208 = pkin(4) - pkin(6);
t197 = cos(t208) / 0.2e1;
t207 = pkin(4) + pkin(6);
t201 = cos(t207);
t194 = t197 + t201 / 0.2e1;
t151 = t189 * t184 + t188 * t194;
t186 = cos(pkin(6));
t196 = sin(t207) / 0.2e1;
t200 = sin(t208);
t193 = t196 - t200 / 0.2e1;
t152 = t186 * t189 - t188 * t193;
t185 = sin(pkin(4));
t216 = t185 * t188;
t228 = -t235 * t151 + t237 * t152 - t234 * t216;
t149 = t184 * t188 - t189 * t194;
t150 = t188 * t186 + t189 * t193;
t215 = t185 * t189;
t227 = -t235 * t149 + t237 * t150 + t234 * t215;
t226 = t232 * t151 - t235 * t152 + t233 * t216;
t225 = t232 * t149 - t235 * t150 - t233 * t215;
t224 = t233 * t151 - t234 * t152 + t236 * t216;
t223 = t233 * t149 - t234 * t150 - t236 * t215;
t159 = t196 + t200 / 0.2e1;
t160 = t197 - t201 / 0.2e1;
t187 = cos(pkin(4));
t222 = t235 * t159 + t237 * t160 - t234 * t187;
t221 = -t232 * t159 - t235 * t160 + t233 * t187;
t220 = -t233 * t159 - t234 * t160 + t236 * t187;
t218 = Icges(2,4) * t188;
t217 = qJ(2) * t187;
t214 = rSges(5,2) * t151 + t230 * t152 + t231 * t216;
t213 = t149 * rSges(5,2) + t230 * t150 - t231 * t215;
t126 = pkin(2) * t150 + qJ(3) * t149;
t161 = t188 * pkin(1) - qJ(2) * t215;
t212 = -t126 - t161;
t127 = pkin(2) * t152 + qJ(3) * t151;
t162 = pkin(1) * t189 + qJ(2) * t216;
t211 = -t127 - t162;
t210 = -rSges(5,2) * t159 + t230 * t160 + t231 * t187;
t209 = qJD(2) * t185;
t206 = V_base(5) * pkin(5) + V_base(1);
t203 = -pkin(5) - t217;
t202 = qJD(2) * t187 + V_base(4) * t161 + V_base(3);
t143 = pkin(2) * t160 - qJ(3) * t159;
t199 = -t143 + t203;
t198 = t188 * t209 + V_base(5) * t217 + t206;
t195 = -qJD(3) * t159 + V_base(4) * t126 + t202;
t192 = t181 * t162 - t189 * t209 + V_base(2);
t191 = qJD(3) * t151 + V_base(5) * t143 + t198;
t190 = qJD(3) * t149 + t181 * t127 + t192;
t182 = Icges(2,4) * t189;
t173 = rSges(2,1) * t189 - t188 * rSges(2,2);
t172 = t188 * rSges(2,1) + rSges(2,2) * t189;
t171 = Icges(2,1) * t189 - t218;
t170 = Icges(2,1) * t188 + t182;
t169 = -Icges(2,2) * t188 + t182;
t168 = Icges(2,2) * t189 + t218;
t165 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t164 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t163 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t146 = V_base(5) * rSges(2,3) - t172 * t181 + t206;
t145 = t173 * t181 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t144 = t172 * V_base(4) - t173 * V_base(5) + V_base(3);
t139 = rSges(3,1) * t160 + rSges(3,2) * t159 + rSges(3,3) * t187;
t138 = rSges(4,1) * t187 - rSges(4,2) * t160 - rSges(4,3) * t159;
t123 = rSges(3,1) * t152 - rSges(3,2) * t151 + rSges(3,3) * t216;
t122 = t150 * rSges(3,1) - t149 * rSges(3,2) - rSges(3,3) * t215;
t121 = -rSges(4,1) * t215 - t150 * rSges(4,2) + t149 * rSges(4,3);
t119 = rSges(4,1) * t216 - rSges(4,2) * t152 + rSges(4,3) * t151;
t99 = t139 * V_base(5) + (-t122 - t161) * t181 + t198;
t98 = t181 * t123 + (-t139 + t203) * V_base(4) + t192;
t97 = t122 * V_base(4) + (-t123 - t162) * V_base(5) + t202;
t96 = t138 * V_base(5) + (-t121 + t212) * t181 + t191;
t95 = t181 * t119 + (-t138 + t199) * V_base(4) + t190;
t94 = t121 * V_base(4) + (-t119 + t211) * V_base(5) + t195;
t93 = qJD(4) * t152 + t210 * V_base(5) + (t212 - t213) * t181 + t191;
t92 = qJD(4) * t150 + t214 * t181 + (t199 - t210) * V_base(4) + t190;
t91 = qJD(4) * t160 + t213 * V_base(4) + (t211 - t214) * V_base(5) + t195;
t1 = m(1) * (t163 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(3) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + (Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (Icges(1,1) * t229 + Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6)) * V_base(4) + ((-t225 * t159 + t227 * t160 + t223 * t187) * V_base(5) + (-t226 * t159 + t228 * t160 + t224 * t187) * V_base(4)) * t219 + (V_base(4) * (Icges(2,5) * t189 - Icges(2,6) * t188) + V_base(5) * (Icges(2,5) * t188 + Icges(2,6) * t189) + (-t221 * t159 + t222 * t160 + t220 * t187 + Icges(2,3)) * t219) * t181 + ((t221 * t151 + t222 * t152 + t220 * t216) * t181 + (t225 * t151 + t227 * t152 - t188 * t168 + t170 * t189 + t223 * t216) * V_base(5) + (t226 * t151 + t228 * t152 - t188 * t169 + t189 * t171 + t224 * t216) * V_base(4)) * t229 + ((t221 * t149 + t222 * t150 - t220 * t215) * t181 + (t225 * t149 + t227 * t150 + t189 * t168 + t188 * t170 - t223 * t215) * V_base(5) + (t226 * t149 + t228 * t150 + t169 * t189 + t188 * t171 - t224 * t215) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
