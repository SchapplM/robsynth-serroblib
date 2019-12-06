% Calculate kinetic energy for
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:13
% EndTime: 2019-12-05 15:26:16
% DurationCPUTime: 3.06s
% Computational Cost: add. (1059->257), mult. (1346->351), div. (0->0), fcn. (1184->8), ass. (0->133)
t269 = Icges(4,4) + Icges(5,6);
t268 = Icges(4,1) + Icges(5,2);
t267 = -Icges(4,2) - Icges(5,3);
t184 = qJ(2) + pkin(8);
t181 = cos(t184);
t266 = t269 * t181;
t180 = sin(t184);
t265 = t269 * t180;
t264 = Icges(5,4) - Icges(4,5);
t263 = Icges(5,5) - Icges(4,6);
t262 = t267 * t180 + t266;
t261 = t268 * t181 - t265;
t185 = sin(pkin(7));
t186 = cos(pkin(7));
t260 = t262 * t185 + t263 * t186;
t259 = -t263 * t185 + t262 * t186;
t258 = t261 * t185 + t264 * t186;
t257 = -t264 * t185 + t261 * t186;
t256 = t267 * t181 - t265;
t255 = t268 * t180 + t266;
t254 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t189 = sin(qJ(2));
t191 = cos(qJ(2));
t253 = Icges(3,5) * t191 - Icges(3,6) * t189 + t263 * t180 - t264 * t181;
t240 = Icges(3,4) * t191;
t211 = -Icges(3,2) * t189 + t240;
t127 = -Icges(3,6) * t186 + t185 * t211;
t128 = Icges(3,6) * t185 + t186 * t211;
t241 = Icges(3,4) * t189;
t213 = Icges(3,1) * t191 - t241;
t129 = -Icges(3,5) * t186 + t185 * t213;
t130 = Icges(3,5) * t185 + t186 * t213;
t172 = Icges(3,2) * t191 + t241;
t173 = Icges(3,1) * t189 + t240;
t175 = -qJD(2) * t186 + V_base(5);
t176 = qJD(2) * t185 + V_base(4);
t250 = (-t172 * t189 + t173 * t191 + t256 * t180 + t255 * t181) * V_base(6) + (-t128 * t189 + t130 * t191 - t259 * t180 + t257 * t181) * t176 + (-t127 * t189 + t129 * t191 - t260 * t180 + t258 * t181) * t175;
t249 = (Icges(3,5) * t189 + Icges(3,6) * t191 - t264 * t180 - t263 * t181) * V_base(6) + (t254 * t185 + t253 * t186) * t176 + (t253 * t185 - t254 * t186) * t175;
t246 = pkin(2) * t189;
t245 = pkin(6) * t180;
t244 = pkin(2) * t191;
t242 = Icges(2,4) * t185;
t235 = t181 * t185;
t234 = t181 * t186;
t188 = sin(qJ(5));
t233 = t185 * t188;
t190 = cos(qJ(5));
t232 = t185 * t190;
t231 = t186 * t188;
t230 = t186 * t190;
t103 = -qJ(3) * t186 + t185 * t244;
t169 = pkin(1) * t185 - pkin(5) * t186;
t229 = -t103 - t169;
t104 = qJ(3) * t185 + t186 * t244;
t214 = pkin(3) * t181 + qJ(4) * t180;
t134 = t214 * t186;
t228 = -t104 - t134;
t227 = qJD(4) * t180;
t226 = qJD(5) * t181;
t225 = V_base(5) * qJ(1) + V_base(1);
t221 = qJD(1) + V_base(3);
t133 = t214 * t185;
t220 = -t133 + t229;
t149 = pkin(3) * t180 - qJ(4) * t181;
t219 = -t149 - t246;
t218 = qJD(3) * t185 + t175 * t246 + t225;
t217 = rSges(3,1) * t191 - rSges(3,2) * t189;
t216 = rSges(4,1) * t181 - rSges(4,2) * t180;
t215 = -rSges(5,2) * t181 + rSges(5,3) * t180;
t170 = pkin(1) * t186 + pkin(5) * t185;
t204 = -V_base(4) * qJ(1) + V_base(6) * t170 + V_base(2);
t203 = t175 * t149 + t186 * t227 + t218;
t202 = V_base(4) * t169 - t170 * V_base(5) + t221;
t201 = t176 * t103 + t202;
t197 = -qJD(3) * t186 + V_base(6) * t104 + t204;
t196 = V_base(6) * t134 + t185 * t227 + t197;
t195 = -qJD(4) * t181 + t176 * t133 + t201;
t182 = Icges(2,4) * t186;
t174 = t189 * rSges(3,1) + rSges(3,2) * t191;
t166 = rSges(2,1) * t186 - rSges(2,2) * t185;
t165 = rSges(2,1) * t185 + rSges(2,2) * t186;
t164 = qJD(5) * t180 + V_base(6);
t163 = Icges(2,1) * t186 - t242;
t162 = Icges(2,1) * t185 + t182;
t161 = -Icges(2,2) * t185 + t182;
t160 = Icges(2,2) * t186 + t242;
t157 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t156 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t155 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t151 = rSges(4,1) * t180 + rSges(4,2) * t181;
t150 = -rSges(5,2) * t180 - rSges(5,3) * t181;
t142 = -pkin(4) * t186 + pkin(6) * t235;
t141 = pkin(4) * t185 + pkin(6) * t234;
t140 = t180 * t233 - t230;
t139 = t180 * t232 + t231;
t138 = t180 * t231 + t232;
t137 = t180 * t230 - t233;
t136 = t186 * t226 + t176;
t135 = t185 * t226 + t175;
t132 = t185 * rSges(3,3) + t186 * t217;
t131 = -t186 * rSges(3,3) + t185 * t217;
t122 = -rSges(5,1) * t186 + t185 * t215;
t121 = rSges(5,1) * t185 + t186 * t215;
t120 = rSges(4,3) * t185 + t186 * t216;
t119 = -rSges(4,3) * t186 + t185 * t216;
t118 = V_base(5) * rSges(2,3) - t165 * V_base(6) + t225;
t117 = t166 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t102 = rSges(6,3) * t180 + (-rSges(6,1) * t188 - rSges(6,2) * t190) * t181;
t101 = Icges(6,5) * t180 + (-Icges(6,1) * t188 - Icges(6,4) * t190) * t181;
t100 = Icges(6,6) * t180 + (-Icges(6,4) * t188 - Icges(6,2) * t190) * t181;
t99 = Icges(6,3) * t180 + (-Icges(6,5) * t188 - Icges(6,6) * t190) * t181;
t96 = t165 * V_base(4) - t166 * V_base(5) + t221;
t94 = rSges(6,1) * t140 + rSges(6,2) * t139 + rSges(6,3) * t235;
t93 = rSges(6,1) * t138 + rSges(6,2) * t137 + rSges(6,3) * t234;
t92 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t235;
t91 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t234;
t90 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t235;
t89 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t234;
t88 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t235;
t87 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t234;
t86 = t174 * t175 + (-t131 - t169) * V_base(6) + t225;
t85 = t132 * V_base(6) - t174 * t176 + t204;
t84 = t131 * t176 - t132 * t175 + t202;
t83 = t151 * t175 + (-t119 + t229) * V_base(6) + t218;
t82 = t120 * V_base(6) + (-t151 - t246) * t176 + t197;
t81 = t119 * t176 + (-t104 - t120) * t175 + t201;
t80 = t150 * t175 + (-t122 + t220) * V_base(6) + t203;
t79 = t121 * V_base(6) + (-t150 + t219) * t176 + t196;
t78 = t122 * t176 + (-t121 + t228) * t175 + t195;
t77 = t175 * t245 + t102 * t135 - t164 * t94 + (-t142 + t220) * V_base(6) + t203;
t76 = -t102 * t136 + t141 * V_base(6) + t164 * t93 + (t219 - t245) * t176 + t196;
t75 = -t135 * t93 + t136 * t94 + t142 * t176 + (-t141 + t228) * t175 + t195;
t1 = m(1) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(2) * (t117 ^ 2 + t118 ^ 2 + t96 ^ 2) / 0.2e1 + m(3) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(4) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + t136 * ((t137 * t89 + t138 * t91 + t87 * t234) * t136 + (t137 * t90 + t138 * t92 + t234 * t88) * t135 + (t100 * t137 + t101 * t138 + t234 * t99) * t164) / 0.2e1 + t135 * ((t139 * t89 + t140 * t91 + t235 * t87) * t136 + (t139 * t90 + t140 * t92 + t88 * t235) * t135 + (t100 * t139 + t101 * t140 + t235 * t99) * t164) / 0.2e1 + t164 * ((t88 * t135 + t87 * t136 + t99 * t164) * t180 + ((-t188 * t91 - t190 * t89) * t136 + (-t188 * t92 - t190 * t90) * t135 + (-t100 * t190 - t101 * t188) * t164) * t181) / 0.2e1 + ((-t160 * t185 + t162 * t186 + Icges(1,4)) * V_base(5) + (-t161 * t185 + t163 * t186 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t160 * t186 + t162 * t185 + Icges(1,2)) * V_base(5) + (t161 * t186 + t163 * t185 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t250 * t185 - t249 * t186) * t175 / 0.2e1 + (t249 * t185 + t250 * t186) * t176 / 0.2e1 + ((t128 * t191 + t189 * t130 + t257 * t180 + t259 * t181) * t176 + (t127 * t191 + t189 * t129 + t258 * t180 + t260 * t181) * t175 + (t172 * t191 + t189 * t173 + t255 * t180 - t256 * t181 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t186 - Icges(2,6) * t185 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t185 + Icges(2,6) * t186 + Icges(1,6));
T = t1;
