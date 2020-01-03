% Calculate kinetic energy for
% S5RPRPP3
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:16
% EndTime: 2019-12-31 18:12:19
% DurationCPUTime: 2.28s
% Computational Cost: add. (929->207), mult. (1105->273), div. (0->0), fcn. (887->6), ass. (0->113)
t254 = Icges(4,4) + Icges(5,6) - Icges(6,6);
t253 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t252 = -Icges(4,2) - Icges(6,2) - Icges(5,3);
t167 = pkin(7) + qJ(3);
t162 = cos(t167);
t251 = t254 * t162;
t161 = sin(t167);
t250 = t254 * t161;
t249 = Icges(5,4) - Icges(4,5) - Icges(6,5);
t248 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t247 = t253 * t162 - t250;
t246 = t252 * t161 + t251;
t245 = rSges(6,3) + qJ(5);
t171 = sin(qJ(1));
t172 = cos(qJ(1));
t244 = t246 * t171 + t248 * t172;
t243 = -t248 * t171 + t246 * t172;
t242 = t247 * t171 + t249 * t172;
t241 = -t249 * t171 + t247 * t172;
t240 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t239 = t253 * t161 + t251;
t238 = t252 * t162 - t250;
t237 = t248 * t161 - t249 * t162;
t236 = rSges(6,1) + pkin(4);
t235 = rSges(6,2) * t161 + t245 * t162;
t157 = -qJD(3) * t172 + V_base(5);
t158 = qJD(3) * t171 + V_base(4);
t163 = V_base(6) + qJD(1);
t234 = (t238 * t161 + t239 * t162) * t163 + (-t243 * t161 + t241 * t162) * t158 + (-t244 * t161 + t242 * t162) * t157;
t233 = (-t249 * t161 - t248 * t162) * t163 + (t240 * t171 + t237 * t172) * t158 + (t237 * t171 - t240 * t172) * t157;
t168 = sin(pkin(7));
t229 = pkin(2) * t168;
t169 = cos(pkin(7));
t228 = pkin(2) * t169;
t153 = t171 * pkin(1) - qJ(2) * t172;
t83 = -pkin(6) * t172 + t228 * t171;
t227 = -t153 - t83;
t226 = Icges(2,4) * t171;
t225 = Icges(3,4) * t168;
t224 = Icges(3,4) * t169;
t215 = t236 * t171 + t235 * t172;
t214 = t235 * t171 - t236 * t172;
t213 = qJD(4) * t161;
t212 = qJD(5) * t162;
t211 = V_base(4) * t153 + V_base(3);
t210 = V_base(5) * pkin(5) + V_base(1);
t199 = pkin(3) * t162 + qJ(4) * t161;
t119 = t199 * t171;
t207 = -t119 + t227;
t206 = qJD(2) * t171 + t210;
t205 = -rSges(6,2) * t162 + t245 * t161;
t204 = V_base(5) * t229 + t206;
t203 = rSges(3,1) * t169 - rSges(3,2) * t168;
t202 = rSges(4,1) * t162 - rSges(4,2) * t161;
t201 = -rSges(5,2) * t162 + rSges(5,3) * t161;
t198 = Icges(3,1) * t169 - t225;
t196 = -Icges(3,2) * t168 + t224;
t192 = Icges(3,5) * t169 - Icges(3,6) * t168;
t155 = pkin(1) * t172 + t171 * qJ(2);
t186 = -qJD(2) * t172 + t163 * t155 + V_base(2);
t133 = pkin(3) * t161 - qJ(4) * t162;
t185 = t157 * t133 + t172 * t213 + t204;
t84 = pkin(6) * t171 + t228 * t172;
t181 = V_base(4) * t83 + (-t155 - t84) * V_base(5) + t211;
t180 = (-Icges(3,3) * t172 + t192 * t171) * V_base(5) + (Icges(3,3) * t171 + t192 * t172) * V_base(4) + (Icges(3,5) * t168 + Icges(3,6) * t169) * t163;
t179 = -qJD(4) * t162 + t158 * t119 + t181;
t178 = t163 * t84 + (-pkin(5) - t229) * V_base(4) + t186;
t120 = t199 * t172;
t177 = t163 * t120 + t171 * t213 + t178;
t113 = -Icges(3,6) * t172 + t196 * t171;
t114 = Icges(3,6) * t171 + t196 * t172;
t115 = -Icges(3,5) * t172 + t198 * t171;
t116 = Icges(3,5) * t171 + t198 * t172;
t142 = Icges(3,2) * t169 + t225;
t143 = Icges(3,1) * t168 + t224;
t173 = (-t114 * t168 + t116 * t169) * V_base(4) + (-t113 * t168 + t115 * t169) * V_base(5) + (-t142 * t168 + t143 * t169) * t163;
t165 = Icges(2,4) * t172;
t156 = rSges(2,1) * t172 - t171 * rSges(2,2);
t154 = t171 * rSges(2,1) + rSges(2,2) * t172;
t150 = Icges(2,1) * t172 - t226;
t149 = Icges(2,1) * t171 + t165;
t148 = -Icges(2,2) * t171 + t165;
t147 = Icges(2,2) * t172 + t226;
t146 = Icges(2,5) * t172 - Icges(2,6) * t171;
t145 = Icges(2,5) * t171 + Icges(2,6) * t172;
t144 = rSges(3,1) * t168 + rSges(3,2) * t169;
t140 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t139 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t138 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t135 = rSges(4,1) * t161 + rSges(4,2) * t162;
t134 = -rSges(5,2) * t161 - rSges(5,3) * t162;
t118 = t171 * rSges(3,3) + t203 * t172;
t117 = -rSges(3,3) * t172 + t203 * t171;
t109 = -rSges(5,1) * t172 + t201 * t171;
t107 = t171 * rSges(5,1) + t201 * t172;
t105 = t171 * rSges(4,3) + t202 * t172;
t104 = -rSges(4,3) * t172 + t202 * t171;
t82 = V_base(5) * rSges(2,3) - t154 * t163 + t210;
t81 = t156 * t163 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t79 = t154 * V_base(4) - t156 * V_base(5) + V_base(3);
t76 = t144 * V_base(5) + (-t117 - t153) * t163 + t206;
t75 = t163 * t118 + (-pkin(5) - t144) * V_base(4) + t186;
t74 = t117 * V_base(4) + (-t118 - t155) * V_base(5) + t211;
t73 = t135 * t157 + (-t104 + t227) * t163 + t204;
t72 = t163 * t105 - t158 * t135 + t178;
t71 = t104 * t158 - t105 * t157 + t181;
t70 = t134 * t157 + (-t109 + t207) * t163 + t185;
t69 = t163 * t107 + (-t133 - t134) * t158 + t177;
t68 = t109 * t158 + (-t107 - t120) * t157 + t179;
t67 = t172 * t212 + t205 * t157 + (t207 - t214) * t163 + t185;
t66 = t171 * t212 + t215 * t163 + (-t133 - t205) * t158 + t177;
t65 = qJD(5) * t161 + t214 * t158 + (-t120 - t215) * t157 + t179;
t1 = m(1) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(2) * (t79 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(3) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + (t234 * t171 - t233 * t172) * t157 / 0.2e1 + (t233 * t171 + t234 * t172) * t158 / 0.2e1 + (t146 * t163 + t180 * t171 + t173 * t172 + (-t171 * t147 + t149 * t172 + Icges(1,4)) * V_base(5) + (-t171 * t148 + t150 * t172 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t145 * t163 + t173 * t171 - t180 * t172 + (t147 * t172 + t171 * t149 + Icges(1,2)) * V_base(5) + (t148 * t172 + t171 * t150 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t113 * t169 + t115 * t168 + t145) * V_base(5) + (t114 * t169 + t116 * t168 + t146) * V_base(4) + (t241 * t161 + t243 * t162) * t158 + (t242 * t161 + t244 * t162) * t157 + (t142 * t169 + t143 * t168 + t239 * t161 - t238 * t162 + Icges(2,3)) * t163) * t163 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
