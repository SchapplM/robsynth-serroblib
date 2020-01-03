% Calculate kinetic energy for
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:37:58
% EndTime: 2019-12-31 17:38:01
% DurationCPUTime: 2.91s
% Computational Cost: add. (914->270), mult. (1961->373), div. (0->0), fcn. (2033->8), ass. (0->126)
t272 = Icges(3,4) - Icges(4,5);
t271 = Icges(3,1) + Icges(4,1);
t270 = Icges(3,2) + Icges(4,3);
t206 = sin(qJ(2));
t269 = t272 * t206;
t208 = cos(qJ(2));
t268 = t272 * t208;
t267 = Icges(4,4) + Icges(3,5);
t266 = Icges(3,6) - Icges(4,6);
t265 = t270 * t206 - t268;
t264 = t271 * t208 - t269;
t203 = sin(pkin(7));
t204 = cos(pkin(7));
t263 = t265 * t203 + t266 * t204;
t262 = -t266 * t203 + t265 * t204;
t261 = t264 * t203 - t267 * t204;
t260 = t267 * t203 + t264 * t204;
t259 = -t270 * t208 - t269;
t258 = t271 * t206 + t268;
t257 = t266 * t206 - t267 * t208;
t256 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t194 = -qJD(2) * t204 + V_base(5);
t195 = qJD(2) * t203 + V_base(4);
t253 = (t259 * t206 + t258 * t208) * V_base(6) + (t262 * t206 + t260 * t208) * t195 + (t263 * t206 + t261 * t208) * t194;
t202 = sin(pkin(8));
t247 = cos(pkin(8));
t164 = -t208 * t202 + t206 * t247;
t153 = t164 * t203;
t163 = t206 * t202 + t208 * t247;
t154 = t163 * t203;
t155 = t164 * t204;
t156 = t163 * t204;
t252 = (Icges(5,5) * t164 - Icges(5,6) * t163 - t267 * t206 - t266 * t208) * V_base(6) + (Icges(5,5) * t156 + Icges(5,6) * t155 - t256 * t203 + t257 * t204) * t195 + (Icges(5,5) * t154 + Icges(5,6) * t153 + t257 * t203 + t256 * t204) * t194;
t249 = pkin(3) * t206;
t248 = pkin(3) * t208;
t246 = Icges(2,4) * t203;
t226 = pkin(2) * t208 + qJ(3) * t206;
t161 = t226 * t203;
t181 = pkin(1) * t203 - pkin(5) * t204;
t240 = -t161 - t181;
t162 = t226 * t204;
t168 = -t203 * qJ(4) + t204 * t248;
t239 = -t162 - t168;
t238 = qJD(3) * t206;
t237 = V_base(5) * qJ(1) + V_base(1);
t233 = qJD(1) + V_base(3);
t167 = t204 * qJ(4) + t203 * t248;
t232 = -t167 + t240;
t189 = t206 * pkin(2) - qJ(3) * t208;
t231 = -t189 - t249;
t229 = t194 * t189 + t204 * t238 + t237;
t228 = rSges(3,1) * t208 - rSges(3,2) * t206;
t227 = rSges(4,1) * t208 + rSges(4,3) * t206;
t182 = pkin(1) * t204 + pkin(5) * t203;
t219 = -V_base(4) * qJ(1) + V_base(6) * t182 + V_base(2);
t218 = V_base(4) * t181 - V_base(5) * t182 + t233;
t217 = V_base(6) * t162 + t203 * t238 + t219;
t216 = -qJD(4) * t203 + t194 * t249 + t229;
t213 = qJD(4) * t204 + V_base(6) * t168 + t217;
t212 = -qJD(3) * t208 + t195 * t161 + t218;
t211 = t195 * t167 + t212;
t207 = cos(qJ(5));
t205 = sin(qJ(5));
t200 = Icges(2,4) * t204;
t191 = t206 * rSges(3,1) + rSges(3,2) * t208;
t190 = t206 * rSges(4,1) - rSges(4,3) * t208;
t180 = rSges(2,1) * t204 - rSges(2,2) * t203;
t179 = rSges(2,1) * t203 + rSges(2,2) * t204;
t178 = Icges(2,1) * t204 - t246;
t177 = Icges(2,1) * t203 + t200;
t176 = -Icges(2,2) * t203 + t200;
t175 = Icges(2,2) * t204 + t246;
t172 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t171 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t170 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t157 = qJD(5) * t163 + V_base(6);
t152 = t203 * rSges(3,3) + t204 * t228;
t151 = t203 * rSges(4,2) + t204 * t227;
t150 = -t204 * rSges(3,3) + t203 * t228;
t149 = -t204 * rSges(4,2) + t203 * t227;
t134 = V_base(5) * rSges(2,3) - t179 * V_base(6) + t237;
t133 = t180 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t132 = t156 * t207 - t203 * t205;
t131 = -t156 * t205 - t203 * t207;
t130 = t154 * t207 + t204 * t205;
t129 = -t154 * t205 + t204 * t207;
t128 = -qJD(5) * t155 + t195;
t127 = -qJD(5) * t153 + t194;
t126 = pkin(4) * t164 + pkin(6) * t163;
t125 = rSges(5,1) * t164 - rSges(5,2) * t163;
t124 = t179 * V_base(4) - t180 * V_base(5) + t233;
t123 = Icges(5,1) * t164 - Icges(5,4) * t163;
t122 = Icges(5,4) * t164 - Icges(5,2) * t163;
t120 = pkin(4) * t156 - pkin(6) * t155;
t119 = pkin(4) * t154 - pkin(6) * t153;
t118 = rSges(5,1) * t156 + rSges(5,2) * t155 - rSges(5,3) * t203;
t117 = rSges(5,1) * t154 + rSges(5,2) * t153 + rSges(5,3) * t204;
t116 = Icges(5,1) * t156 + Icges(5,4) * t155 - Icges(5,5) * t203;
t115 = Icges(5,1) * t154 + Icges(5,4) * t153 + Icges(5,5) * t204;
t114 = Icges(5,4) * t156 + Icges(5,2) * t155 - Icges(5,6) * t203;
t113 = Icges(5,4) * t154 + Icges(5,2) * t153 + Icges(5,6) * t204;
t110 = rSges(6,3) * t163 + (rSges(6,1) * t207 - rSges(6,2) * t205) * t164;
t109 = Icges(6,5) * t163 + (Icges(6,1) * t207 - Icges(6,4) * t205) * t164;
t108 = Icges(6,6) * t163 + (Icges(6,4) * t207 - Icges(6,2) * t205) * t164;
t107 = Icges(6,3) * t163 + (Icges(6,5) * t207 - Icges(6,6) * t205) * t164;
t106 = t191 * t194 + (-t150 - t181) * V_base(6) + t237;
t105 = t152 * V_base(6) - t191 * t195 + t219;
t104 = rSges(6,1) * t132 + rSges(6,2) * t131 - rSges(6,3) * t155;
t103 = rSges(6,1) * t130 + rSges(6,2) * t129 - rSges(6,3) * t153;
t102 = Icges(6,1) * t132 + Icges(6,4) * t131 - Icges(6,5) * t155;
t101 = Icges(6,1) * t130 + Icges(6,4) * t129 - Icges(6,5) * t153;
t100 = Icges(6,4) * t132 + Icges(6,2) * t131 - Icges(6,6) * t155;
t99 = Icges(6,4) * t130 + Icges(6,2) * t129 - Icges(6,6) * t153;
t98 = Icges(6,5) * t132 + Icges(6,6) * t131 - Icges(6,3) * t155;
t97 = Icges(6,5) * t130 + Icges(6,6) * t129 - Icges(6,3) * t153;
t96 = t150 * t195 - t152 * t194 + t218;
t95 = t190 * t194 + (-t149 + t240) * V_base(6) + t229;
t94 = t151 * V_base(6) + (-t189 - t190) * t195 + t217;
t93 = t195 * t149 + (-t151 - t162) * t194 + t212;
t92 = t125 * t194 + (-t117 + t232) * V_base(6) + t216;
t91 = t118 * V_base(6) + (-t125 + t231) * t195 + t213;
t90 = t195 * t117 + (-t118 + t239) * t194 + t211;
t89 = -t103 * t157 + t110 * t127 + t126 * t194 + (-t119 + t232) * V_base(6) + t216;
t88 = t104 * t157 - t110 * t128 + t120 * V_base(6) + (-t126 + t231) * t195 + t213;
t87 = t128 * t103 - t127 * t104 + t195 * t119 + (-t120 + t239) * t194 + t211;
t1 = m(1) * (t170 ^ 2 + t171 ^ 2 + t172 ^ 2) / 0.2e1 + m(2) * (t124 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(3) * (t105 ^ 2 + t106 ^ 2 + t96 ^ 2) / 0.2e1 + m(4) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + t128 * ((t131 * t100 + t132 * t102 - t155 * t98) * t128 + (t101 * t132 + t131 * t99 - t155 * t97) * t127 + (-t107 * t155 + t108 * t131 + t109 * t132) * t157) / 0.2e1 + t127 * ((t100 * t129 + t102 * t130 - t153 * t98) * t128 + (t130 * t101 + t129 * t99 - t153 * t97) * t127 + (-t107 * t153 + t108 * t129 + t109 * t130) * t157) / 0.2e1 + t157 * ((t107 * t157 + t97 * t127 + t98 * t128) * t163 + ((-t100 * t205 + t102 * t207) * t128 + (t101 * t207 - t205 * t99) * t127 + (-t108 * t205 + t109 * t207) * t157) * t164) / 0.2e1 + ((-t175 * t203 + t177 * t204 + Icges(1,4)) * V_base(5) + (-t203 * t176 + t204 * t178 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t204 * t175 + t203 * t177 + Icges(1,2)) * V_base(5) + (t176 * t204 + t178 * t203 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t114 * t153 + t116 * t154) * t195 + (t153 * t113 + t154 * t115) * t194 + (t122 * t153 + t123 * t154) * V_base(6) + t252 * t204 + t253 * t203) * t194 / 0.2e1 + ((t155 * t114 + t156 * t116) * t195 + (t113 * t155 + t115 * t156) * t194 + (t122 * t155 + t123 * t156) * V_base(6) + t253 * t204 - t252 * t203) * t195 / 0.2e1 + ((-t114 * t163 + t116 * t164 + t260 * t206 - t262 * t208) * t195 + (-t113 * t163 + t115 * t164 + t261 * t206 - t263 * t208) * t194 + (-t163 * t122 + t164 * t123 + t258 * t206 - t259 * t208 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t204 - Icges(2,6) * t203 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t203 + Icges(2,6) * t204 + Icges(1,6));
T = t1;
