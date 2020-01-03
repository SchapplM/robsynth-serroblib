% Calculate kinetic energy for
% S5RPRRR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:25
% EndTime: 2019-12-31 19:12:27
% DurationCPUTime: 2.21s
% Computational Cost: add. (911->255), mult. (1174->369), div. (0->0), fcn. (1012->8), ass. (0->132)
t244 = Icges(2,4) + Icges(3,6);
t243 = Icges(2,1) + Icges(3,2);
t242 = -Icges(3,4) + Icges(2,5);
t241 = Icges(3,5) - Icges(2,6);
t240 = Icges(2,2) + Icges(3,3);
t177 = cos(qJ(1));
t239 = t244 * t177;
t174 = sin(qJ(1));
t238 = t244 * t174;
t237 = -t240 * t177 - t238;
t236 = t240 * t174 - t239;
t235 = t243 * t174 + t239;
t234 = t243 * t177 - t238;
t171 = qJ(3) + qJ(4);
t166 = sin(t171);
t167 = cos(t171);
t219 = Icges(5,4) * t167;
t122 = -Icges(5,2) * t166 + t219;
t220 = Icges(5,4) * t166;
t123 = Icges(5,1) * t167 - t220;
t157 = qJD(3) * t174 + V_base(5);
t126 = qJD(4) * t174 + t157;
t158 = qJD(3) * t177 + V_base(4);
t127 = qJD(4) * t177 + t158;
t160 = V_base(6) + qJD(1);
t193 = Icges(5,2) * t167 + t220;
t95 = Icges(5,6) * t177 + t174 * t193;
t96 = Icges(5,6) * t174 - t177 * t193;
t195 = Icges(5,1) * t166 + t219;
t97 = Icges(5,5) * t177 + t174 * t195;
t98 = Icges(5,5) * t174 - t177 * t195;
t231 = (t166 * t97 + t167 * t95) * t127 + (t166 * t98 + t167 * t96) * t126 + (t122 * t167 + t123 * t166) * t160;
t176 = cos(qJ(3));
t173 = sin(qJ(3));
t222 = Icges(4,4) * t173;
t194 = Icges(4,2) * t176 + t222;
t106 = Icges(4,6) * t177 + t174 * t194;
t107 = Icges(4,6) * t174 - t177 * t194;
t221 = Icges(4,4) * t176;
t196 = Icges(4,1) * t173 + t221;
t108 = Icges(4,5) * t177 + t174 * t196;
t109 = Icges(4,5) * t174 - t177 * t196;
t141 = -Icges(4,2) * t173 + t221;
t146 = Icges(4,1) * t176 - t222;
t230 = (t106 * t176 + t108 * t173) * t158 + (t107 * t176 + t109 * t173) * t157 + (t141 * t176 + t146 * t173) * t160;
t228 = pkin(3) * t173;
t227 = pkin(3) * t176;
t226 = t174 * pkin(6);
t225 = t177 * pkin(6);
t216 = t167 * t174;
t215 = t167 * t177;
t172 = sin(qJ(5));
t214 = t172 * t174;
t213 = t172 * t177;
t175 = cos(qJ(5));
t212 = t174 * t175;
t211 = t175 * t177;
t210 = qJD(5) * t167;
t149 = pkin(1) * t174 - qJ(2) * t177;
t209 = V_base(4) * t149 + V_base(3);
t208 = V_base(5) * pkin(5) + V_base(1);
t205 = -t149 - t226;
t204 = qJD(2) * t174 + t208;
t114 = pkin(7) * t174 - t177 * t228;
t203 = -t114 + t205;
t202 = V_base(5) * pkin(2) + t204;
t201 = pkin(4) * t166 - pkin(8) * t167;
t200 = rSges(4,1) * t173 + rSges(4,2) * t176;
t199 = rSges(5,1) * t166 + rSges(5,2) * t167;
t192 = Icges(4,5) * t173 + Icges(4,6) * t176;
t191 = Icges(5,5) * t166 + Icges(5,6) * t167;
t153 = pkin(1) * t177 + qJ(2) * t174;
t186 = -qJD(2) * t177 + t160 * t153 + V_base(2);
t185 = t157 * t227 + t202;
t184 = (Icges(5,5) * t167 - Icges(5,6) * t166) * t160 + (Icges(5,3) * t174 - t177 * t191) * t126 + (Icges(5,3) * t177 + t174 * t191) * t127;
t183 = (Icges(4,3) * t177 + t174 * t192) * t158 + (Icges(4,3) * t174 - t177 * t192) * t157 + (Icges(4,5) * t176 - Icges(4,6) * t173) * t160;
t182 = V_base(4) * t226 + (-t153 - t225) * V_base(5) + t209;
t181 = t160 * t225 + (-pkin(2) - pkin(5)) * V_base(4) + t186;
t115 = pkin(7) * t177 + t174 * t228;
t180 = t158 * t114 - t115 * t157 + t182;
t179 = t160 * t115 - t158 * t227 + t181;
t155 = rSges(2,1) * t177 - rSges(2,2) * t174;
t154 = -rSges(3,2) * t177 + rSges(3,3) * t174;
t152 = rSges(4,1) * t176 - rSges(4,2) * t173;
t151 = rSges(2,1) * t174 + rSges(2,2) * t177;
t150 = -rSges(3,2) * t174 - rSges(3,3) * t177;
t133 = qJD(5) * t166 + t160;
t132 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t131 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t130 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t125 = pkin(4) * t167 + pkin(8) * t166;
t124 = rSges(5,1) * t167 - rSges(5,2) * t166;
t119 = -t166 * t211 + t214;
t118 = t166 * t213 + t212;
t117 = t166 * t212 + t213;
t116 = -t166 * t214 + t211;
t113 = t201 * t177;
t112 = t201 * t174;
t111 = rSges(4,3) * t174 - t177 * t200;
t110 = rSges(4,3) * t177 + t174 * t200;
t103 = -t174 * t210 + t127;
t102 = t177 * t210 + t126;
t100 = rSges(5,3) * t174 - t177 * t199;
t99 = rSges(5,3) * t177 + t174 * t199;
t92 = rSges(6,3) * t166 + (rSges(6,1) * t175 - rSges(6,2) * t172) * t167;
t91 = Icges(6,5) * t166 + (Icges(6,1) * t175 - Icges(6,4) * t172) * t167;
t90 = Icges(6,6) * t166 + (Icges(6,4) * t175 - Icges(6,2) * t172) * t167;
t89 = Icges(6,3) * t166 + (Icges(6,5) * t175 - Icges(6,6) * t172) * t167;
t88 = V_base(5) * rSges(2,3) - t151 * t160 + t208;
t87 = t155 * t160 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t85 = t151 * V_base(4) - t155 * V_base(5) + V_base(3);
t84 = V_base(5) * rSges(3,1) + (-t149 - t150) * t160 + t204;
t83 = t154 * t160 + (-rSges(3,1) - pkin(5)) * V_base(4) + t186;
t82 = rSges(6,1) * t119 + rSges(6,2) * t118 + rSges(6,3) * t215;
t81 = rSges(6,1) * t117 + rSges(6,2) * t116 - rSges(6,3) * t216;
t80 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t215;
t79 = Icges(6,1) * t117 + Icges(6,4) * t116 - Icges(6,5) * t216;
t78 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t215;
t77 = Icges(6,4) * t117 + Icges(6,2) * t116 - Icges(6,6) * t216;
t76 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t215;
t75 = Icges(6,5) * t117 + Icges(6,6) * t116 - Icges(6,3) * t216;
t74 = t150 * V_base(4) + (-t153 - t154) * V_base(5) + t209;
t73 = t152 * t157 + (-t111 + t205) * t160 + t202;
t72 = t110 * t160 - t152 * t158 + t181;
t71 = -t110 * t157 + t111 * t158 + t182;
t70 = t124 * t126 + (-t100 + t203) * t160 + t185;
t69 = -t124 * t127 + t160 * t99 + t179;
t68 = t100 * t127 - t126 * t99 + t180;
t67 = t102 * t92 + t125 * t126 - t133 * t82 + (t113 + t203) * t160 + t185;
t66 = -t103 * t92 + t112 * t160 - t125 * t127 + t133 * t81 + t179;
t65 = -t102 * t81 + t103 * t82 - t112 * t126 - t113 * t127 + t180;
t1 = m(1) * (t130 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(2) * (t85 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(3) * (t74 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(4) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + t158 * (t230 * t174 + t183 * t177) / 0.2e1 + t157 * (t183 * t174 - t230 * t177) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + t127 * (t231 * t174 + t184 * t177) / 0.2e1 + t126 * (t184 * t174 - t231 * t177) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + t103 * ((t116 * t77 + t117 * t79 - t75 * t216) * t103 + (t116 * t78 + t117 * t80 - t216 * t76) * t102 + (t116 * t90 + t117 * t91 - t216 * t89) * t133) / 0.2e1 + t102 * ((t118 * t77 + t119 * t79 + t215 * t75) * t103 + (t118 * t78 + t119 * t80 + t76 * t215) * t102 + (t118 * t90 + t119 * t91 + t215 * t89) * t133) / 0.2e1 + t133 * ((t76 * t102 + t75 * t103 + t89 * t133) * t166 + ((-t172 * t77 + t175 * t79) * t103 + (-t172 * t78 + t175 * t80) * t102 + (-t172 * t90 + t175 * t91) * t133) * t167) / 0.2e1 + ((t174 * t237 + t235 * t177 + Icges(1,4)) * V_base(5) + (t236 * t174 + t234 * t177 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t235 * t174 - t237 * t177 + Icges(1,2)) * V_base(5) + (t174 * t234 - t177 * t236 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t106 * t173 + t108 * t176) * t158 + (-t107 * t173 + t109 * t176) * t157 + (-t166 * t95 + t167 * t97) * t127 + (-t166 * t96 + t167 * t98) * t126 + (-t122 * t166 + t123 * t167 - t141 * t173 + t146 * t176 + Icges(3,1) + Icges(2,3)) * t160) * t160 / 0.2e1 + t160 * V_base(5) * (t242 * t174 - t241 * t177) + t160 * V_base(4) * (t241 * t174 + t242 * t177) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
