% Calculate joint inertia matrix for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:21
% EndTime: 2019-12-31 20:03:26
% DurationCPUTime: 1.87s
% Computational Cost: add. (1762->240), mult. (4140->347), div. (0->0), fcn. (4434->6), ass. (0->124)
t235 = Icges(5,2) + Icges(6,2);
t237 = Icges(3,1) + Icges(4,1);
t232 = Icges(5,1) + Icges(6,1);
t231 = Icges(5,4) + Icges(6,4);
t226 = Icges(3,5) + Icges(4,4);
t236 = Icges(5,5) + Icges(6,5);
t234 = Icges(5,6) + Icges(6,6);
t233 = Icges(5,3) + Icges(6,3);
t120 = sin(qJ(2));
t230 = (Icges(3,4) - Icges(4,5)) * t120;
t121 = sin(qJ(1));
t122 = cos(qJ(4));
t119 = sin(qJ(4));
t123 = cos(qJ(2));
t163 = t119 * t123;
t87 = t120 * t122 - t163;
t79 = t87 * t121;
t229 = t235 * t79;
t124 = cos(qJ(1));
t81 = t87 * t124;
t164 = t119 * t120;
t133 = t122 * t123 + t164;
t82 = t133 * t124;
t228 = -t234 * t121 + t231 * t82 + t235 * t81;
t227 = -t236 * t121 + t231 * t81 + t232 * t82;
t211 = -rSges(6,3) - qJ(5) - pkin(7);
t225 = Icges(4,2) + Icges(3,3);
t80 = t133 * t121;
t180 = Icges(6,4) * t80;
t181 = Icges(5,4) * t80;
t224 = t234 * t124 + t180 + t181 + t229;
t223 = t236 * t124 + t231 * t79 + t232 * t80;
t222 = t226 * t123 + (-Icges(3,6) + Icges(4,6)) * t120;
t221 = t237 * t123 - t230;
t220 = t233 * t124;
t219 = (t233 * t121 - t234 * t81 - t236 * t82) * t121;
t218 = -t121 / 0.2e1;
t217 = t121 / 0.2e1;
t216 = -t124 / 0.2e1;
t215 = t124 / 0.2e1;
t214 = -t234 * t133 + t236 * t87;
t213 = -t235 * t133 + t231 * t87;
t212 = -t231 * t133 + t232 * t87;
t174 = Icges(6,6) * t79;
t175 = Icges(5,6) * t79;
t178 = Icges(6,5) * t80;
t179 = Icges(5,5) * t80;
t192 = (-t227 * t80 - t228 * t79) * t121 + ((0.2e1 * t175 + 0.2e1 * t179 + 0.2e1 * t174 + 0.2e1 * t178 + t220) * t124 + (0.2e1 * t181 + 0.2e1 * t180 + t229) * t79 + t232 * t80 ^ 2 + t219) * t124;
t191 = (-t223 * t82 - t224 * t81) * t124 + (t219 + t227 * t82 + t228 * t81 + (t175 + t179 + t174 + t178 + t220) * t124) * t121;
t199 = -t222 * t121 + t225 * t124;
t198 = t225 * t121 + t222 * t124;
t109 = pkin(4) * t122 + pkin(3);
t157 = pkin(4) * t164;
t161 = t123 * t124;
t197 = t82 * rSges(6,1) + t81 * rSges(6,2) + t109 * t161 + t211 * t121 + t124 * t157;
t147 = -t80 * rSges(6,1) - t79 * rSges(6,2);
t106 = pkin(3) * t161;
t152 = -pkin(7) * t121 + t106;
t188 = -pkin(3) + t109;
t6 = -t121 * (-t147 + (t123 * t188 + t157) * t121 + (-pkin(7) - t211) * t124) - t124 * (-t152 + t197);
t196 = -t191 * t121 - t192 * t124;
t116 = t121 ^ 2;
t117 = t124 ^ 2;
t194 = m(4) / 0.2e1;
t193 = m(6) / 0.2e1;
t98 = rSges(3,1) * t120 + rSges(3,2) * t123;
t190 = m(3) * t98;
t189 = -rSges(5,3) - pkin(7);
t185 = rSges(6,1) * t87 - rSges(6,2) * t133 - pkin(4) * t163 + t120 * t188;
t184 = t82 * rSges(5,1) + t81 * rSges(5,2);
t162 = t120 * t124;
t160 = pkin(2) * t161 + qJ(3) * t162;
t165 = qJ(3) * t120;
t183 = t116 * (pkin(2) * t123 + t165) + t124 * t160;
t96 = pkin(2) * t120 - qJ(3) * t123;
t182 = -rSges(4,1) * t120 + rSges(4,3) * t123 - t96;
t173 = t124 * rSges(4,2);
t172 = t124 * rSges(3,3);
t170 = Icges(3,4) * t123;
t168 = Icges(4,5) * t123;
t159 = t124 * pkin(1) + t121 * pkin(6);
t158 = t116 + t117;
t156 = t226 * t120 / 0.2e1 + (Icges(3,6) / 0.2e1 - Icges(4,6) / 0.2e1) * t123;
t154 = rSges(4,1) * t161 + t121 * rSges(4,2) + rSges(4,3) * t162;
t153 = -pkin(3) * t120 - t96;
t151 = t121 * (t121 * t123 * pkin(3) + pkin(7) * t124) + t124 * t152 + t183;
t150 = t159 + t160;
t55 = rSges(5,1) * t87 - rSges(5,2) * t133;
t149 = t153 - t55;
t148 = -t80 * rSges(5,1) - t79 * rSges(5,2);
t146 = rSges(3,1) * t123 - rSges(3,2) * t120;
t30 = t149 * t121;
t31 = t149 * t124;
t141 = t121 * t30 + t124 * t31;
t20 = t121 * t148 - t124 * t184;
t140 = t153 - t185;
t137 = -Icges(3,2) * t120 + t170;
t134 = Icges(4,3) * t120 + t168;
t132 = rSges(3,1) * t161 - rSges(3,2) * t162 + t121 * rSges(3,3);
t131 = t213 * t79 / 0.2e1 + t212 * t80 / 0.2e1 - t224 * t133 / 0.2e1 + t223 * t87 / 0.2e1 + t214 * t215;
t130 = -t213 * t81 / 0.2e1 - t212 * t82 / 0.2e1 + t228 * t133 / 0.2e1 - t227 * t87 / 0.2e1 + t214 * t217;
t113 = t124 * pkin(6);
t24 = t113 + t189 * t124 + (-t165 - pkin(1) + (-pkin(2) - pkin(3)) * t123) * t121 + t148;
t25 = t121 * t189 + t106 + t150 + t184;
t129 = m(5) * (t121 * t25 + t124 * t24);
t100 = rSges(2,1) * t124 - t121 * rSges(2,2);
t99 = -t121 * rSges(2,1) - rSges(2,2) * t124;
t59 = t182 * t124;
t58 = t182 * t121;
t57 = t132 + t159;
t56 = t172 + t113 + (-pkin(1) - t146) * t121;
t44 = t150 + t154;
t43 = t173 + t113 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t123 + (-rSges(4,3) - qJ(3)) * t120) * t121;
t32 = t124 * t132 + (t121 * t146 - t172) * t121;
t29 = t185 * t124;
t28 = t185 * t121;
t27 = t140 * t124;
t26 = t140 * t121;
t23 = t124 * t154 + (-t173 + (rSges(4,1) * t123 + rSges(4,3) * t120) * t121) * t121 + t183;
t22 = t150 + t197;
t21 = t113 + t211 * t124 + (-pkin(1) + (-pkin(2) - t109) * t123 + (-pkin(4) * t119 - qJ(3)) * t120) * t121 + t147;
t11 = -t20 + t151;
t5 = t151 - t6;
t1 = [Icges(2,3) + t212 * t87 - t213 * t133 + m(6) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2) + m(4) * (t43 ^ 2 + t44 ^ 2) + m(3) * (t56 ^ 2 + t57 ^ 2) + m(2) * (t100 ^ 2 + t99 ^ 2) + ((Icges(3,2) + Icges(4,3)) * t123 + t230) * t123 + (t237 * t120 - t168 + t170) * t120; m(6) * (t21 * t27 + t22 * t26) + m(5) * (t24 * t31 + t25 * t30) + m(4) * (t43 * t59 + t44 * t58) + (-t56 * t190 + t156 * t124 + (Icges(3,6) * t215 + Icges(4,6) * t216 + t134 * t217 + t137 * t218) * t123 - t131) * t124 + (-t57 * t190 + t156 * t121 + (Icges(3,6) * t217 + Icges(4,6) * t218 + t134 * t216 + t137 * t215) * t123 - t130) * t121 + (t221 * t124 * t218 + t226 * t121 * t217 + (t221 * t121 + t226 * t124) * t215) * t120; m(6) * (t26 ^ 2 + t27 ^ 2 + t5 ^ 2) + m(5) * (t11 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(4) * (t23 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(3) * (t158 * t98 ^ 2 + t32 ^ 2) + (t199 * t117 + t192) * t124 + (t198 * t116 + (t199 * t121 + t198 * t124) * t124 + t191) * t121; 0.2e1 * ((t121 * t22 + t124 * t21) * t193 + t129 / 0.2e1 + (t121 * t44 + t124 * t43) * t194) * t120; m(6) * (-t123 * t5 + (t121 * t26 + t124 * t27) * t120) + m(5) * (-t123 * t11 + t120 * t141) + m(4) * (-t123 * t23 + (t121 * t58 + t124 * t59) * t120); 0.2e1 * (t194 + m(5) / 0.2e1 + t193) * (t120 ^ 2 * t158 + t123 ^ 2); t131 * t124 + t130 * t121 + m(6) * (t21 * t29 + t22 * t28) + t55 * t129; m(6) * (t26 * t28 + t27 * t29 + t5 * t6) + m(5) * (t20 * t11 + t141 * t55) + t196; m(5) * (t120 * t158 * t55 - t123 * t20) + m(6) * (-t6 * t123 + (t121 * t28 + t124 * t29) * t120); m(5) * (t158 * t55 ^ 2 + t20 ^ 2) + m(6) * (t28 ^ 2 + t29 ^ 2 + t6 ^ 2) - t196; m(6) * (-t121 * t21 + t124 * t22); m(6) * (-t121 * t27 + t124 * t26); 0; m(6) * (-t121 * t29 + t124 * t28); m(6) * t158;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
