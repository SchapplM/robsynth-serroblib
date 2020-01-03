% Calculate joint inertia matrix for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:08
% EndTime: 2019-12-31 21:19:17
% DurationCPUTime: 3.36s
% Computational Cost: add. (3866->294), mult. (5033->438), div. (0->0), fcn. (5198->8), ass. (0->152)
t247 = Icges(4,4) + Icges(5,6);
t246 = Icges(4,1) + Icges(5,2);
t245 = -Icges(4,2) - Icges(5,3);
t150 = qJ(2) + qJ(3);
t142 = cos(t150);
t244 = t247 * t142;
t141 = sin(t150);
t243 = t247 * t141;
t242 = Icges(5,4) - Icges(4,5);
t241 = Icges(5,5) - Icges(4,6);
t240 = t245 * t141 + t244;
t239 = t246 * t142 - t243;
t153 = sin(qJ(1));
t156 = cos(qJ(1));
t238 = t240 * t153 + t241 * t156;
t237 = -t241 * t153 + t240 * t156;
t236 = t239 * t153 + t242 * t156;
t235 = -t242 * t153 + t239 * t156;
t234 = Icges(5,1) + Icges(4,3);
t233 = t241 * t141 - t242 * t142;
t232 = t245 * t142 - t243;
t231 = t246 * t141 + t244;
t230 = -t233 * t153 + t234 * t156;
t229 = t234 * t153 + t233 * t156;
t228 = t238 * t141 - t236 * t142;
t227 = t237 * t141 - t235 * t142;
t226 = t153 * pkin(6);
t200 = t142 * t156;
t154 = cos(qJ(5));
t196 = t154 * t156;
t151 = sin(qJ(5));
t199 = t151 * t153;
t108 = t141 * t196 - t199;
t197 = t153 * t154;
t198 = t151 * t156;
t109 = t141 * t198 + t197;
t61 = t109 * rSges(6,1) + t108 * rSges(6,2) + rSges(6,3) * t200;
t225 = t153 * pkin(4) + pkin(8) * t200 + t61;
t224 = -t242 * t141 - t241 * t142;
t223 = t232 * t141 + t231 * t142;
t149 = t156 ^ 2;
t110 = t141 * t197 + t198;
t111 = t141 * t199 - t196;
t201 = t142 * t153;
t55 = Icges(6,5) * t109 + Icges(6,6) * t108 + Icges(6,3) * t200;
t57 = Icges(6,4) * t109 + Icges(6,2) * t108 + Icges(6,6) * t200;
t59 = Icges(6,1) * t109 + Icges(6,4) * t108 + Icges(6,5) * t200;
t16 = t110 * t57 + t111 * t59 + t55 * t201;
t56 = Icges(6,5) * t111 + Icges(6,6) * t110 + Icges(6,3) * t201;
t58 = Icges(6,4) * t111 + Icges(6,2) * t110 + Icges(6,6) * t201;
t60 = Icges(6,1) * t111 + Icges(6,4) * t110 + Icges(6,5) * t201;
t17 = t110 * t58 + t111 * t60 + t56 * t201;
t9 = t153 * t16 - t156 * t17;
t222 = -t9 + t230 * t149 + (t227 * t153 + (-t228 + t229) * t156) * t153;
t148 = t153 ^ 2;
t221 = m(5) / 0.2e1;
t220 = m(6) / 0.2e1;
t219 = t153 / 0.2e1;
t218 = -t156 / 0.2e1;
t152 = sin(qJ(2));
t217 = pkin(2) * t152;
t155 = cos(qJ(2));
t139 = pkin(2) * t155 + pkin(1);
t132 = t156 * t139;
t147 = t156 * pkin(6);
t216 = t153 * (t147 + (-pkin(1) + t139) * t153) + t156 * (-t156 * pkin(1) + t132 - t226);
t202 = t141 * t156;
t163 = rSges(4,1) * t200 - rSges(4,2) * t202 + t153 * rSges(4,3);
t183 = rSges(4,1) * t142 - rSges(4,2) * t141;
t52 = t153 * (-t156 * rSges(4,3) + t183 * t153) + t156 * t163;
t194 = pkin(3) * t200 + qJ(4) * t202;
t203 = qJ(4) * t141;
t215 = t148 * (pkin(3) * t142 + t203) + t156 * t194;
t214 = rSges(3,1) * t155;
t213 = rSges(3,2) * t152;
t212 = t156 * rSges(3,3);
t23 = t141 * t55 + (-t151 * t59 - t154 * t57) * t142;
t211 = t23 * t153;
t24 = t141 * t56 + (-t151 * t60 - t154 * t58) * t142;
t210 = t24 * t156;
t209 = Icges(3,4) * t152;
t208 = Icges(3,4) * t155;
t121 = pkin(3) * t141 - qJ(4) * t142;
t195 = rSges(5,2) * t141 + rSges(5,3) * t142 - t121;
t193 = t153 * rSges(3,3) + t156 * t214;
t191 = t148 + t149;
t14 = t108 * t57 + t109 * t59 + t55 * t200;
t15 = t108 * t58 + t109 * t60 + t56 * t200;
t8 = t14 * t153 - t15 * t156;
t190 = (t8 + t229 * t148 + ((-t227 + t230) * t153 + t228 * t156) * t156) * t153;
t123 = rSges(4,1) * t141 + rSges(4,2) * t142;
t189 = -t123 - t217;
t162 = t153 * rSges(5,1) - rSges(5,2) * t200 + rSges(5,3) * t202;
t33 = t153 * (-t156 * rSges(5,1) + (-rSges(5,2) * t142 + rSges(5,3) * t141) * t153) + t156 * t162 + t215;
t157 = -pkin(7) - pkin(6);
t188 = -t153 * t157 + t132;
t79 = Icges(6,3) * t141 + (-Icges(6,5) * t151 - Icges(6,6) * t154) * t142;
t80 = Icges(6,6) * t141 + (-Icges(6,4) * t151 - Icges(6,2) * t154) * t142;
t81 = Icges(6,5) * t141 + (-Icges(6,1) * t151 - Icges(6,4) * t154) * t142;
t28 = t108 * t80 + t109 * t81 + t79 * t200;
t3 = t141 * t28 + (t14 * t156 + t15 * t153) * t142;
t29 = t110 * t80 + t111 * t81 + t79 * t201;
t4 = t141 * t29 + (t153 * t17 + t156 * t16) * t142;
t187 = t9 * t201 / 0.2e1 + t4 * t218 + t3 * t219 + t141 * (-t210 + t211) / 0.2e1 + t8 * t200 / 0.2e1;
t82 = rSges(6,3) * t141 + (-rSges(6,1) * t151 - rSges(6,2) * t154) * t142;
t186 = -pkin(8) * t141 - t121 - t82;
t185 = t195 - t217;
t184 = -t213 + t214;
t182 = -t111 * rSges(6,1) - t110 * rSges(6,2);
t177 = -t151 * t81 - t154 * t80;
t176 = Icges(3,1) * t155 - t209;
t174 = -Icges(3,2) * t152 + t208;
t171 = Icges(3,5) * t155 - Icges(3,6) * t152;
t62 = rSges(6,3) * t201 - t182;
t20 = t215 + t225 * t156 + (-t156 * pkin(4) + pkin(8) * t201 + t62) * t153;
t161 = t188 + t194;
t160 = t186 - t217;
t159 = t222 * t156 + t190;
t158 = -t210 / 0.2e1 + t211 / 0.2e1 + (t235 * t141 + t237 * t142 + t224 * t153 + t223 * t156 + t28) * t219 + (t236 * t141 + t238 * t142 + t223 * t153 - t224 * t156 + t29) * t218;
t130 = rSges(2,1) * t156 - rSges(2,2) * t153;
t129 = -rSges(2,1) * t153 - rSges(2,2) * t156;
t128 = rSges(3,1) * t152 + rSges(3,2) * t155;
t101 = Icges(3,3) * t153 + t171 * t156;
t100 = -Icges(3,3) * t156 + t171 * t153;
t84 = t189 * t156;
t83 = t189 * t153;
t72 = t226 + (pkin(1) - t213) * t156 + t193;
t71 = t212 + t147 + (-pkin(1) - t184) * t153;
t70 = t141 * t79;
t69 = t195 * t156;
t68 = t195 * t153;
t67 = t163 + t188;
t66 = (rSges(4,3) - t157) * t156 + (-t139 - t183) * t153;
t65 = t185 * t156;
t64 = t185 * t153;
t63 = t156 * (-t156 * t213 + t193) + (t184 * t153 - t212) * t153;
t51 = t186 * t156;
t50 = t186 * t153;
t49 = t161 + t162;
t48 = (rSges(5,1) - t157) * t156 + (-t139 + (rSges(5,2) - pkin(3)) * t142 + (-rSges(5,3) - qJ(4)) * t141) * t153;
t47 = t160 * t156;
t46 = t160 * t153;
t37 = t141 * t61 - t82 * t200;
t36 = -t141 * t62 + t82 * t201;
t35 = t161 + t225;
t34 = (pkin(4) - t157) * t156 + (-t203 - t139 + (-rSges(6,3) - pkin(3) - pkin(8)) * t142) * t153 + t182;
t32 = (t177 * t142 + t70) * t141;
t31 = t52 + t216;
t30 = (-t153 * t61 + t156 * t62) * t142;
t25 = t33 + t216;
t11 = t20 + t216;
t1 = [t155 * (Icges(3,2) * t155 + t209) + t152 * (Icges(3,1) * t152 + t208) + Icges(2,3) + t70 + t231 * t141 + (t177 - t232) * t142 + m(6) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t48 ^ 2 + t49 ^ 2) + m(4) * (t66 ^ 2 + t67 ^ 2) + m(3) * (t71 ^ 2 + t72 ^ 2) + m(2) * (t129 ^ 2 + t130 ^ 2); ((-Icges(3,6) * t156 + t174 * t153) * t155 + (-Icges(3,5) * t156 + t176 * t153) * t152) * t218 + ((Icges(3,6) * t153 + t174 * t156) * t155 + (Icges(3,5) * t153 + t176 * t156) * t152) * t219 + t158 + (t148 / 0.2e1 + t149 / 0.2e1) * (Icges(3,5) * t152 + Icges(3,6) * t155) + m(3) * (-t153 * t72 - t156 * t71) * t128 + m(6) * (t34 * t47 + t35 * t46) + m(5) * (t48 * t65 + t49 * t64) + m(4) * (t66 * t84 + t67 * t83); m(6) * (t11 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(5) * (t25 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(4) * (t31 ^ 2 + t83 ^ 2 + t84 ^ 2) + t153 * t148 * t101 + m(3) * (t191 * t128 ^ 2 + t63 ^ 2) + t190 + (-t149 * t100 + (-t153 * t100 + t156 * t101) * t153 + t222) * t156; m(4) * (-t153 * t67 - t156 * t66) * t123 + t158 + m(6) * (t34 * t51 + t35 * t50) + m(5) * (t48 * t69 + t49 * t68); m(6) * (t11 * t20 + t46 * t50 + t47 * t51) + m(5) * (t25 * t33 + t64 * t68 + t65 * t69) + m(4) * (t31 * t52 + (-t153 * t83 - t156 * t84) * t123) + t159; m(6) * (t20 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t33 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t191 * t123 ^ 2 + t52 ^ 2) + t159; 0.2e1 * ((t153 * t35 + t156 * t34) * t220 + (t153 * t49 + t156 * t48) * t221) * t141; m(6) * (-t11 * t142 + (t153 * t46 + t156 * t47) * t141) + m(5) * (-t142 * t25 + (t153 * t64 + t156 * t65) * t141); m(6) * (-t142 * t20 + (t153 * t50 + t156 * t51) * t141) + m(5) * (-t142 * t33 + (t153 * t68 + t156 * t69) * t141); 0.2e1 * (t221 + t220) * (t191 * t141 ^ 2 + t142 ^ 2); t32 + m(6) * (t34 * t36 + t35 * t37) + ((t28 / 0.2e1 + t23 / 0.2e1) * t156 + (t29 / 0.2e1 + t24 / 0.2e1) * t153) * t142; m(6) * (t11 * t30 + t36 * t47 + t37 * t46) + t187; m(6) * (t20 * t30 + t36 * t51 + t37 * t50) + t187; m(6) * (-t142 * t30 + (t153 * t37 + t156 * t36) * t141); m(6) * (t30 ^ 2 + t36 ^ 2 + t37 ^ 2) + t141 * t32 + (t156 * t3 + t153 * t4 + t141 * (t153 * t24 + t156 * t23)) * t142;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
