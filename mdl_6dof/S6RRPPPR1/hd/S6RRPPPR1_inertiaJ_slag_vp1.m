% Calculate joint inertia matrix for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:06:05
% EndTime: 2019-03-09 08:06:12
% DurationCPUTime: 3.15s
% Computational Cost: add. (4786->381), mult. (6704->558), div. (0->0), fcn. (7615->10), ass. (0->178)
t248 = Icges(5,1) + Icges(6,1);
t247 = Icges(5,4) - Icges(6,5);
t246 = Icges(6,4) + Icges(5,5);
t245 = -Icges(5,2) - Icges(6,3);
t244 = Icges(5,6) - Icges(6,6);
t243 = Icges(3,3) + Icges(4,3);
t146 = qJ(2) + pkin(9);
t141 = sin(t146);
t142 = cos(t146);
t153 = sin(qJ(2));
t156 = cos(qJ(2));
t242 = Icges(3,5) * t156 + Icges(4,5) * t142 - Icges(3,6) * t153 - Icges(4,6) * t141;
t154 = sin(qJ(1));
t241 = -t154 / 0.2e1;
t223 = t154 / 0.2e1;
t157 = cos(qJ(1));
t240 = t157 / 0.2e1;
t149 = sin(pkin(10));
t201 = t154 * t149;
t150 = cos(pkin(10));
t203 = t150 * t157;
t113 = t142 * t201 + t203;
t200 = t154 * t150;
t204 = t149 * t157;
t114 = t142 * t200 - t204;
t207 = t141 * t154;
t51 = Icges(5,5) * t114 - Icges(5,6) * t113 + Icges(5,3) * t207;
t53 = Icges(6,4) * t114 + Icges(6,2) * t207 + Icges(6,6) * t113;
t239 = t51 + t53;
t115 = t142 * t204 - t200;
t116 = t142 * t203 + t201;
t206 = t141 * t157;
t52 = Icges(5,5) * t116 - Icges(5,6) * t115 + Icges(5,3) * t206;
t54 = Icges(6,4) * t116 + Icges(6,2) * t206 + Icges(6,6) * t115;
t238 = t52 + t54;
t233 = t113 * t245 + t114 * t247 + t207 * t244;
t232 = t115 * t245 + t116 * t247 + t206 * t244;
t231 = t247 * t113 - t114 * t248 - t246 * t207;
t230 = -t247 * t115 + t116 * t248 + t246 * t206;
t224 = m(7) / 0.2e1;
t225 = m(6) / 0.2e1;
t196 = t225 + t224;
t237 = 0.2e1 * t196;
t236 = t142 / 0.2e1;
t235 = t153 / 0.2e1;
t234 = t156 / 0.2e1;
t229 = -t242 * t154 + t243 * t157;
t228 = t243 * t154 + t242 * t157;
t147 = t154 ^ 2;
t148 = t157 ^ 2;
t227 = 0.2e1 * t141;
t226 = m(5) / 0.2e1;
t222 = -t157 / 0.2e1;
t221 = -rSges(7,3) - pkin(8);
t220 = pkin(2) * t153;
t219 = pkin(3) * t142;
t152 = sin(qJ(6));
t155 = cos(qJ(6));
t70 = t115 * t155 - t116 * t152;
t71 = t115 * t152 + t116 * t155;
t218 = t71 * rSges(7,1) + t70 * rSges(7,2);
t139 = pkin(2) * t156 + pkin(1);
t135 = t157 * t139;
t145 = t157 * pkin(7);
t151 = -qJ(3) - pkin(7);
t202 = t151 * t157;
t217 = t154 * (t202 + t145 + (-pkin(1) + t139) * t154) + t157 * (-pkin(1) * t157 + t135 + (-pkin(7) - t151) * t154);
t216 = rSges(3,1) * t156;
t215 = rSges(3,2) * t153;
t214 = t113 * rSges(6,3);
t213 = t157 * rSges(3,3);
t212 = Icges(3,4) * t153;
t211 = Icges(3,4) * t156;
t210 = Icges(4,4) * t141;
t209 = Icges(4,4) * t142;
t208 = t141 * t149;
t205 = t142 * t157;
t199 = pkin(3) * t205 + qJ(4) * t206;
t198 = t154 * rSges(3,3) + t157 * t216;
t197 = t147 + t148;
t94 = (t149 * t155 - t150 * t152) * t141;
t95 = (t149 * t152 + t150 * t155) * t141;
t43 = Icges(7,5) * t95 + Icges(7,6) * t94 + Icges(7,3) * t142;
t44 = Icges(7,4) * t95 + Icges(7,2) * t94 + Icges(7,6) * t142;
t45 = Icges(7,1) * t95 + Icges(7,4) * t94 + Icges(7,5) * t142;
t195 = t142 * t43 + t94 * t44 + t95 * t45;
t67 = t113 * t155 - t114 * t152;
t68 = t113 * t152 + t114 * t155;
t27 = Icges(7,5) * t68 + Icges(7,6) * t67 - Icges(7,3) * t207;
t29 = Icges(7,4) * t68 + Icges(7,2) * t67 - Icges(7,6) * t207;
t31 = Icges(7,1) * t68 + Icges(7,4) * t67 - Icges(7,5) * t207;
t10 = t142 * t27 + t29 * t94 + t31 * t95;
t13 = -t43 * t207 + t44 * t67 + t45 * t68;
t194 = -t10 / 0.2e1 - t13 / 0.2e1;
t28 = Icges(7,5) * t71 + Icges(7,6) * t70 - Icges(7,3) * t206;
t30 = Icges(7,4) * t71 + Icges(7,2) * t70 - Icges(7,6) * t206;
t32 = Icges(7,1) * t71 + Icges(7,4) * t70 - Icges(7,5) * t206;
t11 = t142 * t28 + t30 * t94 + t32 * t95;
t14 = -t43 * t206 + t70 * t44 + t71 * t45;
t193 = -t11 / 0.2e1 - t14 / 0.2e1;
t78 = -Icges(6,6) * t142 + (Icges(6,5) * t150 + Icges(6,3) * t149) * t141;
t81 = -Icges(5,6) * t142 + (Icges(5,4) * t150 - Icges(5,2) * t149) * t141;
t192 = t78 / 0.2e1 - t81 / 0.2e1;
t82 = -Icges(6,4) * t142 + (Icges(6,1) * t150 + Icges(6,5) * t149) * t141;
t83 = -Icges(5,5) * t142 + (Icges(5,1) * t150 - Icges(5,4) * t149) * t141;
t191 = t82 / 0.2e1 + t83 / 0.2e1;
t190 = t116 * rSges(6,1) + rSges(6,2) * t206 + t115 * rSges(6,3);
t189 = t116 * rSges(5,1) - t115 * rSges(5,2) + rSges(5,3) * t206;
t188 = Icges(4,5) * t141 / 0.2e1 + Icges(4,6) * t236 + Icges(3,5) * t235 + Icges(3,6) * t234;
t187 = -pkin(3) * t141 + qJ(4) * t142 - t220;
t186 = -rSges(4,1) * t141 - rSges(4,2) * t142 - t220;
t185 = -t139 - t219;
t184 = t147 * (qJ(4) * t141 + t219) + t157 * t199 + t217;
t99 = t113 * qJ(5);
t183 = -t99 - t202;
t182 = t116 * pkin(4) + t115 * qJ(5);
t181 = -t154 * t151 + t135;
t180 = t226 + t196;
t179 = t187 + rSges(5,3) * t142 - (rSges(5,1) * t150 - rSges(5,2) * t149) * t141;
t178 = -(pkin(4) * t150 + qJ(5) * t149) * t141 + t187;
t177 = -t68 * rSges(7,1) - t67 * rSges(7,2);
t176 = -t215 + t216;
t175 = rSges(4,1) * t142 - rSges(4,2) * t141;
t174 = -t114 * rSges(5,1) + t113 * rSges(5,2);
t171 = t154 * (t114 * pkin(4) + t99) + t157 * t182 + t184;
t170 = Icges(3,1) * t156 - t212;
t169 = Icges(4,1) * t142 - t210;
t168 = -Icges(3,2) * t153 + t211;
t167 = -Icges(4,2) * t141 + t209;
t162 = t178 + rSges(6,2) * t142 - (rSges(6,1) * t150 + rSges(6,3) * t149) * t141;
t161 = rSges(4,1) * t205 - rSges(4,2) * t206 + t154 * rSges(4,3);
t160 = t181 + t199;
t46 = rSges(7,1) * t95 + rSges(7,2) * t94 + rSges(7,3) * t142;
t159 = -pkin(5) * t141 * t150 - pkin(8) * t142 + t178 - t46;
t158 = t160 + t182;
t140 = t141 ^ 2;
t129 = rSges(2,1) * t157 - t154 * rSges(2,2);
t128 = -t154 * rSges(2,1) - rSges(2,2) * t157;
t127 = rSges(3,1) * t153 + rSges(3,2) * t156;
t105 = t116 * pkin(5);
t87 = t186 * t157;
t86 = t186 * t154;
t77 = t154 * pkin(7) + (pkin(1) - t215) * t157 + t198;
t76 = t213 + t145 + (-pkin(1) - t176) * t154;
t73 = t161 + t181;
t72 = (rSges(4,3) - t151) * t157 + (-t139 - t175) * t154;
t62 = t157 * (-t157 * t215 + t198) + (t176 * t154 - t213) * t154;
t48 = t179 * t157;
t47 = t179 * t154;
t40 = t162 * t157;
t39 = t162 * t154;
t36 = t160 + t189;
t35 = -t202 + ((-rSges(5,3) - qJ(4)) * t141 + t185) * t154 + t174;
t34 = -rSges(7,3) * t206 + t218;
t33 = -rSges(7,3) * t207 - t177;
t26 = t157 * t161 + (-t157 * rSges(4,3) + t175 * t154) * t154 + t217;
t25 = t158 + t190;
t24 = -t214 + (-rSges(6,1) - pkin(4)) * t114 + ((-rSges(6,2) - qJ(4)) * t141 + t185) * t154 + t183;
t23 = t159 * t157;
t22 = t159 * t154;
t21 = t142 * t34 + t46 * t206;
t20 = -t142 * t33 - t46 * t207;
t19 = t221 * t206 + t105 + t158 + t218;
t18 = (-pkin(4) - pkin(5)) * t114 + ((-qJ(4) - t221) * t141 + t185) * t154 + t177 + t183;
t17 = t154 * (rSges(5,3) * t207 - t174) + t157 * t189 + t184;
t16 = (t154 * t34 - t157 * t33) * t141;
t15 = t195 * t142;
t12 = t154 * (t114 * rSges(6,1) + rSges(6,2) * t207 + t214) + t157 * t190 + t171;
t9 = -t28 * t206 + t70 * t30 + t71 * t32;
t8 = -t27 * t206 + t70 * t29 + t71 * t31;
t7 = -t28 * t207 + t30 * t67 + t32 * t68;
t6 = -t27 * t207 + t29 * t67 + t31 * t68;
t5 = (-pkin(8) * t206 + t105 + t34) * t157 + (t114 * pkin(5) - pkin(8) * t207 + t33) * t154 + t171;
t4 = t9 * t154 - t157 * t8;
t3 = t7 * t154 - t157 * t6;
t2 = t14 * t142 + (-t154 * t8 - t157 * t9) * t141;
t1 = t13 * t142 + (-t154 * t6 - t157 * t7) * t141;
t37 = [t156 * (Icges(3,2) * t156 + t212) + t153 * (Icges(3,1) * t153 + t211) + Icges(2,3) + (t210 + (t149 * t244 - t150 * t246) * t141 + (Icges(6,2) + Icges(4,2) + Icges(5,3)) * t142) * t142 + (Icges(4,1) * t141 + t209 + (t82 + t83) * t150 + (t78 - t81) * t149) * t141 + m(7) * (t18 ^ 2 + t19 ^ 2) + m(6) * (t24 ^ 2 + t25 ^ 2) + m(4) * (t72 ^ 2 + t73 ^ 2) + m(5) * (t35 ^ 2 + t36 ^ 2) + m(3) * (t76 ^ 2 + t77 ^ 2) + m(2) * (t128 ^ 2 + t129 ^ 2) + t195; (-t156 * (-Icges(3,6) * t157 + t168 * t154) / 0.2e1 - t153 * (-Icges(3,5) * t157 + t170 * t154) / 0.2e1 + t188 * t157 - t191 * t114 - t192 * t113 + t194) * t157 + ((Icges(3,6) * t154 + t168 * t157) * t234 + (Icges(3,5) * t154 + t170 * t157) * t235 + t192 * t115 + t191 * t116 + t188 * t154 - t193) * t154 + m(5) * (t35 * t48 + t36 * t47) + m(6) * (t24 * t40 + t25 * t39) + m(7) * (t18 * t23 + t19 * t22) + m(4) * (t72 * t87 + t73 * t86) + m(3) * (-t154 * t77 - t157 * t76) * t127 + ((Icges(4,6) * t240 + t167 * t241 + t51 / 0.2e1 + t53 / 0.2e1) * t157 + (Icges(4,6) * t223 + t167 * t240 - t54 / 0.2e1 - t52 / 0.2e1) * t154) * t142 + ((Icges(4,5) * t154 - t232 * t149 + t230 * t150 + t169 * t157) * t223 + (-Icges(4,5) * t157 - t233 * t149 - t231 * t150 + t169 * t154) * t222) * t141; m(7) * (t22 ^ 2 + t23 ^ 2 + t5 ^ 2) + m(6) * (t12 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t17 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(4) * (t26 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(3) * (t197 * t127 ^ 2 + t62 ^ 2) + (-t3 + (-t233 * t113 - t231 * t114 + t239 * t207) * t157 + t229 * t148) * t157 + (t4 + (-t232 * t115 + t230 * t116 + t238 * t206) * t154 + t228 * t147 + (t232 * t113 - t230 * t114 + t233 * t115 + t231 * t116 + t229 * t154 + t228 * t157 - t239 * t206 - t238 * t207) * t157) * t154; m(7) * (t154 * t18 - t157 * t19) + m(6) * (t154 * t24 - t157 * t25) + m(4) * (t154 * t72 - t157 * t73) + m(5) * (t154 * t35 - t157 * t36); m(7) * (t154 * t23 - t157 * t22) + m(6) * (t154 * t40 - t157 * t39) + m(5) * (t154 * t48 - t157 * t47) + m(4) * (t154 * t87 - t157 * t86); 0.2e1 * (m(4) / 0.2e1 + t180) * t197; ((t154 * t19 + t157 * t18) * t224 + (t154 * t25 + t157 * t24) * t225 + (t154 * t36 + t157 * t35) * t226) * t227; m(7) * (-t142 * t5 + (t154 * t22 + t157 * t23) * t141) + m(6) * (-t142 * t12 + (t154 * t39 + t157 * t40) * t141) + m(5) * (-t142 * t17 + (t154 * t47 + t157 * t48) * t141); 0; 0.2e1 * t180 * (t197 * t140 + t142 ^ 2); m(7) * (t113 * t19 + t115 * t18) + m(6) * (t113 * t25 + t115 * t24); m(7) * (t113 * t22 + t115 * t23 + t5 * t208) + m(6) * (t113 * t39 + t115 * t40 + t12 * t208); (-t113 * t157 + t115 * t154) * t237; t196 * (t113 * t154 + t115 * t157 - t142 * t149) * t227; (t140 * t149 ^ 2 + t113 ^ 2 + t115 ^ 2) * t237; m(7) * (t18 * t20 + t19 * t21) + t15 + (t194 * t154 + t193 * t157) * t141; m(7) * (t16 * t5 + t20 * t23 + t21 * t22) + (-t10 * t157 + t11 * t154) * t236 + t2 * t223 + t1 * t222 + (t4 * t222 + t3 * t241) * t141; m(7) * (t20 * t154 - t157 * t21); m(7) * (-t16 * t142 + (t154 * t21 + t157 * t20) * t141); m(7) * (t113 * t21 + t115 * t20 + t16 * t208); t142 * t15 + m(7) * (t16 ^ 2 + t20 ^ 2 + t21 ^ 2) + (-t157 * t2 - t154 * t1 + t142 * (-t10 * t154 - t11 * t157)) * t141;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t37(1) t37(2) t37(4) t37(7) t37(11) t37(16); t37(2) t37(3) t37(5) t37(8) t37(12) t37(17); t37(4) t37(5) t37(6) t37(9) t37(13) t37(18); t37(7) t37(8) t37(9) t37(10) t37(14) t37(19); t37(11) t37(12) t37(13) t37(14) t37(15) t37(20); t37(16) t37(17) t37(18) t37(19) t37(20) t37(21);];
Mq  = res;
