% Calculate joint inertia matrix for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:14
% EndTime: 2019-03-09 03:11:18
% DurationCPUTime: 2.21s
% Computational Cost: add. (4468->332), mult. (5348->482), div. (0->0), fcn. (5675->8), ass. (0->158)
t229 = Icges(5,1) + Icges(4,3);
t146 = sin(qJ(3));
t149 = cos(qJ(3));
t228 = (-Icges(5,4) + Icges(4,5)) * t149 + (Icges(5,5) - Icges(4,6)) * t146;
t221 = rSges(7,3) + qJ(6);
t222 = rSges(7,1) + pkin(5);
t142 = qJ(1) + pkin(9);
t139 = sin(t142);
t140 = cos(t142);
t145 = sin(qJ(5));
t148 = cos(qJ(5));
t186 = t146 * t148;
t96 = t139 * t186 + t140 * t145;
t187 = t145 * t146;
t98 = t139 * t187 - t140 * t148;
t227 = t221 * t96 - t222 * t98;
t226 = -t139 / 0.2e1;
t225 = t139 / 0.2e1;
t224 = -t140 / 0.2e1;
t223 = t140 / 0.2e1;
t100 = Icges(6,3) * t146 + (-Icges(6,5) * t145 - Icges(6,6) * t148) * t149;
t101 = Icges(7,2) * t146 + (-Icges(7,4) * t145 + Icges(7,6) * t148) * t149;
t184 = t148 * t149;
t99 = Icges(7,6) * t146 + (-Icges(7,5) * t145 + Icges(7,3) * t148) * t149;
t220 = t99 * t184 + (t100 + t101) * t146;
t102 = Icges(6,6) * t146 + (-Icges(6,4) * t145 - Icges(6,2) * t148) * t149;
t103 = Icges(7,4) * t146 + (-Icges(7,1) * t145 + Icges(7,5) * t148) * t149;
t104 = Icges(6,5) * t146 + (-Icges(6,1) * t145 - Icges(6,4) * t148) * t149;
t219 = -t148 * t102 + (-t103 - t104) * t145;
t188 = t140 * t149;
t94 = t139 * t145 - t140 * t186;
t95 = t139 * t148 + t140 * t187;
t45 = Icges(7,5) * t95 + Icges(7,6) * t188 + Icges(7,3) * t94;
t49 = Icges(7,4) * t95 + Icges(7,2) * t188 + Icges(7,6) * t94;
t53 = Icges(7,1) * t95 + Icges(7,4) * t188 + Icges(7,5) * t94;
t12 = t49 * t188 + t45 * t94 + t53 * t95;
t190 = t139 * t149;
t46 = Icges(7,5) * t98 + Icges(7,6) * t190 - Icges(7,3) * t96;
t50 = Icges(7,4) * t98 + Icges(7,2) * t190 - Icges(7,6) * t96;
t54 = Icges(7,1) * t98 + Icges(7,4) * t190 - Icges(7,5) * t96;
t13 = t50 * t188 + t46 * t94 + t54 * t95;
t47 = Icges(6,5) * t95 - Icges(6,6) * t94 + Icges(6,3) * t188;
t51 = Icges(6,4) * t95 - Icges(6,2) * t94 + Icges(6,6) * t188;
t55 = Icges(6,1) * t95 - Icges(6,4) * t94 + Icges(6,5) * t188;
t14 = t47 * t188 - t51 * t94 + t55 * t95;
t48 = Icges(6,5) * t98 + Icges(6,6) * t96 + Icges(6,3) * t190;
t52 = Icges(6,4) * t98 + Icges(6,2) * t96 + Icges(6,6) * t190;
t56 = Icges(6,1) * t98 + Icges(6,4) * t96 + Icges(6,5) * t190;
t15 = t48 * t188 - t52 * t94 + t56 * t95;
t30 = t101 * t188 + t103 * t95 + t94 * t99;
t31 = t100 * t188 - t102 * t94 + t104 * t95;
t218 = ((t12 + t14) * t140 + (t13 + t15) * t139) * t149 + (t30 + t31) * t146;
t16 = t49 * t190 - t45 * t96 + t53 * t98;
t17 = t50 * t190 - t46 * t96 + t54 * t98;
t18 = t47 * t190 + t51 * t96 + t55 * t98;
t19 = t48 * t190 + t52 * t96 + t56 * t98;
t32 = t101 * t190 + t103 * t98 - t96 * t99;
t33 = t100 * t190 + t102 * t96 + t104 * t98;
t217 = ((t16 + t18) * t140 + (t17 + t19) * t139) * t149 + (t32 + t33) * t146;
t216 = t146 / 0.2e1;
t20 = t146 * t49 + (-t145 * t53 + t148 * t45) * t149;
t22 = t146 * t47 + (-t145 * t55 - t148 * t51) * t149;
t215 = t20 + t22;
t21 = t146 * t50 + (-t145 * t54 + t148 * t46) * t149;
t23 = t146 * t48 + (-t145 * t56 - t148 * t52) * t149;
t214 = t21 + t23;
t213 = -t228 * t139 + t140 * t229;
t212 = t229 * t139 + t228 * t140;
t137 = t139 ^ 2;
t138 = t140 ^ 2;
t211 = m(5) / 0.2e1;
t210 = m(6) / 0.2e1;
t209 = m(7) / 0.2e1;
t208 = -pkin(3) - pkin(8);
t120 = rSges(4,1) * t146 + rSges(4,2) * t149;
t204 = m(4) * t120;
t147 = sin(qJ(1));
t203 = pkin(1) * t147;
t202 = (t149 * t219 + t220) * t146;
t201 = rSges(7,2) * t188 + t221 * t94 + t222 * t95;
t200 = rSges(7,2) * t190 - t227;
t189 = t140 * t146;
t181 = pkin(3) * t188 + qJ(4) * t189;
t191 = qJ(4) * t146;
t198 = t137 * (pkin(3) * t149 + t191) + t140 * t181;
t197 = t140 * rSges(5,1);
t196 = t140 * rSges(4,3);
t195 = Icges(4,4) * t146;
t194 = Icges(4,4) * t149;
t193 = Icges(5,6) * t146;
t192 = Icges(5,6) * t149;
t183 = rSges(7,2) * t146 + (-t222 * t145 + t221 * t148) * t149;
t118 = pkin(3) * t146 - qJ(4) * t149;
t182 = rSges(5,2) * t146 + rSges(5,3) * t149 - t118;
t180 = t139 * pkin(4) + pkin(8) * t188;
t179 = t137 + t138;
t178 = -m(5) - m(6) - m(7);
t58 = t95 * rSges(6,1) - t94 * rSges(6,2) + rSges(6,3) * t188;
t150 = cos(qJ(1));
t141 = t150 * pkin(1);
t177 = t140 * pkin(2) + t139 * pkin(7) + t141;
t176 = -Icges(5,4) * t146 / 0.2e1 + Icges(4,5) * t216 + (-Icges(5,5) / 0.2e1 + Icges(4,6) / 0.2e1) * t149;
t175 = t140 * pkin(7) - t203;
t174 = -pkin(2) - t191;
t173 = -pkin(8) * t146 - t118;
t135 = t140 * pkin(4);
t172 = t139 * (pkin(8) * t190 - t135) + t140 * t180 + t198;
t171 = t135 + t175;
t106 = rSges(6,3) * t146 + (-rSges(6,1) * t145 - rSges(6,2) * t148) * t149;
t170 = -t106 + t173;
t169 = -rSges(6,1) * t98 - rSges(6,2) * t96;
t168 = rSges(4,1) * t149 - rSges(4,2) * t146;
t163 = t177 + t181;
t162 = Icges(4,1) * t149 - t195;
t161 = -Icges(4,2) * t146 + t194;
t158 = -Icges(5,2) * t149 + t193;
t157 = Icges(5,3) * t146 - t192;
t156 = t173 - t183;
t155 = rSges(4,1) * t188 - rSges(4,2) * t189 + t139 * rSges(4,3);
t154 = t139 * rSges(5,1) - rSges(5,2) * t188 + rSges(5,3) * t189;
t153 = t23 / 0.2e1 + t21 / 0.2e1 + t33 / 0.2e1 + t32 / 0.2e1;
t152 = t31 / 0.2e1 + t30 / 0.2e1 + t22 / 0.2e1 + t20 / 0.2e1;
t151 = t163 + t180;
t144 = t149 ^ 2;
t122 = rSges(2,1) * t150 - t147 * rSges(2,2);
t121 = -t147 * rSges(2,1) - rSges(2,2) * t150;
t108 = rSges(3,1) * t140 - rSges(3,2) * t139 + t141;
t107 = -rSges(3,1) * t139 - rSges(3,2) * t140 - t203;
t68 = t182 * t140;
t67 = t182 * t139;
t66 = t155 + t177;
t65 = t196 + (-pkin(2) - t168) * t139 + t175;
t62 = t170 * t140;
t61 = t170 * t139;
t60 = rSges(6,3) * t190 - t169;
t44 = t154 + t163;
t43 = t197 + (-pkin(2) + (rSges(5,2) - pkin(3)) * t149 + (-rSges(5,3) - qJ(4)) * t146) * t139 + t175;
t42 = t140 * t155 + (t168 * t139 - t196) * t139;
t41 = t156 * t140;
t40 = t156 * t139;
t37 = -t106 * t188 + t146 * t58;
t36 = t106 * t190 - t146 * t60;
t35 = t151 + t58;
t34 = ((-rSges(6,3) + t208) * t149 + t174) * t139 + t169 + t171;
t29 = t140 * t154 + (-t197 + (-rSges(5,2) * t149 + rSges(5,3) * t146) * t139) * t139 + t198;
t28 = (-t139 * t58 + t140 * t60) * t149;
t27 = t151 + t201;
t26 = ((-rSges(7,2) + t208) * t149 + t174) * t139 + t171 + t227;
t25 = t201 * t146 - t183 * t188;
t24 = -t200 * t146 + t183 * t190;
t11 = t139 * t60 + t140 * t58 + t172;
t10 = (-t201 * t139 + t200 * t140) * t149;
t9 = t200 * t139 + t201 * t140 + t172;
t8 = t139 * t18 - t140 * t19;
t7 = t139 * t16 - t140 * t17;
t6 = t139 * t14 - t140 * t15;
t5 = t12 * t139 - t13 * t140;
t1 = [Icges(2,3) + Icges(3,3) + m(7) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t43 ^ 2 + t44 ^ 2) + m(4) * (t65 ^ 2 + t66 ^ 2) + m(3) * (t107 ^ 2 + t108 ^ 2) + m(2) * (t121 ^ 2 + t122 ^ 2) + (t193 + t195 + (Icges(5,3) + Icges(4,2)) * t149 + t219) * t149 + (t192 + t194 + (Icges(4,1) + Icges(5,2)) * t146) * t146 + t220; 0; m(3) + m(4) - t178; m(7) * (t26 * t41 + t27 * t40) + m(6) * (t34 * t62 + t35 * t61) + m(5) * (t43 * t68 + t44 * t67) + (-t65 * t204 + t176 * t140 + (Icges(5,5) * t224 + Icges(4,6) * t223 + t157 * t225 + t161 * t226) * t149 + (Icges(5,4) * t224 + Icges(4,5) * t223 + t158 * t225 + t162 * t226) * t146 - t153) * t140 + (-t66 * t204 + t176 * t139 + (Icges(5,5) * t226 + Icges(4,6) * t225 + t157 * t224 + t161 * t223) * t149 + (Icges(5,4) * t226 + Icges(4,5) * t225 + t158 * t224 + t162 * t223) * t146 + t152) * t139; m(4) * t42 + m(5) * t29 + m(6) * t11 + m(7) * t9; m(7) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(6) * (t11 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t29 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(4) * (t179 * t120 ^ 2 + t42 ^ 2) + (t213 * t138 - t7 - t8) * t140 + (t5 + t6 + t212 * t137 + (t213 * t139 + t212 * t140) * t140) * t139; 0.2e1 * ((t139 * t27 + t140 * t26) * t209 + (t139 * t35 + t140 * t34) * t210 + (t139 * t44 + t140 * t43) * t211) * t146; t178 * t149; m(7) * (-t149 * t9 + (t139 * t40 + t140 * t41) * t146) + m(6) * (-t11 * t149 + (t139 * t61 + t140 * t62) * t146) + m(5) * (-t149 * t29 + (t139 * t67 + t140 * t68) * t146); 0.2e1 * (t211 + t210 + t209) * (t179 * t146 ^ 2 + t144); m(7) * (t24 * t26 + t25 * t27) + m(6) * (t34 * t36 + t35 * t37) + (t153 * t139 + t152 * t140) * t149 + t202; m(6) * t28 + m(7) * t10; m(7) * (t10 * t9 + t24 * t41 + t25 * t40) + m(6) * (t11 * t28 + t36 * t62 + t37 * t61) + ((t6 / 0.2e1 + t5 / 0.2e1) * t140 + (t8 / 0.2e1 + t7 / 0.2e1) * t139) * t149 + t218 * t225 + t217 * t224 + (t215 * t139 - t214 * t140) * t216; m(6) * (-t149 * t28 + (t139 * t37 + t140 * t36) * t146) + m(7) * (-t10 * t149 + (t139 * t25 + t140 * t24) * t146); t202 * t146 + m(7) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t28 ^ 2 + t36 ^ 2 + t37 ^ 2) + ((t215 * t146 + t218) * t140 + (t214 * t146 + t217) * t139) * t149; m(7) * (t26 * t94 - t27 * t96); m(7) * t184; m(7) * (t9 * t184 - t40 * t96 + t41 * t94); m(7) * (-t144 * t148 + (-t139 * t96 + t140 * t94) * t146); m(7) * (t10 * t184 + t24 * t94 - t25 * t96); m(7) * (t144 * t148 ^ 2 + t94 ^ 2 + t96 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
