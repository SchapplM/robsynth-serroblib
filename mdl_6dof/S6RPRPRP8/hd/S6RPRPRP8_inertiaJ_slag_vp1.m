% Calculate joint inertia matrix for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:27
% EndTime: 2019-03-09 03:24:32
% DurationCPUTime: 2.12s
% Computational Cost: add. (3886->325), mult. (5260->474), div. (0->0), fcn. (5550->8), ass. (0->150)
t230 = rSges(7,1) + pkin(5);
t229 = rSges(7,3) + qJ(6);
t234 = Icges(4,3) + Icges(5,3);
t146 = qJ(3) + pkin(9);
t139 = sin(t146);
t140 = cos(t146);
t151 = sin(qJ(3));
t154 = cos(qJ(3));
t233 = Icges(4,5) * t151 + Icges(5,5) * t139 + Icges(4,6) * t154 + Icges(5,6) * t140;
t152 = sin(qJ(1));
t153 = cos(qJ(5));
t190 = t152 * t153;
t150 = sin(qJ(5));
t155 = cos(qJ(1));
t193 = t150 * t155;
t108 = t139 * t193 + t190;
t188 = t155 * t153;
t191 = t152 * t150;
t110 = -t139 * t188 + t191;
t232 = -t108 * t229 + t110 * t230;
t231 = Icges(4,5) * t154 + Icges(5,5) * t140 - Icges(4,6) * t151 - Icges(5,6) * t139;
t197 = t140 * t150;
t72 = Icges(7,6) * t139 + (Icges(7,5) * t153 + Icges(7,3) * t150) * t140;
t73 = Icges(6,3) * t139 + (Icges(6,5) * t153 - Icges(6,6) * t150) * t140;
t74 = Icges(7,2) * t139 + (Icges(7,4) * t153 + Icges(7,6) * t150) * t140;
t76 = Icges(7,4) * t139 + (Icges(7,1) * t153 + Icges(7,5) * t150) * t140;
t77 = Icges(6,5) * t139 + (Icges(6,1) * t153 - Icges(6,4) * t150) * t140;
t228 = t72 * t197 + (t76 + t77) * t140 * t153 + (t73 + t74) * t139;
t106 = t139 * t191 - t188;
t107 = t139 * t190 + t193;
t196 = t140 * t152;
t44 = Icges(7,5) * t107 - Icges(7,6) * t196 + Icges(7,3) * t106;
t48 = Icges(7,4) * t107 - Icges(7,2) * t196 + Icges(7,6) * t106;
t52 = Icges(7,1) * t107 - Icges(7,4) * t196 + Icges(7,5) * t106;
t11 = t106 * t44 + t107 * t52 - t196 * t48;
t194 = t140 * t155;
t45 = Icges(7,5) * t110 + Icges(7,6) * t194 - Icges(7,3) * t108;
t49 = Icges(7,4) * t110 + Icges(7,2) * t194 - Icges(7,6) * t108;
t53 = Icges(7,1) * t110 + Icges(7,4) * t194 - Icges(7,5) * t108;
t12 = t106 * t45 + t107 * t53 - t196 * t49;
t46 = Icges(6,5) * t107 - Icges(6,6) * t106 - Icges(6,3) * t196;
t50 = Icges(6,4) * t107 - Icges(6,2) * t106 - Icges(6,6) * t196;
t54 = Icges(6,1) * t107 - Icges(6,4) * t106 - Icges(6,5) * t196;
t13 = -t106 * t50 + t107 * t54 - t196 * t46;
t47 = Icges(6,5) * t110 + Icges(6,6) * t108 + Icges(6,3) * t194;
t51 = Icges(6,4) * t110 + Icges(6,2) * t108 + Icges(6,6) * t194;
t55 = Icges(6,1) * t110 + Icges(6,4) * t108 + Icges(6,5) * t194;
t14 = -t106 * t51 + t107 * t55 - t196 * t47;
t28 = t106 * t72 + t107 * t76 - t196 * t74;
t75 = Icges(6,6) * t139 + (Icges(6,4) * t153 - Icges(6,2) * t150) * t140;
t29 = -t106 * t75 + t107 * t77 - t196 * t73;
t227 = ((t12 + t14) * t155 + (-t11 - t13) * t152) * t140 + (t28 + t29) * t139;
t15 = -t108 * t44 + t110 * t52 + t194 * t48;
t16 = -t108 * t45 + t110 * t53 + t194 * t49;
t17 = t108 * t50 + t110 * t54 + t194 * t46;
t18 = t108 * t51 + t110 * t55 + t194 * t47;
t30 = -t108 * t72 + t110 * t76 + t194 * t74;
t31 = t108 * t75 + t110 * t77 + t194 * t73;
t226 = ((t16 + t18) * t155 + (-t15 - t17) * t152) * t140 + (t30 + t31) * t139;
t20 = t139 * t48 + (t150 * t44 + t153 * t52) * t140;
t22 = t139 * t46 + (-t150 * t50 + t153 * t54) * t140;
t225 = t20 + t22;
t21 = t139 * t49 + (t150 * t45 + t153 * t53) * t140;
t23 = t139 * t47 + (-t150 * t51 + t153 * t55) * t140;
t224 = t21 + t23;
t223 = t152 * t233 + t155 * t234;
t222 = t152 * t234 - t155 * t233;
t221 = t106 * t229 + t107 * t230;
t220 = (rSges(4,1) * t151 + rSges(4,2) * t154) * t155;
t147 = t152 ^ 2;
t148 = t155 ^ 2;
t125 = rSges(4,1) * t154 - rSges(4,2) * t151;
t212 = m(4) * t125;
t211 = pkin(3) * t154;
t210 = (-t197 * t75 + t228) * t139;
t209 = -rSges(7,2) * t196 + t221;
t208 = rSges(7,2) * t194 + t232;
t206 = rSges(7,2) * t139 + (t150 * t229 + t153 * t230) * t140;
t130 = t155 * t139 * pkin(4);
t149 = -qJ(4) - pkin(7);
t185 = t155 * t151 * pkin(3) + t152 * t149;
t89 = t155 * (-t152 * pkin(7) - t185);
t205 = t155 * (pkin(8) * t194 - t130) + t89;
t143 = t155 * rSges(5,3);
t203 = t107 * rSges(6,1) - t106 * rSges(6,2);
t198 = t139 * t152;
t192 = t151 * t152;
t189 = t152 * t154;
t134 = pkin(3) * t192;
t105 = t134 + (-pkin(7) - t149) * t155;
t129 = pkin(4) * t198;
t187 = pkin(8) * t196 - t105 - t129;
t116 = pkin(4) * t140 + pkin(8) * t139;
t135 = pkin(3) * t189;
t186 = t152 * t116 + t135;
t184 = t155 * pkin(1) + t152 * qJ(2);
t131 = t147 + t148;
t182 = (m(5) + m(6) + m(7)) * t131;
t181 = -rSges(5,1) * t198 - rSges(5,2) * t196 - t143;
t180 = rSges(4,1) * t192 + rSges(4,2) * t189 + t155 * rSges(4,3);
t142 = t155 * qJ(2);
t179 = t142 + t185;
t178 = t140 * (-rSges(7,2) - pkin(8));
t177 = t140 * (-rSges(6,3) - pkin(8));
t175 = t206 * t152;
t174 = -t116 - t211;
t172 = rSges(5,1) * t139 + rSges(5,2) * t140;
t171 = -t110 * rSges(6,1) - t108 * rSges(6,2);
t160 = -t149 * t155 + t134 + t184;
t159 = -t22 / 0.2e1 - t20 / 0.2e1 - t29 / 0.2e1 - t28 / 0.2e1;
t158 = t23 / 0.2e1 + t21 / 0.2e1 + t31 / 0.2e1 + t30 / 0.2e1;
t157 = -t152 * pkin(1) + t130 + t179;
t156 = t129 + t160;
t126 = rSges(2,1) * t155 - t152 * rSges(2,2);
t124 = -t152 * rSges(2,1) - rSges(2,2) * t155;
t115 = rSges(5,1) * t140 - rSges(5,2) * t139;
t104 = -rSges(3,2) * t155 + t152 * rSges(3,3) + t184;
t103 = rSges(3,3) * t155 + t142 + (rSges(3,2) - pkin(1)) * t152;
t81 = (-t115 - t211) * t155;
t80 = t115 * t152 + t135;
t79 = rSges(6,3) * t139 + (rSges(6,1) * t153 - rSges(6,2) * t150) * t140;
t69 = pkin(7) * t155 + t180 + t184;
t68 = t142 + t220 + (-rSges(4,3) - pkin(1) - pkin(7)) * t152;
t62 = t160 - t181;
t61 = t172 * t155 + (-rSges(5,3) - pkin(1)) * t152 + t179;
t60 = -t152 * t180 + (t152 * rSges(4,3) - t220) * t155;
t59 = rSges(6,3) * t194 - t171;
t57 = -rSges(6,3) * t196 + t203;
t43 = (t174 - t79) * t155;
t42 = t152 * t79 + t186;
t41 = (t174 - t206) * t155;
t40 = t175 + t186;
t39 = t152 * t177 + t156 + t203;
t38 = t155 * t177 + t157 + t171;
t37 = -t139 * t59 + t194 * t79;
t36 = t139 * t57 + t196 * t79;
t35 = t89 - t172 * t148 + (-t105 + t181 + t143) * t152;
t32 = (-t152 * t59 - t155 * t57) * t140;
t27 = t152 * t178 + t156 + t221;
t26 = t155 * t178 + t157 - t232;
t25 = -t139 * t208 + t194 * t206;
t24 = t139 * t209 + t140 * t175;
t19 = t155 * t59 + (-t57 + t187) * t152 + t205;
t10 = (-t152 * t208 - t155 * t209) * t140;
t9 = t208 * t155 + (t187 - t209) * t152 + t205;
t8 = t18 * t152 + t155 * t17;
t7 = t15 * t155 + t16 * t152;
t6 = t13 * t155 + t14 * t152;
t5 = t11 * t155 + t12 * t152;
t1 = [Icges(4,1) * t154 ^ 2 + Icges(3,1) + Icges(2,3) + (Icges(5,1) * t140 - t150 * t75) * t140 + m(7) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t38 ^ 2 + t39 ^ 2) + m(5) * (t61 ^ 2 + t62 ^ 2) + m(4) * (t68 ^ 2 + t69 ^ 2) + m(3) * (t103 ^ 2 + t104 ^ 2) + m(2) * (t124 ^ 2 + t126 ^ 2) + t228 + (-0.2e1 * Icges(4,4) * t154 + Icges(4,2) * t151) * t151 + (-0.2e1 * Icges(5,4) * t140 + Icges(5,2) * t139) * t139; m(7) * (t152 * t26 - t155 * t27) + m(6) * (t152 * t38 - t155 * t39) + m(5) * (t152 * t61 - t155 * t62) + m(4) * (t152 * t68 - t155 * t69) + m(3) * (t152 * t103 - t104 * t155); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t131 + t182; m(7) * (t26 * t40 + t27 * t41) + m(6) * (t38 * t42 + t39 * t43) + m(5) * (t61 * t80 + t62 * t81) + (t155 * t231 - t69 * t212 - t159) * t155 + (t152 * t231 + t68 * t212 + t158) * t152; m(5) * (t80 * t152 - t155 * t81) + m(6) * (t42 * t152 - t155 * t43) + m(7) * (t40 * t152 - t155 * t41) + t131 * t212; m(7) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(6) * (t19 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t35 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t125 ^ 2 * t131 + t60 ^ 2) + (t148 * t223 + t5 + t6) * t155 + (t7 + t8 + t222 * t147 + (t152 * t223 + t155 * t222) * t155) * t152; m(7) * (t152 * t27 + t155 * t26) + m(6) * (t152 * t39 + t155 * t38) + m(5) * (t152 * t62 + t155 * t61); 0; m(7) * (t152 * t41 + t155 * t40) + m(6) * (t152 * t43 + t155 * t42) + m(5) * (t152 * t81 + t155 * t80); t182; m(7) * (t24 * t27 + t25 * t26) + m(6) * (t36 * t39 + t37 * t38) + (t152 * t159 + t155 * t158) * t140 + t210; m(6) * (t37 * t152 - t155 * t36) + m(7) * (t25 * t152 - t155 * t24); m(7) * (t10 * t9 + t24 * t41 + t25 * t40) + m(6) * (t19 * t32 + t36 * t43 + t37 * t42) + ((t8 / 0.2e1 + t7 / 0.2e1) * t155 + (-t6 / 0.2e1 - t5 / 0.2e1) * t152) * t140 + (t152 * t224 + t155 * t225) * t139 / 0.2e1 + t226 * t152 / 0.2e1 + t227 * t155 / 0.2e1; m(6) * (t36 * t152 + t155 * t37) + m(7) * (t24 * t152 + t155 * t25); t210 * t139 + m(7) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t32 ^ 2 + t36 ^ 2 + t37 ^ 2) + (t226 * t155 - t227 * t152 + (-t152 * t225 + t155 * t224) * t139) * t140; m(7) * (t106 * t26 - t108 * t27); m(7) * (t106 * t152 + t108 * t155); m(7) * (t106 * t40 - t108 * t41 + t197 * t9); m(7) * (t106 * t155 - t108 * t152); m(7) * (t10 * t197 + t106 * t25 - t108 * t24); m(7) * (t140 ^ 2 * t150 ^ 2 + t106 ^ 2 + t108 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
