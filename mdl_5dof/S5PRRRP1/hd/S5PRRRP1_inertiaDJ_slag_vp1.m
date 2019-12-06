% Calculate time derivative of joint inertia matrix for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:48
% EndTime: 2019-12-05 16:39:56
% DurationCPUTime: 3.56s
% Computational Cost: add. (4920->268), mult. (3752->371), div. (0->0), fcn. (2756->6), ass. (0->159)
t147 = cos(qJ(4));
t146 = sin(qJ(4));
t218 = Icges(6,4) * t146;
t220 = Icges(5,4) * t146;
t254 = t218 + t220 + (Icges(5,2) + Icges(6,2)) * t147;
t217 = Icges(6,4) * t147;
t219 = Icges(5,4) * t147;
t253 = t217 + t219 + (Icges(5,1) + Icges(6,1)) * t146;
t170 = -Icges(6,2) * t146 + t217;
t171 = -Icges(5,2) * t146 + t219;
t172 = Icges(6,1) * t147 - t218;
t173 = Icges(5,1) * t147 - t220;
t252 = (t170 + t171) * t147 + (t172 + t173) * t146;
t249 = t254 * t146 - t253 * t147;
t143 = pkin(8) + qJ(2);
t142 = qJ(3) + t143;
t138 = cos(t142);
t144 = qJD(2) + qJD(3);
t213 = t138 * t144;
t248 = -rSges(6,1) - pkin(4);
t118 = Icges(6,5) * t146 + Icges(6,6) * t147;
t119 = Icges(5,5) * t146 + Icges(5,6) * t147;
t247 = t118 + t119;
t137 = sin(t142);
t203 = qJD(4) * t147;
t208 = t144 * t146;
t244 = -t137 * t203 - t138 * t208;
t168 = Icges(6,5) * t147 - Icges(6,6) * t146;
t169 = Icges(5,5) * t147 - Icges(5,6) * t146;
t243 = t249 * t144 + (t168 + t169) * qJD(4);
t128 = t137 * rSges(6,3);
t139 = pkin(4) * t147 + pkin(3);
t145 = -qJ(5) - pkin(7);
t210 = t138 * t147;
t211 = t138 * t146;
t52 = rSges(6,1) * t210 - rSges(6,2) * t211 - t137 * t145 + t138 * t139 + t128;
t162 = t171 * t138;
t68 = Icges(5,6) * t137 + t162;
t164 = t173 * t138;
t72 = Icges(5,5) * t137 + t164;
t176 = t146 * t68 - t147 * t72;
t242 = t137 * t176;
t161 = t170 * t138;
t66 = Icges(6,6) * t137 + t161;
t163 = t172 * t138;
t70 = Icges(6,5) * t137 + t163;
t180 = t146 * t66 - t147 * t70;
t241 = t137 * t180;
t67 = -Icges(5,6) * t138 + t137 * t171;
t71 = -Icges(5,5) * t138 + t137 * t173;
t178 = t146 * t67 - t147 * t71;
t240 = t138 * t178;
t65 = -Icges(6,6) * t138 + t137 * t170;
t69 = -Icges(6,5) * t138 + t137 * t172;
t182 = t146 * t65 - t147 * t69;
t239 = t138 * t182;
t212 = t138 * t145;
t215 = t137 * t146;
t238 = rSges(6,2) * t215 + t138 * rSges(6,3) - t212;
t204 = qJD(4) * t146;
t151 = t252 * qJD(4) + t253 * t203 - t254 * t204;
t235 = 2 * m(4);
t234 = 2 * m(5);
t233 = 2 * m(6);
t225 = rSges(5,1) * t147;
t104 = (-rSges(5,2) * t146 + t225) * qJD(4);
t230 = m(5) * t104;
t125 = t146 * rSges(5,1) + rSges(5,2) * t147;
t229 = m(5) * t125;
t140 = sin(t143);
t228 = pkin(2) * t140;
t134 = t138 * pkin(7);
t214 = t137 * t147;
t227 = t134 + t137 * (-pkin(3) + t139) + rSges(6,1) * t214 - t238;
t205 = -t138 * pkin(3) - t137 * pkin(7);
t226 = t205 + t52;
t224 = rSges(6,1) * t147;
t223 = pkin(2) * qJD(2);
t129 = t137 * rSges(5,3);
t222 = t147 * rSges(6,2);
t198 = t137 * t208;
t221 = rSges(5,2) * t198 + rSges(5,3) * t213;
t216 = t137 * t144;
t209 = t144 * t145;
t206 = rSges(5,2) * t215 + t138 * rSges(5,3);
t202 = t140 * t223;
t141 = cos(t143);
t201 = t141 * t223;
t196 = t137 * t204;
t200 = -rSges(5,1) * t196 + t244 * rSges(5,2);
t199 = rSges(6,2) * t198 + rSges(6,3) * t213 + qJD(5) * t137;
t194 = t138 * t203;
t191 = -pkin(3) - t225;
t124 = t146 * rSges(6,1) + t222;
t190 = -pkin(4) * t146 - t124;
t189 = -t139 - t224;
t84 = t138 * rSges(4,1) - rSges(4,2) * t137;
t80 = -rSges(4,1) * t213 + rSges(4,2) * t216;
t83 = -rSges(4,1) * t137 - rSges(4,2) * t138;
t183 = t146 * t69 + t147 * t65;
t181 = t146 * t70 + t147 * t66;
t179 = t146 * t71 + t147 * t67;
t177 = t146 * t72 + t147 * t68;
t175 = t189 * t137;
t76 = rSges(5,1) * t210 - rSges(5,2) * t211 + t129;
t165 = t244 * rSges(6,2) - qJD(5) * t138 - t137 * t209 + t248 * t196;
t79 = t83 * t144;
t160 = t169 * t138;
t159 = t168 * t138;
t56 = t76 - t205;
t158 = (t248 * t146 - t222) * qJD(4);
t55 = t137 * t191 + t134 + t206;
t153 = Icges(5,3) * t144 - qJD(4) * t119;
t152 = Icges(6,3) * t144 - qJD(4) * t118;
t51 = t175 + t238;
t150 = (t243 * t137 + (-t176 - t180) * qJD(4) - t252 * t216) * t137 / 0.2e1 - (-t243 * t138 + (-t178 - t182) * qJD(4)) * t138 / 0.2e1 + (-t137 * t249 - t247 * t138 + t179 + t183) * t216 / 0.2e1 + (t247 * t137 - t138 * t249 + t177 + t181) * t213 / 0.2e1 - ((t161 + t162) * t147 + (t163 + t164) * t146) * t213 / 0.2e1;
t26 = (t191 * t138 + (-rSges(5,3) - pkin(7)) * t137) * t144 - t200;
t22 = (t138 * t189 - t128) * t144 - t165;
t113 = pkin(7) * t213;
t25 = -rSges(5,2) * t194 - pkin(3) * t216 + t113 + (-t138 * t204 - t144 * t214) * rSges(5,1) + t221;
t21 = t144 * t175 + (t158 - t209) * t138 + t199;
t136 = pkin(2) * t141;
t103 = (-rSges(6,2) * t146 + t224) * qJD(4);
t82 = t136 + t84;
t81 = t83 - t228;
t78 = t190 * t138;
t77 = t190 * t137;
t74 = rSges(5,1) * t214 - t206;
t64 = Icges(5,3) * t137 + t160;
t63 = -Icges(5,3) * t138 + t137 * t169;
t62 = Icges(6,3) * t137 + t159;
t61 = -Icges(6,3) * t138 + t137 * t168;
t60 = t80 - t201;
t59 = t79 - t202;
t54 = t136 + t56;
t53 = t55 - t228;
t50 = t136 + t52;
t49 = t51 - t228;
t40 = t137 * t153 + t144 * t160;
t39 = t138 * t153 - t169 * t216;
t38 = t137 * t152 + t144 * t159;
t37 = t138 * t152 - t168 * t216;
t36 = t244 * pkin(4) - t137 * t103 - t124 * t213;
t35 = t124 * t216 - t138 * t103 + (-t194 + t198) * pkin(4);
t24 = t26 - t201;
t23 = t25 - t202;
t20 = t137 * t64 - t176 * t138;
t19 = t137 * t63 - t240;
t18 = t137 * t62 - t180 * t138;
t17 = t137 * t61 - t239;
t16 = -t138 * t64 - t242;
t15 = -t178 * t137 - t138 * t63;
t14 = -t138 * t62 - t241;
t13 = -t182 * t137 - t138 * t61;
t12 = t22 - t201;
t11 = t21 - t202;
t2 = ((-t76 + t129) * t144 + t200) * t137 + (-qJD(4) * t125 * t138 + t144 * t74 + t221) * t138;
t1 = (((rSges(6,3) - pkin(7)) * t137 - t226) * t144 + t165) * t137 + (-t113 + (-t212 + t227) * t144 + t138 * t158 + t199) * t138;
t3 = [0; 0; (t11 * t50 + t12 * t49) * t233 + (t23 * t54 + t24 * t53) * t234 + (t59 * t82 + t60 * t81) * t235 + t151; 0; m(6) * (t11 * t52 + t12 * t51 + t21 * t50 + t22 * t49) + m(5) * (t23 * t56 + t24 * t55 + t25 * t54 + t26 * t53) + m(4) * (t59 * t84 + t60 * t83 + t79 * t82 + t80 * t81) + t151; (t21 * t52 + t22 * t51) * t233 + (t25 * t56 + t26 * t55) * t234 + (t79 * t84 + t80 * t83) * t235 + t151; m(5) * t2 + m(6) * t1; ((-t144 * t54 - t24) * t138 + (t144 * t53 - t23) * t137) * t229 + m(6) * (t11 * t77 + t12 * t78 + t35 * t49 + t36 * t50) + (-t137 * t54 - t138 * t53) * t230 + t150; m(6) * (t21 * t77 + t22 * t78 + t35 * t51 + t36 * t52) + ((-t144 * t56 - t26) * t138 + (t144 * t55 - t25) * t137) * t229 + (-t137 * t56 - t138 * t55) * t230 + t150; ((t137 * t74 + t138 * t76) * t2 + (t137 ^ 2 + t138 ^ 2) * t125 * t104) * t234 + t137 * ((t137 * t39 + (t19 + t242) * t144) * t137 + (t20 * t144 + (t203 * t67 + t204 * t71) * t138 + (-t177 * qJD(4) - t144 * t178 - t40) * t137) * t138) - t138 * ((t138 * t40 + (t16 + t240) * t144) * t138 + (t15 * t144 + (-t203 * t68 - t204 * t72) * t137 + (t179 * qJD(4) - t144 * t176 - t39) * t138) * t137) + (t78 * t35 + t77 * t36 + (t137 * t227 + t138 * t226) * t1) * t233 + t137 * ((t137 * t37 + (t17 + t241) * t144) * t137 + (t18 * t144 + (t203 * t65 + t204 * t69) * t138 + (-t181 * qJD(4) - t144 * t182 - t38) * t137) * t138) - t138 * ((t138 * t38 + (t14 + t239) * t144) * t138 + (t13 * t144 + (-t203 * t66 - t204 * t70) * t137 + (t183 * qJD(4) - t144 * t180 - t37) * t138) * t137) + ((-t13 - t15) * t138 + (t14 + t16) * t137) * t216 + ((-t17 - t19) * t138 + (t18 + t20) * t137) * t213; 0; m(6) * ((t144 * t49 - t11) * t138 + (t144 * t50 + t12) * t137); m(6) * ((t144 * t51 - t21) * t138 + (t144 * t52 + t22) * t137); m(6) * ((t144 * t78 - t36) * t138 + (t144 * t77 + t35) * t137); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
