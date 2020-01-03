% Calculate time derivative of joint inertia matrix for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:01
% EndTime: 2019-12-31 17:14:07
% DurationCPUTime: 3.71s
% Computational Cost: add. (3563->259), mult. (4058->370), div. (0->0), fcn. (3016->6), ass. (0->155)
t145 = cos(qJ(3));
t143 = sin(qJ(3));
t212 = Icges(5,5) * t143;
t214 = Icges(4,4) * t143;
t252 = t212 - t214 + (-Icges(4,2) - Icges(5,3)) * t145;
t211 = Icges(5,5) * t145;
t213 = Icges(4,4) * t145;
t251 = -t211 + t213 + (Icges(4,1) + Icges(5,1)) * t143;
t168 = Icges(5,3) * t143 + t211;
t171 = -Icges(4,2) * t143 + t213;
t172 = Icges(5,1) * t145 + t212;
t173 = Icges(4,1) * t145 - t214;
t250 = (t172 + t173) * t143 - (t168 - t171) * t145;
t239 = rSges(5,3) + qJ(4);
t245 = rSges(5,1) + pkin(3);
t241 = t143 * t239 + t145 * t245;
t242 = t252 * t143 + t251 * t145;
t142 = qJ(1) + qJ(2);
t138 = sin(t142);
t200 = qJD(3) * t143;
t193 = t138 * t200;
t198 = qJD(4) * t143;
t199 = qJD(3) * t145;
t247 = -(t239 * t199 + t198) * t138 + t245 * t193;
t139 = cos(t142);
t141 = qJD(1) + qJD(2);
t208 = t139 * t141;
t109 = Icges(4,5) * t143 + Icges(4,6) * t145;
t110 = Icges(5,4) * t143 - Icges(5,6) * t145;
t244 = t109 + t110;
t169 = Icges(4,5) * t145 - Icges(4,6) * t143;
t170 = Icges(5,4) * t145 + Icges(5,6) * t143;
t240 = t242 * t141 + (-t169 - t170) * qJD(3);
t163 = t171 * t139;
t68 = Icges(4,6) * t138 + t163;
t165 = t173 * t139;
t72 = Icges(4,5) * t138 + t165;
t174 = t143 * t68 - t145 * t72;
t238 = t138 * t174;
t160 = t168 * t139;
t62 = Icges(5,6) * t138 + t160;
t164 = t172 * t139;
t70 = Icges(5,4) * t138 + t164;
t178 = t143 * t62 + t145 * t70;
t237 = t138 * t178;
t67 = -Icges(4,6) * t139 + t138 * t171;
t71 = -Icges(4,5) * t139 + t138 * t173;
t176 = t143 * t67 - t145 * t71;
t236 = t139 * t176;
t61 = -Icges(5,6) * t139 + t138 * t168;
t69 = -Icges(5,4) * t139 + t138 * t172;
t180 = t143 * t61 + t145 * t69;
t235 = t139 * t180;
t232 = 2 * m(3);
t231 = 2 * m(4);
t230 = 2 * m(5);
t220 = rSges(4,2) * t143;
t221 = rSges(4,1) * t145;
t101 = (-t220 + t221) * qJD(3);
t226 = m(4) * t101;
t120 = rSges(4,1) * t143 + rSges(4,2) * t145;
t225 = m(4) * t120;
t144 = sin(qJ(1));
t224 = pkin(1) * t144;
t131 = t139 * rSges(5,2);
t223 = t241 * t138 - t131;
t129 = t138 * rSges(5,2);
t206 = t139 * t145;
t207 = t139 * t143;
t222 = t245 * t206 + t239 * t207 + t129;
t219 = pkin(1) * qJD(1);
t128 = t138 * rSges(4,3);
t204 = t245 * t143 - t239 * t145;
t58 = t204 * t139;
t218 = t141 * t58;
t216 = -t241 * qJD(3) + qJD(4) * t145;
t123 = t138 * t220;
t215 = rSges(4,3) * t208 + t141 * t123;
t210 = t138 * t141;
t209 = t138 * t145;
t203 = t139 * rSges(4,3) + t123;
t202 = t139 * pkin(2) + t138 * pkin(6);
t201 = t138 ^ 2 + t139 ^ 2;
t197 = rSges(4,2) * t207;
t196 = t144 * t219;
t146 = cos(qJ(1));
t195 = t146 * t219;
t194 = -t138 * rSges(4,2) * t199 - rSges(4,1) * t193 - t141 * t197;
t192 = t139 * t200;
t191 = t139 * t199;
t188 = -pkin(2) - t221;
t91 = t139 * rSges(3,1) - rSges(3,2) * t138;
t78 = -rSges(3,1) * t208 + rSges(3,2) * t210;
t90 = -rSges(3,1) * t138 - rSges(3,2) * t139;
t181 = t143 * t69 - t145 * t61;
t179 = t143 * t70 - t145 * t62;
t177 = t143 * t71 + t145 * t67;
t175 = t143 * t72 + t145 * t68;
t76 = rSges(4,1) * t206 + t128 - t197;
t77 = t90 * t141;
t162 = t170 * t139;
t161 = t169 * t139;
t52 = t202 + t222;
t56 = t76 + t202;
t134 = t139 * pkin(6);
t55 = t138 * t188 + t134 + t203;
t159 = -pkin(2) - t241;
t158 = rSges(5,2) * t208 + t139 * t198 + t239 * t191 - t192 * t245;
t157 = t159 * t138;
t154 = Icges(5,2) * t141 - qJD(3) * t110;
t151 = Icges(4,3) * t141 - qJD(3) * t109;
t150 = t250 * qJD(3) + t251 * t199 + t252 * t200;
t51 = t131 + t134 + t157;
t149 = (-t240 * t138 + (-t174 + t178) * qJD(3) - t250 * t210) * t138 / 0.2e1 - (t240 * t139 + (-t176 + t180) * qJD(3)) * t139 / 0.2e1 + (t242 * t138 - t244 * t139 + t177 + t181) * t210 / 0.2e1 + (t244 * t138 + t242 * t139 + t175 + t179) * t208 / 0.2e1 - ((-t160 + t163) * t145 + (t164 + t165) * t143) * t208 / 0.2e1;
t28 = (t188 * t139 + (-rSges(4,3) - pkin(6)) * t138) * t141 - t194;
t117 = pkin(6) * t208;
t27 = -rSges(4,2) * t191 - pkin(2) * t210 + t117 + (-t141 * t209 - t192) * rSges(4,1) + t215;
t12 = t141 * t157 + t117 + t158;
t13 = ((-rSges(5,2) - pkin(6)) * t138 + t159 * t139) * t141 + t247;
t140 = t146 * pkin(1);
t82 = t140 + t91;
t81 = t90 - t224;
t74 = rSges(4,1) * t209 - t203;
t66 = Icges(5,2) * t138 + t162;
t65 = -Icges(5,2) * t139 + t138 * t170;
t64 = Icges(4,3) * t138 + t161;
t63 = -Icges(4,3) * t139 + t138 * t169;
t60 = t78 - t195;
t59 = t77 - t196;
t57 = t204 * t138;
t54 = t140 + t56;
t53 = t55 - t224;
t44 = t138 * t154 + t141 * t162;
t43 = t139 * t154 - t170 * t210;
t42 = t138 * t151 + t141 * t161;
t41 = t139 * t151 - t169 * t210;
t38 = t140 + t52;
t37 = t51 - t224;
t26 = t138 * t216 - t218;
t25 = t139 * t216 + t204 * t210;
t24 = t28 - t195;
t23 = t27 - t196;
t22 = t138 * t64 - t174 * t139;
t21 = t138 * t63 - t236;
t20 = t138 * t66 + t178 * t139;
t19 = t138 * t65 + t235;
t18 = -t139 * t64 - t238;
t17 = -t176 * t138 - t139 * t63;
t16 = -t139 * t66 + t237;
t15 = t180 * t138 - t139 * t65;
t14 = t138 * t223 + t139 * t222;
t11 = t13 - t195;
t10 = t12 - t196;
t1 = (t141 * t223 + t158) * t139 + ((t129 - t222) * t141 - t247) * t138;
t2 = [(t59 * t82 + t60 * t81) * t232 + (t23 * t54 + t24 * t53) * t231 + (t10 * t38 + t11 * t37) * t230 + t150; m(3) * (t59 * t91 + t60 * t90 + t77 * t82 + t78 * t81) + m(4) * (t23 * t56 + t24 * t55 + t27 * t54 + t28 * t53) + m(5) * (t10 * t52 + t11 * t51 + t12 * t38 + t13 * t37) + t150; (t77 * t91 + t78 * t90) * t232 + (t27 * t56 + t28 * t55) * t231 + (t12 * t52 + t13 * t51) * t230 + t150; t149 + m(5) * (-t10 * t57 - t11 * t58 + t25 * t37 + t26 * t38) + ((-t141 * t54 - t24) * t139 + (t141 * t53 - t23) * t138) * t225 + (-t138 * t54 - t139 * t53) * t226; t149 + m(5) * (-t12 * t57 - t13 * t58 + t25 * t51 + t26 * t52) + ((-t141 * t56 - t28) * t139 + (t141 * t55 - t27) * t138) * t225 + (-t138 * t56 - t139 * t55) * t226; ((t138 * t74 + t139 * t76) * (((-t76 + t128) * t141 + t194) * t138 + (-qJD(3) * t120 * t139 + t141 * t74 + t215) * t139) + t201 * t120 * t101) * t231 + t138 * ((t138 * t41 + (t21 + t238) * t141) * t138 + (t22 * t141 + (t199 * t67 + t200 * t71) * t139 + (-t175 * qJD(3) - t141 * t176 - t42) * t138) * t139) - t139 * ((t139 * t42 + (t18 + t236) * t141) * t139 + (t17 * t141 + (-t199 * t68 - t200 * t72) * t138 + (t177 * qJD(3) - t141 * t174 - t41) * t139) * t138) + (t1 * t14 - t25 * t58 - t26 * t57) * t230 + t138 * ((t138 * t43 + (t19 - t237) * t141) * t138 + (t20 * t141 + (-t199 * t61 + t200 * t69) * t139 + (-t179 * qJD(3) + t141 * t180 - t44) * t138) * t139) - t139 * ((t139 * t44 + (t16 - t235) * t141) * t139 + (t15 * t141 + (t199 * t62 - t200 * t70) * t138 + (t181 * qJD(3) + t141 * t178 - t43) * t139) * t138) + ((-t15 - t17) * t139 + (t16 + t18) * t138) * t210 + ((-t19 - t21) * t139 + (t20 + t22) * t138) * t208; m(5) * ((t138 * t38 + t139 * t37) * t199 + ((t141 * t38 + t11) * t139 + (-t141 * t37 + t10) * t138) * t143); m(5) * ((t138 * t52 + t139 * t51) * t199 + ((t141 * t52 + t13) * t139 + (-t141 * t51 + t12) * t138) * t143); m(5) * ((-t1 + (-t138 * t57 - t139 * t58) * qJD(3)) * t145 + (qJD(3) * t14 + (-t141 * t57 + t25) * t139 + (t26 + t218) * t138) * t143); (-0.1e1 + t201) * t143 * t199 * t230;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
