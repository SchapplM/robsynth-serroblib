% Calculate time derivative of joint inertia matrix for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:28
% EndTime: 2020-01-03 11:27:42
% DurationCPUTime: 3.32s
% Computational Cost: add. (6403->302), mult. (4524->436), div. (0->0), fcn. (3474->10), ass. (0->183)
t122 = qJ(1) + pkin(8);
t112 = sin(t122);
t114 = cos(t122);
t120 = pkin(9) + qJ(4);
t115 = qJ(5) + t120;
t105 = sin(t115);
t217 = rSges(6,2) * t105;
t106 = cos(t115);
t219 = rSges(6,1) * t106;
t158 = -t217 + t219;
t67 = -rSges(6,3) * t114 + t158 * t112;
t141 = Icges(6,5) * t106 - Icges(6,6) * t105;
t233 = Icges(6,3) * t112 + t141 * t114;
t243 = t233 * qJD(1);
t118 = cos(qJ(1)) * pkin(1);
t108 = qJD(1) * t118;
t111 = sin(t120);
t113 = cos(t120);
t202 = Icges(5,4) * t111;
t149 = Icges(5,1) * t113 - t202;
t228 = Icges(5,5) * t112 + t149 * t114;
t201 = Icges(5,4) * t113;
t146 = -Icges(5,2) * t111 + t201;
t230 = Icges(5,6) * t112 + t146 * t114;
t151 = -t111 * t230 + t113 * t228;
t241 = t151 * t112;
t200 = Icges(6,4) * t105;
t147 = Icges(6,1) * t106 - t200;
t229 = Icges(6,5) * t112 + t147 * t114;
t199 = Icges(6,4) * t106;
t144 = -Icges(6,2) * t105 + t199;
t231 = Icges(6,6) * t112 + t144 * t114;
t156 = -t105 * t231 + t106 * t229;
t240 = t156 * t112;
t218 = rSges(5,2) * t111;
t221 = rSges(5,1) * t113;
t159 = -t218 + t221;
t239 = t159 * t114;
t125 = -pkin(6) - qJ(3);
t119 = -pkin(7) + t125;
t184 = t119 - t125;
t238 = t184 * t112;
t100 = t112 * t221;
t124 = cos(pkin(9));
t107 = t124 * pkin(3) + pkin(2);
t163 = t107 - t218;
t117 = sin(qJ(1)) * pkin(1);
t206 = rSges(5,3) - t125;
t234 = t206 * t114 - t117;
t50 = t163 * t112 + t100 - t234;
t99 = t114 * t107;
t51 = t206 * t112 + t118 + t239 + t99;
t236 = -t112 * t51 + t114 * t50;
t235 = rSges(5,3) * t114 - t100;
t143 = Icges(5,5) * t113 - Icges(5,6) * t111;
t232 = Icges(5,3) * t112 + t143 * t114;
t121 = qJD(4) + qJD(5);
t78 = t144 * t121;
t79 = t147 * t121;
t87 = Icges(6,5) * t105 + Icges(6,6) * t106;
t88 = Icges(6,2) * t106 + t200;
t89 = Icges(6,1) * t105 + t199;
t227 = (t121 * t89 + t78) * t105 + (t121 * t88 - t79) * t106 - qJD(1) * t87;
t140 = rSges(4,1) * t124 - rSges(4,2) * sin(pkin(9)) + pkin(2);
t207 = rSges(4,3) + qJ(3);
t54 = t140 * t112 - t207 * t114 + t117;
t226 = 2 * m(5);
t225 = 2 * m(6);
t109 = t112 ^ 2;
t110 = t114 ^ 2;
t95 = rSges(5,1) * t111 + rSges(5,2) * t113;
t224 = m(5) * t95;
t223 = -t112 / 0.2e1;
t222 = -t114 / 0.2e1;
t220 = rSges(6,1) * t105;
t73 = -Icges(5,5) * t114 + t149 * t112;
t214 = t111 * t73;
t213 = t111 * t228;
t80 = t158 * t121;
t211 = t112 * t80;
t98 = pkin(4) * t113 + t107;
t82 = t112 * t98;
t71 = -Icges(5,6) * t114 + t146 * t112;
t210 = t113 * t71;
t209 = t113 * t230;
t205 = rSges(6,3) - t119;
t204 = t114 * t119 + t82;
t182 = qJD(1) * t114;
t183 = qJD(1) * t112;
t203 = rSges(6,3) * t183 + t182 * t219;
t61 = -Icges(6,3) * t114 + t141 * t112;
t192 = qJD(1) * t61;
t69 = -Icges(5,3) * t114 + t143 * t112;
t191 = qJD(1) * t69;
t189 = t105 * t121;
t188 = t106 * t121;
t187 = t112 * t121;
t186 = t114 * t121;
t185 = t109 + t110;
t181 = qJD(4) * t111;
t180 = qJD(4) * t112;
t179 = qJD(4) * t113;
t178 = qJD(4) * t114;
t68 = -rSges(6,3) * t112 - t158 * t114;
t177 = t112 * (-t187 * t220 + (-t105 * t182 - t106 * t187) * rSges(6,2) + t203) + t68 * t183 + t67 * t182;
t176 = pkin(4) * t181;
t90 = rSges(6,2) * t106 + t220;
t175 = t90 * t183;
t63 = -Icges(6,6) * t114 + t144 * t112;
t65 = -Icges(6,5) * t114 + t147 * t112;
t157 = t105 * t63 - t106 * t65;
t13 = -t157 * t112 - t61 * t114;
t14 = t114 * t233 - t240;
t137 = t157 * t114;
t15 = -t112 * t61 + t137;
t16 = t112 * t233 + t156 * t114;
t37 = t65 * qJD(1) + t89 * t186;
t167 = t121 * t231 + t37;
t35 = t63 * qJD(1) + t88 * t186;
t169 = -t121 * t229 + t35;
t33 = t87 * t186 + t192;
t34 = -t87 * t187 + t243;
t36 = t231 * qJD(1) - t88 * t187;
t38 = t229 * qJD(1) - t89 * t187;
t174 = -t112 * ((t112 * t33 + (t15 + t240) * qJD(1)) * t112 + (-t16 * qJD(1) + (-t105 * t36 + t106 * t38 - t63 * t188 - t65 * t189 + t192) * t114 + (t34 + t167 * t106 - t169 * t105 + (t157 - t233) * qJD(1)) * t112) * t114) + (-t112 * t14 - t114 * t13) * t183;
t173 = pkin(4) * t111 + t90;
t166 = t121 * t63 - t38;
t168 = t121 * t65 + t36;
t2 = (t114 * t34 + (-t14 + t137) * qJD(1)) * t114 + (t13 * qJD(1) + (t105 * t35 - t106 * t37 - t188 * t231 - t189 * t229 + t243) * t112 + (t33 + t166 * t106 + t168 * t105 + (t156 - t61) * qJD(1)) * t114) * t112;
t5 = -t112 * t16 - t114 * t15;
t172 = -qJD(1) * t5 - t2;
t162 = -qJD(3) * t114 + t108;
t155 = t105 * t88 - t106 * t89;
t152 = t111 * t71 - t113 * t73;
t40 = t117 + t67 + t204;
t139 = t158 + t98;
t41 = t205 * t112 + t139 * t114 + t118;
t150 = t112 * t41 - t114 * t40;
t148 = Icges(5,1) * t111 + t201;
t145 = Icges(5,2) * t113 + t202;
t142 = Icges(5,5) * t111 + Icges(5,6) * t113;
t131 = -t155 * qJD(1) - t141 * t121;
t138 = (t167 * t105 + t169 * t106 + t131 * t112 + t227 * t114) * t223 + (-t166 * t105 + t168 * t106 - t227 * t112 + t131 * t114) * t222 + (t105 * t65 + t106 * t63 - t155 * t112 - t114 * t87) * t183 / 0.2e1 - (-t105 * t229 - t106 * t231 - t112 * t87 + t155 * t114) * t182 / 0.2e1;
t136 = t152 * t114;
t135 = t95 * qJD(4);
t134 = qJD(4) * t148;
t133 = qJD(4) * t145;
t132 = t114 * t135;
t130 = rSges(5,3) * t183 - t112 * t135 + t182 * t221;
t129 = -t90 * t121 - t176;
t128 = t207 * t112 + t140 * t114;
t104 = qJD(3) * t112;
t86 = t159 * qJD(4);
t81 = t98 * t182;
t76 = -rSges(5,3) * t112 - t239;
t75 = -t112 * t218 - t235;
t60 = t173 * t114;
t59 = t173 * t112;
t58 = t112 * t67;
t55 = t118 + t128;
t53 = -t114 * t98 + t238 + t99;
t52 = -t107 * t112 - t114 * t125 + t204;
t49 = t128 * qJD(1) + t162;
t48 = -t54 * qJD(1) + t104;
t43 = t232 * qJD(1) - t142 * t180;
t42 = t142 * t178 + t191;
t39 = t67 * qJD(1) + t90 * t186;
t31 = -t114 * t68 + t58;
t30 = -t90 * t182 - t211 + (-t111 * t182 - t112 * t179) * pkin(4);
t29 = -t175 + t114 * t80 + (-t111 * t183 + t113 * t178) * pkin(4);
t24 = (-t112 * t125 + t163 * t114) * qJD(1) + t130 + t162;
t23 = t104 - t132 + ((-t107 - t159) * t112 + t234) * qJD(1);
t22 = t112 * t232 + t151 * t114;
t21 = -t112 * t69 + t136;
t20 = t114 * t232 - t241;
t19 = -t152 * t112 - t69 * t114;
t18 = t108 + t81 + (-qJD(1) * t217 - qJD(3)) * t114 + (-qJD(1) * t119 + t129) * t112 + t203;
t17 = t104 + t129 * t114 + (-t139 * t112 + t205 * t114 - t117) * qJD(1);
t12 = t112 * t52 + t58 + (-t53 - t68) * t114;
t11 = (qJD(1) * t76 + t130) * t112 + (-t132 + (t75 + t235) * qJD(1)) * t114;
t10 = -t114 * t39 + t177;
t3 = (-t112 * t176 + t81 + (t53 - t238) * qJD(1)) * t112 + (-t114 * t176 - t39 + (-t184 * t114 + t52 - t82) * qJD(1)) * t114 + t177;
t1 = [0.2e1 * m(4) * (t48 * t55 + t49 * t54) + t113 * t134 + t149 * t181 - t111 * t133 + t146 * t179 + (t23 * t51 + t24 * t50) * t226 + t89 * t188 + t105 * t79 - t88 * t189 + t106 * t78 + (t17 * t41 + t18 * t40) * t225; 0; 0; m(4) * (-t112 * t49 - t114 * t48 + (t112 * t55 - t114 * t54) * qJD(1)) + m(5) * (-t236 * qJD(1) - t112 * t24 - t114 * t23) + m(6) * (t150 * qJD(1) - t112 * t18 - t114 * t17); 0; 0; (-t152 * qJD(4) + t111 * (t228 * qJD(1) - t148 * t180) + t113 * (t230 * qJD(1) - t145 * t180)) * t222 + (-t151 * qJD(4) + t111 * (t73 * qJD(1) + t114 * t134) + t113 * (t71 * qJD(1) + t114 * t133)) * t223 + m(5) * ((-t112 * t23 + t114 * t24) * t95 + t236 * t86) + m(6) * (-t59 * t17 + t18 * t60 + t29 * t40 + t30 * t41) + (t110 / 0.2e1 + t109 / 0.2e1) * t143 * qJD(4) + ((t214 / 0.2e1 + t210 / 0.2e1 - t50 * t224) * t112 + (t213 / 0.2e1 + t209 / 0.2e1 - t51 * t224) * t114) * qJD(1) + t138; m(5) * t11 + m(6) * t3; m(6) * (-t112 * t29 - t114 * t30 + (-t112 * t59 - t114 * t60) * qJD(1)); ((t112 * t75 - t114 * t76) * t11 + t185 * t95 * t86) * t226 + (-t112 * t20 - t114 * t19) * t183 - t114 * ((t114 * t43 + (-t20 + t136) * qJD(1)) * t114 + (t19 * qJD(1) + (-t179 * t230 - t181 * t228) * t112 + (t42 + (t210 + t214) * qJD(4) + (t151 - t69) * qJD(1)) * t114) * t112) - t112 * ((t112 * t42 + (t21 + t241) * qJD(1)) * t112 + (-t22 * qJD(1) + (-t71 * t179 - t73 * t181 + t191) * t114 + (t43 + (t209 + t213) * qJD(4) + t152 * qJD(1)) * t112) * t114) + (t12 * t3 + t29 * t60 - t59 * t30) * t225 - t114 * t2 + t174 + (t112 * t22 + t114 * t21 - t5) * t182; m(6) * (-t150 * t80 + (-t112 * t17 + t114 * t18 + (-t112 * t40 - t114 * t41) * qJD(1)) * t90) + t138; m(6) * t10; 0; m(6) * (-t112 * t30 * t90 + t10 * t12 - t60 * t175 + t59 * t211 + t31 * t3) + (m(6) * (t60 * t80 + (qJD(1) * t59 + t29) * t90) + t172) * t114 + t174; (t185 * t90 * t80 + t31 * t10) * t225 + t172 * t114 + t174;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
