% Calculate time derivative of joint inertia matrix for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:37
% EndTime: 2019-12-05 15:12:43
% DurationCPUTime: 2.62s
% Computational Cost: add. (10580->258), mult. (10347->429), div. (0->0), fcn. (9681->8), ass. (0->154)
t142 = qJD(3) + qJD(4);
t143 = sin(pkin(8));
t139 = t143 ^ 2;
t144 = cos(pkin(8));
t140 = t144 ^ 2;
t238 = t139 + t140;
t233 = t142 * t238;
t141 = pkin(9) + qJ(3);
t134 = sin(t141);
t135 = cos(t141);
t239 = t143 * qJD(3) * (-Icges(4,5) * t134 - Icges(4,6) * t135);
t136 = qJ(4) + t141;
t131 = sin(t136);
t132 = cos(t136);
t237 = -Icges(5,5) * t131 - Icges(5,6) * t132;
t227 = t237 * t233;
t147 = cos(qJ(5));
t199 = t144 * t147;
t146 = sin(qJ(5));
t202 = t143 * t146;
t122 = -t132 * t202 - t199;
t200 = t144 * t146;
t201 = t143 * t147;
t123 = t132 * t201 - t200;
t205 = t132 * t142;
t204 = t142 * t143;
t191 = t132 * t204;
t208 = t131 * t142;
t193 = t146 * t208;
t82 = -qJD(5) * t123 + t143 * t193;
t192 = t147 * t208;
t83 = qJD(5) * t122 - t143 * t192;
t46 = Icges(6,5) * t83 + Icges(6,6) * t82 + Icges(6,3) * t191;
t207 = t131 * t143;
t67 = Icges(6,5) * t123 + Icges(6,6) * t122 + Icges(6,3) * t207;
t159 = t131 * t46 + t67 * t205;
t48 = Icges(6,4) * t83 + Icges(6,2) * t82 + Icges(6,6) * t191;
t50 = Icges(6,1) * t83 + Icges(6,4) * t82 + Icges(6,5) * t191;
t69 = Icges(6,4) * t123 + Icges(6,2) * t122 + Icges(6,6) * t207;
t71 = Icges(6,1) * t123 + Icges(6,4) * t122 + Icges(6,5) * t207;
t13 = t122 * t48 + t123 * t50 + t143 * t159 + t69 * t82 + t71 * t83;
t203 = t142 * t144;
t190 = t132 * t203;
t125 = t132 * t199 + t202;
t84 = -qJD(5) * t125 + t144 * t193;
t124 = -t132 * t200 + t201;
t85 = qJD(5) * t124 - t144 * t192;
t47 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t190;
t206 = t131 * t144;
t68 = Icges(6,5) * t125 + Icges(6,6) * t124 + Icges(6,3) * t206;
t158 = t131 * t47 + t68 * t205;
t49 = Icges(6,4) * t85 + Icges(6,2) * t84 + Icges(6,6) * t190;
t51 = Icges(6,1) * t85 + Icges(6,4) * t84 + Icges(6,5) * t190;
t70 = Icges(6,4) * t125 + Icges(6,2) * t124 + Icges(6,6) * t206;
t72 = Icges(6,1) * t125 + Icges(6,4) * t124 + Icges(6,5) * t206;
t14 = t122 * t49 + t123 * t51 + t143 * t158 + t70 * t82 + t72 * t83;
t9 = -t13 * t144 + t14 * t143;
t152 = t142 * t237;
t97 = t143 * t152;
t98 = t144 * t152;
t236 = -t140 * t97 - (-t144 * t98 + t227) * t143 - t9;
t235 = t144 * t239;
t230 = qJD(3) * (rSges(4,1) * t134 + rSges(4,2) * t135);
t226 = 2 * m(5);
t225 = 2 * m(6);
t15 = t124 * t48 + t125 * t50 + t144 * t159 + t69 * t84 + t71 * t85;
t16 = t124 * t49 + t125 * t51 + t144 * t158 + t70 * t84 + t72 * t85;
t10 = t143 * t16 - t144 * t15;
t224 = (t139 * t98 + t10 + (-t143 * t97 + t227) * t144) * t143;
t223 = pkin(3) * t134;
t222 = t238 * pkin(3) * t135;
t127 = rSges(5,1) * t131 + rSges(5,2) * t132;
t57 = t127 * t233;
t181 = rSges(5,1) * t132 - rSges(5,2) * t131;
t64 = t238 * t181;
t221 = pkin(3) * qJD(3);
t166 = Icges(6,5) * t147 - Icges(6,6) * t146;
t220 = t132 * (t166 * t205 + (Icges(6,3) * t142 + (-Icges(6,5) * t146 - Icges(6,6) * t147) * qJD(5)) * t131);
t94 = -Icges(6,3) * t132 + t166 * t131;
t219 = t132 * t94;
t30 = t124 * t69 + t125 * t71 + t206 * t67;
t218 = t143 * t30;
t29 = t122 * t70 + t123 * t72 + t207 * t68;
t217 = t144 * t29;
t183 = pkin(4) * t132 + pkin(7) * t131;
t180 = rSges(6,1) * t147 - rSges(6,2) * t146;
t61 = t180 * t205 + (rSges(6,3) * t142 + (-rSges(6,1) * t146 - rSges(6,2) * t147) * qJD(5)) * t131;
t216 = -t183 * t142 - t61;
t211 = Icges(6,4) * t146;
t210 = Icges(6,4) * t147;
t103 = -t132 * rSges(6,3) + t180 * t131;
t128 = pkin(4) * t131 - pkin(7) * t132;
t198 = -t103 - t128;
t194 = t135 * t221;
t188 = -t127 - t223;
t54 = rSges(6,1) * t83 + rSges(6,2) * t82 + rSges(6,3) * t191;
t55 = rSges(6,1) * t85 + rSges(6,2) * t84 + rSges(6,3) * t190;
t26 = -t128 * t233 + t143 * t54 + t144 * t55;
t73 = rSges(6,1) * t123 + rSges(6,2) * t122 + rSges(6,3) * t207;
t74 = rSges(6,1) * t125 + rSges(6,2) * t124 + rSges(6,3) * t206;
t34 = t143 * t73 + t144 * t74 + t183 * t238;
t185 = t198 - t223;
t114 = t181 * t142;
t184 = -t114 - t194;
t177 = -t146 * t69 + t147 * t71;
t32 = t131 * t177 - t132 * t67;
t176 = -t146 * t70 + t147 * t72;
t33 = t131 * t176 - t132 * t68;
t179 = t32 * t143 + t33 * t144;
t178 = -t143 * t74 + t144 * t73;
t167 = -Icges(6,2) * t146 + t210;
t95 = -Icges(6,6) * t132 + t167 * t131;
t170 = Icges(6,1) * t147 - t211;
t96 = -Icges(6,5) * t132 + t170 * t131;
t175 = t146 * t95 - t147 * t96;
t162 = t236 * t144 + t224;
t161 = -t194 + t216;
t160 = t238 * t134 * t221;
t11 = (t142 * t177 - t46) * t132 + (t142 * t67 - t146 * t48 + t147 * t50 + (-t146 * t71 - t147 * t69) * qJD(5)) * t131;
t12 = (t142 * t176 - t47) * t132 + (t142 * t68 - t146 * t49 + t147 * t51 + (-t146 * t72 - t147 * t70) * qJD(5)) * t131;
t28 = t122 * t69 + t123 * t71 + t207 * t67;
t36 = t122 * t95 + t123 * t96 + t207 * t94;
t59 = t167 * t205 + (Icges(6,6) * t142 + (-Icges(6,2) * t147 - t211) * qJD(5)) * t131;
t60 = t170 * t205 + (Icges(6,5) * t142 + (-Icges(6,1) * t146 - t210) * qJD(5)) * t131;
t3 = (-t122 * t59 - t123 * t60 - t82 * t95 - t83 * t96 + (t217 + (t28 - t219) * t143) * t142) * t132 + (t14 * t144 + t36 * t142 + (t13 - t220) * t143) * t131;
t31 = t124 * t70 + t125 * t72 + t206 * t68;
t37 = t124 * t95 + t125 * t96 + t206 * t94;
t4 = (-t124 * t59 - t125 * t60 - t84 * t95 - t85 * t96 + (t218 + (t31 - t219) * t144) * t142) * t132 + (t37 * t142 + t15 * t143 + (t16 - t220) * t144) * t131;
t156 = t143 * t4 / 0.2e1 - t144 * t3 / 0.2e1 + (t143 * t33 - t144 * t32) * t208 / 0.2e1 - t132 * (-t11 * t144 + t12 * t143) / 0.2e1 + t9 * t207 / 0.2e1 + t10 * t206 / 0.2e1 + (t143 * (t143 * t29 - t144 * t28) + t144 * (t143 * t31 - t144 * t30)) * t205 / 0.2e1;
t93 = t188 * t144;
t92 = t188 * t143;
t81 = t184 * t144;
t80 = t184 * t143;
t77 = t238 * t230;
t76 = t198 * t144;
t75 = t198 * t143;
t63 = t185 * t144;
t62 = t185 * t143;
t56 = -t160 - t57;
t53 = t216 * t144;
t52 = t216 * t143;
t43 = -t103 * t206 - t132 * t74;
t42 = t103 * t207 + t132 * t73;
t41 = t161 * t144;
t40 = t161 * t143;
t39 = -t131 * t175 - t219;
t38 = t178 * t131;
t35 = t64 + t222;
t27 = t34 + t222;
t25 = (-t103 * t203 - t55) * t132 + (t142 * t74 - t144 * t61) * t131;
t24 = (t103 * t204 + t54) * t132 + (-t142 * t73 + t143 * t61) * t131;
t23 = -t160 + t26;
t22 = t178 * t205 + (-t143 * t55 + t144 * t54) * t131;
t1 = [0; 0; 0; -m(4) * t77 + m(5) * t56 + m(6) * t23; m(5) * (t143 * t81 - t144 * t80) + m(6) * (t143 * t41 - t144 * t40); (t23 * t27 + t40 * t62 + t41 * t63) * t225 + t139 * t235 + (t35 * t56 + t80 * t92 + t81 * t93) * t226 + t224 + 0.2e1 * m(4) * (t230 - t77) * t238 * (rSges(4,1) * t135 - rSges(4,2) * t134) + (t144 * t235 - t238 * t239 + t236) * t144; -m(5) * t57 + m(6) * t26; m(6) * (t143 * t53 - t144 * t52); m(6) * (t34 * t23 + t26 * t27 + t40 * t75 + t41 * t76 + t52 * t62 + t53 * t63) + m(5) * (-t57 * t35 + t64 * t56 + (-t143 * t80 - t144 * t81) * t127 + (-t143 * t92 - t144 * t93) * t114) + t162; (t114 * t127 * t238 - t64 * t57) * t226 + (t34 * t26 + t52 * t75 + t53 * t76) * t225 + t162; m(6) * t22; m(6) * (t143 * t24 - t144 * t25); m(6) * (t22 * t27 + t23 * t38 + t24 * t63 + t25 * t62 + t40 * t43 + t41 * t42) + t156; m(6) * (t22 * t34 + t24 * t76 + t25 * t75 + t38 * t26 + t42 * t53 + t43 * t52) + t156; (t22 * t38 + t24 * t42 + t25 * t43) * t225 + (-t132 * t37 + (t144 * t31 + t218) * t131) * t190 + t4 * t206 + (-t132 * t36 + (t143 * t28 + t217) * t131) * t191 + t3 * t207 + (t131 * t179 - t39 * t132) * t208 - t132 * ((t220 + (t132 * t175 + t179) * t142) * t132 + (t12 * t144 + t11 * t143 - (t142 * t94 - t146 * t59 + t147 * t60 + (-t146 * t96 - t147 * t95) * qJD(5)) * t132 + t39 * t142) * t131);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
