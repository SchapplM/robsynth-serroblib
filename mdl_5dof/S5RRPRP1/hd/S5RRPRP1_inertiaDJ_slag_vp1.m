% Calculate time derivative of joint inertia matrix for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:34
% EndTime: 2022-01-20 10:19:43
% DurationCPUTime: 3.89s
% Computational Cost: add. (5142->285), mult. (3946->387), div. (0->0), fcn. (2872->8), ass. (0->174)
t159 = cos(qJ(4));
t157 = sin(qJ(4));
t232 = Icges(6,4) * t157;
t234 = Icges(5,4) * t157;
t272 = t232 + t234 + (Icges(5,2) + Icges(6,2)) * t159;
t231 = Icges(6,4) * t159;
t233 = Icges(5,4) * t159;
t271 = t231 + t233 + (Icges(5,1) + Icges(6,1)) * t157;
t186 = -Icges(6,2) * t157 + t231;
t187 = -Icges(5,2) * t157 + t233;
t188 = Icges(6,1) * t159 - t232;
t189 = Icges(5,1) * t159 - t234;
t270 = (t186 + t187) * t159 + (t188 + t189) * t157;
t267 = t272 * t157 - t271 * t159;
t155 = qJ(1) + qJ(2);
t150 = pkin(8) + t155;
t147 = cos(t150);
t154 = qJD(1) + qJD(2);
t227 = t147 * t154;
t266 = -rSges(6,1) - pkin(4);
t126 = Icges(6,5) * t157 + Icges(6,6) * t159;
t127 = Icges(5,5) * t157 + Icges(5,6) * t159;
t265 = t126 + t127;
t151 = sin(t155);
t152 = cos(t155);
t102 = t152 * rSges(3,1) - rSges(3,2) * t151;
t146 = sin(t150);
t217 = qJD(4) * t159;
t222 = t154 * t157;
t262 = -t146 * t217 - t147 * t222;
t184 = Icges(6,5) * t159 - Icges(6,6) * t157;
t185 = Icges(5,5) * t159 - Icges(5,6) * t157;
t261 = t267 * t154 + (t184 + t185) * qJD(4);
t137 = t146 * rSges(6,3);
t149 = pkin(4) * t159 + pkin(3);
t156 = -qJ(5) - pkin(7);
t224 = t147 * t159;
t225 = t147 * t157;
t260 = rSges(6,1) * t224 - rSges(6,2) * t225 - t146 * t156 + t147 * t149 + t137;
t178 = t187 * t147;
t70 = Icges(5,6) * t146 + t178;
t180 = t189 * t147;
t74 = Icges(5,5) * t146 + t180;
t191 = t157 * t70 - t159 * t74;
t259 = t146 * t191;
t177 = t186 * t147;
t68 = Icges(6,6) * t146 + t177;
t179 = t188 * t147;
t72 = Icges(6,5) * t146 + t179;
t195 = t157 * t68 - t159 * t72;
t258 = t146 * t195;
t69 = -Icges(5,6) * t147 + t146 * t187;
t73 = -Icges(5,5) * t147 + t146 * t189;
t193 = t157 * t69 - t159 * t73;
t257 = t147 * t193;
t67 = -Icges(6,6) * t147 + t146 * t186;
t71 = -Icges(6,5) * t147 + t146 * t188;
t197 = t157 * t67 - t159 * t71;
t256 = t147 * t197;
t148 = pkin(2) * t152;
t255 = -t147 * rSges(4,1) - t148;
t226 = t147 * t156;
t229 = t146 * t157;
t254 = rSges(6,2) * t229 + t147 * rSges(6,3) - t226;
t218 = qJD(4) * t157;
t164 = t270 * qJD(4) + t271 * t217 - t272 * t218;
t251 = 2 * m(3);
t250 = 2 * m(4);
t249 = 2 * m(5);
t248 = 2 * m(6);
t239 = rSges(5,1) * t159;
t112 = (-rSges(5,2) * t157 + t239) * qJD(4);
t245 = m(5) * t112;
t134 = rSges(5,1) * t157 + rSges(5,2) * t159;
t244 = m(5) * t134;
t158 = sin(qJ(1));
t243 = pkin(1) * t158;
t242 = pkin(2) * t151;
t143 = t147 * pkin(7);
t228 = t146 * t159;
t241 = t143 + t146 * (-pkin(3) + t149) + rSges(6,1) * t228 - t254;
t219 = -t147 * pkin(3) - t146 * pkin(7);
t240 = t219 + t260;
t238 = rSges(6,1) * t159;
t236 = rSges(6,2) * t159;
t235 = pkin(1) * qJD(1);
t138 = t146 * rSges(5,3);
t230 = t146 * t154;
t223 = t154 * t156;
t220 = rSges(5,2) * t229 + t147 * rSges(5,3);
t216 = t158 * t235;
t160 = cos(qJ(1));
t215 = t160 * t235;
t213 = t146 * t222;
t214 = rSges(6,2) * t213 + rSges(6,3) * t227 + qJD(5) * t146;
t210 = t146 * t218;
t211 = -rSges(5,1) * t210 + t262 * rSges(5,2);
t206 = -pkin(3) - t239;
t133 = rSges(6,1) * t157 + t236;
t205 = -pkin(4) * t157 - t133;
t204 = -t149 - t238;
t88 = t102 * t154;
t86 = -rSges(4,2) * t146 - t255;
t101 = -rSges(3,1) * t151 - rSges(3,2) * t152;
t198 = t157 * t71 + t159 * t67;
t196 = t157 * t72 + t159 * t68;
t194 = t157 * t73 + t159 * t69;
t192 = t157 * t74 + t159 * t70;
t78 = rSges(5,1) * t224 - rSges(5,2) * t225 + t138;
t181 = t262 * rSges(6,2) - qJD(5) * t147 - t146 * t223 + t266 * t210;
t87 = t101 * t154;
t176 = t185 * t147;
t175 = t184 * t147;
t85 = -rSges(4,1) * t146 - rSges(4,2) * t147 - t242;
t174 = t146 * t206 - t242;
t173 = (t266 * t157 - t236) * qJD(4);
t62 = rSges(4,2) * t230 + t255 * t154;
t172 = t146 * t204 - t242;
t56 = t148 + t78 - t219;
t61 = t85 * t154;
t167 = Icges(5,3) * t154 - qJD(4) * t127;
t166 = Icges(6,3) * t154 - qJD(4) * t126;
t165 = -qJD(4) * t134 * t147 + rSges(5,2) * t213 + rSges(5,3) * t227;
t52 = t148 + t260;
t55 = t143 + t174 + t220;
t51 = t172 + t254;
t163 = (t261 * t146 + (-t191 - t195) * qJD(4) - t270 * t230) * t146 / 0.2e1 - (-t261 * t147 + (-t193 - t197) * qJD(4)) * t147 / 0.2e1 + (-t146 * t267 - t265 * t147 + t194 + t198) * t230 / 0.2e1 + (t265 * t146 - t147 * t267 + t192 + t196) * t227 / 0.2e1 - ((t177 + t178) * t159 + (t179 + t180) * t157) * t227 / 0.2e1;
t26 = (-t148 + t206 * t147 + (-rSges(5,3) - pkin(7)) * t146) * t154 - t211;
t14 = (t147 * t204 - t137 - t148) * t154 - t181;
t121 = pkin(7) * t227;
t25 = t154 * t174 + t121 + t165;
t13 = t172 * t154 + (t173 - t223) * t147 + t214;
t153 = t160 * pkin(1);
t111 = (-rSges(6,2) * t157 + t238) * qJD(4);
t90 = t102 + t153;
t89 = t101 - t243;
t84 = -t88 - t215;
t83 = t87 - t216;
t82 = t153 + t86;
t81 = t85 - t243;
t80 = t205 * t147;
t79 = t205 * t146;
t76 = rSges(5,1) * t228 - t220;
t66 = Icges(5,3) * t146 + t176;
t65 = -Icges(5,3) * t147 + t146 * t185;
t64 = Icges(6,3) * t146 + t175;
t63 = -Icges(6,3) * t147 + t146 * t184;
t58 = t62 - t215;
t57 = t61 - t216;
t54 = t153 + t56;
t53 = t55 - t243;
t50 = t153 + t52;
t49 = t51 - t243;
t40 = t146 * t167 + t154 * t176;
t39 = t147 * t167 - t185 * t230;
t38 = t146 * t166 + t154 * t175;
t37 = t147 * t166 - t184 * t230;
t36 = t262 * pkin(4) - t111 * t146 - t133 * t227;
t35 = t133 * t230 - t111 * t147 + (-t147 * t217 + t213) * pkin(4);
t24 = t26 - t215;
t23 = t25 - t216;
t22 = t146 * t66 - t191 * t147;
t21 = t146 * t65 - t257;
t20 = t146 * t64 - t195 * t147;
t19 = t146 * t63 - t256;
t18 = -t147 * t66 - t259;
t17 = -t193 * t146 - t147 * t65;
t16 = -t147 * t64 - t258;
t15 = -t197 * t146 - t147 * t63;
t12 = t14 - t215;
t11 = t13 - t216;
t2 = ((-t78 + t138) * t154 + t211) * t146 + (t154 * t76 + t165) * t147;
t1 = (((rSges(6,3) - pkin(7)) * t146 - t240) * t154 + t181) * t146 + (-t121 + (-t226 + t241) * t154 + t147 * t173 + t214) * t147;
t3 = [(t11 * t50 + t12 * t49) * t248 + (t23 * t54 + t24 * t53) * t249 + (t57 * t82 + t58 * t81) * t250 + (t83 * t90 + t84 * t89) * t251 + t164; m(6) * (t11 * t52 + t12 * t51 + t13 * t50 + t14 * t49) + m(5) * (t23 * t56 + t24 * t55 + t25 * t54 + t26 * t53) + m(4) * (t57 * t86 + t58 * t85 + t61 * t82 + t62 * t81) + m(3) * (t101 * t84 + t102 * t83 + t87 * t90 - t88 * t89) + t164; (t13 * t52 + t14 * t51) * t248 + (t25 * t56 + t26 * t55) * t249 + (t61 * t86 + t62 * t85) * t250 + (-t101 * t88 + t102 * t87) * t251 + t164; 0; 0; 0; ((-t154 * t54 - t24) * t147 + (t154 * t53 - t23) * t146) * t244 + m(6) * (t11 * t79 + t12 * t80 + t35 * t49 + t36 * t50) + (-t146 * t54 - t147 * t53) * t245 + t163; m(6) * (t13 * t79 + t14 * t80 + t35 * t51 + t36 * t52) + ((-t154 * t56 - t26) * t147 + (t154 * t55 - t25) * t146) * t244 + (-t146 * t56 - t147 * t55) * t245 + t163; m(5) * t2 + m(6) * t1; ((t146 * t76 + t147 * t78) * t2 + (t146 ^ 2 + t147 ^ 2) * t134 * t112) * t249 + t146 * ((t146 * t39 + (t21 + t259) * t154) * t146 + (t22 * t154 + (t217 * t69 + t218 * t73) * t147 + (-t192 * qJD(4) - t154 * t193 - t40) * t146) * t147) - t147 * ((t147 * t40 + (t18 + t257) * t154) * t147 + (t17 * t154 + (-t217 * t70 - t218 * t74) * t146 + (t194 * qJD(4) - t154 * t191 - t39) * t147) * t146) + (t80 * t35 + t79 * t36 + (t146 * t241 + t147 * t240) * t1) * t248 + t146 * ((t146 * t37 + (t19 + t258) * t154) * t146 + (t20 * t154 + (t217 * t67 + t218 * t71) * t147 + (-t196 * qJD(4) - t154 * t197 - t38) * t146) * t147) - t147 * ((t147 * t38 + (t16 + t256) * t154) * t147 + (t15 * t154 + (-t217 * t68 - t218 * t72) * t146 + (t198 * qJD(4) - t154 * t195 - t37) * t147) * t146) + ((-t15 - t17) * t147 + (t16 + t18) * t146) * t230 + ((-t19 - t21) * t147 + (t20 + t22) * t146) * t227; m(6) * ((t154 * t49 - t11) * t147 + (t154 * t50 + t12) * t146); m(6) * ((t154 * t51 - t13) * t147 + (t154 * t52 + t14) * t146); 0; m(6) * ((t154 * t80 - t36) * t147 + (t154 * t79 + t35) * t146); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
