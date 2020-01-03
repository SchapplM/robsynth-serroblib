% Calculate time derivative of joint inertia matrix for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:23
% EndTime: 2019-12-31 16:57:35
% DurationCPUTime: 7.30s
% Computational Cost: add. (2660->349), mult. (4196->497), div. (0->0), fcn. (3212->6), ass. (0->191)
t139 = qJ(2) + pkin(6);
t131 = sin(t139);
t132 = cos(t139);
t242 = Icges(5,5) * t132;
t247 = Icges(4,4) * t132;
t284 = -t242 + t247 + (Icges(4,1) + Icges(5,1)) * t131;
t251 = rSges(5,3) + qJ(4);
t269 = rSges(5,1) + pkin(3);
t278 = t131 * t251 + t132 * t269;
t275 = t132 * t251;
t276 = t131 * t269;
t229 = -t275 + t276;
t283 = t284 * qJD(2);
t144 = sin(qJ(1));
t282 = t144 / 0.2e1;
t146 = cos(qJ(1));
t281 = -t146 / 0.2e1;
t280 = -qJD(1) / 0.2e1;
t279 = qJD(1) / 0.2e1;
t145 = cos(qJ(2));
t127 = pkin(2) * t145 + pkin(1);
t265 = pkin(1) - t127;
t274 = t144 * t265;
t143 = sin(qJ(2));
t249 = Icges(3,4) * t145;
t178 = -Icges(3,2) * t143 + t249;
t80 = Icges(3,6) * t144 + t146 * t178;
t250 = Icges(3,4) * t143;
t184 = Icges(3,1) * t145 - t250;
t82 = Icges(3,5) * t144 + t146 * t184;
t187 = t143 * t80 - t145 * t82;
t161 = t187 * t144;
t79 = -Icges(3,6) * t146 + t144 * t178;
t81 = -Icges(3,5) * t146 + t144 * t184;
t188 = t143 * t79 - t145 * t81;
t162 = t188 * t146;
t176 = -Icges(4,2) * t131 + t247;
t67 = Icges(4,6) * t144 + t146 * t176;
t248 = Icges(4,4) * t131;
t182 = Icges(4,1) * t132 - t248;
t71 = Icges(4,5) * t144 + t146 * t182;
t191 = t131 * t67 - t132 * t71;
t163 = t191 * t144;
t66 = -Icges(4,6) * t146 + t144 * t176;
t70 = -Icges(4,5) * t146 + t144 * t182;
t192 = t131 * t66 - t132 * t70;
t164 = t192 * t146;
t171 = Icges(5,3) * t131 + t242;
t61 = Icges(5,6) * t144 + t146 * t171;
t243 = Icges(5,5) * t131;
t180 = Icges(5,1) * t132 + t243;
t69 = Icges(5,4) * t144 + t146 * t180;
t193 = t131 * t61 + t132 * t69;
t165 = t193 * t144;
t60 = -Icges(5,6) * t146 + t144 * t171;
t68 = -Icges(5,4) * t146 + t144 * t180;
t194 = t131 * t60 + t132 * t68;
t166 = t194 * t146;
t135 = t144 * rSges(4,3);
t232 = t131 * t146;
t273 = -rSges(4,2) * t232 + t135;
t172 = Icges(4,5) * t132 - Icges(4,6) * t131;
t62 = -Icges(4,3) * t146 + t144 * t172;
t173 = Icges(3,5) * t145 - Icges(3,6) * t143;
t77 = -Icges(3,3) * t146 + t144 * t173;
t174 = Icges(5,4) * t132 + Icges(5,6) * t131;
t64 = -Icges(5,2) * t146 + t144 * t174;
t272 = 2 * m(3);
t271 = 2 * m(4);
t270 = 2 * m(5);
t140 = t144 ^ 2;
t141 = t146 ^ 2;
t115 = rSges(3,1) * t143 + rSges(3,2) * t145;
t268 = m(3) * t115;
t267 = pkin(2) * t143;
t266 = t144 * pkin(5);
t138 = t146 * pkin(5);
t142 = -qJ(3) - pkin(5);
t264 = -pkin(5) - t142;
t230 = t142 * t146;
t58 = t138 + t230 - t274;
t121 = t146 * t127;
t59 = -pkin(1) * t146 + t144 * t264 + t121;
t263 = t144 * t58 + t146 * t59;
t262 = -rSges(5,2) * t146 + t144 * t278;
t137 = t144 * rSges(5,2);
t231 = t132 * t146;
t261 = t269 * t231 + t251 * t232 + t137;
t260 = rSges(3,1) * t145;
t259 = rSges(4,1) * t132;
t258 = rSges(3,2) * t143;
t257 = rSges(4,2) * t131;
t256 = rSges(3,3) * t146;
t136 = t144 * rSges(3,3);
t63 = Icges(4,3) * t144 + t146 * t172;
t235 = qJD(1) * t63;
t65 = Icges(5,2) * t144 + t146 * t174;
t234 = qJD(1) * t65;
t78 = Icges(3,3) * t144 + t146 * t173;
t233 = qJD(1) * t78;
t224 = qJD(1) * t146;
t225 = qJD(1) * t144;
t228 = rSges(4,3) * t224 + t225 * t257;
t227 = t146 * t260 + t136;
t226 = t140 + t141;
t223 = qJD(2) * t131;
t222 = qJD(2) * t132;
t221 = qJD(2) * t143;
t220 = qJD(2) * t144;
t219 = qJD(2) * t145;
t218 = qJD(2) * t146;
t217 = qJD(4) * t131;
t133 = qJD(3) * t144;
t213 = pkin(2) * t221;
t210 = qJD(3) * t146 + t142 * t225 + t144 * t213;
t216 = t144 * ((-t146 * t265 - t266) * qJD(1) - t210) + t146 * (-t146 * t213 + t133 + (t146 * t264 + t274) * qJD(1)) + t58 * t224;
t215 = t146 * t258;
t212 = pkin(2) * t219;
t211 = (t172 + t174) * qJD(2) / 0.2e1;
t209 = t143 * t225;
t207 = t269 * qJD(2);
t103 = rSges(4,1) * t131 + rSges(4,2) * t132;
t206 = -t103 - t267;
t205 = -t144 * t142 + t121;
t204 = rSges(5,2) * t224 + t146 * t217 + t218 * t275;
t203 = -t229 - t267;
t202 = -t258 + t260;
t201 = -t257 + t259;
t149 = -t127 - t278;
t147 = t149 * t144;
t21 = (rSges(5,2) - t142) * t146 + t147;
t22 = t205 + t261;
t186 = t144 * t22 + t146 * t21;
t41 = t203 * t144;
t42 = t203 * t146;
t185 = t144 * t41 + t146 * t42;
t175 = Icges(4,2) * t132 + t248;
t169 = -qJD(2) * t278 + qJD(4) * t132 - t212;
t76 = rSges(4,1) * t231 + t273;
t168 = -pkin(1) - t202;
t167 = -t127 - t201;
t160 = qJD(2) * t115;
t159 = qJD(2) * t103;
t154 = qJD(2) * t175;
t153 = qJD(2) * (-Icges(5,4) * t131 + Icges(5,6) * t132);
t152 = qJD(2) * (-Icges(3,5) * t143 - Icges(3,6) * t145);
t151 = qJD(2) * (-Icges(4,5) * t131 - Icges(4,6) * t132);
t148 = rSges(3,2) * t209 + rSges(3,3) * t224 - t146 * t160;
t124 = pkin(2) * t209;
t107 = t202 * qJD(2);
t94 = t201 * qJD(2);
t86 = -t215 + t227;
t85 = t144 * t202 - t256;
t74 = -rSges(4,3) * t146 + t144 * t201;
t57 = t206 * t146;
t56 = t206 * t144;
t55 = t266 + (pkin(1) - t258) * t146 + t227;
t54 = t144 * t168 + t138 + t256;
t50 = t205 + t76;
t49 = (rSges(4,3) - t142) * t146 + t167 * t144;
t44 = t144 * t152 + t233;
t43 = -qJD(1) * t77 + t146 * t152;
t32 = t144 * t153 + t234;
t31 = -qJD(1) * t64 + t146 * t153;
t30 = t144 * t151 + t235;
t29 = -qJD(1) * t62 + t146 * t151;
t26 = t115 * t220 + ((-rSges(3,3) - pkin(5)) * t144 + t168 * t146) * qJD(1);
t25 = (t138 + (-pkin(1) - t260) * t144) * qJD(1) + t148;
t24 = -t103 * t224 - t144 * t94 + (-t143 * t224 - t144 * t219) * pkin(2);
t23 = t103 * t225 + t124 + (-t94 - t212) * t146;
t20 = t144 * t78 - t187 * t146;
t19 = t144 * t77 - t162;
t18 = -t146 * t78 - t161;
t17 = -t144 * t188 - t146 * t77;
t16 = t103 * t220 + (t146 * t167 - t135) * qJD(1) + t210;
t15 = t133 + (-t127 - t259) * t225 + (-qJD(1) * t142 + qJD(2) * t206) * t146 + t228;
t14 = t144 * t63 - t191 * t146;
t13 = t144 * t62 - t164;
t12 = t144 * t65 + t193 * t146;
t11 = t144 * t64 + t166;
t10 = -t146 * t63 - t163;
t9 = -t144 * t192 - t146 * t62;
t8 = -t146 * t65 + t165;
t7 = t144 * t194 - t146 * t64;
t6 = qJD(1) * t42 + t144 * t169;
t5 = t146 * t169 + t225 * t229 + t124;
t4 = (qJD(2) * t229 - t217) * t144 + (t146 * t149 - t137) * qJD(1) + t210;
t3 = t133 + (-t267 - t276) * t218 + (t147 - t230) * qJD(1) + t204;
t2 = t144 * t262 + t146 * t261 + t263;
t1 = (qJD(1) * t262 - t207 * t232 + t204) * t146 + (t220 * t275 + (-t59 + t137 - t261) * qJD(1) + (qJD(4) - t207) * t144 * t131) * t144 + t216;
t27 = [(t25 * t55 + t26 * t54) * t272 + (t15 * t50 + t16 * t49) * t271 + (t21 * t4 + t22 * t3) * t270 + (-Icges(3,2) * t145 + t184 - t250) * t221 + (Icges(3,1) * t143 + t178 + t249) * t219 + (-Icges(5,3) * t132 - t175 + t180 + t182 + t243) * t223 + (t176 - t171 + t284) * t222; m(4) * (t15 * t56 + t16 * t57 + t23 * t49 + t24 * t50) + m(5) * (t21 * t5 + t22 * t6 + t3 * t41 + t4 * t42) + (m(3) * (-t107 * t54 - t115 * t26) + (t154 * t282 + t61 * t279 + t67 * t280) * t132 + (t283 * t282 + (t69 + t71) * t280) * t131 + t211 * t146) * t146 + (m(3) * (-t107 * t55 - t115 * t25) + (t154 * t281 + t60 * t279 + t66 * t280) * t132 + (t283 * t281 + (t68 + t70) * t280) * t131 + t211 * t144) * t144 + ((t140 / 0.2e1 + t141 / 0.2e1) * t173 + t162 / 0.2e1 - t166 / 0.2e1 + t164 / 0.2e1 - t161 / 0.2e1 + t165 / 0.2e1 - t163 / 0.2e1) * qJD(2) + ((-t55 * t268 + (t67 / 0.2e1 - t61 / 0.2e1) * t132 + (t71 / 0.2e1 + t69 / 0.2e1) * t131) * t146 + (t54 * t268 + (t66 / 0.2e1 - t60 / 0.2e1) * t132 + (t70 / 0.2e1 + t68 / 0.2e1) * t131) * t144) * qJD(1); (t1 * t2 + t41 * t6 + t42 * t5) * t270 + (t57 * t23 + t56 * t24 + (t144 * t74 + t146 * t76 + t263) * ((qJD(1) * t74 - t146 * t159 + t228) * t146 + (-t144 * t159 + (-t59 - t76 + t273) * qJD(1)) * t144 + t216)) * t271 - t146 * ((t146 * t30 + (t10 + t164) * qJD(1)) * t146 + (t9 * qJD(1) + (-t222 * t67 - t223 * t71 + t235) * t144 + (-t29 + (t131 * t70 + t132 * t66) * qJD(2) - t191 * qJD(1)) * t146) * t144) - t146 * ((t146 * t32 + (t8 - t166) * qJD(1)) * t146 + (t7 * qJD(1) + (t222 * t61 - t223 * t69 + t234) * t144 + (-t31 + (t131 * t68 - t132 * t60) * qJD(2) + t193 * qJD(1)) * t146) * t144) + t144 * ((t144 * t43 + (t19 + t161) * qJD(1)) * t144 + (t20 * qJD(1) + (t219 * t79 + t221 * t81) * t146 + (-t44 + (-t143 * t82 - t145 * t80) * qJD(2) + (-t188 + t78) * qJD(1)) * t144) * t146) + t144 * ((t144 * t31 + (t11 - t165) * qJD(1)) * t144 + (t12 * qJD(1) + (-t222 * t60 + t223 * t68) * t146 + (-t32 + (-t131 * t69 + t132 * t61) * qJD(2) + (t194 + t65) * qJD(1)) * t144) * t146) - t146 * ((t146 * t44 + (t18 + t162) * qJD(1)) * t146 + (t17 * qJD(1) + (-t219 * t80 - t221 * t82 + t233) * t144 + (-t43 + (t143 * t81 + t145 * t79) * qJD(2) - t187 * qJD(1)) * t146) * t144) + ((t144 * t85 + t146 * t86) * ((qJD(1) * t85 + t148) * t146 + (-t144 * t160 + (-t215 - t86 + t136) * qJD(1)) * t144) + t226 * t115 * t107) * t272 + t144 * ((t144 * t29 + (t13 + t163) * qJD(1)) * t144 + (t14 * qJD(1) + (t222 * t66 + t223 * t70) * t146 + (-t30 + (-t131 * t71 - t132 * t67) * qJD(2) + (-t192 + t63) * qJD(1)) * t144) * t146) + ((-t17 - t7 - t9) * t146 + (t10 + t18 + t8) * t144) * t225 + ((-t11 - t13 - t19) * t146 + (t12 + t14 + t20) * t144) * t224; m(4) * (t144 * t16 - t146 * t15 + (t144 * t50 + t146 * t49) * qJD(1)) + m(5) * (qJD(1) * t186 + t144 * t4 - t146 * t3); m(5) * (qJD(1) * t185 + t144 * t5 - t146 * t6) + m(4) * (t144 * t23 - t146 * t24 + (t144 * t56 + t146 * t57) * qJD(1)); 0; m(5) * (t186 * t222 + (t144 * t3 + t146 * t4 + (-t144 * t21 + t146 * t22) * qJD(1)) * t131); m(5) * ((qJD(2) * t185 - t1) * t132 + (qJD(2) * t2 + t144 * t6 + t146 * t5 + (-t144 * t42 + t146 * t41) * qJD(1)) * t131); 0; (-0.1e1 + t226) * t131 * t222 * t270;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t27(1), t27(2), t27(4), t27(7); t27(2), t27(3), t27(5), t27(8); t27(4), t27(5), t27(6), t27(9); t27(7), t27(8), t27(9), t27(10);];
Mq = res;
