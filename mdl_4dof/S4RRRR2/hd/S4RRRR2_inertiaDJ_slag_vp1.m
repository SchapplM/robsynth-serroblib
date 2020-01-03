% Calculate time derivative of joint inertia matrix for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR2_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR2_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:08
% EndTime: 2019-12-31 17:23:13
% DurationCPUTime: 2.79s
% Computational Cost: add. (6497->330), mult. (5482->478), div. (0->0), fcn. (4148->8), ass. (0->201)
t156 = qJ(1) + qJ(2);
t149 = sin(t156);
t151 = cos(t156);
t159 = cos(qJ(3));
t220 = qJD(3) * t159;
t154 = qJD(1) + qJD(2);
t157 = sin(qJ(3));
t230 = t154 * t157;
t266 = -t149 * t220 - t151 * t230;
t247 = rSges(4,2) * t157;
t249 = rSges(4,1) * t159;
t265 = -t247 + t249;
t241 = Icges(4,4) * t159;
t184 = -Icges(4,2) * t157 + t241;
t175 = t184 * t151;
t87 = Icges(4,6) * t149 + t175;
t242 = Icges(4,4) * t157;
t186 = Icges(4,1) * t159 - t242;
t177 = t186 * t151;
t89 = Icges(4,5) * t149 + t177;
t188 = t157 * t87 - t159 * t89;
t264 = t188 * t149;
t86 = -Icges(4,6) * t151 + t184 * t149;
t88 = -Icges(4,5) * t151 + t186 * t149;
t190 = t157 * t86 - t159 * t88;
t263 = t190 * t151;
t155 = qJ(3) + qJ(4);
t148 = sin(t155);
t150 = cos(t155);
t239 = Icges(5,4) * t150;
t183 = -Icges(5,2) * t148 + t239;
t174 = t183 * t151;
t73 = Icges(5,6) * t149 + t174;
t240 = Icges(5,4) * t148;
t185 = Icges(5,1) * t150 - t240;
t176 = t185 * t151;
t75 = Icges(5,5) * t149 + t176;
t194 = t148 * t73 - t150 * t75;
t262 = t194 * t149;
t72 = -Icges(5,6) * t151 + t183 * t149;
t74 = -Icges(5,5) * t151 + t185 * t149;
t195 = t148 * t72 - t150 * t74;
t261 = t195 * t151;
t137 = t149 * rSges(5,3);
t248 = rSges(5,1) * t150;
t260 = t151 * t248 + t137;
t153 = qJD(3) + qJD(4);
t109 = Icges(5,2) * t150 + t240;
t110 = Icges(5,1) * t148 + t239;
t180 = t109 * t148 - t110 * t150;
t181 = Icges(5,5) * t150 - Icges(5,6) * t148;
t259 = t181 * t153 + t180 * t154;
t128 = Icges(4,2) * t159 + t242;
t129 = Icges(4,1) * t157 + t241;
t179 = t128 * t157 - t129 * t159;
t182 = Icges(4,5) * t159 - Icges(4,6) * t157;
t258 = t182 * qJD(3) + t179 * t154;
t257 = 2 * m(3);
t256 = 2 * m(4);
t255 = 2 * m(5);
t254 = t149 / 0.2e1;
t253 = -t151 / 0.2e1;
t119 = t265 * qJD(3);
t252 = m(4) * t119;
t134 = rSges(4,1) * t157 + rSges(4,2) * t159;
t251 = m(4) * t134;
t142 = t149 * pkin(6);
t158 = sin(qJ(1));
t250 = t158 * pkin(1);
t246 = rSges(5,2) * t148;
t122 = t149 * t246;
t225 = t151 * rSges(5,3) + t122;
t76 = t149 * t248 - t225;
t215 = t151 * t246;
t77 = -t215 + t260;
t46 = t149 * t76 + t151 * t77;
t245 = pkin(1) * qJD(1);
t138 = t149 * rSges(4,3);
t238 = t109 * t153;
t237 = t110 * t153;
t236 = t148 * t153;
t235 = t149 * t154;
t234 = t149 * t159;
t233 = t150 * t153;
t232 = t151 * t154;
t161 = -pkin(7) - pkin(6);
t231 = t151 * t161;
t229 = t154 * t161;
t228 = rSges(5,3) * t232 + t154 * t122;
t214 = t149 * t230;
t227 = rSges(4,2) * t214 + rSges(4,3) * t232;
t221 = qJD(3) * t157;
t210 = t149 * t221;
t226 = -pkin(3) * t210 - t149 * t229;
t224 = t151 * rSges(4,3) + t149 * t247;
t223 = -t151 * pkin(2) - t142;
t222 = t149 ^ 2 + t151 ^ 2;
t218 = rSges(5,2) * t233;
t212 = -t154 * t215 + (-t236 * rSges(5,1) - t218) * t149;
t219 = t149 * (t260 * t154 + t212) + t151 * (-t151 * t218 + (-t150 * t235 - t151 * t236) * rSges(5,1) + t228) + t76 * t232;
t160 = cos(qJ(1));
t217 = t160 * t245;
t216 = t158 * t245;
t211 = -rSges(4,1) * t210 + t266 * rSges(4,2);
t208 = t151 * t221;
t207 = t151 * t220;
t206 = t235 / 0.2e1;
t205 = t232 / 0.2e1;
t204 = -pkin(2) - t249;
t111 = rSges(5,1) * t148 + rSges(5,2) * t150;
t203 = -pkin(3) * t157 - t111;
t170 = Icges(5,6) * t154 - t238;
t42 = t170 * t151 - t183 * t235;
t202 = t153 * t75 + t42;
t43 = t170 * t149 + t154 * t174;
t201 = t153 * t74 + t43;
t171 = Icges(5,5) * t154 - t237;
t44 = t171 * t151 - t185 * t235;
t200 = -t153 * t73 + t44;
t45 = t171 * t149 + t154 * t176;
t199 = -t153 * t72 + t45;
t145 = pkin(3) * t159 + pkin(2);
t198 = -t145 - t248;
t113 = t151 * rSges(3,1) - rSges(3,2) * t149;
t197 = t151 * t145 - t149 * t161;
t98 = -rSges(3,1) * t232 + rSges(3,2) * t235;
t70 = -Icges(5,3) * t151 + t181 * t149;
t16 = -t195 * t149 - t151 * t70;
t172 = t181 * t151;
t71 = Icges(5,3) * t149 + t172;
t17 = -t151 * t71 - t262;
t18 = t149 * t70 - t261;
t19 = t149 * t71 - t194 * t151;
t108 = Icges(5,5) * t148 + Icges(5,6) * t150;
t169 = Icges(5,3) * t154 - t108 * t153;
t40 = t169 * t151 - t181 * t235;
t41 = t169 * t149 + t154 * t172;
t196 = -t151 * ((t151 * t41 + (t17 + t261) * t154) * t151 + (t16 * t154 + (-t148 * t42 + t150 * t44 - t73 * t233 - t75 * t236) * t149 + (-t40 + (t154 * t75 - t199) * t150 + (-t154 * t73 + t201) * t148) * t151) * t149) + t149 * ((t149 * t40 + (t18 + t262) * t154) * t149 + (t19 * t154 + (t148 * t43 - t150 * t45 + t72 * t233 + t74 * t236) * t151 + (-t41 + (t154 * t74 + t200) * t150 + (-t154 * t72 - t202) * t148) * t149) * t151) + (t149 * t17 - t151 * t16) * t235 + (t149 * t19 - t151 * t18) * t232;
t112 = -rSges(3,1) * t149 - rSges(3,2) * t151;
t191 = t157 * t88 + t159 * t86;
t189 = t157 * t89 + t159 * t87;
t187 = t198 * t149;
t127 = Icges(4,5) * t157 + Icges(4,6) * t159;
t91 = t265 * t151 + t138;
t94 = t183 * t153;
t95 = t185 * t153;
t163 = t108 * t154 + (t95 - t238) * t150 + (-t94 - t237) * t148;
t178 = (t200 * t148 + t259 * t149 + t202 * t150 + t163 * t151) * t254 + (t199 * t148 + t163 * t149 + t201 * t150 - t259 * t151) * t253 + (-t108 * t151 + t148 * t74 - t180 * t149 + t150 * t72) * t206 + (t108 * t149 + t148 * t75 + t150 * t73 - t180 * t151) * t205;
t97 = t112 * t154;
t173 = t182 * t151;
t65 = t91 - t223;
t143 = t151 * pkin(6);
t64 = t204 * t149 + t143 + t224;
t60 = t197 + t77;
t168 = Icges(4,5) * t154 - t129 * qJD(3);
t167 = Icges(4,6) * t154 - t128 * qJD(3);
t166 = Icges(4,3) * t154 - t127 * qJD(3);
t59 = t187 + t225 - t231;
t115 = t184 * qJD(3);
t116 = t186 * qJD(3);
t165 = -t109 * t236 + t110 * t233 + t159 * t115 + t157 * t116 - t128 * t221 + t129 * t220 + t148 * t95 + t150 * t94;
t162 = -t115 * t157 + t116 * t159 + t127 * t154 + (-t128 * t159 - t129 * t157) * qJD(3);
t164 = t178 + (-t188 * qJD(3) + t258 * t149 + t162 * t151 + t157 * (t168 * t151 - t186 * t235) + t159 * (t167 * t151 - t184 * t235)) * t254 + (-t190 * qJD(3) + t162 * t149 - t258 * t151 + t157 * (t168 * t149 + t154 * t177) + t159 * (t167 * t149 + t154 * t175)) * t253 + (-t127 * t151 - t179 * t149 + t191) * t206 + (t127 * t149 - t179 * t151 + t189) * t205;
t35 = (t204 * t151 + (-rSges(4,3) - pkin(6)) * t149) * t154 - t211;
t23 = (t198 * t151 - t137) * t154 - t212 - t226;
t133 = pkin(6) * t232;
t34 = -rSges(4,2) * t207 - pkin(2) * t235 + t133 + (-t154 * t234 - t208) * rSges(4,1) + t227;
t22 = t154 * t187 + (-pkin(3) * t221 - t111 * t153 - t229) * t151 + t228;
t152 = t160 * pkin(1);
t100 = t113 + t152;
t99 = t112 - t250;
t96 = (-t246 + t248) * t153;
t90 = rSges(4,1) * t234 - t224;
t85 = Icges(4,3) * t149 + t173;
t84 = -Icges(4,3) * t151 + t182 * t149;
t81 = t98 - t217;
t80 = t97 - t216;
t79 = t203 * t151;
t78 = t203 * t149;
t69 = t197 + t223;
t68 = t231 + t143 + t149 * (-pkin(2) + t145);
t62 = t152 + t65;
t61 = t64 - t250;
t58 = t152 + t60;
t57 = t59 - t250;
t52 = t166 * t149 + t154 * t173;
t51 = t166 * t151 - t182 * t235;
t39 = t266 * pkin(3) - t111 * t232 - t149 * t96;
t38 = t111 * t235 - t151 * t96 + (-t207 + t214) * pkin(3);
t29 = t35 - t217;
t28 = t34 - t216;
t27 = t149 * t85 - t188 * t151;
t26 = t149 * t84 - t263;
t25 = -t151 * t85 - t264;
t24 = -t190 * t149 - t151 * t84;
t21 = t23 - t217;
t20 = t22 - t216;
t15 = t149 * t68 + t151 * t69 + t46;
t10 = -t77 * t235 + t219;
t3 = t149 * t226 + t151 * (-pkin(3) * t208 - t133) + ((t68 - t231) * t151 + (-t69 - t77 - t142) * t149) * t154 + t219;
t1 = [(t100 * t80 + t81 * t99) * t257 + (t28 * t62 + t29 * t61) * t256 + (t20 * t58 + t21 * t57) * t255 + t165; m(3) * (t100 * t97 + t112 * t81 + t113 * t80 + t98 * t99) + m(4) * (t28 * t65 + t29 * t64 + t34 * t62 + t35 * t61) + m(5) * (t20 * t60 + t59 * t21 + t22 * t58 + t23 * t57) + t165; (t22 * t60 + t59 * t23) * t255 + (t34 * t65 + t35 * t64) * t256 + (t112 * t98 + t113 * t97) * t257 + t165; t164 + (-t149 * t62 - t151 * t61) * t252 + ((-t154 * t62 - t29) * t151 + (t154 * t61 - t28) * t149) * t251 + m(5) * (t20 * t78 + t21 * t79 + t38 * t57 + t39 * t58); t164 + m(5) * (t22 * t78 + t23 * t79 + t38 * t59 + t39 * t60) + ((-t154 * t65 - t35) * t151 + (t154 * t64 - t34) * t149) * t251 + (-t149 * t65 - t151 * t64) * t252; ((t149 * t90 + t151 * t91) * (((-t91 + t138) * t154 + t211) * t149 + (-t134 * t151 * qJD(3) + t154 * t90 + t227) * t151) + t222 * t134 * t119) * t256 + (t27 * t149 - t151 * t26) * t232 + t149 * ((t149 * t51 + (t26 + t264) * t154) * t149 + (t27 * t154 + (t86 * t220 + t88 * t221) * t151 + (-t189 * qJD(3) - t190 * t154 - t52) * t149) * t151) + (t149 * t25 - t24 * t151) * t235 - t151 * ((t151 * t52 + (t25 + t263) * t154) * t151 + (t24 * t154 + (-t87 * t220 - t89 * t221) * t149 + (t191 * qJD(3) - t188 * t154 - t51) * t151) * t149) + (t15 * t3 + t38 * t79 + t39 * t78) * t255 + t196; m(5) * ((-t149 * t58 - t151 * t57) * t96 + ((-t154 * t58 - t21) * t151 + (t154 * t57 - t20) * t149) * t111) + t178; m(5) * ((-t149 * t60 - t151 * t59) * t96 + ((-t154 * t60 - t23) * t151 + (t154 * t59 - t22) * t149) * t111) + t178; m(5) * (t10 * t15 + t46 * t3 + (-t149 * t78 - t151 * t79) * t96 + ((-t154 * t78 - t38) * t151 + (t154 * t79 - t39) * t149) * t111) + t196; (t222 * t96 * t111 + t46 * t10) * t255 + t196;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
