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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:58:54
% EndTime: 2020-01-03 11:59:02
% DurationCPUTime: 4.14s
% Computational Cost: add. (5142->294), mult. (3946->395), div. (0->0), fcn. (2872->8), ass. (0->171)
t165 = cos(qJ(4));
t163 = sin(qJ(4));
t242 = Icges(6,4) * t163;
t244 = Icges(5,4) * t163;
t282 = t242 + t244 + (Icges(5,2) + Icges(6,2)) * t165;
t241 = Icges(6,4) * t165;
t243 = Icges(5,4) * t165;
t281 = t241 + t243 + (Icges(5,1) + Icges(6,1)) * t163;
t192 = -Icges(6,2) * t163 + t241;
t193 = -Icges(5,2) * t163 + t243;
t194 = Icges(6,1) * t165 - t242;
t195 = Icges(5,1) * t165 - t244;
t280 = (t192 + t193) * t165 + (t194 + t195) * t163;
t276 = t163 * t282 - t165 * t281;
t160 = qJD(1) + qJD(2);
t190 = Icges(6,5) * t165 - Icges(6,6) * t163;
t191 = Icges(5,5) * t165 - Icges(5,6) * t163;
t279 = (-t190 - t191) * qJD(4) + (-t276 + t280) * t160;
t126 = Icges(6,5) * t163 + Icges(6,6) * t165;
t127 = Icges(5,5) * t163 + Icges(5,6) * t165;
t278 = -t126 - t127;
t161 = qJ(1) + qJ(2);
t155 = pkin(8) + t161;
t149 = sin(t155);
t150 = cos(t155);
t153 = pkin(4) * t165 + pkin(3);
t228 = t150 * t165;
t229 = t150 * t163;
t162 = -qJ(5) - pkin(7);
t233 = t149 * t162;
t274 = rSges(6,1) * t228 - rSges(6,2) * t229 + rSges(6,3) * t149 + t150 * t153 - t233;
t134 = t150 * t162;
t231 = t149 * t165;
t232 = t149 * t163;
t273 = rSges(6,1) * t231 - rSges(6,2) * t232 - rSges(6,3) * t150 + t149 * t153 + t134;
t70 = -Icges(5,6) * t149 - t150 * t193;
t74 = -Icges(5,5) * t149 - t150 * t195;
t197 = t163 * t70 - t165 * t74;
t272 = t149 * t197;
t68 = -Icges(6,6) * t149 - t150 * t192;
t72 = -Icges(6,5) * t149 - t150 * t194;
t201 = t163 * t68 - t165 * t72;
t271 = t149 * t201;
t156 = sin(t161);
t157 = cos(t161);
t270 = -rSges(3,1) * t156 - rSges(3,2) * t157;
t151 = pkin(2) * t156;
t85 = rSges(4,1) * t149 + rSges(4,2) * t150 + t151;
t269 = -Icges(6,3) * t160 + qJD(4) * t126;
t268 = -Icges(5,3) * t160 + qJD(4) * t127;
t220 = qJD(4) * t165;
t221 = qJD(4) * t163;
t169 = qJD(4) * t280 + t220 * t281 - t221 * t282;
t261 = 2 * m(3);
t260 = 2 * m(4);
t259 = 2 * m(5);
t258 = 2 * m(6);
t248 = rSges(5,2) * t163;
t251 = rSges(5,1) * t165;
t113 = (-t248 + t251) * qJD(4);
t255 = m(5) * t113;
t136 = rSges(5,1) * t163 + rSges(5,2) * t165;
t254 = m(5) * t136;
t144 = t149 * pkin(3);
t253 = pkin(7) * t150 - t144 + t273;
t222 = pkin(3) * t150 + pkin(7) * t149;
t252 = -t222 + t274;
t250 = rSges(6,1) * t165;
t249 = rSges(3,2) * t156;
t247 = rSges(6,2) * t163;
t246 = rSges(6,2) * t165;
t245 = pkin(1) * qJD(1);
t234 = t149 * t160;
t230 = t150 * t160;
t227 = t157 * t160;
t226 = t160 * t163;
t215 = t160 * t228;
t225 = rSges(5,1) * t215 + rSges(5,3) * t234;
t223 = pkin(3) * t230 + pkin(7) * t234;
t164 = sin(qJ(1));
t219 = t164 * t245;
t218 = rSges(6,1) * t215 + rSges(6,3) * t234 + t153 * t230;
t216 = t149 * t226;
t217 = rSges(6,2) * t216 + rSges(6,3) * t230 + qJD(5) * t149;
t135 = rSges(6,1) * t163 + t246;
t212 = pkin(4) * t163 + t135;
t211 = -t153 - t250;
t103 = rSges(3,1) * t157 - t249;
t88 = rSges(3,1) * t227 - t160 * t249;
t209 = rSges(5,1) * t231 - rSges(5,2) * t232;
t152 = pkin(2) * t157;
t86 = rSges(4,1) * t150 - rSges(4,2) * t149 + t152;
t67 = -Icges(6,6) * t150 + t149 * t192;
t71 = -Icges(6,5) * t150 + t149 * t194;
t204 = t163 * t71 + t165 * t67;
t203 = t163 * t67 - t165 * t71;
t202 = t163 * t72 + t165 * t68;
t69 = -Icges(5,6) * t150 + t149 * t193;
t73 = -Icges(5,5) * t150 + t149 * t195;
t200 = t163 * t73 + t165 * t69;
t199 = t163 * t69 - t165 * t73;
t198 = t163 * t74 + t165 * t70;
t133 = pkin(2) * t227;
t62 = rSges(4,1) * t230 - rSges(4,2) * t234 + t133;
t78 = -rSges(5,1) * t228 + rSges(5,2) * t229 - rSges(5,3) * t149;
t87 = t270 * t160;
t185 = t203 * t150;
t184 = t199 * t150;
t183 = qJD(4) * t136;
t178 = t191 * t160;
t177 = t190 * t160;
t176 = -t149 * t220 - t150 * t226;
t173 = (-t246 + (-rSges(6,1) - pkin(4)) * t163) * qJD(4);
t56 = t152 - t78 + t222;
t61 = t85 * t160;
t51 = t151 + t273;
t170 = rSges(5,2) * t216 + rSges(5,3) * t230 - t150 * t183;
t52 = t152 + t274;
t55 = t144 + t151 + (-rSges(5,3) - pkin(7)) * t150 + t209;
t168 = -t160 * t162 + t173;
t167 = -((-t197 - t201) * qJD(4) + t279 * t149) * t149 / 0.2e1 - ((-t199 - t203) * qJD(4) + t279 * t150) * t150 / 0.2e1 + (-t149 * t276 + t150 * t278 + t200 + t204) * t234 / 0.2e1 - (t149 * t278 + t150 * t276 + t198 + t202) * t230 / 0.2e1;
t26 = -rSges(5,1) * t149 * t221 + rSges(5,2) * t176 + t133 + t223 + t225;
t120 = pkin(7) * t230;
t25 = t120 + (-t151 + (-pkin(3) - t251) * t149) * t160 + t170;
t14 = t133 + (-rSges(6,2) * t226 - qJD(5)) * t150 + t168 * t149 + t218;
t13 = (t149 * t211 - t151) * t160 + t168 * t150 + t217;
t166 = cos(qJ(1));
t159 = t166 * pkin(1);
t158 = t164 * pkin(1);
t154 = t166 * t245;
t112 = (-t247 + t250) * qJD(4);
t90 = t103 + t159;
t89 = t158 - t270;
t84 = t154 + t88;
t83 = t87 - t219;
t82 = t159 + t86;
t81 = t158 + t85;
t80 = t212 * t150;
t79 = t212 * t149;
t76 = -rSges(5,3) * t150 + t209;
t66 = -Icges(5,3) * t149 - t150 * t191;
t65 = -Icges(5,3) * t150 + t149 * t191;
t64 = -Icges(6,3) * t149 - t150 * t190;
t63 = -Icges(6,3) * t150 + t149 * t190;
t58 = t154 + t62;
t57 = -t61 - t219;
t54 = t159 + t56;
t53 = t158 + t55;
t50 = t159 + t52;
t49 = t158 + t51;
t40 = -t149 * t268 + t150 * t178;
t39 = t149 * t178 + t150 * t268;
t38 = -t149 * t269 + t150 * t177;
t37 = t149 * t177 + t150 * t269;
t36 = pkin(4) * t176 - t112 * t149 - t135 * t230;
t35 = -t135 * t234 + t112 * t150 + (t150 * t220 - t216) * pkin(4);
t24 = t154 + t26;
t23 = t25 - t219;
t22 = -t149 * t66 + t150 * t197;
t21 = -t149 * t65 + t184;
t20 = -t149 * t64 + t150 * t201;
t19 = -t149 * t63 + t185;
t18 = -t150 * t66 - t272;
t17 = -t149 * t199 - t150 * t65;
t16 = -t150 * t64 - t271;
t15 = -t149 * t203 - t150 * t63;
t12 = t154 + t14;
t11 = t13 - t219;
t2 = (t160 * t76 + t170) * t150 + (-t149 * t183 + (t78 + (-t248 - t251) * t150) * t160 + t225) * t149;
t1 = (-t120 + (-t134 + t253) * t160 + t150 * t173 + t217) * t150 + (-qJD(5) * t150 + t149 * t173 + (-t233 + (pkin(3) + t211 - t247) * t150 - t252) * t160 + t218 - t223) * t149;
t3 = [(t83 * t90 + t84 * t89) * t261 + (t57 * t82 + t58 * t81) * t260 + (t23 * t54 + t24 * t53) * t259 + (t11 * t50 + t12 * t49) * t258 + t169; m(3) * (t103 * t83 - t270 * t84 + t87 * t90 + t88 * t89) + m(4) * (t57 * t86 + t58 * t85 - t61 * t82 + t62 * t81) + m(5) * (t23 * t56 + t24 * t55 + t25 * t54 + t26 * t53) + m(6) * (t11 * t52 + t12 * t51 + t13 * t50 + t14 * t49) + t169; (t13 * t52 + t14 * t51) * t258 + (t25 * t56 + t26 * t55) * t259 + (-t61 * t86 + t62 * t85) * t260 + (t103 * t87 - t270 * t88) * t261 + t169; 0; 0; 0; ((-t160 * t54 + t24) * t150 + (-t160 * t53 - t23) * t149) * t254 + (-t149 * t54 + t150 * t53) * t255 + m(6) * (-t11 * t79 + t12 * t80 + t35 * t49 + t36 * t50) + t167; ((-t160 * t56 + t26) * t150 + (-t160 * t55 - t25) * t149) * t254 + (-t149 * t56 + t150 * t55) * t255 + m(6) * (-t13 * t79 + t14 * t80 + t35 * t51 + t36 * t52) + t167; m(5) * t2 + m(6) * t1; ((t149 * t76 - t150 * t78) * t2 + (t149 ^ 2 + t150 ^ 2) * t136 * t113) * t259 - t150 * ((t150 * t40 + (-t18 + t184) * t160) * t150 + (t17 * t160 + (t220 * t70 + t221 * t74) * t149 + (t200 * qJD(4) + t160 * t197 + t39) * t150) * t149) - t149 * ((t149 * t39 + (t21 + t272) * t160) * t149 + (-t22 * t160 + (-t220 * t69 - t221 * t73) * t150 + (-t198 * qJD(4) + t160 * t199 + t40) * t149) * t150) + ((t149 * t253 + t150 * t252) * t1 - t79 * t36 + t80 * t35) * t258 - t150 * ((t150 * t38 + (-t16 + t185) * t160) * t150 + (t15 * t160 + (t220 * t68 + t221 * t72) * t149 + (t204 * qJD(4) + t160 * t201 + t37) * t150) * t149) - t149 * ((t149 * t37 + (t19 + t271) * t160) * t149 + (-t20 * t160 + (-t220 * t67 - t221 * t71) * t150 + (-t202 * qJD(4) + t160 * t203 + t38) * t149) * t150) + ((-t15 - t17) * t150 + (-t16 - t18) * t149) * t234 + ((t19 + t21) * t150 + (t20 + t22) * t149) * t230; m(6) * ((-t160 * t49 - t11) * t150 + (t160 * t50 - t12) * t149); m(6) * ((-t160 * t51 - t13) * t150 + (t160 * t52 - t14) * t149); 0; m(6) * ((-t160 * t80 - t36) * t150 + (-t160 * t79 - t35) * t149); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
