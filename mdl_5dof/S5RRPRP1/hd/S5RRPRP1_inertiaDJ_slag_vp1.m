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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:22:00
% EndTime: 2019-12-05 18:22:09
% DurationCPUTime: 3.81s
% Computational Cost: add. (5142->298), mult. (3946->399), div. (0->0), fcn. (2872->8), ass. (0->181)
t161 = sin(qJ(4));
t163 = cos(qJ(4));
t243 = Icges(6,4) * t161;
t191 = Icges(6,1) * t163 - t243;
t245 = Icges(5,4) * t161;
t192 = Icges(5,1) * t163 - t245;
t293 = (t191 + t192) * t161;
t292 = t243 + t245 + (Icges(5,2) + Icges(6,2)) * t163;
t242 = Icges(6,4) * t163;
t244 = Icges(5,4) * t163;
t291 = t242 + t244 + (Icges(5,1) + Icges(6,1)) * t161;
t158 = qJD(1) + qJD(2);
t289 = t158 * t163;
t285 = t292 * t161 - t291 * t163;
t187 = Icges(6,5) * t163 - Icges(6,6) * t161;
t188 = Icges(5,5) * t163 - Icges(5,6) * t161;
t190 = -Icges(5,2) * t161 + t244;
t288 = -t190 * t289 + (t187 + t188) * qJD(4) + (t285 - t293) * t158;
t132 = Icges(6,5) * t161 + Icges(6,6) * t163;
t133 = Icges(5,5) * t161 + Icges(5,6) * t163;
t287 = t132 + t133;
t159 = qJ(1) + qJ(2);
t155 = pkin(8) + t159;
t152 = cos(t155);
t222 = qJD(4) * t161;
t213 = t152 * t222;
t151 = sin(t155);
t233 = t151 * t163;
t284 = -t158 * t233 - t213;
t221 = qJD(4) * t163;
t225 = t158 * t161;
t283 = t151 * t221 + t152 * t225;
t69 = Icges(5,6) * t152 - t151 * t190;
t73 = Icges(5,5) * t152 - t151 * t192;
t196 = t161 * t69 - t163 * t73;
t281 = t152 * t196;
t189 = -Icges(6,2) * t161 + t242;
t175 = t189 * t151;
t67 = Icges(6,6) * t152 - t175;
t71 = Icges(6,5) * t152 - t151 * t191;
t200 = t161 * t67 - t163 * t71;
t280 = t152 * t200;
t160 = -qJ(5) - pkin(7);
t230 = t152 * t161;
t279 = rSges(6,2) * t230 + t151 * t160;
t234 = t151 * t161;
t278 = rSges(6,2) * t234 + (rSges(6,3) - t160) * t152;
t277 = -Icges(6,3) * t158 + qJD(4) * t132;
t276 = -Icges(5,3) * t158 + qJD(4) * t133;
t167 = -t292 * t222 + t291 * t221 + ((t189 + t190) * t163 + t293) * qJD(4);
t269 = 2 * m(3);
t268 = 2 * m(4);
t267 = 2 * m(5);
t266 = 2 * m(6);
t263 = -rSges(5,3) - pkin(7);
t262 = -pkin(7) + rSges(6,3);
t250 = rSges(5,2) * t161;
t252 = rSges(5,1) * t163;
t113 = (-t250 + t252) * qJD(4);
t261 = m(5) * t113;
t144 = rSges(5,1) * t161 + rSges(5,2) * t163;
t260 = m(5) * t144;
t162 = sin(qJ(1));
t259 = pkin(1) * t162;
t164 = cos(qJ(1));
t258 = pkin(1) * t164;
t156 = sin(t159);
t257 = pkin(2) * t156;
t157 = cos(t159);
t256 = pkin(2) * t157;
t153 = pkin(4) * t163 + pkin(3);
t255 = pkin(3) - t153;
t149 = t152 * pkin(7);
t254 = rSges(6,1) * t233 - t151 * t255 + t149 - t278;
t229 = t152 * t163;
t246 = t151 * rSges(6,3);
t253 = rSges(6,1) * t229 - t151 * pkin(7) - t152 * t255 + t246 - t279;
t251 = rSges(6,1) * t163;
t249 = rSges(6,2) * t161;
t248 = pkin(1) * qJD(1);
t247 = t151 * rSges(5,3);
t148 = t152 * rSges(5,3);
t235 = t151 * t158;
t232 = t152 * t158;
t228 = t156 * t158;
t227 = t157 * t158;
t226 = t158 * t160;
t129 = rSges(5,2) * t234;
t223 = t129 + t148;
t87 = rSges(3,1) * t228 + rSges(3,2) * t227;
t145 = qJD(5) * t152;
t220 = t164 * t248;
t212 = t152 * t221;
t219 = t284 * rSges(5,1) - rSges(5,2) * t212;
t215 = t151 * t222;
t216 = rSges(5,1) * t215 + t283 * rSges(5,2);
t141 = pkin(2) * t228;
t61 = rSges(4,1) * t235 + rSges(4,2) * t232 + t141;
t209 = -pkin(3) - t252;
t143 = rSges(6,1) * t161 + rSges(6,2) * t163;
t208 = pkin(4) * t161 + t143;
t103 = -t157 * rSges(3,1) + t156 * rSges(3,2);
t207 = -t153 - t251;
t88 = -rSges(3,1) * t227 + rSges(3,2) * t228;
t206 = -t152 * rSges(4,1) - t256;
t102 = -rSges(3,1) * t156 - rSges(3,2) * t157;
t201 = t161 * t71 + t163 * t67;
t68 = Icges(6,6) * t151 + t152 * t189;
t72 = Icges(6,5) * t151 + t152 * t191;
t199 = t161 * t72 + t163 * t68;
t198 = t161 * t68 - t163 * t72;
t197 = t161 * t73 + t163 * t69;
t70 = Icges(5,6) * t151 + t152 * t190;
t74 = Icges(5,5) * t151 + t152 * t192;
t195 = t161 * t74 + t163 * t70;
t194 = t161 * t70 - t163 * t74;
t182 = t284 * rSges(6,1) - rSges(6,2) * t212 - pkin(4) * t213 - t152 * t226 - t153 * t235;
t86 = t151 * rSges(4,2) + t206;
t181 = t151 * t226 + t145 + (rSges(6,1) + pkin(4)) * t215 + t283 * rSges(6,2);
t180 = t198 * t151;
t179 = t194 * t151;
t174 = t188 * t158;
t173 = t187 * t158;
t85 = -rSges(4,1) * t151 - rSges(4,2) * t152 - t257;
t62 = rSges(4,2) * t235 + t158 * t206;
t168 = t152 * t207 - t246 - t256;
t55 = t151 * t209 + t149 + t223 - t257;
t166 = t151 * t263 + t152 * t209 - t256;
t131 = rSges(5,2) * t230;
t56 = t131 + t166;
t52 = t168 + t279;
t51 = t151 * t207 - t257 + t278;
t165 = (-t175 * t289 + (-t194 - t198) * qJD(4) + t288 * t151) * t151 / 0.2e1 + (-t189 * t232 * t163 + (-t196 - t200) * qJD(4) + t288 * t152) * t152 / 0.2e1 - (t285 * t151 + t287 * t152 + t197 + t201) * t235 / 0.2e1 + (t287 * t151 - t285 * t152 + t195 + t199) * t232 / 0.2e1;
t127 = pkin(3) * t235;
t25 = t127 + t141 + (t152 * t263 - t129) * t158 - t219;
t26 = t158 * t166 + t216;
t13 = -rSges(6,3) * t232 + t141 + (-rSges(6,2) * t225 - qJD(5)) * t151 - t182;
t14 = t158 * t168 + t181;
t154 = t162 * t248;
t112 = (-t249 + t251) * qJD(4);
t90 = t103 - t258;
t89 = t102 - t259;
t84 = t88 - t220;
t83 = t154 + t87;
t82 = t86 - t258;
t81 = t85 - t259;
t80 = t208 * t152;
t79 = t208 * t151;
t78 = rSges(5,1) * t229 - t131 + t247;
t76 = -rSges(5,1) * t233 + t223;
t66 = Icges(5,3) * t151 + t152 * t188;
t65 = Icges(5,3) * t152 - t151 * t188;
t64 = Icges(6,3) * t151 + t152 * t187;
t63 = Icges(6,3) * t152 - t151 * t187;
t58 = t62 - t220;
t57 = t154 + t61;
t54 = t56 - t258;
t53 = t55 - t259;
t50 = t52 - t258;
t49 = t51 - t259;
t40 = t276 * t151 - t152 * t174;
t39 = -t151 * t174 - t276 * t152;
t38 = t277 * t151 - t152 * t173;
t37 = -t151 * t173 - t277 * t152;
t36 = t283 * pkin(4) + t112 * t151 + t143 * t232;
t35 = t143 * t235 - t112 * t152 + (t151 * t225 - t212) * pkin(4);
t24 = t26 - t220;
t23 = t154 + t25;
t22 = t151 * t66 - t194 * t152;
t21 = t151 * t65 - t281;
t20 = t151 * t64 - t198 * t152;
t19 = t151 * t63 - t280;
t18 = t152 * t66 + t179;
t17 = t196 * t151 + t152 * t65;
t16 = t152 * t64 + t180;
t15 = t200 * t151 + t152 * t63;
t12 = t14 - t220;
t11 = t154 + t13;
t2 = t152 * t219 - t151 * t216 + ((-t76 + t148) * t152 + (t247 - t78 + (t250 + t252) * t152) * t151) * t158;
t1 = (t127 + t182) * t152 + (-t181 + t145) * t151 + ((t262 * t151 - t253) * t151 + (t262 * t152 + (-pkin(3) - t207 + t249) * t151 + t254) * t152) * t158;
t3 = [(t83 * t90 + t84 * t89) * t269 + (t57 * t82 + t58 * t81) * t268 + (t23 * t54 + t24 * t53) * t267 + (t11 * t50 + t12 * t49) * t266 + t167; m(3) * (t102 * t84 + t103 * t83 + t87 * t90 + t88 * t89) + m(4) * (t57 * t86 + t58 * t85 + t61 * t82 + t62 * t81) + m(5) * (t23 * t56 + t24 * t55 + t25 * t54 + t26 * t53) + m(6) * (t11 * t52 + t12 * t51 + t13 * t50 + t14 * t49) + t167; (t13 * t52 + t14 * t51) * t266 + (t25 * t56 + t26 * t55) * t267 + (t61 * t86 + t62 * t85) * t268 + (t102 * t88 + t103 * t87) * t269 + t167; 0; 0; 0; ((t158 * t54 - t24) * t152 + (t158 * t53 + t23) * t151) * t260 + (t151 * t54 - t152 * t53) * t261 + t165 + m(6) * (t11 * t79 - t12 * t80 + t35 * t49 + t36 * t50); ((t158 * t56 - t26) * t152 + (t158 * t55 + t25) * t151) * t260 + m(6) * (t13 * t79 - t14 * t80 + t35 * t51 + t36 * t52) + (t151 * t56 - t152 * t55) * t261 + t165; m(5) * t2 + m(6) * t1; ((-t151 * t76 + t152 * t78) * t2 + (t151 ^ 2 + t152 ^ 2) * t144 * t113) * t267 + t152 * ((t152 * t40 + (t18 + t281) * t158) * t152 + (-t17 * t158 + (t221 * t70 + t222 * t74) * t151 + (t197 * qJD(4) + t158 * t194 + t39) * t152) * t151) + t151 * ((t151 * t39 + (-t21 + t179) * t158) * t151 + (t22 * t158 + (-t221 * t69 - t222 * t73) * t152 + (-t195 * qJD(4) + t158 * t196 + t40) * t151) * t152) + ((t151 * t254 + t152 * t253) * t1 + t79 * t36 - t80 * t35) * t266 + t152 * ((t152 * t38 + (t16 + t280) * t158) * t152 + (-t15 * t158 + (t221 * t68 + t222 * t72) * t151 + (t201 * qJD(4) + t158 * t198 + t37) * t152) * t151) + t151 * ((t151 * t37 + (-t19 + t180) * t158) * t151 + (t20 * t158 + (-t221 * t67 - t222 * t71) * t152 + (-t199 * qJD(4) + t158 * t200 + t38) * t151) * t152) + ((-t15 - t17) * t152 + (-t16 - t18) * t151) * t235 + ((t19 + t21) * t152 + (t20 + t22) * t151) * t232; m(6) * ((t158 * t49 + t11) * t152 + (-t158 * t50 + t12) * t151); m(6) * ((t158 * t51 + t13) * t152 + (-t158 * t52 + t14) * t151); 0; m(6) * ((-t158 * t80 + t36) * t152 + (-t158 * t79 + t35) * t151); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
