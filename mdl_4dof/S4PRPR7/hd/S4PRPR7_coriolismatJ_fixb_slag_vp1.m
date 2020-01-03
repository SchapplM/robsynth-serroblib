% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPR7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR7_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:38
% EndTime: 2019-12-31 16:25:43
% DurationCPUTime: 4.16s
% Computational Cost: add. (4669->334), mult. (12430->537), div. (0->0), fcn. (13268->6), ass. (0->202)
t175 = sin(qJ(2));
t290 = t175 ^ 2;
t177 = cos(qJ(2));
t228 = t175 * t177;
t287 = Icges(4,2) - Icges(4,3);
t288 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t228 + (0.2e1 * t177 ^ 2 - 0.2e1 * t290) * Icges(3,4);
t172 = sin(pkin(6));
t196 = -Icges(3,5) * t175 - Icges(3,6) * t177;
t144 = t196 * t172;
t173 = cos(pkin(6));
t145 = t196 * t173;
t198 = Icges(4,4) * t175 + Icges(4,5) * t177;
t286 = t198 * t172 + t144;
t285 = t198 * t173 + t145;
t170 = t172 ^ 2;
t171 = t173 ^ 2;
t218 = t170 + t171;
t234 = t172 * t177;
t239 = Icges(4,6) * t175;
t282 = -0.2e1 * Icges(4,6) * t177 - t287 * t175;
t284 = (Icges(4,5) * t173 - t282 * t172) * t177 + (Icges(4,4) * t173 - 0.2e1 * t172 * t239 + t287 * t234) * t175 + t145 + t288 * t172;
t233 = t173 * t177;
t283 = t144 - t288 * t173 + (Icges(4,5) * t172 + t282 * t173) * t177 + (Icges(4,4) * t172 + 0.2e1 * t173 * t239 - t287 * t233) * t175;
t174 = sin(qJ(4));
t176 = cos(qJ(4));
t230 = (-Icges(5,5) * t176 + Icges(5,6) * t174) * t228;
t215 = m(4) / 0.4e1 + m(5) / 0.4e1;
t219 = t218 * t228;
t279 = t215 * (t219 - t228);
t278 = t218 * t175;
t242 = Icges(5,4) * t174;
t197 = Icges(5,2) * t176 + t242;
t181 = -Icges(5,6) * t175 + t197 * t177;
t225 = -t181 - (-Icges(5,1) * t176 + t242) * t177;
t241 = Icges(5,4) * t176;
t201 = Icges(5,1) * t174 + t241;
t182 = -Icges(5,5) * t175 + t201 * t177;
t224 = -t182 + (Icges(5,2) * t174 - t241) * t177;
t195 = Icges(5,5) * t174 + Icges(5,6) * t176;
t180 = -Icges(5,3) * t175 + t195 * t177;
t208 = rSges(5,1) * t174 + rSges(5,2) * t176;
t183 = -t175 * rSges(5,3) + t208 * t177;
t274 = 2 * qJD(2);
t273 = m(4) / 0.2e1;
t272 = m(5) / 0.2e1;
t157 = t175 * pkin(2) - qJ(3) * t177;
t207 = t175 * rSges(4,2) + rSges(4,3) * t177;
t221 = -t157 + t207;
t107 = t221 * t172;
t109 = t221 * t173;
t253 = t107 * t234 + t109 * t233;
t161 = -rSges(4,2) * t177 + t175 * rSges(4,3);
t160 = pkin(2) * t177 + t175 * qJ(3);
t222 = t218 * t160;
t58 = t218 * t161 + t222;
t271 = m(4) * (t278 * t58 + t253);
t229 = t175 * t176;
t138 = t172 * t229 + t173 * t174;
t232 = t174 * t175;
t139 = -t172 * t232 + t173 * t176;
t86 = -t139 * rSges(5,1) + t138 * rSges(5,2) + rSges(5,3) * t234;
t247 = t175 * t86;
t66 = -t183 * t234 - t247;
t136 = -t172 * t174 + t173 * t229;
t137 = t172 * t176 + t173 * t232;
t85 = t137 * rSges(5,1) + t136 * rSges(5,2) + rSges(5,3) * t233;
t248 = t175 * t85;
t67 = t183 * t233 + t248;
t259 = t66 * t233 + t67 * t234;
t105 = t183 * t172;
t106 = t183 * t173;
t45 = (t105 * t177 - t247) * t173 + (-t106 * t177 + t248) * t172;
t126 = rSges(5,3) * t177 + t208 * t175;
t191 = t126 * t177 + t175 * t183;
t51 = -t175 * t105 + t191 * t172 - t177 * t86;
t52 = t175 * t106 - t191 * t173 + t177 * t85;
t55 = (-t172 * t85 + t173 * t86) * t177;
t270 = m(5) * (-t177 * t45 + (t172 * t52 + t173 * t51 + t55) * t175 + t259);
t269 = m(5) * (t45 * t55 + t51 * t66 + t52 * t67);
t210 = -pkin(5) * t175 - t157 + t183;
t81 = t210 * t172;
t83 = t210 * t173;
t258 = t83 * t233 + t81 * t234;
t260 = pkin(5) * t177;
t47 = t172 * t86 + t173 * t85 + t218 * t260 + t222;
t268 = m(5) * (t278 * t47 + t258);
t267 = m(5) * (t278 * t55 + t259);
t155 = (-rSges(5,1) * t176 + rSges(5,2) * t174) * t177;
t94 = rSges(5,1) * t136 - rSges(5,2) * t137;
t95 = rSges(5,1) * t138 + rSges(5,2) * t139;
t61 = t172 * t95 + t173 * t94;
t266 = m(5) * (-t155 * t278 - t177 * t61);
t265 = -t172 / 0.2e1;
t264 = t172 / 0.2e1;
t263 = -t173 / 0.2e1;
t262 = t173 / 0.2e1;
t261 = t175 / 0.2e1;
t244 = Icges(5,4) * t137;
t77 = Icges(5,2) * t136 + Icges(5,6) * t233 + t244;
t257 = -Icges(5,1) * t136 + t244 + t77;
t243 = Icges(5,4) * t139;
t78 = Icges(5,2) * t138 + Icges(5,6) * t234 - t243;
t256 = -Icges(5,1) * t138 - t243 + t78;
t134 = Icges(5,4) * t136;
t79 = Icges(5,1) * t137 + Icges(5,5) * t233 + t134;
t255 = -Icges(5,2) * t137 + t134 + t79;
t135 = Icges(5,4) * t138;
t80 = -Icges(5,1) * t139 + Icges(5,5) * t234 + t135;
t254 = Icges(5,2) * t139 + t135 + t80;
t252 = m(5) * qJD(4);
t75 = Icges(5,5) * t137 + Icges(5,6) * t136 + Icges(5,3) * t233;
t250 = t175 * t75;
t76 = -Icges(5,5) * t139 + Icges(5,6) * t138 + Icges(5,3) * t234;
t249 = t175 * t76;
t231 = t175 * t180;
t34 = 0.2e1 * (t45 / 0.4e1 - t61 / 0.4e1) * m(5);
t227 = t34 * qJD(1);
t209 = 0.2e1 * t215 * t278;
t216 = t273 + t272;
t74 = -t216 * t175 + t209;
t226 = t74 * qJD(1);
t223 = t218 * t157;
t220 = -t160 - t161;
t217 = qJD(4) * t177;
t88 = Icges(5,5) * t136 - Icges(5,6) * t137;
t30 = t255 * t136 - t257 * t137 + t88 * t233;
t89 = Icges(5,5) * t138 + Icges(5,6) * t139;
t31 = t254 * t136 - t256 * t137 + t89 * t233;
t16 = t172 * t30 - t173 * t31;
t122 = Icges(5,6) * t177 + t197 * t175;
t124 = Icges(5,5) * t177 + t201 * t175;
t192 = t174 * t182 + t176 * t181;
t188 = (Icges(5,3) * t177 + t195 * t175 - t192) * t175;
t102 = t181 * t173;
t104 = t182 * t173;
t205 = -t174 * t79 - t176 * t77;
t189 = t180 * t173 - t205;
t178 = t189 * t177 - t250;
t25 = t136 * t102 + t137 * t104 + t178 * t173;
t101 = t181 * t172;
t103 = t182 * t172;
t204 = -t174 * t80 - t176 * t78;
t190 = t180 * t172 - t204;
t179 = t190 * t177 - t249;
t26 = t136 * t101 + t137 * t103 + t179 * t173;
t41 = t136 * t77 + t137 * t79 + t75 * t233;
t42 = t136 * t78 + t137 * t80 + t76 * t233;
t56 = -t136 * t181 - t137 * t182 - t180 * t233;
t4 = (t136 * t122 + t137 * t124) * t175 + t56 * t177 + (-t175 * t42 + t177 * t26) * t172 + ((-t41 + t231) * t175 + (t25 + t188) * t177) * t173;
t214 = t16 / 0.2e1 - t4 / 0.2e1;
t32 = t255 * t138 + t257 * t139 + t88 * t234;
t33 = t254 * t138 + t256 * t139 + t89 * t234;
t17 = t172 * t32 - t173 * t33;
t23 = t138 * t102 - t139 * t104 + t178 * t172;
t24 = t138 * t101 - t139 * t103 + t179 * t172;
t43 = t138 * t77 - t139 * t79 + t75 * t234;
t44 = t138 * t78 - t139 * t80 + t76 * t234;
t57 = -t138 * t181 + t139 * t182 - t180 * t234;
t3 = (t138 * t122 - t139 * t124) * t175 + t57 * t177 + (-t175 * t43 + t177 * t23) * t173 + ((-t44 + t231) * t175 + (t24 + t188) * t177) * t172;
t213 = t17 / 0.2e1 - t3 / 0.2e1;
t211 = -t126 - t160 - t260;
t159 = t175 * rSges(3,1) + rSges(3,2) * t177;
t48 = t205 * t177 + t250;
t49 = t204 * t177 + t249;
t206 = t49 * t172 + t48 * t173;
t110 = t220 * t173;
t108 = t220 * t172;
t87 = t218 * t159;
t84 = t211 * t173;
t82 = t211 * t172;
t73 = t209 + (m(4) + m(5)) * t261;
t70 = 0.4e1 * t279;
t69 = -t155 * t233 + t175 * t94;
t68 = t155 * t234 - t175 * t95;
t65 = t192 * t177 - t231;
t64 = t218 * t207 - t223;
t59 = (-t172 * t94 + t173 * t95) * t177;
t54 = -pkin(5) * t278 + t105 * t172 + t106 * t173 - t223;
t50 = t266 / 0.2e1;
t39 = t175 * t89 + (t256 * t174 - t254 * t176) * t177;
t38 = t175 * t88 + (t257 * t174 - t255 * t176) * t177;
t37 = (-t101 * t176 - t103 * t174 + t76) * t177 + t190 * t175;
t36 = (-t102 * t176 - t104 * t174 + t75) * t177 + t189 * t175;
t35 = (t45 + t61) * t272;
t28 = t267 / 0.2e1;
t22 = t47 * t61 + (-t172 * t81 - t173 * t83) * t155;
t21 = t65 * t175 + t206 * t177;
t20 = t268 + t271;
t19 = t57 * t175 + (t172 * t44 + t173 * t43) * t177;
t18 = t56 * t175 + (t172 * t42 + t173 * t41) * t177;
t15 = t172 * t25 - t173 * t26;
t14 = t172 * t23 - t173 * t24;
t11 = t270 / 0.2e1;
t10 = (t224 * t138 + t225 * t139) * t175 + (t32 * t173 + (t33 + t230) * t172) * t177;
t9 = (t224 * t136 - t225 * t137) * t175 + (t31 * t172 + (t30 + t230) * t173) * t177;
t8 = (t37 * t172 + t36 * t173 + t65) * t177 + (t188 + (-t122 * t176 - t124 * t174 - t180) * t177 - t206) * t175;
t7 = t28 + t11 - t266 / 0.2e1;
t6 = t50 + t28 - t270 / 0.2e1;
t5 = t50 + t11 - t267 / 0.2e1;
t2 = m(5) * t22 + t16 * t264 + t17 * t263;
t1 = t269 + (t4 * t262 + t3 * t264 + t21 / 0.2e1) * t177 + (t18 * t263 + t19 * t265 + t8 / 0.2e1) * t175;
t12 = [0, t73 * qJD(3) + t35 * qJD(4) + (-m(3) * t87 / 0.2e1 + t64 * t273 + t54 * t272) * t274, t73 * qJD(2), t35 * qJD(2) + t59 * t252; qJD(3) * t74 - qJD(4) * t34, t20 * qJD(3) + t2 * qJD(4) + (m(4) * (t107 * t108 + t109 * t110 + t58 * t64) + m(5) * (t47 * t54 + t81 * t82 + t83 * t84) + m(3) * (t159 - t87) * t218 * (rSges(3,1) * t177 - t175 * rSges(3,2)) + (t15 + t285 * t170 + (t284 * t173 + (t283 - t286) * t172) * t173) * t264 + (t14 + t286 * t171 + (t283 * t172 + (t284 - t285) * t173) * t172) * t263) * qJD(2), t20 * qJD(2) + t6 * qJD(4) + t226 + (-0.4e1 * t279 + 0.2e1 * t216 * (-t177 * t278 + t219)) * qJD(3), -t227 + t2 * qJD(2) + t6 * qJD(3) + (-t21 / 0.2e1 + t214 * t173 + t213 * t172) * t217 + (t9 * t264 + t10 * t263 + m(5) * (t59 * t47 + t55 * t61 + t68 * t83 + t69 * t81 + (-t172 * t67 - t173 * t66) * t155) - t269 + (-t8 / 0.2e1 + (-t39 / 0.2e1 + t18 / 0.2e1) * t173 + (t38 / 0.2e1 + t19 / 0.2e1) * t172) * t175) * qJD(4); -t74 * qJD(2), -t226 + t70 * qJD(3) + t5 * qJD(4) + 0.4e1 * (-t271 / 0.4e1 - t268 / 0.4e1) * qJD(2) + ((-t177 * t64 + t253) * t273 + (-t177 * t54 + t258) * t272 + ((t108 * t172 + t110 * t173 + t58) * t273 + (t172 * t82 + t173 * t84 + t47) * t272) * t175) * t274, t70 * qJD(2), t5 * qJD(2) + (-t177 * t59 + (t172 * t69 + t173 * t68) * t175) * t252; t34 * qJD(2), t227 + (((-t49 / 0.2e1 + t15 / 0.2e1) * t177 + t213) * t173 + ((t48 / 0.2e1 + t14 / 0.2e1) * t177 - t214) * t172 + (t36 * t264 + (t172 * t43 - t173 * t44) * t265 + (t172 * t41 - t173 * t42 + t37) * t263) * t175 + (t45 * t47 + t51 * t83 + t52 * t81 + t54 * t55 + t66 * t84 + t67 * t82 - t22) * m(5)) * qJD(2) + t7 * qJD(3) + t1 * qJD(4), t7 * qJD(2), t1 * qJD(2) + (m(5) * (t55 * t59 + t66 * t68 + t67 * t69) + t290 * t230 / 0.2e1) * qJD(4) + (t9 * t262 + t10 * t264 + (t38 * t173 + t39 * t172 + (t225 * t174 - t224 * t176) * t175) * t261) * t217;];
Cq = t12;
