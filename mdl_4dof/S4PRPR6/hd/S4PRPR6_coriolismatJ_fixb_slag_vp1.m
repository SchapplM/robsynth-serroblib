% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:21
% EndTime: 2019-12-31 16:24:27
% DurationCPUTime: 4.99s
% Computational Cost: add. (8865->371), mult. (13030->590), div. (0->0), fcn. (13928->8), ass. (0->214)
t188 = cos(qJ(2));
t298 = t188 ^ 2;
t187 = sin(qJ(2));
t243 = t187 * t188;
t296 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t243 + (-0.2e1 * t187 ^ 2 + 0.2e1 * t298) * Icges(3,4);
t183 = sin(pkin(6));
t211 = -Icges(3,5) * t187 - Icges(3,6) * t188;
t160 = t211 * t183;
t185 = cos(pkin(6));
t161 = t211 * t185;
t179 = t183 ^ 2;
t180 = t185 ^ 2;
t231 = t179 + t180;
t181 = pkin(7) + qJ(4);
t177 = sin(t181);
t178 = cos(t181);
t250 = (-Icges(5,5) * t177 - Icges(5,6) * t178) * t243;
t228 = m(4) / 0.4e1 + m(5) / 0.4e1;
t232 = t231 * t243;
t293 = t228 * (t232 - t243);
t166 = t231 * t187;
t255 = Icges(5,4) * t178;
t212 = -Icges(5,2) * t177 + t255;
t128 = -Icges(5,6) * t188 + t212 * t187;
t238 = -t128 + (-Icges(5,1) * t177 - t255) * t187;
t256 = Icges(5,4) * t177;
t215 = Icges(5,1) * t178 - t256;
t130 = -Icges(5,5) * t188 + t215 * t187;
t237 = t130 + (-Icges(5,2) * t178 - t256) * t187;
t210 = Icges(5,5) * t178 - Icges(5,6) * t177;
t126 = -Icges(5,3) * t188 + t210 * t187;
t221 = rSges(5,1) * t178 - rSges(5,2) * t177;
t132 = -rSges(5,3) * t188 + t221 * t187;
t182 = sin(pkin(7));
t184 = cos(pkin(7));
t222 = rSges(4,1) * t184 - rSges(4,2) * t182;
t291 = rSges(4,3) * t188 - t222 * t187;
t247 = t183 * t188;
t156 = -t182 * t247 - t185 * t184;
t157 = -t185 * t182 + t184 * t247;
t245 = t185 * t188;
t158 = -t182 * t245 + t183 * t184;
t159 = t183 * t182 + t184 * t245;
t246 = t185 * t187;
t248 = t183 * t187;
t290 = (t183 * (Icges(4,5) * t159 + Icges(4,6) * t158 + Icges(4,3) * t246) - t185 * (Icges(4,5) * t157 + Icges(4,6) * t156 + Icges(4,3) * t248)) * t188 + (-(t182 * (Icges(4,4) * t157 + Icges(4,2) * t156 + Icges(4,6) * t248) - t184 * (Icges(4,1) * t157 + Icges(4,4) * t156 + Icges(4,5) * t248)) * t185 + (t182 * (Icges(4,4) * t159 + Icges(4,2) * t158 + Icges(4,6) * t246) - t184 * (Icges(4,1) * t159 + Icges(4,4) * t158 + Icges(4,5) * t246)) * t183) * t187;
t287 = 2 * qJD(2);
t286 = m(4) / 0.2e1;
t285 = m(5) / 0.2e1;
t167 = t187 * pkin(2) - qJ(3) * t188;
t236 = -t167 + t291;
t104 = t236 * t183;
t106 = t236 * t185;
t239 = t104 * t247 + t106 * t245;
t169 = pkin(2) * t188 + t187 * qJ(3);
t233 = t231 * t169;
t55 = t183 * (rSges(4,1) * t157 + rSges(4,2) * t156 + rSges(4,3) * t248) + t185 * (rSges(4,1) * t159 + rSges(4,2) * t158 + rSges(4,3) * t246) + t233;
t284 = m(4) * (t166 * t55 + t239);
t140 = -t177 * t247 - t185 * t178;
t141 = -t185 * t177 + t178 * t247;
t83 = rSges(5,1) * t141 + rSges(5,2) * t140 + rSges(5,3) * t248;
t262 = t188 * t83;
t66 = t132 * t248 + t262;
t142 = -t177 * t245 + t183 * t178;
t143 = t183 * t177 + t178 * t245;
t84 = rSges(5,1) * t143 + rSges(5,2) * t142 + rSges(5,3) * t246;
t261 = t188 * t84;
t67 = -t132 * t246 - t261;
t273 = t66 * t245 + t67 * t247;
t115 = t132 * t183;
t116 = t132 * t185;
t47 = (-t115 * t187 + t262) * t185 + (t116 * t187 - t261) * t183;
t133 = t187 * rSges(5,3) + t221 * t188;
t208 = t132 * t188 + t133 * t187;
t53 = -t115 * t188 + t208 * t183 - t187 * t83;
t54 = t116 * t188 - t208 * t185 + t187 * t84;
t58 = (-t183 * t84 + t185 * t83) * t187;
t283 = m(5) * (-t188 * t47 + (t183 * t54 + t185 * t53 + t58) * t187 + t273);
t282 = m(5) * (t47 * t58 + t53 * t66 + t54 * t67);
t186 = -pkin(5) - qJ(3);
t240 = qJ(3) + t186;
t174 = pkin(3) * t184 + pkin(2);
t274 = -pkin(2) + t174;
t225 = -t274 * t187 - t240 * t188 - t132 - t167;
t72 = t225 * t183;
t74 = t225 * t185;
t272 = t74 * t245 + t72 * t247;
t124 = -t240 * t187 + t274 * t188;
t40 = (t124 * t185 + t84) * t185 + (t124 * t183 + t83) * t183 + t233;
t281 = m(5) * (t166 * t40 + t272);
t280 = m(5) * (t166 * t58 + t273);
t149 = (-rSges(5,1) * t177 - rSges(5,2) * t178) * t187;
t100 = rSges(5,1) * t142 - rSges(5,2) * t143;
t99 = rSges(5,1) * t140 - rSges(5,2) * t141;
t62 = t100 * t185 + t183 * t99;
t279 = m(5) * (-t149 * t166 - t188 * t62);
t278 = t183 / 0.2e1;
t277 = -t185 / 0.2e1;
t276 = t185 / 0.2e1;
t258 = Icges(5,4) * t141;
t79 = Icges(5,2) * t140 + Icges(5,6) * t248 + t258;
t271 = Icges(5,1) * t140 - t258 - t79;
t257 = Icges(5,4) * t143;
t80 = Icges(5,2) * t142 + Icges(5,6) * t246 + t257;
t270 = Icges(5,1) * t142 - t257 - t80;
t134 = Icges(5,4) * t140;
t81 = Icges(5,1) * t141 + Icges(5,5) * t248 + t134;
t269 = -Icges(5,2) * t141 + t134 + t81;
t135 = Icges(5,4) * t142;
t82 = Icges(5,1) * t143 + Icges(5,5) * t246 + t135;
t268 = -Icges(5,2) * t143 + t135 + t82;
t267 = m(5) * qJD(4);
t77 = Icges(5,5) * t141 + Icges(5,6) * t140 + Icges(5,3) * t248;
t264 = t188 * t77;
t78 = Icges(5,5) * t143 + Icges(5,6) * t142 + Icges(5,3) * t246;
t263 = t188 * t78;
t251 = t126 * t188;
t244 = t187 * t126;
t34 = 0.2e1 * (t47 / 0.4e1 - t62 / 0.4e1) * m(5);
t242 = t34 * qJD(1);
t223 = 0.2e1 * t228 * t166;
t229 = t286 + t285;
t86 = -t229 * t187 + t223;
t241 = t86 * qJD(1);
t235 = -t187 * rSges(4,3) - t222 * t188 - t169;
t234 = t231 * t167;
t230 = qJD(4) * t187;
t87 = Icges(5,5) * t140 - Icges(5,6) * t141;
t28 = t269 * t140 + t271 * t141 + t87 * t248;
t88 = Icges(5,5) * t142 - Icges(5,6) * t143;
t29 = t268 * t140 + t270 * t141 + t88 * t248;
t16 = t183 * t29 - t185 * t28;
t129 = Icges(5,6) * t187 + t212 * t188;
t131 = Icges(5,5) * t187 + t215 * t188;
t209 = -t128 * t177 + t130 * t178;
t202 = (Icges(5,3) * t187 + t210 * t188 - t209) * t188;
t111 = t128 * t183;
t113 = t130 * t183;
t220 = -t177 * t79 + t178 * t81;
t207 = -t126 * t183 - t220;
t190 = t207 * t187 + t264;
t24 = -t140 * t111 - t141 * t113 + t190 * t183;
t112 = t128 * t185;
t114 = t130 * t185;
t219 = -t177 * t80 + t178 * t82;
t206 = -t126 * t185 - t219;
t189 = t206 * t187 + t263;
t25 = -t140 * t112 - t141 * t114 + t189 * t183;
t42 = t140 * t79 + t141 * t81 + t77 * t248;
t43 = t140 * t80 + t141 * t82 + t78 * t248;
t56 = t128 * t140 + t130 * t141 + t183 * t244;
t3 = -(t140 * t129 + t141 * t131) * t188 + t56 * t187 + (t187 * t25 + t188 * t43) * t185 + ((t42 - t251) * t188 + (t24 - t202) * t187) * t183;
t227 = t16 / 0.2e1 - t3 / 0.2e1;
t30 = t269 * t142 + t271 * t143 + t87 * t246;
t31 = t268 * t142 + t270 * t143 + t88 * t246;
t17 = t183 * t31 - t185 * t30;
t26 = -t142 * t111 - t143 * t113 + t190 * t185;
t27 = -t142 * t112 - t143 * t114 + t189 * t185;
t44 = t142 * t79 + t143 * t81 + t77 * t246;
t45 = t142 * t80 + t143 * t82 + t78 * t246;
t57 = t128 * t142 + t130 * t143 + t185 * t244;
t4 = -(t142 * t129 + t143 * t131) * t188 + t57 * t187 + (t187 * t26 + t188 * t44) * t183 + ((t45 - t251) * t188 + (t27 - t202) * t187) * t185;
t226 = t4 / 0.2e1 - t17 / 0.2e1;
t224 = -t124 - t133 - t169;
t168 = t187 * rSges(3,1) + rSges(3,2) * t188;
t49 = t220 * t187 - t264;
t50 = t219 * t187 - t263;
t218 = t49 * t183 + t50 * t185;
t201 = t296 * t183 + t161;
t200 = -t296 * t185 + t160;
t197 = Icges(4,5) * t188 + (-Icges(4,1) * t184 + Icges(4,4) * t182) * t187;
t195 = Icges(4,6) * t188 + (-Icges(4,4) * t184 + Icges(4,2) * t182) * t187;
t191 = -t187 * t174 - t186 * t188 + t167;
t122 = t197 * t185;
t121 = t197 * t183;
t120 = t195 * t185;
t119 = t195 * t183;
t107 = t235 * t185;
t105 = t235 * t183;
t103 = t231 * t168;
t85 = t223 + (m(4) + m(5)) * t187 / 0.2e1;
t76 = 0.4e1 * t293;
t75 = t224 * t185;
t73 = t224 * t183;
t71 = -t100 * t188 - t149 * t246;
t70 = t149 * t248 + t188 * t99;
t63 = t209 * t187 - t251;
t60 = (-t100 * t183 + t185 * t99) * t187;
t59 = t231 * t291 - t234;
t51 = t279 / 0.2e1;
t48 = (t191 * t185 - t116) * t185 + (t191 * t183 - t115) * t183 - t234;
t39 = -t188 * t88 + (-t268 * t177 + t270 * t178) * t187;
t38 = -t188 * t87 + (-t269 * t177 + t271 * t178) * t187;
t37 = -t206 * t188 + (t112 * t177 - t114 * t178 + t78) * t187;
t36 = -t207 * t188 + (t111 * t177 - t113 * t178 + t77) * t187;
t35 = (t47 + t62) * t285;
t32 = t280 / 0.2e1;
t22 = t218 * t187 - t188 * t63;
t21 = t40 * t62 + (-t183 * t72 - t185 * t74) * t149;
t20 = -t188 * t57 + (t183 * t44 + t185 * t45) * t187;
t19 = -t188 * t56 + (t183 * t42 + t185 * t43) * t187;
t18 = t281 + t284;
t15 = t183 * t27 - t185 * t26;
t14 = t183 * t25 - t185 * t24;
t11 = t283 / 0.2e1;
t10 = -(t237 * t142 + t238 * t143) * t188 + (t30 * t183 + (t31 - t250) * t185) * t187;
t9 = -(t237 * t140 + t238 * t141) * t188 + (t29 * t185 + (t28 - t250) * t183) * t187;
t8 = (t202 + t218) * t188 + (t37 * t185 + t36 * t183 - (-t129 * t177 + t131 * t178 + t126) * t188 + t63) * t187;
t7 = t32 + t11 - t279 / 0.2e1;
t6 = t51 + t32 - t283 / 0.2e1;
t5 = t51 + t11 - t280 / 0.2e1;
t2 = m(5) * t21 + t16 * t277 + t17 * t278;
t1 = t282 + (t20 * t276 + t19 * t278 - t8 / 0.2e1) * t188 + (t4 * t276 + t3 * t278 + t22 / 0.2e1) * t187;
t12 = [0, t85 * qJD(3) + t35 * qJD(4) + (-m(3) * t103 / 0.2e1 + t59 * t286 + t48 * t285) * t287, t85 * qJD(2), t35 * qJD(2) + t60 * t267; qJD(3) * t86 - qJD(4) * t34, t18 * qJD(3) + t2 * qJD(4) + (m(5) * (t40 * t48 + t72 * t73 + t74 * t75) + m(4) * (t104 * t105 + t106 * t107 + t55 * t59) + m(3) * (-t103 + t168) * t231 * (rSges(3,1) * t188 - t187 * rSges(3,2)) + (t15 + t179 * t161 + (t158 * t120 + t159 * t122) * t183 + ((-t160 + t200) * t183 - t158 * t119 - t159 * t121 + t290 + t201 * t185) * t185) * t278 + (t14 + t180 * t160 - (t156 * t119 + t157 * t121) * t185 + ((-t161 + t201) * t185 + t156 * t120 + t157 * t122 + t290 + t200 * t183) * t183) * t277) * qJD(2), t18 * qJD(2) + t6 * qJD(4) + t241 + (-0.4e1 * t293 + 0.2e1 * t229 * (-t166 * t188 + t232)) * qJD(3), -t242 + t2 * qJD(2) + t6 * qJD(3) + (-t22 / 0.2e1 - t226 * t185 + t227 * t183) * t230 + (t10 * t278 + t9 * t277 + m(5) * (t40 * t60 + t58 * t62 + t70 * t74 + t71 * t72 + (-t183 * t67 - t185 * t66) * t149) - t282 + (t8 / 0.2e1 + (t38 / 0.2e1 - t20 / 0.2e1) * t185 + (-t39 / 0.2e1 - t19 / 0.2e1) * t183) * t188) * qJD(4); -t86 * qJD(2), -t241 + t76 * qJD(3) + t5 * qJD(4) + 0.4e1 * (-t284 / 0.4e1 - t281 / 0.4e1) * qJD(2) + ((-t188 * t59 + t239) * t286 + (-t188 * t48 + t272) * t285 + ((t105 * t183 + t107 * t185 + t55) * t286 + (t183 * t73 + t185 * t75 + t40) * t285) * t187) * t287, t76 * qJD(2), t5 * qJD(2) + (-t188 * t60 + (t183 * t71 + t185 * t70) * t187) * t267; t34 * qJD(2), t242 + (t227 * t185 + t226 * t183 + (-t183 * t37 / 0.2e1 + (t183 * t43 - t185 * t42) * t278 + (t183 * t45 - t185 * t44 + t36) * t276) * t188 + ((-t49 / 0.2e1 + t15 / 0.2e1) * t185 + (t50 / 0.2e1 + t14 / 0.2e1) * t183) * t187 + (t40 * t47 + t48 * t58 + t53 * t74 + t54 * t72 + t66 * t75 + t67 * t73 - t21) * m(5)) * qJD(2) + t7 * qJD(3) + t1 * qJD(4), t7 * qJD(2), t1 * qJD(2) + (m(5) * (t58 * t60 + t66 * t70 + t67 * t71) - t298 * t250 / 0.2e1) * qJD(4) + (t10 * t276 + t9 * t278 - t188 * (t39 * t185 + t38 * t183 - (-t237 * t177 + t238 * t178) * t188) / 0.2e1) * t230;];
Cq = t12;
