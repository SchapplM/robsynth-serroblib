% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:36
% EndTime: 2022-01-23 09:20:43
% DurationCPUTime: 4.52s
% Computational Cost: add. (10877->355), mult. (9167->497), div. (0->0), fcn. (8673->10), ass. (0->186)
t174 = sin(pkin(9));
t175 = cos(pkin(9));
t176 = sin(qJ(5));
t178 = cos(qJ(5));
t173 = qJ(1) + pkin(8);
t170 = qJ(3) + t173;
t166 = sin(t170);
t167 = cos(t170);
t250 = t175 * t176;
t113 = t166 * t250 + t167 * t178;
t249 = t175 * t178;
t114 = t166 * t249 - t167 * t176;
t256 = t166 * t174;
t68 = Icges(6,5) * t114 - Icges(6,6) * t113 + Icges(6,3) * t256;
t107 = Icges(6,4) * t114;
t71 = -Icges(6,2) * t113 + Icges(6,6) * t256 + t107;
t106 = Icges(6,4) * t113;
t75 = -Icges(6,1) * t114 - Icges(6,5) * t256 + t106;
t30 = -(t176 * t71 + t178 * t75) * t174 - t175 * t68;
t203 = (-rSges(6,1) * t176 - rSges(6,2) * t178) * t174;
t129 = qJD(5) * t203;
t239 = qJD(5) * t174;
t291 = t129 * t239;
t290 = -t113 * t71 - t114 * t75;
t115 = t166 * t178 - t167 * t250;
t116 = t166 * t176 + t167 * t249;
t289 = t115 * t71 - t116 * t75;
t180 = qJD(1) ^ 2;
t254 = t167 * t174;
t286 = t68 * t254;
t172 = qJD(1) + qJD(3);
t155 = -qJD(5) * t175 + t172;
t125 = -rSges(6,3) * t175 + (rSges(6,1) * t178 - rSges(6,2) * t176) * t174;
t233 = t125 * t239;
t217 = rSges(6,1) * t114 - rSges(6,2) * t113;
t77 = rSges(6,3) * t256 + t217;
t284 = -t155 * t77 + t166 * t233;
t26 = t286 + t289;
t282 = t26 - t286;
t132 = t167 * pkin(3) + t166 * qJ(4);
t253 = t167 * t175;
t280 = -rSges(5,1) * t253 - t166 * rSges(5,3);
t196 = -rSges(5,2) * t254 + t132 - t280;
t281 = t172 * t196;
t279 = pkin(2) * cos(t173) + cos(qJ(1)) * pkin(1);
t70 = Icges(6,5) * t116 + Icges(6,6) * t115 + Icges(6,3) * t254;
t260 = Icges(6,4) * t116;
t73 = Icges(6,2) * t115 + Icges(6,6) * t254 + t260;
t108 = Icges(6,4) * t115;
t76 = Icges(6,1) * t116 + Icges(6,5) * t254 + t108;
t25 = -t113 * t73 + t114 * t76 + t70 * t256;
t271 = pkin(4) * t175;
t220 = pkin(7) * t174 + t271;
t119 = t220 * t166;
t277 = -t172 * t119 + t284;
t248 = pkin(4) * t253 + pkin(7) * t254 + t132;
t79 = t116 * rSges(6,1) + t115 * rSges(6,2) + rSges(6,3) * t254;
t275 = -t155 * t79 + t167 * t233 - t172 * t248;
t274 = t166 * (-Icges(6,2) * t114 - t106 - t75) + t167 * (-Icges(6,2) * t116 + t108 + t76);
t273 = -t167 / 0.2e1;
t82 = -qJD(5) * t116 + t113 * t172;
t83 = qJD(5) * t115 - t114 * t172;
t268 = t83 * rSges(6,1) + t82 * rSges(6,2);
t267 = rSges(5,1) * t175;
t131 = rSges(4,1) * t166 + rSges(4,2) * t167;
t162 = t167 * rSges(4,1);
t133 = -rSges(4,2) * t166 + t162;
t205 = t279 * qJD(1);
t98 = t133 * t172 + t205;
t265 = t131 * t98;
t121 = -Icges(6,3) * t175 + (Icges(6,5) * t178 - Icges(6,6) * t176) * t174;
t258 = Icges(6,4) * t178;
t122 = -Icges(6,6) * t175 + (-Icges(6,2) * t176 + t258) * t174;
t259 = Icges(6,4) * t176;
t123 = -Icges(6,5) * t175 + (Icges(6,1) * t178 - t259) * t174;
t50 = -t113 * t122 + t114 * t123 + t121 * t256;
t264 = t155 * t50;
t24 = t256 * t68 + t290;
t262 = t166 * t24;
t157 = qJD(4) * t167;
t100 = t132 * t172 - t157;
t187 = -t172 * t220 + t233;
t209 = t279 * t180;
t240 = qJD(4) * t172;
t195 = t167 * t240 - t209;
t84 = -qJD(5) * t114 + t115 * t172;
t85 = -qJD(5) * t113 + t116 * t172;
t223 = -rSges(6,1) * t85 - rSges(6,2) * t84;
t251 = t172 * t174;
t235 = t167 * t251;
t49 = rSges(6,3) * t235 - t223;
t22 = t166 * t291 - t155 * t49 + (t167 * t187 - t100) * t172 + t195;
t261 = t22 * t166;
t257 = t166 * t172;
t255 = t167 * t172;
t252 = t172 * t131;
t136 = (-Icges(6,1) * t176 - t258) * t174;
t247 = -t122 + t136;
t135 = (-Icges(6,2) * t178 - t259) * t174;
t246 = t123 + t135;
t236 = t166 * t251;
t245 = rSges(5,2) * t236 + rSges(5,3) * t255;
t156 = qJD(4) * t166;
t244 = qJ(4) * t255 + t156;
t243 = rSges(5,2) * t256 + t167 * rSges(5,3);
t159 = t167 * qJ(4);
t130 = pkin(3) * t166 - t159;
t124 = t172 * t130;
t242 = t156 - t124;
t27 = t115 * t73 + t116 * t76 + t70 * t254;
t238 = t166 * t267;
t231 = -pkin(3) - t267;
t230 = -t239 / 0.2e1;
t229 = t239 / 0.2e1;
t227 = t166 * t230;
t226 = t166 * t229;
t225 = t167 * t230;
t224 = t167 * t229;
t222 = -pkin(2) * sin(t173) - sin(qJ(1)) * pkin(1);
t206 = t222 * qJD(1);
t192 = t156 + t206;
t34 = (-t119 - t130) * t172 + t192 + t284;
t191 = t205 - t157;
t35 = t191 - t275;
t216 = t166 * t34 - t167 * t35;
t215 = -t166 * t79 + t167 * t77;
t214 = t166 * (-Icges(6,5) * t113 - Icges(6,6) * t114) + t167 * (Icges(6,5) * t115 - Icges(6,6) * t116);
t213 = t172 * t227;
t212 = t172 * t224;
t211 = t159 - t217;
t210 = t222 * t180;
t204 = t174 * t214;
t202 = (t167 * t25 + t262) * t174;
t201 = (t166 * t26 + t167 * t27) * t174;
t134 = (-Icges(6,5) * t176 - Icges(6,6) * t178) * t174;
t199 = -t271 - pkin(3) + (-rSges(6,3) - pkin(7)) * t174;
t198 = t79 + t248;
t194 = (Icges(6,1) * t115 - t260 - t73) * t167 + (-Icges(6,1) * t113 - t107 - t71) * t166;
t190 = t166 * t231 + t159 + t243;
t189 = t166 * t240 + t172 * (-pkin(3) * t257 + t244) + t210;
t188 = -t124 + t192;
t97 = t206 - t252;
t51 = t115 * t122 + t116 * t123 + t121 * t254;
t41 = t51 * t155;
t10 = qJD(5) * t201 + t41;
t43 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t235;
t45 = Icges(6,4) * t85 + Icges(6,2) * t84 + Icges(6,6) * t235;
t47 = Icges(6,1) * t85 + Icges(6,4) * t84 + Icges(6,5) * t235;
t13 = -t175 * t43 + (-t176 * t45 + t178 * t47 + (t176 * t75 - t178 * t71) * qJD(5)) * t174;
t42 = Icges(6,5) * t83 + Icges(6,6) * t82 - Icges(6,3) * t236;
t44 = Icges(6,4) * t83 + Icges(6,2) * t82 - Icges(6,6) * t236;
t46 = Icges(6,1) * t83 + Icges(6,4) * t82 - Icges(6,5) * t236;
t14 = -t175 * t42 + (-t176 * t44 + t178 * t46 + (-t176 * t76 - t178 * t73) * qJD(5)) * t174;
t126 = qJD(5) * t134;
t127 = qJD(5) * t135;
t128 = qJD(5) * t136;
t19 = t115 * t127 + t116 * t128 + t122 * t82 + t123 * t83 + (-t121 * t257 + t126 * t167) * t174;
t20 = -t113 * t127 + t114 * t128 + t122 * t84 + t123 * t85 + (t121 * t255 + t126 * t166) * t174;
t31 = -t175 * t70 + (-t176 * t73 + t178 * t76) * t174;
t40 = -t126 * t175 + (-t127 * t176 + t128 * t178 + (-t122 * t178 - t123 * t176) * qJD(5)) * t174;
t37 = t40 * t155;
t9 = qJD(5) * t202 + t264;
t186 = (t41 + ((t24 + t27 - t290) * t167 + t282 * t166) * t239) * t227 + t37 + (t14 + t19) * t224 + (t31 + t51) * t213 + (t30 + t50) * t212 + (t13 + t20 + t10) * t226 + (-t264 + ((-t25 + t282 - t289) * t167 - t262) * t239 + t9) * t225;
t185 = ((t172 * t24 - t113 * t44 + t114 * t46 + t73 * t84 + t76 * t85 + (t166 * t42 + t255 * t70) * t174) * t167 + (-t172 * t25 - t113 * t45 + t114 * t47 + t71 * t84 - t75 * t85 + (t166 * t43 + t255 * t68) * t174) * t166) * t174;
t184 = ((t172 * t26 + t115 * t44 + t116 * t46 + t73 * t82 + t76 * t83 + (t167 * t42 - t257 * t70) * t174) * t167 + (-t172 * t27 + t115 * t45 + t116 * t47 + t71 * t82 - t75 * t83 + (t167 * t43 - t257 * t68) * t174) * t166) * t174;
t183 = ((t172 * t30 + t14) * t167 + (-t172 * t31 + t13) * t166) * t174;
t139 = rSges(5,2) * t235;
t103 = t238 - t243;
t62 = (-t103 - t130) * t172 + t192;
t63 = t191 + t281;
t182 = t62 * (t139 + t157) + t63 * (t244 + t245) + (t62 * t231 * t167 + (t62 * (-rSges(5,3) - qJ(4)) + t63 * t231) * t166) * t172;
t181 = t34 * (t157 + t223) + t35 * (t244 + t268) + (t34 * t199 * t167 + (-t34 * qJ(4) + t35 * (-rSges(6,3) * t174 - pkin(3) - t220)) * t166) * t172 + t199 * t261;
t145 = rSges(4,2) * t257;
t112 = rSges(4,1) * t255 - t145;
t99 = t172 * t103;
t95 = -t112 * t172 - t209;
t94 = -t172 * t252 + t210;
t93 = rSges(6,1) * t115 - rSges(6,2) * t116;
t92 = -rSges(6,1) * t113 - rSges(6,2) * t114;
t53 = (t280 * t172 - t100 + t139) * t172 + t195;
t52 = t172 * (-t172 * t238 + t245) + t189;
t48 = -rSges(6,3) * t236 + t268;
t36 = t215 * t239 + qJD(2);
t21 = t155 * t48 - t167 * t291 + t187 * t257 + t189;
t16 = ((-t172 * t79 + t49) * t167 + (-t172 * t77 - t48) * t166) * t239;
t1 = [m(4) * (t95 * (-t131 + t222) + t97 * t145 + t94 * (t133 + t279) + (-t162 * t97 - t265) * t172 + (t222 * t98 - t279 * t97) * qJD(1)) + t186 + (t22 * (t211 + t222) + t21 * (t198 + t279) + (t222 * t35 - t279 * t34) * qJD(1) + t181 - (-t34 + t188 + t277) * t35) * m(6) + (-(-t62 - t99 + t188) * t63 + t53 * (t190 + t222) + t52 * (t196 + t279) + (t222 * t63 - t279 * t62) * qJD(1) + t182) * m(5); m(6) * t16; t186 + (t21 * t198 + t22 * t211 + t181 - t34 * (t157 + t275) - t35 * (t242 + t277)) * m(6) + (t53 * t190 + t52 * t196 + t182 - t62 * (t157 - t281) - t63 * (-t99 + t242)) * m(5) + (-(-t133 * t97 - t265) * t172 - t252 * t98 - t112 * t97 - t131 * t95 + t133 * t94) * m(4); 0.2e1 * (t261 / 0.2e1 + t21 * t273) * m(6) + 0.2e1 * (t53 * t166 / 0.2e1 + t52 * t273) * m(5); -t10 * t236 / 0.2e1 + (-t175 * t51 + t201) * t213 + (-t175 * t19 + t184) * t224 + (qJD(5) * t185 + t155 * t20) * t256 / 0.2e1 + (-t175 * t50 + t202) * t212 + (-t175 * t20 + t185) * t226 - t175 * (qJD(5) * t183 + t37) / 0.2e1 + t155 * (-t175 * t40 + t183) / 0.2e1 + ((t115 * t246 + t116 * t247 + t134 * t254) * t155 + (t274 * t115 + t194 * t116 + t167 * t204) * t239) * t225 + ((-t113 * t246 + t114 * t247 + t134 * t256) * t155 + (-t113 * t274 + t114 * t194 + t166 * t204) * t239) * t227 - t155 * (-t175 * t134 * t155 + ((-t176 * t246 + t178 * t247) * t155 + ((-t176 * t274 + t194 * t178) * t174 - t214 * t175) * qJD(5)) * t174) / 0.2e1 + (qJD(5) * t184 + t155 * t19 + t172 * t9) * t254 / 0.2e1 + ((-t21 * t79 + t22 * t77 + t34 * t49 - t35 * t48) * t175 + (t16 * t215 + t36 * (-t166 * t48 + t167 * t49 - t255 * t79 - t257 * t77) + t216 * t129 + ((t172 * t34 - t21) * t167 + (t172 * t35 + t22) * t166) * t125) * t174 - (-t34 * t92 + t35 * t93) * t155 - (t36 * (-t166 * t93 + t167 * t92) + t216 * t203) * t239) * m(6);];
tauc = t1(:);
