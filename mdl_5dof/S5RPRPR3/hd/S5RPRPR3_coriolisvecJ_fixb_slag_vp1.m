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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:51:14
% EndTime: 2019-12-05 17:51:25
% DurationCPUTime: 5.75s
% Computational Cost: add. (10877->375), mult. (9167->521), div. (0->0), fcn. (8673->10), ass. (0->194)
t172 = qJD(1) ^ 2;
t165 = qJ(1) + pkin(8);
t162 = qJ(3) + t165;
t157 = sin(t162);
t158 = cos(t162);
t170 = cos(qJ(5));
t167 = cos(pkin(9));
t168 = sin(qJ(5));
t243 = t167 * t168;
t110 = t157 * t243 + t158 * t170;
t242 = t167 * t170;
t111 = t157 * t242 - t158 * t168;
t113 = t157 * t168 + t158 * t242;
t105 = Icges(6,4) * t113;
t112 = -t157 * t170 + t158 * t243;
t166 = sin(pkin(9));
t246 = t158 * t166;
t69 = -Icges(6,2) * t112 + Icges(6,6) * t246 + t105;
t104 = Icges(6,4) * t112;
t73 = -Icges(6,1) * t113 - Icges(6,5) * t246 + t104;
t266 = t110 * t69 + t111 * t73;
t205 = -t112 * t69 - t113 * t73;
t66 = Icges(6,5) * t113 - Icges(6,6) * t112 + Icges(6,3) * t246;
t31 = -t167 * t66 + (-t168 * t69 - t170 * t73) * t166;
t152 = qJD(4) * t158;
t164 = qJD(1) + qJD(3);
t267 = pkin(4) * t167;
t213 = pkin(7) * t166 + t267;
t185 = -rSges(6,3) * t166 - pkin(3) - t213;
t191 = (-rSges(6,1) * t168 - rSges(6,2) * t170) * t166;
t127 = qJD(5) * t191;
t150 = -qJD(5) * t167 + t164;
t233 = qJD(5) * t166;
t226 = t158 * t233;
t123 = -rSges(6,3) * t167 + (rSges(6,1) * t170 - rSges(6,2) * t168) * t166;
t228 = t123 * t233;
t169 = sin(qJ(1));
t269 = pkin(1) * t172;
t159 = t169 * t269;
t160 = sin(t165);
t268 = pkin(2) * t160;
t234 = t172 * t268 + t159;
t249 = t157 * t164;
t230 = t167 * t249;
t244 = t164 * t166;
t231 = t157 * t244;
t237 = pkin(4) * t230 + pkin(7) * t231;
t78 = -qJD(5) * t113 + t110 * t164;
t79 = -qJD(5) * t112 - t111 * t164;
t47 = rSges(6,1) * t79 + rSges(6,2) * t78 - rSges(6,3) * t231;
t151 = qJD(4) * t157;
t236 = -pkin(3) * t249 + t151;
t247 = t158 * t164;
t97 = qJ(4) * t247 + t236;
t21 = t127 * t226 - t150 * t47 + (-t97 + (-qJD(4) - t228) * t157 + t237) * t164 + t234;
t161 = cos(t165);
t171 = cos(qJ(1));
t270 = pkin(1) * t171;
t214 = pkin(2) * t161 + t270;
t198 = t214 * t172;
t128 = pkin(3) * t158 + qJ(4) * t157;
t235 = -t164 * t128 + t152;
t183 = -t198 + (t152 + t235) * t164;
t197 = t213 * t164;
t227 = t157 * t233;
t229 = t158 * t244;
t80 = qJD(5) * t111 + t112 * t164;
t81 = qJD(5) * t110 - t113 * t164;
t260 = t81 * rSges(6,1) + t80 * rSges(6,2);
t48 = -rSges(6,3) * t229 + t260;
t22 = t127 * t227 + t150 * t48 + (-t197 + t228) * t247 + t183;
t271 = pkin(1) * t169;
t215 = t268 + t271;
t194 = t215 * qJD(1);
t238 = -t164 * (-pkin(3) * t157 + qJ(4) * t158) - t151;
t181 = t194 + t238;
t241 = -t111 * rSges(6,1) + t110 * rSges(6,2);
t248 = t157 * t166;
t200 = -rSges(6,3) * t248 + t241;
t273 = -t123 * t227 - t150 * t200 + t157 * t197;
t33 = t181 + t273;
t116 = t213 * t158;
t274 = t214 * qJD(1);
t184 = t152 - t274;
t206 = rSges(6,1) * t113 - rSges(6,2) * t112;
t232 = rSges(6,3) * t246;
t74 = t206 + t232;
t278 = t123 * t226 - t150 * t74;
t34 = (-t116 - t128) * t164 + t184 + t278;
t279 = t34 * (-t47 - t236 + t237) + (t22 * (-t267 - pkin(3) + (-rSges(6,3) - pkin(7)) * t166) + (t164 * t33 - t21) * qJ(4)) * t157 + (t21 * t185 + t22 * qJ(4) + (-t34 * qJ(4) - t185 * t33) * t164) * t158 - t33 * (t152 + t260);
t208 = rSges(4,1) * t157 + rSges(4,2) * t158;
t258 = rSges(4,1) * t158;
t129 = -t157 * rSges(4,2) + t258;
t245 = t164 * t129;
t94 = -t274 - t245;
t277 = t208 * t94;
t108 = t208 * t164;
t156 = t160 * rSges(3,2);
t259 = rSges(3,1) * t161;
t212 = -t259 - t270;
t276 = t156 + t212;
t65 = -Icges(6,5) * t111 + Icges(6,6) * t110 - Icges(6,3) * t248;
t252 = Icges(6,4) * t111;
t68 = Icges(6,2) * t110 - Icges(6,6) * t248 - t252;
t103 = Icges(6,4) * t110;
t71 = -Icges(6,1) * t111 - Icges(6,5) * t248 + t103;
t26 = -t112 * t68 + t113 * t71 + t65 * t246;
t272 = -t164 * t116 + t278;
t186 = t157 * (Icges(6,2) * t111 + t103 + t71) - t158 * (-Icges(6,2) * t113 - t104 - t73);
t187 = t157 * (-Icges(6,1) * t110 - t252 + t68) - t158 * (Icges(6,1) * t112 + t105 + t69);
t265 = -t110 * t68 + t111 * t71;
t257 = rSges(5,1) * t167;
t256 = rSges(5,2) * t166;
t255 = rSges(5,3) * t158;
t117 = -Icges(6,3) * t167 + (Icges(6,5) * t170 - Icges(6,6) * t168) * t166;
t250 = Icges(6,4) * t170;
t118 = -Icges(6,6) * t167 + (-Icges(6,2) * t168 + t250) * t166;
t251 = Icges(6,4) * t168;
t119 = -Icges(6,5) * t167 + (Icges(6,1) * t170 - t251) * t166;
t50 = -t112 * t118 + t113 * t119 + t117 * t246;
t254 = t150 * t50;
t253 = rSges(5,3) + qJ(4);
t133 = (-Icges(6,1) * t168 - t250) * t166;
t240 = t118 - t133;
t132 = (-Icges(6,2) * t170 - t251) * t166;
t239 = t119 + t132;
t224 = -pkin(3) - t257;
t223 = -t233 / 0.2e1;
t222 = t233 / 0.2e1;
t220 = t157 * t223;
t219 = t157 * t222;
t218 = t158 * t223;
t217 = t158 * t222;
t216 = t164 * t223;
t209 = rSges(3,1) * t160 + rSges(3,2) * t161;
t207 = t256 - t257;
t204 = t157 * (Icges(6,5) * t110 + Icges(6,6) * t111) - t158 * (-Icges(6,5) * t112 - Icges(6,6) * t113);
t203 = t248 * t65 + t265;
t202 = t157 * t216;
t201 = t158 * t216;
t199 = -rSges(5,3) * t157 - t158 * t257;
t25 = -t248 * t66 + t266;
t190 = (t157 * t203 + t158 * t25) * t166;
t27 = t246 * t66 + t205;
t189 = (-t157 * t26 + t158 * t27) * t166;
t131 = (-Icges(6,5) * t168 - Icges(6,6) * t170) * t166;
t182 = t274 - t235;
t10 = qJD(5) * t189 + t254;
t42 = Icges(6,5) * t81 + Icges(6,6) * t80 - Icges(6,3) * t229;
t44 = Icges(6,4) * t81 + Icges(6,2) * t80 - Icges(6,6) * t229;
t46 = Icges(6,1) * t81 + Icges(6,4) * t80 - Icges(6,5) * t229;
t13 = -t167 * t42 + (-t168 * t44 + t170 * t46 + (-t168 * t71 - t170 * t68) * qJD(5)) * t166;
t41 = Icges(6,5) * t79 + Icges(6,6) * t78 - Icges(6,3) * t231;
t43 = Icges(6,4) * t79 + Icges(6,2) * t78 - Icges(6,6) * t231;
t45 = Icges(6,1) * t79 + Icges(6,4) * t78 - Icges(6,5) * t231;
t14 = -t167 * t41 + (-t168 * t43 + t170 * t45 + (t168 * t73 - t170 * t69) * qJD(5)) * t166;
t124 = qJD(5) * t131;
t125 = qJD(5) * t132;
t126 = qJD(5) * t133;
t19 = -t112 * t125 + t113 * t126 + t118 * t78 + t119 * t79 + (-t117 * t249 + t124 * t158) * t166;
t20 = t110 * t125 - t111 * t126 + t118 * t80 + t119 * t81 + (-t117 * t247 - t124 * t157) * t166;
t30 = -t167 * t65 + (-t168 * t68 + t170 * t71) * t166;
t39 = -t124 * t167 + (-t125 * t168 + t126 * t170 + (-t118 * t170 - t119 * t168) * qJD(5)) * t166;
t36 = t39 * t150;
t49 = t110 * t118 - t111 * t119 - t117 * t248;
t40 = t49 * t150;
t9 = qJD(5) * t190 + t40;
t180 = (t40 + (t266 * t158 + (t203 + t205 - t27) * t157) * t233) * t218 + t36 + (-t254 + (-(-t26 - t266) * t157 + t203 * t158 + (-t205 - t265) * t158 - t25 * t157 + (-t66 * t157 ^ 2 + (-t157 * t65 - t158 * t66) * t158) * t166) * t233 + t10) * t219 + (t31 + t50) * t202 + (t30 + t49) * t201 + (t13 + t20) * t220 + (t14 + t19 + t9) * t217;
t179 = t166 * (-t158 * t241 + (-t74 + t232) * t157);
t178 = ((t164 * t203 + t110 * t43 - t111 * t45 + t69 * t80 - t73 * t81 + (-t157 * t41 - t247 * t66) * t166) * t158 + (-t164 * t25 - t110 * t44 + t111 * t46 - t68 * t80 - t71 * t81 - (-t157 * t42 - t247 * t65) * t166) * t157) * t166;
t177 = ((-t164 * t26 - t112 * t43 + t113 * t45 + t69 * t78 - t73 * t79 + (t158 * t41 - t249 * t66) * t166) * t158 + (-t164 * t27 + t112 * t44 - t113 * t46 - t68 * t78 - t71 * t79 - (t158 * t42 - t249 * t65) * t166) * t157) * t166;
t176 = ((-t164 * t30 + t14) * t158 + (-t164 * t31 - t13) * t157) * t166;
t175 = t166 * ((-t164 * t74 - t48) * t158 + (t164 * t200 - t47) * t157);
t136 = rSges(5,1) * t230;
t137 = rSges(5,2) * t229;
t51 = (-t151 + t136 - t97 + (-rSges(5,2) * t248 - t255) * t164) * t164 + t234;
t52 = (t164 * t199 + t137) * t164 + t183;
t95 = t164 * (t157 * t207 + t255);
t60 = -t95 + t181;
t145 = rSges(5,2) * t246;
t100 = -t145 - t199;
t61 = (-t100 - t128) * t164 + t184;
t174 = t61 * (t136 - t236) + (-t51 * t253 + t52 * (-pkin(3) + t207)) * t157 + (t51 * t224 + t52 * t253) * t158 + ((t253 * t60 - t256 * t61) * t157 + (-t224 * t60 - t253 * t61) * t158) * t164 - t60 * (t137 + t152);
t143 = rSges(4,2) * t249;
t109 = -rSges(4,1) * t247 + t143;
t96 = t164 * t100;
t93 = t194 + t108;
t91 = t109 * t164 - t198;
t90 = t108 * t164 + t234;
t89 = -rSges(6,1) * t112 - rSges(6,2) * t113;
t88 = rSges(6,1) * t110 + rSges(6,2) * t111;
t35 = qJD(5) * t179 + qJD(2);
t16 = qJD(5) * t175;
t1 = [m(3) * ((t172 * t209 + t159) * t276 + (t171 * t269 + (-0.2e1 * t156 - t212 + t259 + t276) * t172) * (t209 + t271)) + t180 + (-(t34 + t182 - t272) * t33 + t21 * (-t206 - t214) + t22 * (-t215 + t241) + (t214 * t33 + t215 * t34) * qJD(1) + t279) * m(6) + (-(t61 + t96 + t182) * t60 + t51 * (t145 - t214) - t52 * t215 + (t214 * t60 + t215 * t61) * qJD(1) + t174) * m(5) + (t90 * (-t129 - t214) + t91 * (-t208 - t215) + t277 * t164 + t94 * t194 + (t258 * t164 - t143 - t245 - t94) * t93) * m(4); m(6) * t16; t180 + (-t34 * (t238 + t273) + t33 * (t235 + t272) - t21 * t206 + t22 * t241 + t279) * m(6) + (-t61 * (-t95 + t238) + t60 * (-t96 + t235) + t51 * t145 + t174) * m(5) + (-(t93 * t129 + t277) * t164 + t94 * t108 - t93 * t109 - t90 * t129 - t91 * t208) * m(4); m(5) * (t157 * t52 + t51 * t158) + m(6) * (t157 * t22 + t158 * t21); -t167 * (qJD(5) * t176 + t36) / 0.2e1 + t150 * (-t167 * t39 + t176) / 0.2e1 - (qJD(5) * t178 + t150 * t20) * t248 / 0.2e1 + (-t167 * t49 + t190) * t201 + (-t167 * t20 + t178) * t220 + (qJD(5) * t177 + t150 * t19) * t246 / 0.2e1 + (-t167 * t50 + t189) * t202 + (-t167 * t19 + t177) * t217 - t150 * (-t167 * t131 * t150 + ((-t168 * t239 - t170 * t240) * t150 + ((t168 * t186 + t170 * t187) * t166 + t204 * t167) * qJD(5)) * t166) / 0.2e1 + ((t110 * t239 + t111 * t240 - t131 * t248) * t150 + (-t110 * t186 - t187 * t111 + t204 * t248) * t233) * t219 + ((-t112 * t239 - t113 * t240 + t131 * t246) * t150 + (t112 * t186 + t113 * t187 - t204 * t246) * t233) * t218 - (t10 * t157 + t158 * t9) * t244 / 0.2e1 + (t16 * t179 + t35 * t175 + t21 * (t123 * t246 + t167 * t74) + t34 * (t167 * t47 + (-t123 * t249 + t127 * t158) * t166) + t22 * (t123 * t248 - t167 * t200) - t33 * (-t167 * t48 + (t123 * t247 + t127 * t157) * t166) - (-t33 * t88 - t34 * t89) * t150 - (t35 * (-t157 * t89 - t158 * t88) + (-t157 * t33 + t158 * t34) * t191) * t233) * m(6);];
tauc = t1(:);
