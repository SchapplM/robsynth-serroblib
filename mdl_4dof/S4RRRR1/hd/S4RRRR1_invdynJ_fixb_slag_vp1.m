% Calculate vector of inverse dynamics joint torques for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:08
% EndTime: 2019-12-31 17:22:12
% DurationCPUTime: 3.03s
% Computational Cost: add. (8006->324), mult. (5734->419), div. (0->0), fcn. (4348->8), ass. (0->201)
t164 = sin(qJ(1));
t251 = pkin(1) * qJD(1);
t218 = t164 * t251;
t162 = qJ(1) + qJ(2);
t154 = sin(t162);
t161 = qJD(1) + qJD(2);
t230 = t154 * t161;
t221 = pkin(2) * t230;
t157 = qJ(3) + t162;
t148 = sin(t157);
t149 = cos(t157);
t100 = rSges(4,1) * t148 + rSges(4,2) * t149;
t153 = qJD(3) + t161;
t91 = t153 * t100;
t61 = -t218 - t221 - t91;
t143 = t149 * pkin(7);
t102 = pkin(3) * t148 - t143;
t233 = t149 * t153;
t116 = pkin(7) * t233;
t163 = sin(qJ(4));
t235 = t148 * t163;
t122 = rSges(5,2) * t235;
t165 = cos(qJ(4));
t252 = rSges(5,2) * t165;
t216 = qJD(4) * t252;
t193 = rSges(5,3) * t233 + t153 * t122 - t149 * t216;
t222 = qJD(4) * t163;
t213 = t149 * t222;
t228 = t149 * rSges(5,3) + t122;
t234 = t148 * t165;
t73 = rSges(5,1) * t234 - t228;
t65 = t153 * t73;
t279 = -rSges(5,1) * t213 + t153 * t102 + t116 + t193 + t65;
t155 = cos(t162);
t106 = rSges(3,1) * t154 + rSges(3,2) * t155;
t240 = t106 * t161;
t79 = -t218 - t240;
t253 = rSges(5,1) * t165;
t134 = -rSges(5,2) * t163 + t253;
t113 = t134 * qJD(4);
t132 = rSges(5,1) * t163 + t252;
t160 = qJDD(1) + qJDD(2);
t152 = qJDD(3) + t160;
t159 = t161 ^ 2;
t166 = cos(qJ(1));
t167 = qJD(1) ^ 2;
t189 = (-qJDD(1) * t164 - t166 * t167) * pkin(1);
t174 = t189 + (-t154 * t160 - t155 * t159) * pkin(2);
t224 = qJD(4) * t149;
t245 = -t102 - t73;
t270 = -t149 * pkin(3) - t148 * pkin(7);
t232 = t149 * t163;
t219 = rSges(5,2) * t232;
t215 = t153 * t219 + (rSges(5,1) * t222 + t216) * t148;
t231 = t149 * t165;
t271 = rSges(5,1) * t231 + t148 * rSges(5,3);
t44 = t271 * t153 - t215;
t223 = qJD(4) * t153;
t84 = -qJDD(4) * t149 + t148 * t223;
t11 = -t113 * t224 + t84 * t132 + t245 * t152 + (t270 * t153 - t44) * t153 + t174;
t278 = -g(1) + t11;
t236 = t148 * t153;
t78 = rSges(4,1) * t233 - rSges(4,2) * t236;
t277 = -t100 * t152 - t153 * t78 - g(1) + t174;
t229 = t155 * t161;
t96 = rSges(3,1) * t229 - rSges(3,2) * t230;
t276 = -t106 * t160 - t161 * t96 - g(1) + t189;
t147 = pkin(2) * t155;
t158 = t166 * pkin(1);
t261 = pkin(1) * t164;
t205 = qJDD(1) * t158 - t167 * t261;
t260 = pkin(2) * t154;
t183 = t160 * t147 - t159 * t260 + t205;
t225 = qJD(4) * t148;
t43 = (-t153 * t234 - t213) * rSges(5,1) + t193;
t74 = -t219 + t271;
t57 = t74 - t270;
t83 = qJDD(4) * t148 + t149 * t223;
t12 = -t113 * t225 - t83 * t132 + (-pkin(3) * t236 + t116 + t43) * t153 + t57 * t152 + t183;
t275 = -g(2) + t12;
t101 = t149 * rSges(4,1) - rSges(4,2) * t148;
t274 = t101 * t152 - t153 * t91 - g(2) + t183;
t107 = t155 * rSges(3,1) - rSges(3,2) * t154;
t273 = t107 * t160 - t161 * t240 - g(2) + t205;
t156 = Icges(5,4) * t165;
t196 = -Icges(5,2) * t163 + t156;
t129 = Icges(5,1) * t163 + t156;
t269 = -t221 + t279;
t204 = t132 * t225 - t153 * t57;
t126 = Icges(5,5) * t165 - Icges(5,6) * t163;
t125 = Icges(5,5) * t163 + Icges(5,6) * t165;
t179 = Icges(5,3) * t153 - qJD(4) * t125;
t191 = t196 * t149;
t70 = Icges(5,6) * t148 + t191;
t247 = t163 * t70;
t243 = Icges(5,4) * t163;
t130 = Icges(5,1) * t165 - t243;
t192 = t130 * t149;
t72 = Icges(5,5) * t148 + t192;
t198 = -t165 * t72 + t247;
t268 = -t126 * t236 + t149 * t179 + t153 * t198;
t190 = t126 * t149;
t69 = Icges(5,4) * t234 - Icges(5,2) * t235 - Icges(5,6) * t149;
t248 = t163 * t69;
t121 = Icges(5,4) * t235;
t71 = Icges(5,1) * t234 - Icges(5,5) * t149 - t121;
t199 = -t165 * t71 + t248;
t267 = t148 * t179 + (t190 + t199) * t153;
t127 = Icges(5,2) * t165 + t243;
t194 = t163 * t127 - t165 * t129;
t266 = t126 * qJD(4) + t153 * t194;
t67 = Icges(5,5) * t234 - Icges(5,6) * t235 - Icges(5,3) * t149;
t23 = -t148 * t199 - t149 * t67;
t212 = -pkin(3) - t253;
t214 = t132 * t224;
t186 = -t214 - t221;
t177 = t186 - t218;
t29 = t153 * t245 + t177;
t250 = t149 * t29;
t217 = t166 * t251;
t220 = pkin(2) * t229;
t187 = t217 + t220;
t30 = t187 - t204;
t168 = (t212 * t250 + (t29 * (-rSges(5,3) - pkin(7)) + t30 * t212) * t148) * t153;
t265 = t168 + (-t204 + t215) * t29;
t255 = -Icges(5,2) * t234 - t121 + t71;
t257 = t129 * t148 + t69;
t264 = -t163 * t255 - t165 * t257;
t263 = t83 / 0.2e1;
t262 = t84 / 0.2e1;
t259 = -t148 * t67 - t71 * t231;
t68 = Icges(5,3) * t148 + t190;
t258 = t148 * t68 + t72 * t231;
t256 = -t129 * t149 - t70;
t254 = -t127 * t149 + t72;
t238 = t125 * t149;
t48 = -t148 * t194 - t238;
t246 = t48 * t153;
t241 = t101 * t153;
t239 = t125 * t148;
t237 = t126 * t153;
t227 = -t127 + t130;
t226 = t129 + t196;
t211 = -t225 / 0.2e1;
t210 = t225 / 0.2e1;
t209 = -t224 / 0.2e1;
t208 = t224 / 0.2e1;
t58 = t72 * t234;
t207 = t149 * t68 - t58;
t206 = -t67 + t247;
t82 = t101 + t147;
t135 = rSges(2,1) * t166 - rSges(2,2) * t164;
t133 = rSges(2,1) * t164 + rSges(2,2) * t166;
t24 = -t235 * t70 - t207;
t203 = t148 * t24 - t149 * t23;
t25 = -t232 * t69 - t259;
t26 = -t232 * t70 + t258;
t202 = t148 * t26 - t149 * t25;
t201 = -t148 * t30 - t250;
t200 = t148 * t73 + t149 * t74;
t45 = t163 * t71 + t165 * t69;
t46 = t163 * t72 + t165 * t70;
t195 = t127 * t165 + t163 * t129;
t81 = -t100 - t260;
t185 = -t163 * t254 + t165 * t256;
t184 = -t78 - t220;
t56 = t148 * t212 + t143 + t228;
t55 = t147 + t57;
t182 = (-t163 * t226 + t165 * t227) * t153;
t181 = Icges(5,5) * t153 - qJD(4) * t129;
t180 = Icges(5,6) * t153 - qJD(4) * t127;
t54 = t56 - t260;
t49 = -t149 * t194 + t239;
t47 = t49 * t153;
t10 = qJD(4) * t202 + t47;
t109 = t196 * qJD(4);
t110 = t130 * qJD(4);
t39 = t148 * t180 + t153 * t191;
t41 = t148 * t181 + t153 * t192;
t15 = -qJD(4) * t199 + t163 * t41 + t165 * t39;
t38 = t149 * t180 - t196 * t236;
t40 = -t130 * t236 + t149 * t181;
t16 = -qJD(4) * t198 + t163 * t40 + t165 * t38;
t170 = -qJD(4) * t195 - t109 * t163 + t110 * t165 + t125 * t153;
t19 = t266 * t148 + t170 * t149;
t20 = t170 * t148 - t266 * t149;
t9 = qJD(4) * t203 + t246;
t175 = (t47 + ((t24 - t58 + (t68 + t248) * t149 + t259) * t149 + t258 * t148) * qJD(4)) * t208 + (-qJD(4) * t194 + t109 * t165 + t110 * t163) * t153 + (t46 + t49) * t263 + (t45 + t48) * t262 + (-t246 + ((t149 * t206 - t258 + t26) * t149 + (t148 * t206 + t207 + t25) * t148) * qJD(4) + t9) * t211 + (t16 + t19) * t210 + (Icges(4,3) + t195) * t152 + (t15 + t20 + t10) * t209;
t173 = Icges(3,3) * t160 + t175;
t172 = -qJD(4) * t45 + t153 * t67 - t163 * t39 + t165 * t41;
t171 = -qJD(4) * t46 + t153 * t68 - t163 * t38 + t165 * t40;
t93 = t132 * t149;
t92 = t132 * t148;
t80 = t107 * t161 + t217;
t62 = t187 + t241;
t35 = t200 * qJD(4);
t6 = t171 * t148 - t268 * t149;
t5 = t172 * t148 - t267 * t149;
t4 = t268 * t148 + t171 * t149;
t3 = t267 * t148 + t172 * t149;
t1 = [Icges(2,3) * qJDD(1) + t173 + (t273 * (t107 + t158) + t276 * (-t106 - t261) + (-t96 - t217 + t80) * t79) * m(3) + (g(1) * t133 - g(2) * t135 + (t133 ^ 2 + t135 ^ 2) * qJDD(1)) * m(2) + (t29 * (-t187 + t215) + t168 + t275 * (t158 + t55) + t278 * (t54 - t261) + (-t177 + t29 - t218 + t269) * t30) * m(5) + (t274 * (t158 + t82) + t277 * (t81 - t261) + (t62 + t184 - t217) * t61) * m(4); t173 + (t275 * t55 + t278 * t54 + (-t186 + t269) * t30 + t265) * m(5) + (t274 * t82 + t277 * t81 + (t184 + t220 + t241) * t61) * m(4) + (-t79 * t96 - t80 * t240 + (t79 * t161 + t273) * t107 + (t80 * t161 - t276) * t106) * m(3); t175 + (t275 * t57 + t278 * t56 + (t214 + t279) * t30 + t265) * m(5) + (-t61 * t78 - t62 * t91 + (t153 * t61 + t274) * t101 + (t153 * t62 - t277) * t100) * m(4); t10 * t233 / 0.2e1 + t148 * (t152 * t49 + t153 * t19 + t25 * t84 + t26 * t83 + (t148 * t4 - t149 * t3) * qJD(4)) / 0.2e1 + t202 * t263 + ((t153 * t26 - t3) * t149 + (t153 * t25 + t4) * t148) * t210 + t9 * t236 / 0.2e1 - t149 * (t152 * t48 + t153 * t20 + t23 * t84 + t24 * t83 + (t148 * t6 - t149 * t5) * qJD(4)) / 0.2e1 + t203 * t262 + ((t153 * t24 - t5) * t149 + (t153 * t23 + t6) * t148) * t209 + t152 * (t148 * t46 - t149 * t45) / 0.2e1 + t153 * ((t153 * t46 - t15) * t149 + (t153 * t45 + t16) * t148) / 0.2e1 + ((-t225 * t238 + t237) * t148 + (t182 + (-t264 * t149 + (t239 + t185) * t148) * qJD(4)) * t149) * t211 + ((-t224 * t239 - t237) * t149 + (t182 + (t185 * t148 + (-t264 + t238) * t149) * qJD(4)) * t148) * t208 - t153 * ((t163 * t227 + t165 * t226) * t153 + ((t148 * t254 - t149 * t255) * t165 + (t148 * t256 + t149 * t257) * t163) * qJD(4)) / 0.2e1 + ((t73 * t83 - t74 * t84 + (t148 * t44 + t149 * t43) * qJD(4)) * t200 + t35 * ((t43 + t65) * t149 + (-t153 * t74 + t44) * t148) + t201 * t113 + ((-t153 * t30 - t11) * t149 + (t153 * t29 - t12) * t148) * t132 - (t29 * t92 - t30 * t93) * t153 - (t35 * (-t148 * t92 - t149 * t93) + t201 * t134) * qJD(4) + g(1) * t93 + g(2) * t92 - g(3) * t134) * m(5);];
tau = t1;
