% Calculate vector of inverse dynamics joint torques for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:46
% EndTime: 2019-12-31 17:36:52
% DurationCPUTime: 4.30s
% Computational Cost: add. (6158->403), mult. (8219->517), div. (0->0), fcn. (7938->6), ass. (0->183)
t175 = pkin(7) + qJ(2);
t173 = sin(t175);
t177 = cos(pkin(8));
t238 = t173 * t177;
t174 = cos(t175);
t254 = pkin(6) * t174;
t127 = pkin(4) * t238 + t254;
t176 = sin(pkin(8));
t178 = sin(qJ(5));
t256 = cos(qJ(5));
t142 = t176 * t256 - t177 * t178;
t117 = t142 * t173;
t186 = t176 * t178 + t177 * t256;
t118 = t186 * t173;
t68 = -t118 * rSges(6,1) - t117 * rSges(6,2) - rSges(6,3) * t174;
t275 = -t127 + t68;
t119 = t142 * t174;
t120 = t186 * t174;
t109 = Icges(6,4) * t118;
t61 = Icges(6,2) * t117 + Icges(6,6) * t174 + t109;
t108 = Icges(6,4) * t117;
t65 = -Icges(6,1) * t118 - Icges(6,5) * t174 - t108;
t252 = t119 * t61 - t120 * t65;
t58 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t174;
t16 = -t173 * t58 + t252;
t273 = t117 * t61 - t118 * t65;
t28 = -t142 * t65 - t186 * t61;
t14 = t174 * t58 + t273;
t243 = Icges(6,4) * t120;
t63 = Icges(6,2) * t119 - Icges(6,6) * t173 + t243;
t110 = Icges(6,4) * t119;
t66 = Icges(6,1) * t120 - Icges(6,5) * t173 + t110;
t251 = t119 * t63 + t120 * t66;
t60 = Icges(6,5) * t120 + Icges(6,6) * t119 - Icges(6,3) * t173;
t201 = t173 * t60 - t251;
t271 = t201 - t14;
t270 = -m(2) - m(3);
t269 = t117 * t63 + t118 * t66;
t236 = t174 * t177;
t154 = pkin(4) * t236;
t268 = pkin(6) * t173 - t154;
t138 = t174 * pkin(2) + t173 * qJ(3);
t237 = t174 * t176;
t105 = rSges(5,1) * t236 + t173 * rSges(5,2) + rSges(5,3) * t237;
t124 = pkin(3) * t236 + qJ(4) * t237;
t202 = t124 + t138;
t266 = t202 + t105;
t265 = qJD(5) * t142;
t106 = rSges(4,1) * t236 - rSges(4,2) * t237 + t173 * rSges(4,3);
t264 = t173 * (-Icges(6,2) * t120 + t110 + t66) - t174 * (-Icges(6,2) * t118 + t108 - t65);
t262 = -m(5) - m(6);
t218 = qJD(2) * qJD(5);
t129 = -qJDD(5) * t173 - t174 * t218;
t261 = t129 / 0.2e1;
t130 = qJDD(5) * t174 - t173 * t218;
t260 = t130 / 0.2e1;
t259 = t173 / 0.2e1;
t258 = -t174 / 0.2e1;
t257 = -rSges(6,3) - pkin(6);
t253 = g(1) * t173;
t132 = t186 * qJD(5);
t225 = qJD(2) * t173;
t74 = -t132 * t174 - t142 * t225;
t75 = t174 * t265 - t186 * t225;
t248 = t75 * rSges(6,1) + t74 * rSges(6,2);
t242 = Icges(6,4) * t142;
t98 = -Icges(6,2) * t186 + t242;
t247 = -Icges(6,1) * t186 - t242 - t98;
t135 = Icges(6,4) * t186;
t100 = Icges(6,1) * t142 - t135;
t245 = -Icges(6,2) * t142 + t100 - t135;
t234 = t120 * rSges(6,1) + t119 * rSges(6,2);
t69 = -rSges(6,3) * t173 + t234;
t244 = t69 - t268;
t241 = qJ(4) * t176;
t96 = Icges(6,5) * t142 - Icges(6,6) * t186;
t26 = t100 * t118 + t117 * t98 + t174 * t96;
t240 = qJD(2) * t26;
t239 = t173 * t176;
t87 = t106 + t138;
t195 = pkin(3) * t177 + t241;
t123 = t195 * t173;
t166 = t174 * qJ(3);
t136 = pkin(2) * t173 - t166;
t233 = -t123 - t136;
t216 = rSges(4,1) * t238;
t150 = rSges(4,2) * t239;
t229 = t174 * rSges(4,3) + t150;
t104 = t216 - t229;
t232 = -t136 - t104;
t224 = qJD(2) * t174;
t231 = rSges(4,3) * t224 + qJD(2) * t150;
t223 = qJD(4) * t176;
t148 = t174 * t223;
t163 = qJD(3) * t173;
t230 = t148 + t163;
t219 = qJD(2) * qJD(3);
t228 = qJDD(3) * t173 + t174 * t219;
t227 = qJ(3) * t224 + t163;
t226 = -qJD(2) * t136 + t163;
t222 = qJD(5) * t173;
t221 = qJD(5) * t174;
t220 = -m(4) + t262;
t217 = qJDD(4) * t176;
t15 = t174 * t60 + t269;
t170 = t174 * rSges(5,2);
t196 = rSges(5,1) * t177 + rSges(5,3) * t176;
t103 = t173 * t196 - t170;
t214 = -t103 + t233;
t213 = t174 * t217 + t228;
t212 = t148 + t227;
t210 = t173 * t223;
t209 = -rSges(4,1) * t177 - pkin(2);
t208 = -t222 / 0.2e1;
t207 = t222 / 0.2e1;
t206 = -t221 / 0.2e1;
t205 = t221 / 0.2e1;
t204 = t233 + t275;
t203 = -qJD(2) * t123 + t148 + t226;
t161 = -qJDD(4) * t177 + qJDD(1);
t76 = qJD(2) * t119 - t132 * t173;
t77 = qJD(2) * t120 + t173 * t265;
t200 = -rSges(6,1) * t77 - rSges(6,2) * t76;
t33 = Icges(6,5) * t75 + Icges(6,6) * t74 - Icges(6,3) * t224;
t34 = Icges(6,5) * t77 + Icges(6,6) * t76 - Icges(6,3) * t225;
t35 = Icges(6,4) * t75 + Icges(6,2) * t74 - Icges(6,6) * t224;
t36 = Icges(6,4) * t77 + Icges(6,2) * t76 - Icges(6,6) * t225;
t37 = Icges(6,1) * t75 + Icges(6,4) * t74 - Icges(6,5) * t224;
t38 = Icges(6,1) * t77 + Icges(6,4) * t76 - Icges(6,5) * t225;
t199 = (t119 * t36 + t120 * t38 - t173 * t34 - t224 * t58 + t61 * t74 - t65 * t75) * t174 - t173 * (t119 * t35 + t120 * t37 - t173 * t33 - t224 * t60 + t63 * t74 + t66 * t75);
t198 = -t173 * (t117 * t35 + t118 * t37 + t174 * t33 - t225 * t60 + t63 * t76 + t66 * t77) + t174 * (t117 * t36 + t118 * t38 + t174 * t34 - t225 * t58 + t61 * t76 - t65 * t77);
t164 = qJD(3) * t174;
t197 = t164 - t210;
t139 = rSges(3,1) * t174 - rSges(3,2) * t173;
t137 = rSges(3,1) * t173 + rSges(3,2) * t174;
t194 = t14 * t174 - t15 * t173;
t193 = t16 * t174 + t173 * t201;
t192 = t173 * (Icges(6,5) * t119 - Icges(6,6) * t120) - t174 * (Icges(6,5) * t117 - Icges(6,6) * t118);
t189 = -pkin(2) - t195;
t116 = qJD(2) * t138 - t164;
t188 = -t195 * t224 - t116 - 0.2e1 * t210;
t187 = -qJDD(3) * t174 + qJD(2) * (-pkin(2) * t225 + t227) + qJDD(2) * t138 + t173 * t219;
t185 = -(Icges(6,1) * t119 - t243 - t63) * t173 + (Icges(6,1) * t117 - t109 - t61) * t174;
t183 = -pkin(4) * t177 + t189;
t181 = -pkin(2) + (-rSges(5,1) - pkin(3)) * t177 + (-rSges(5,3) - qJ(4)) * t176;
t180 = qJDD(2) * t124 + t173 * t217 + t187 + (-t195 * t225 + 0.2e1 * t148) * qJD(2);
t160 = rSges(5,2) * t224;
t102 = rSges(6,1) * t142 - rSges(6,2) * t186;
t101 = -rSges(6,1) * t186 - rSges(6,2) * t142;
t95 = -Icges(6,5) * t186 - Icges(6,6) * t142;
t92 = -rSges(6,1) * t132 - rSges(6,2) * t265;
t91 = -Icges(6,1) * t132 - Icges(6,4) * t265;
t90 = -Icges(6,4) * t132 - Icges(6,2) * t265;
t89 = -Icges(6,5) * t132 - Icges(6,6) * t265;
t88 = t102 * t221;
t85 = rSges(6,1) * t119 - rSges(6,2) * t120;
t84 = rSges(6,1) * t117 - rSges(6,2) * t118;
t73 = qJD(2) * t87 - t164;
t72 = qJD(2) * t232 + t163;
t46 = qJD(2) * t266 - t197;
t45 = qJD(2) * t214 + t230;
t40 = -rSges(6,3) * t225 - t200;
t39 = -rSges(6,3) * t224 + t248;
t32 = qJDD(2) * t106 + qJD(2) * (-qJD(2) * t216 + t231) + t187;
t31 = t232 * qJDD(2) + (-qJD(2) * t106 - t116) * qJD(2) + t228;
t30 = -qJD(4) * t177 + qJD(1) + (t173 * t68 - t174 * t69) * qJD(5);
t29 = t142 * t66 - t186 * t63;
t27 = t100 * t120 + t119 * t98 - t173 * t96;
t25 = t27 * qJD(2);
t24 = -t164 + (qJD(5) * t102 + t223) * t173 + (t202 + t244) * qJD(2);
t23 = qJD(2) * t204 + t230 + t88;
t19 = qJDD(2) * t105 + (-t196 * t225 + t160) * qJD(2) + t180;
t18 = t214 * qJDD(2) + (-qJD(2) * t105 + t188) * qJD(2) + t213;
t13 = t100 * t77 + t117 * t90 + t118 * t91 + t174 * t89 - t225 * t96 + t76 * t98;
t12 = t100 * t75 + t119 * t90 + t120 * t91 - t173 * t89 - t224 * t96 + t74 * t98;
t11 = -t129 * t68 - t130 * t69 + (-t173 * t40 - t174 * t39) * qJD(5) + t161;
t10 = -t132 * t66 + t142 * t37 - t186 * t35 - t265 * t63;
t9 = t132 * t65 + t142 * t38 - t186 * t36 - t265 * t61;
t8 = t92 * t222 - t129 * t102 + t244 * qJDD(2) + (-t127 * qJD(2) + t39) * qJD(2) + t180;
t7 = t92 * t221 + t130 * t102 + t204 * qJDD(2) + (qJD(2) * t268 + t188 - t40) * qJD(2) + t213;
t6 = qJD(5) * t193 + t25;
t5 = qJD(5) * t194 + t240;
t1 = [m(5) * t161 + m(6) * t11 + (m(4) - t270) * qJDD(1) + (t220 + t270) * g(3); (t25 + (t252 * t174 + (t271 + t273) * t173) * qJD(5)) * t206 - m(3) * (-g(1) * t137 + g(2) * t139) + (-t100 * t132 + t142 * t91 - t186 * t90 - t265 * t98) * qJD(2) + (t29 + t27) * t261 + (t28 + t26) * t260 + (t5 - t240 + ((t251 + t271) * t174 + t269 * t173) * qJD(5)) * t207 + (t12 + t10) * t208 + (t13 + t9 + t6) * t205 + (t23 * (t164 + t200) + t24 * (t212 + t248) + (t7 * t183 - t23 * t223) * t173 + ((t183 * t23 + t24 * t257) * t174 + (t23 * (-qJ(3) - t257) + t24 * t183) * t173) * qJD(2) - (-t241 - pkin(2) + (-pkin(3) - pkin(4)) * t177) * t253 - (t275 * qJD(2) + t203 - t23 + t88) * t24 + (t8 - g(2)) * (t173 * t257 + t154 + t202 + t234) + (t7 - g(1)) * (t166 + t68 - t254)) * m(6) + (-(-qJD(2) * t103 + t203 - t45) * t46 + t45 * t197 + t46 * (t160 + t212) + (t45 * t181 * t174 + (t45 * (-rSges(5,2) - qJ(3)) + t46 * (t189 - t196)) * t173) * qJD(2) + (-g(2) + t19) * t266 + (-g(1) + t18) * (t173 * t181 + t166 + t170)) * m(5) + (-(-qJD(2) * t104 + t226 - t72) * t73 + t72 * t164 + t73 * (t227 + t231) + (t72 * (rSges(4,2) * t176 + t209) * t174 + (t72 * (-rSges(4,3) - qJ(3)) + t73 * t209) * t173) * qJD(2) + (-g(2) + t32) * t87 + (-g(1) + t31) * (t209 * t173 + t166 + t229)) * m(4) + (m(3) * (t137 ^ 2 + t139 ^ 2) + t100 * t142 - t186 * t98 + Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t177 ^ 2 + ((Icges(4,1) + Icges(5,1)) * t176 + 0.2e1 * (Icges(4,4) - Icges(5,5)) * t177) * t176) * qJDD(2); t220 * (-g(2) * t174 + t253) + 0.2e1 * (t258 * t8 + t259 * t7) * m(6) + 0.2e1 * (t18 * t259 + t19 * t258) * m(5) + 0.2e1 * (t258 * t32 + t259 * t31) * m(4); t262 * (-g(3) * t177 + (g(1) * t174 + g(2) * t173) * t176) + m(5) * (-t161 * t177 + t18 * t237 + t19 * t239) + m(6) * (-t11 * t177 + t237 * t7 + t239 * t8); -t6 * t224 / 0.2e1 - t173 * (t12 * qJD(2) + qJD(5) * t199 + t27 * qJDD(2) - t129 * t201 + t16 * t130) / 0.2e1 + t193 * t261 + ((-t16 * t173 + t174 * t201) * qJD(2) + t199) * t208 - t5 * t225 / 0.2e1 + t174 * (t13 * qJD(2) + qJD(5) * t198 + t26 * qJDD(2) + t15 * t129 + t14 * t130) / 0.2e1 + t194 * t260 + ((-t14 * t173 - t15 * t174) * qJD(2) + t198) * t205 + qJDD(2) * (-t29 * t173 + t28 * t174) / 0.2e1 + qJD(2) * (-t10 * t173 + t9 * t174 + (-t28 * t173 - t174 * t29) * qJD(2)) / 0.2e1 + ((t119 * t245 + t120 * t247 - t173 * t95) * qJD(2) + (-t119 * t264 + t120 * t185 + t173 * t192) * qJD(5)) * t207 + ((t117 * t245 + t118 * t247 + t174 * t95) * qJD(2) + (-t117 * t264 + t185 * t118 - t192 * t174) * qJD(5)) * t206 - qJD(2) * ((t247 * t142 - t186 * t245) * qJD(2) + (t142 * t185 + t186 * t264) * qJD(5)) / 0.2e1 + ((t23 * t92 - t11 * t69 + t30 * (qJD(2) * t68 - t39)) * t174 + (t24 * t92 + t11 * t68 + t30 * (qJD(2) * t69 - t40)) * t173 + (t173 * t8 + t7 * t174 + (-t23 * t173 + t24 * t174) * qJD(2)) * t102 - (-t23 * t84 + t24 * t85) * qJD(2) - (t30 * (-t173 * t84 - t174 * t85) + (t24 * t173 + t23 * t174) * t101) * qJD(5) - g(1) * t85 - g(2) * t84 - g(3) * t101) * m(6);];
tau = t1;
