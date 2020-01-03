% Calculate vector of inverse dynamics joint torques for
% S4RPPR7
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:36
% EndTime: 2019-12-31 16:41:40
% DurationCPUTime: 2.88s
% Computational Cost: add. (3098->342), mult. (4531->446), div. (0->0), fcn. (3366->6), ass. (0->187)
t152 = pkin(6) + qJ(4);
t133 = sin(t152);
t134 = cos(t152);
t184 = rSges(5,1) * t133 + rSges(5,2) * t134;
t156 = sin(qJ(1));
t226 = t134 * t156;
t113 = Icges(5,4) * t226;
t157 = cos(qJ(1));
t227 = t133 * t156;
t233 = Icges(5,5) * t157;
t64 = Icges(5,1) * t227 + t113 + t233;
t234 = Icges(5,4) * t134;
t173 = Icges(5,1) * t133 + t234;
t65 = -Icges(5,5) * t156 + t157 * t173;
t93 = -Icges(5,2) * t133 + t234;
t73 = t93 * t157;
t163 = t156 * (t65 + t73) - t157 * (-Icges(5,2) * t227 + t113 + t64);
t235 = Icges(5,4) * t133;
t172 = Icges(5,2) * t134 + t235;
t62 = Icges(5,6) * t157 + t156 * t172;
t63 = -Icges(5,6) * t156 + t157 * t172;
t95 = Icges(5,1) * t134 - t235;
t74 = t95 * t156;
t75 = t95 * t157;
t164 = t156 * (t63 - t75) - t157 * (t62 - t74);
t274 = -t164 * t133 + t163 * t134;
t245 = t93 + t173;
t246 = -t172 + t95;
t273 = (t133 * t245 - t134 * t246) * qJD(1);
t272 = -qJ(3) * qJD(1) ^ 2 + qJDD(3);
t110 = t157 * pkin(1) + t156 * qJ(2);
t231 = qJ(3) * t157;
t191 = t110 + t231;
t153 = sin(pkin(6));
t225 = t153 * t156;
t154 = cos(pkin(6));
t241 = rSges(4,2) * t154;
t68 = rSges(4,1) * t225 + t157 * rSges(4,3) + t156 * t241;
t271 = t191 + t68;
t181 = t133 * t65 + t134 * t63;
t268 = t181 * t157;
t111 = -rSges(3,2) * t157 + t156 * rSges(3,3);
t124 = pkin(3) * t225;
t155 = -pkin(5) - qJ(3);
t267 = t155 * t157 - t124;
t171 = Icges(5,5) * t133 + Icges(5,6) * t134;
t61 = -Icges(5,3) * t156 + t157 * t171;
t228 = qJD(1) * t61;
t30 = t133 * t63 - t134 * t65;
t37 = qJD(1) * t62 - qJD(4) * t73;
t39 = -qJD(4) * t75 + (t156 * t173 + t233) * qJD(1);
t266 = t30 * qJD(4) + t133 * t39 + t134 * t37 + t228;
t179 = t133 * t93 - t134 * t95;
t86 = t172 * qJD(4);
t87 = t173 * qJD(4);
t91 = Icges(5,5) * t134 - Icges(5,6) * t133;
t265 = qJD(1) * t91 + qJD(4) * t179 + t133 * t87 + t134 * t86;
t182 = t133 * t62 - t134 * t64;
t60 = Icges(5,3) * t157 + t156 * t171;
t229 = qJD(1) * t60;
t209 = t156 * qJD(4);
t38 = qJD(1) * t63 + t209 * t93;
t40 = qJD(1) * t65 + qJD(4) * t74;
t264 = qJD(4) * t182 - t133 * t40 - t134 * t38 + t229;
t262 = -m(4) - m(5);
t206 = qJD(1) * qJD(4);
t100 = qJDD(4) * t156 + t157 * t206;
t261 = t100 / 0.2e1;
t101 = qJDD(4) * t157 - t156 * t206;
t260 = t101 / 0.2e1;
t259 = t156 / 0.2e1;
t258 = -t157 / 0.2e1;
t257 = rSges(3,2) - pkin(1);
t256 = pkin(3) * t153;
t207 = qJD(1) * qJD(3);
t208 = qJD(1) * qJD(2);
t217 = qJDD(2) * t156 + t157 * t208;
t165 = -0.2e1 * t156 * t207 + t272 * t157 + t217;
t142 = t157 * qJ(2);
t107 = pkin(1) * t156 - t142;
t221 = qJ(3) + t155;
t232 = qJ(3) * t156;
t146 = t156 * rSges(5,3);
t67 = t157 * t184 - t146;
t189 = t156 * t221 + t157 * t256 - t232 + t67;
t170 = -t107 + t189;
t149 = t157 * rSges(5,3);
t240 = rSges(5,2) * t133;
t243 = rSges(5,1) * t134;
t97 = -t240 + t243;
t77 = t97 * t157;
t43 = -qJD(4) * t77 + (t156 * t184 + t149) * qJD(1);
t140 = qJD(2) * t157;
t88 = qJD(1) * t110 - t140;
t89 = t184 * qJD(4);
t7 = -t89 * t209 + t100 * t97 + t170 * qJDD(1) + (-t88 - t43 + (t231 + t267) * qJD(1)) * qJD(1) + t165;
t255 = t7 * t156;
t212 = qJD(1) * t156;
t211 = qJD(1) * t157;
t130 = qJ(2) * t211;
t139 = qJD(2) * t156;
t216 = t130 + t139;
t203 = qJDD(1) * t110 + t156 * t208 + qJD(1) * (-pkin(1) * t212 + t216);
t159 = qJDD(1) * t231 + t272 * t156 + 0.2e1 * t157 * t207 + t203;
t197 = t153 * t211;
t218 = pkin(3) * t197 + t155 * t212;
t66 = rSges(5,1) * t227 + rSges(5,2) * t226 + t149;
t247 = -t157 * t221 + t124 + t66;
t202 = qJD(4) * t243;
t200 = t156 * t202 + t184 * t211;
t201 = qJD(4) * t240;
t44 = (-rSges(5,3) * qJD(1) - t201) * t156 + t200;
t8 = -t101 * t97 + (qJD(4) * t89 - qJDD(2)) * t157 + t247 * qJDD(1) + (qJ(3) * t212 + t218 + t44) * qJD(1) + t159;
t254 = t8 * t157;
t253 = -pkin(1) - qJ(3);
t252 = -pkin(1) + t155;
t238 = rSges(3,3) * t157;
t70 = t156 * t91;
t237 = t157 * t91;
t138 = qJD(3) * t157;
t214 = t138 + t139;
t79 = t97 * t209;
t22 = qJD(1) * t170 + t214 + t79;
t236 = t22 * t157;
t180 = t133 * t95 + t134 * t93;
t28 = t157 * t180 - t70;
t223 = t28 * qJD(1);
t222 = t171 * qJD(1);
t108 = rSges(3,2) * t156 + t238;
t220 = -t107 + t108;
t83 = t110 + t111;
t219 = rSges(4,1) * t197 + t211 * t241;
t215 = rSges(3,2) * t212 + rSges(3,3) * t211;
t213 = -qJD(1) * t107 + t139;
t210 = qJD(4) * t157;
t205 = qJDD(2) * t157;
t204 = -rSges(4,3) + t253;
t18 = t157 * t60 + t62 * t226 + t64 * t227;
t19 = -t157 * t61 - t63 * t226 - t65 * t227;
t199 = t130 + t214;
t198 = t138 + t213;
t196 = -t210 / 0.2e1;
t195 = t210 / 0.2e1;
t194 = -t209 / 0.2e1;
t193 = t209 / 0.2e1;
t185 = rSges(4,1) * t153 + t241;
t69 = -t156 * rSges(4,3) + t157 * t185;
t192 = t69 - t232;
t190 = qJD(3) * t156 - t140;
t188 = -t107 + t192;
t183 = t133 * t64 + t134 * t62;
t161 = qJD(1) * t183 + qJD(4) * t70 + t228;
t162 = -qJD(1) * t181 - qJD(4) * t237 + t229;
t187 = (t161 * t156 + t157 * t264) * t157 + t156 * (t162 * t156 - t157 * t266);
t186 = t156 * (t156 * t266 + t162 * t157) + t157 * (-t156 * t264 + t161 * t157);
t112 = rSges(2,1) * t157 - rSges(2,2) * t156;
t109 = rSges(2,1) * t156 + rSges(2,2) * t157;
t178 = t156 * t19 + t157 * t18;
t57 = t156 * t60;
t20 = -t183 * t157 + t57;
t21 = -t156 * t61 + t268;
t177 = t156 * t21 + t157 * t20;
t23 = -t97 * t210 + (t191 + t247) * qJD(1) + t190;
t176 = t22 * t156 - t157 * t23;
t175 = -t156 * t44 + t157 * t43;
t174 = -t156 * t66 - t157 * t67;
t169 = t184 + t256;
t160 = t180 * qJD(1) - t171 * qJD(4);
t76 = t97 * t156;
t56 = qJD(1) * t83 - t140;
t55 = qJD(1) * t220 + t139;
t42 = qJD(1) * t271 + t190;
t41 = qJD(1) * t188 + t214;
t31 = t174 * qJD(4);
t27 = t156 * t180 + t237;
t26 = t27 * qJD(1);
t25 = qJD(1) * t215 + qJDD(1) * t111 + t203 - t205;
t24 = t220 * qJDD(1) + (-qJD(1) * t111 - t88) * qJD(1) + t217;
t14 = -t205 + qJDD(1) * t68 + qJD(1) * (-rSges(4,3) * t212 + t219) + t159;
t13 = t188 * qJDD(1) + (-qJD(1) * t68 - t88) * qJD(1) + t165;
t12 = t181 * qJD(4) - t133 * t37 + t134 * t39;
t11 = -qJD(4) * t183 - t133 * t38 + t134 * t40;
t10 = -t156 * t265 + t160 * t157;
t9 = t160 * t156 + t157 * t265;
t6 = qJD(4) * t177 - t223;
t5 = qJD(4) * t178 + t26;
t1 = [(-qJD(4) * t180 + t133 * t86 - t134 * t87) * qJD(1) - m(2) * (-g(1) * t109 + g(2) * t112) - t100 * t28 / 0.2e1 + (t26 + ((-t20 + t57 + t19) * t156 + (t21 - t268 + (-t183 + t61) * t156 + t18) * t157) * qJD(4)) * t194 + t30 * t261 + (-t182 + t27) * t260 + (t6 + t223 + (t156 ^ 2 * t61 + (-t57 + t19 + (t183 + t61) * t157) * t157) * qJD(4)) * t196 + (t11 + t10) * t195 + (t12 + t9 + t5) * t193 + (-(qJD(1) * t189 + t198 - t22 + t79) * t23 + t22 * (-t190 + (-t201 + t202) * t157) + t23 * (-t156 * t201 + t199 + t200 + t218) + ((-rSges(5,3) + t252) * t236 + (t22 * (-qJ(2) - t169) + t23 * (-rSges(5,3) - pkin(1))) * t156) * qJD(1) + (-g(2) + t8) * (t110 + t66 - t267) + (-g(1) + t7) * (t156 * t252 + t157 * t169 + t142 - t146)) * m(5) + (-(qJD(1) * t192 + t198 - t41) * t42 - t41 * t190 + t42 * (t199 + t219) + (t41 * t204 * t157 + (t41 * (-qJ(2) - t185) + t42 * t204) * t156) * qJD(1) + (-g(2) + t14) * t271 + (-g(1) + t13) * (t156 * t253 + t142 + t69)) * m(4) + (-(qJD(1) * t108 + t213 - t55) * t56 + t55 * t140 + t56 * (t215 + t216) + (t55 * t257 * t157 + (t55 * (-rSges(3,3) - qJ(2)) - t56 * pkin(1)) * t156) * qJD(1) + (-g(2) + t25) * t83 + (-g(1) + t24) * (t156 * t257 + t142 + t238)) * m(3) + (-t179 + m(2) * (t109 ^ 2 + t112 ^ 2) + Icges(4,1) * t154 ^ 2 + (-0.2e1 * Icges(4,4) * t154 + Icges(4,2) * t153) * t153 + Icges(2,3) + Icges(3,1)) * qJDD(1); (-m(3) + t262) * (g(1) * t156 - g(2) * t157) + 0.2e1 * (t255 / 0.2e1 - t254 / 0.2e1) * m(5) + 0.2e1 * (t13 * t259 + t14 * t258) * m(4) + 0.2e1 * (t24 * t259 + t25 * t258) * m(3); t262 * (g(1) * t157 + g(2) * t156) + m(4) * (t13 * t157 + t14 * t156) + m(5) * (t156 * t8 + t157 * t7); -t5 * t212 / 0.2e1 + t157 * (qJD(1) * t10 + t186 * qJD(4) + qJDD(1) * t27 + t100 * t19 + t101 * t18) / 0.2e1 + t178 * t260 + ((-t18 * t156 + t19 * t157) * qJD(1) + t186) * t195 + t6 * t211 / 0.2e1 + (qJD(1) * t9 + qJD(4) * t187 - qJDD(1) * t28 + t100 * t21 + t101 * t20) * t259 + t177 * t261 + ((-t20 * t156 + t21 * t157) * qJD(1) + t187) * t193 + qJDD(1) * (t156 * t30 - t157 * t182) / 0.2e1 + qJD(1) * (t11 * t157 + t12 * t156 + (t156 * t182 + t30 * t157) * qJD(1)) / 0.2e1 + ((t70 * t210 - t222) * t157 + (-t273 + (-t157 * t237 - t274) * qJD(4)) * t156) * t196 + ((-t209 * t237 - t222) * t156 + (t273 + (t156 * t70 + t274) * qJD(4)) * t157) * t194 - qJD(1) * ((-t246 * t133 - t245 * t134) * qJD(1) + (t133 * t163 + t164 * t134) * qJD(4)) / 0.2e1 + ((qJD(4) * t175 - t100 * t66 - t101 * t67) * t174 + t31 * ((t156 * t67 - t157 * t66) * qJD(1) + t175) - t176 * t89 + (t255 - t254 + (t23 * t156 + t236) * qJD(1)) * t97 - (t22 * t77 + t23 * t76) * qJD(1) - (t31 * (-t156 * t76 - t157 * t77) - t176 * t184) * qJD(4) - g(1) * t76 + g(2) * t77 + g(3) * t184) * m(5);];
tau = t1;
