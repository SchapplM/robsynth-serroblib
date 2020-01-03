% Calculate vector of inverse dynamics joint torques for
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:16
% EndTime: 2019-12-31 17:32:22
% DurationCPUTime: 3.81s
% Computational Cost: add. (4804->342), mult. (8149->439), div. (0->0), fcn. (8652->8), ass. (0->176)
t149 = sin(pkin(7));
t151 = cos(pkin(7));
t248 = sin(qJ(3));
t249 = cos(qJ(3));
t136 = -t149 * t249 + t151 * t248;
t122 = t136 * rSges(5,3);
t158 = t149 * t248 + t151 * t249;
t148 = sin(pkin(8));
t150 = cos(pkin(8));
t263 = -rSges(5,1) * t150 + rSges(5,2) * t148;
t165 = t158 * t263 - t122;
t120 = t136 * qJ(4);
t258 = pkin(3) * t158 + t120;
t239 = t165 - t258;
t147 = pkin(8) + qJ(5);
t145 = sin(t147);
t146 = cos(t147);
t226 = Icges(6,4) * t146;
t172 = -Icges(6,2) * t145 + t226;
t55 = Icges(6,6) * t136 + t158 * t172;
t227 = Icges(6,4) * t145;
t174 = Icges(6,1) * t146 - t227;
t58 = Icges(6,5) * t136 + t158 * t174;
t175 = t145 * t55 - t146 * t58;
t170 = Icges(6,5) * t146 - Icges(6,6) * t145;
t52 = Icges(6,3) * t136 + t158 * t170;
t19 = t136 * t52 - t158 * t175;
t51 = Icges(6,3) * t158 - t136 * t170;
t268 = t136 * t51;
t171 = Icges(6,2) * t146 + t227;
t173 = Icges(6,1) * t145 + t226;
t167 = -t145 * t171 + t146 * t173;
t169 = Icges(6,5) * t145 + Icges(6,6) * t146;
t66 = t169 * t136;
t36 = -t158 * t167 - t66;
t266 = qJD(3) * t36;
t176 = t145 * t58 + t146 * t55;
t67 = t169 * t158;
t84 = -rSges(4,1) * t158 + rSges(4,2) * t136;
t264 = t84 * qJD(3);
t134 = -rSges(6,1) * t146 + rSges(6,2) * t145;
t110 = t134 * qJD(5);
t184 = rSges(6,1) * t145 + rSges(6,2) * t146;
t114 = t158 * qJD(3);
t143 = qJDD(2) * t149;
t191 = t136 * pkin(3) - qJ(4) * t158;
t113 = t136 * qJD(3);
t193 = -t114 * pkin(3) - qJ(4) * t113;
t255 = t158 * qJD(4);
t187 = qJD(4) * t114 + qJDD(4) * t136 + t143 + qJD(3) * (t193 + t255) - qJDD(3) * t191;
t208 = qJD(5) * t158;
t142 = pkin(4) * t150 + pkin(3);
t152 = -pkin(6) - qJ(4);
t238 = t113 * t152 - t114 * t142;
t215 = t136 * t146;
t216 = t136 * t145;
t60 = -rSges(6,1) * t215 + rSges(6,2) * t216 + rSges(6,3) * t158;
t259 = -t136 * t142 - t152 * t158 + t60;
t244 = t191 + t259;
t206 = t136 * qJD(5);
t222 = t114 * t145;
t160 = t146 * t206 + t222;
t205 = t145 * qJD(5);
t198 = t136 * t205;
t221 = t114 * t146;
t231 = -rSges(6,1) * t221 - t113 * rSges(6,3);
t34 = rSges(6,1) * t198 + rSges(6,2) * t160 + t231;
t80 = -qJD(5) * t113 + qJDD(5) * t158;
t7 = t110 * t208 - t80 * t184 + t244 * qJDD(3) + (-t193 + t34 + t238) * qJD(3) + t187;
t262 = -g(1) + t7;
t166 = -t113 * rSges(5,3) + t114 * t263;
t63 = t158 * rSges(5,3) + t136 * t263;
t14 = qJD(3) * t166 + qJDD(3) * t63 + t187;
t261 = t14 - g(1);
t112 = t136 * t152;
t121 = t136 * rSges(6,3);
t218 = t158 * t146;
t219 = t158 * t145;
t59 = -rSges(6,1) * t218 + rSges(6,2) * t219 - t121;
t260 = -t142 * t158 + t112 + t59;
t209 = qJD(4) * t136;
t83 = qJD(3) * t191;
t256 = qJD(3) * t63 + t209 - t83;
t254 = t244 * qJD(3) - t83;
t253 = qJD(3) * t258;
t57 = Icges(6,5) * t158 - t136 * t174;
t241 = t171 * t136 + t57;
t54 = Icges(6,6) * t158 - t136 * t172;
t243 = -t173 * t136 + t54;
t252 = t145 * t241 + t146 * t243;
t247 = pkin(3) - t142;
t246 = -t158 * t51 + t57 * t215;
t245 = t57 * t218 + t268;
t242 = t158 * t173 + t55;
t240 = -t158 * t171 + t58;
t210 = qJD(2) * t151;
t190 = t255 - t210;
t61 = -t134 * t158 + t121;
t202 = t158 * t247 + t112 + t120 - t258 - t61;
t21 = qJD(3) * t202 + t184 * t206 + t190;
t235 = t113 * t21;
t234 = t114 * rSges(5,3);
t233 = t114 * rSges(6,3);
t232 = t145 * t54;
t220 = t114 * t152;
t35 = t136 * t167 - t67;
t214 = t35 * qJD(3);
t103 = t114 * qJ(4);
t213 = -t103 - t209;
t212 = -t171 + t174;
t211 = -t172 - t173;
t207 = t170 * qJD(3);
t204 = m(3) + m(4) + m(5);
t203 = qJDD(2) * t151;
t201 = -m(6) - t204;
t197 = -t208 / 0.2e1;
t196 = -t206 / 0.2e1;
t195 = t206 / 0.2e1;
t178 = t145 * t57 + t146 * t54;
t159 = t198 - t221;
t30 = Icges(6,4) * t159 + Icges(6,2) * t160 - Icges(6,6) * t113;
t32 = Icges(6,1) * t159 + Icges(6,4) * t160 - Icges(6,5) * t113;
t154 = qJD(5) * t178 + t145 * t30 - t146 * t32;
t161 = -t113 * t146 - t158 * t205;
t162 = t113 * t145 - t146 * t208;
t29 = Icges(6,4) * t161 + Icges(6,2) * t162 + Icges(6,6) * t114;
t31 = Icges(6,1) * t161 + Icges(6,4) * t162 + Icges(6,5) * t114;
t155 = qJD(5) * t176 + t145 * t29 - t146 * t31;
t177 = -t146 * t57 + t232;
t27 = Icges(6,5) * t161 + Icges(6,6) * t162 + Icges(6,3) * t114;
t28 = Icges(6,5) * t159 + Icges(6,6) * t160 - Icges(6,3) * t113;
t189 = (t113 * t177 + t114 * t51 + t136 * t28 - t154 * t158) * t158 + t136 * (t113 * t175 + t114 * t52 + t136 * t27 - t155 * t158);
t188 = t158 * (-t113 * t51 + t114 * t177 + t136 * t154 + t158 * t28) + t136 * (-t113 * t52 + t114 * t175 + t136 * t155 + t158 * t27);
t186 = -qJD(4) * t113 + qJDD(4) * t158 - t203;
t16 = t216 * t54 - t246;
t17 = t158 * t52 - t58 * t215 + t216 * t55;
t183 = t17 * t136 + t158 * t16;
t18 = -t219 * t54 + t245;
t182 = t19 * t136 + t158 * t18;
t144 = qJD(2) * t149;
t20 = -t184 * t208 + t144 + t209 + t254;
t181 = -t136 * t21 + t158 * t20;
t33 = rSges(6,1) * t161 + rSges(6,2) * t162 + t233;
t180 = t136 * t34 - t158 * t33;
t179 = t136 * t60 - t158 * t61;
t168 = -t145 * t173 - t146 * t171;
t164 = pkin(3) - t263;
t157 = t145 * t240 + t146 * t242;
t156 = (t145 * t211 + t146 * t212) * qJD(3);
t108 = t172 * qJD(5);
t109 = t174 * qJD(5);
t153 = qJD(5) * t168 - t108 * t145 + t109 * t146;
t107 = t170 * qJD(5);
t86 = -rSges(4,1) * t136 - rSges(4,2) * t158;
t79 = qJD(5) * t114 + qJDD(5) * t136;
t78 = -rSges(4,1) * t114 + rSges(4,2) * t113;
t77 = -rSges(4,1) * t113 - rSges(4,2) * t114;
t76 = -t210 + t264;
t73 = t184 * t158;
t72 = t184 * t136;
t65 = -t113 * pkin(3) - t213;
t41 = -qJD(3) * t77 + qJDD(3) * t84 - t203;
t40 = qJD(3) * t78 + qJDD(3) * t86 + t143;
t26 = qJD(3) * t239 + t190;
t25 = t144 + t256;
t22 = qJD(5) * t179 + qJD(1);
t15 = t239 * qJDD(3) + (-t113 * t263 - t234 - t65) * qJD(3) + t186;
t13 = -t107 * t158 + t113 * t169 + t114 * t167 + t136 * t153;
t12 = -t107 * t136 + t113 * t167 - t114 * t169 - t153 * t158;
t11 = qJD(5) * t175 - t145 * t31 - t146 * t29;
t10 = qJD(5) * t177 - t145 * t32 - t146 * t30;
t9 = qJD(5) * t180 + t60 * t79 - t61 * t80 + qJDD(1);
t8 = -t110 * t206 + t79 * t184 + t202 * qJDD(3) + (-t113 * t247 + t103 + t220 - t33 - t65) * qJD(3) + t186;
t6 = qJD(5) * t182 - t266;
t5 = qJD(5) * t183 - t214;
t1 = [m(6) * t9 + (m(2) + t204) * qJDD(1) + (-m(2) + t201) * g(3); t201 * (g(1) * t149 - g(2) * t151) + m(4) * (t149 * t40 - t151 * t41) + m(5) * (t14 * t149 - t15 * t151) + m(6) * (t149 * t7 - t151 * t8) + m(3) * (t149 ^ 2 + t151 ^ 2) * qJDD(2); -qJD(3) * (-qJD(5) * t167 - t108 * t146 - t109 * t145) + t5 * t195 - (t36 - t176) * t79 / 0.2e1 - (t35 - t178) * t80 / 0.2e1 + (-t214 + ((t17 - t18 + t245) * t136 - t246 * t158) * qJD(5) + t11 + t12) * t196 + (t150 ^ 2 * Icges(5,2) - (-Icges(5,1) * t148 - 0.2e1 * Icges(5,4) * t150) * t148 - t168 + Icges(4,3)) * qJDD(3) + (-g(2) * t260 + t262 * t259 + (-t158 * t8 + t235) * (-t134 + t142) + (-(t258 + t260) * qJD(3) + rSges(6,2) * t222 + t231 + t238 + t253) * t20 + (t8 * (-rSges(6,3) + t152) - t22 * (t59 + t61) * qJD(5)) * t136 + (t220 - t233 + t254) * t21) * m(6) + (-g(2) * t239 + t261 * (-t191 + t63) + (-t158 * t164 - t120 - t122) * t15 + (t164 * t113 + t213 - t234 + t256) * t26 + (-qJD(3) * t165 + t166 + t193 + t253) * t25) * m(5) + (-t76 * t77 + (t78 - t264) * (qJD(3) * t86 + t144) + (t76 * qJD(3) - g(1) + t40) * t86 + (-g(2) + t41) * t84) * m(4) + (((-t16 + (-t52 + t232) * t136 - t246) * t136 - (-t17 - (t177 - t52) * t158 + t268) * t158) * qJD(5) + t10 + t13 + t6 + t266) * t197; (-t113 * t26 + t114 * t25 + (qJD(3) * t26 + t261) * t136 - (t25 * qJD(3) + g(2) - t15) * t158) * m(5) + (-qJD(3) * t181 + t114 * t20 - t235 + t262 * t136 - (g(2) - t8) * t158) * m(6); t114 * t6 / 0.2e1 + t136 * (-t12 * qJD(3) + qJD(5) * t189 - t36 * qJDD(3) + t18 * t80 + t19 * t79) / 0.2e1 + t79 * t182 / 0.2e1 + (-t113 * t18 + t114 * t19 + t189) * t195 - t113 * t5 / 0.2e1 + t158 * (-t13 * qJD(3) + qJD(5) * t188 - t35 * qJDD(3) + t16 * t80 + t17 * t79) / 0.2e1 + t80 * t183 / 0.2e1 + (-t113 * t16 + t114 * t17 + t188) * t208 / 0.2e1 - qJDD(3) * (-t136 * t176 - t158 * t178) / 0.2e1 - qJD(3) * (t10 * t158 + t11 * t136 + t113 * t178 - t114 * t176) / 0.2e1 + ((-t67 * t206 + t207) * t136 - (-t156 + (t252 * t158 + (-t66 + t157) * t136) * qJD(5)) * t158) * t196 + (-(-t66 * t208 - t207) * t158 + (-t156 + (t157 * t136 - (-t252 + t67) * t158) * qJD(5)) * t136) * t197 + qJD(3) * (-(t145 * t212 - t146 * t211) * qJD(3) + ((-t136 * t240 - t158 * t241) * t146 + (t136 * t242 + t158 * t243) * t145) * qJD(5)) / 0.2e1 + (t9 * t179 + t22 * (t113 * t61 + t114 * t60 + t180) + t181 * t110 - (-t113 * t20 - t114 * t21 - t136 * t8 + t158 * t7) * t184 - (t20 * t72 + t21 * t73) * qJD(3) - (t22 * (t136 * t72 + t158 * t73) + t181 * t134) * qJD(5) + g(1) * t73 - g(2) * t72 - g(3) * t134) * m(6);];
tau = t1;
