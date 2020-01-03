% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:46
% EndTime: 2019-12-31 17:36:51
% DurationCPUTime: 3.67s
% Computational Cost: add. (5468->338), mult. (7504->454), div. (0->0), fcn. (7272->6), ass. (0->162)
t154 = pkin(7) + qJ(2);
t152 = sin(t154);
t155 = sin(pkin(8));
t156 = cos(pkin(8));
t157 = sin(qJ(5));
t233 = cos(qJ(5));
t125 = t155 * t233 - t156 * t157;
t153 = cos(t154);
t109 = t125 * t153;
t169 = t155 * t157 + t156 * t233;
t110 = t169 * t153;
t108 = t169 * t152;
t102 = Icges(6,4) * t108;
t107 = t125 * t152;
t58 = Icges(6,2) * t107 + Icges(6,6) * t153 + t102;
t101 = Icges(6,4) * t107;
t62 = -Icges(6,1) * t108 - Icges(6,5) * t153 - t101;
t230 = t109 * t58 - t110 * t62;
t55 = Icges(6,5) * t108 + Icges(6,6) * t107 + Icges(6,3) * t153;
t16 = -t152 * t55 + t230;
t252 = t107 * t58 - t108 * t62;
t251 = -t125 * t62 - t169 * t58;
t218 = t152 * t156;
t232 = pkin(6) * t153;
t115 = pkin(4) * t218 + t232;
t64 = rSges(6,1) * t108 + rSges(6,2) * t107 + rSges(6,3) * t153;
t249 = -t115 - t64;
t14 = t153 * t55 + t252;
t220 = Icges(6,4) * t110;
t60 = Icges(6,2) * t109 - Icges(6,6) * t152 + t220;
t103 = Icges(6,4) * t109;
t63 = Icges(6,1) * t110 - Icges(6,5) * t152 + t103;
t229 = t109 * t60 + t110 * t63;
t57 = Icges(6,5) * t110 + Icges(6,6) * t109 - Icges(6,3) * t152;
t183 = t152 * t57 - t229;
t248 = t183 - t14;
t91 = Icges(6,5) * t125 - Icges(6,6) * t169;
t219 = Icges(6,4) * t125;
t93 = -Icges(6,2) * t169 + t219;
t119 = Icges(6,4) * t169;
t95 = Icges(6,1) * t125 - t119;
t165 = -t107 * t93 - t108 * t95 - t153 * t91;
t247 = t165 * qJD(2);
t246 = 2 * qJD(5);
t245 = t107 * t60 + t108 * t63;
t122 = t153 * pkin(2) + t152 * qJ(3);
t216 = t153 * t156;
t217 = t153 * t155;
t171 = rSges(4,1) * t216 - rSges(4,2) * t217 + t152 * rSges(4,3);
t244 = t122 + t171;
t135 = pkin(4) * t216;
t184 = pkin(3) * t216 + qJ(4) * t217 + t122;
t243 = t135 + t184;
t197 = rSges(5,1) * t216 + t152 * rSges(5,2) + rSges(5,3) * t217;
t242 = t184 + t197;
t241 = qJD(5) * t125;
t240 = t152 * (-Icges(6,2) * t110 + t103 + t63) - t153 * (-Icges(6,2) * t108 + t101 - t62);
t205 = qJD(2) * t153;
t116 = t169 * qJD(5);
t206 = qJD(2) * t152;
t71 = -t116 * t153 - t125 * t206;
t72 = t153 * t241 - t169 * t206;
t226 = t72 * rSges(6,1) + t71 * rSges(6,2);
t35 = -rSges(6,3) * t205 + t226;
t73 = qJD(2) * t109 - t116 * t152;
t74 = qJD(2) * t110 + t152 * t241;
t179 = -rSges(6,1) * t74 - rSges(6,2) * t73;
t36 = -rSges(6,3) * t206 - t179;
t213 = t110 * rSges(6,1) + t109 * rSges(6,2);
t66 = -rSges(6,3) * t152 + t213;
t9 = (-t152 * t36 - t153 * t35 + (t152 * t66 - t153 * t64) * qJD(2)) * qJD(5);
t239 = m(6) * t9;
t237 = t152 / 0.2e1;
t236 = -t153 / 0.2e1;
t234 = -rSges(6,3) - pkin(6);
t231 = -qJD(2) / 0.2e1;
t225 = -Icges(6,2) * t125 - t119 + t95;
t224 = -Icges(6,1) * t169 - t219 - t93;
t223 = rSges(4,2) * t155;
t222 = -rSges(5,3) - qJ(4);
t143 = qJD(3) * t153;
t106 = qJD(2) * t122 - t143;
t176 = pkin(3) * t156 + qJ(4) * t155;
t204 = qJD(4) * t155;
t195 = t152 * t204;
t221 = -t176 * t205 - t106 - t195;
t201 = qJD(2) * qJD(3);
t142 = qJD(3) * t152;
t208 = qJ(3) * t205 + t142;
t214 = qJD(2) * (-pkin(2) * t206 + t208) + t152 * t201;
t112 = t176 * t152;
t145 = t153 * qJ(3);
t120 = pkin(2) * t152 - t145;
t212 = -t112 - t120;
t131 = t152 * t223;
t211 = rSges(4,3) * t205 + qJD(2) * t131;
t129 = t153 * t204;
t210 = t129 + t142;
t209 = t153 * rSges(4,3) + t131;
t207 = -qJD(2) * t120 + t142;
t203 = qJD(5) * t152;
t202 = qJD(5) * t153;
t15 = t153 * t57 + t245;
t200 = rSges(4,1) * t218;
t97 = rSges(6,1) * t125 - rSges(6,2) * t169;
t83 = t97 * t202;
t198 = t129 + t208;
t194 = -rSges(4,1) * t156 - pkin(2);
t190 = t203 / 0.2e1;
t189 = -t202 / 0.2e1;
t186 = t214 + (-t176 * t206 + 0.2e1 * t129) * qJD(2);
t185 = -qJD(2) * t112 + t129 + t207;
t180 = -pkin(2) + (-rSges(5,1) - pkin(3)) * t156;
t177 = rSges(5,1) * t156 + rSges(5,3) * t155;
t175 = t152 * (Icges(6,5) * t109 - Icges(6,6) * t110) - t153 * (Icges(6,5) * t107 - Icges(6,6) * t108);
t172 = qJD(5) * t97 + t204;
t170 = -pkin(2) - t176;
t167 = (t14 * t153 - t15 * t152) * qJD(5);
t166 = (t152 * t183 + t153 * t16) * qJD(5);
t164 = -(Icges(6,1) * t109 - t220 - t60) * t152 + (Icges(6,1) * t107 - t102 - t58) * t153;
t162 = -pkin(4) * t156 + t170;
t149 = t153 * rSges(5,2);
t141 = rSges(5,2) * t205;
t138 = t153 * t201;
t99 = t200 - t209;
t98 = t152 * t177 - t149;
t90 = -Icges(6,5) * t169 - Icges(6,6) * t125;
t87 = -rSges(6,1) * t116 - rSges(6,2) * t241;
t86 = -Icges(6,1) * t116 - Icges(6,4) * t241;
t85 = -Icges(6,4) * t116 - Icges(6,2) * t241;
t84 = -Icges(6,5) * t116 - Icges(6,6) * t241;
t82 = rSges(6,1) * t109 - rSges(6,2) * t110;
t81 = rSges(6,1) * t107 - rSges(6,2) * t108;
t70 = qJD(2) * t244 - t143;
t69 = t142 + (-t120 - t99) * qJD(2);
t51 = t138 + (-qJD(2) * t171 - t106) * qJD(2);
t50 = qJD(2) * (-qJD(2) * t200 + t211) + t214;
t43 = qJD(2) * t242 - t143 + t195;
t42 = (-t98 + t212) * qJD(2) + t210;
t38 = t138 + (-qJD(2) * t197 - t195 + t221) * qJD(2);
t37 = (-t177 * t206 + t141) * qJD(2) + t186;
t34 = Icges(6,1) * t74 + Icges(6,4) * t73 - Icges(6,5) * t206;
t33 = Icges(6,1) * t72 + Icges(6,4) * t71 - Icges(6,5) * t205;
t32 = Icges(6,4) * t74 + Icges(6,2) * t73 - Icges(6,6) * t206;
t31 = Icges(6,4) * t72 + Icges(6,2) * t71 - Icges(6,6) * t205;
t30 = Icges(6,5) * t74 + Icges(6,6) * t73 - Icges(6,3) * t206;
t29 = Icges(6,5) * t72 + Icges(6,6) * t71 - Icges(6,3) * t205;
t28 = -qJD(4) * t156 + qJD(1) + (-t152 * t64 - t153 * t66) * qJD(5);
t27 = t125 * t63 - t169 * t60;
t25 = t109 * t93 + t110 * t95 - t152 * t91;
t23 = t25 * qJD(2);
t22 = -t143 + t172 * t152 + (-pkin(6) * t152 + t243 + t66) * qJD(2);
t21 = t83 + (t212 + t249) * qJD(2) + t210;
t13 = t87 * t202 + t138 + (-qJD(2) * t135 - t36 + (pkin(6) * qJD(2) - t172) * t152 + t221) * qJD(2);
t12 = t87 * t203 + (-t115 * qJD(2) + t35 + t83) * qJD(2) + t186;
t11 = t107 * t85 + t108 * t86 + t153 * t84 - t206 * t91 + t73 * t93 + t74 * t95;
t10 = t109 * t85 + t110 * t86 - t152 * t84 - t205 * t91 + t71 * t93 + t72 * t95;
t8 = -t116 * t63 + t125 * t33 - t169 * t31 - t241 * t60;
t7 = t116 * t62 + t125 * t34 - t169 * t32 - t241 * t58;
t6 = t23 + t166;
t5 = t167 - t247;
t1 = [t239; (-t116 * t95 + t125 * t86 - t169 * t85 - t241 * t93) * qJD(2) + (t23 + (t230 * t153 + (t248 + t252) * t152) * qJD(5)) * t189 - (t8 + t10) * t203 / 0.2e1 + (t5 + t247 + ((t229 + t248) * t153 + t245 * t152) * qJD(5)) * t190 + (t13 * (t145 - t64 - t232) + t21 * (t143 + t179) + t12 * (t213 + t243) + t22 * (t198 + t226) + (t12 * t234 + t13 * t162 - t21 * t204) * t152 + ((t162 * t21 + t22 * t234) * t153 + (t21 * (-qJ(3) - t234) + t22 * t162) * t152) * qJD(2) - (qJD(2) * t249 + t185 - t21 + t83) * t22) * m(6) + (t38 * (t145 + t149) + t42 * t143 + t37 * t242 + t43 * (t141 + t198) + (t38 * t180 + (-t42 * qJD(4) + t222 * t38) * t155) * t152 + (t42 * (t155 * t222 + t180) * t153 + (t42 * (-rSges(5,2) - qJ(3)) + t43 * (t170 - t177)) * t152) * qJD(2) - (-qJD(2) * t98 + t185 - t42) * t43) * m(5) + (t51 * (t152 * t194 + t145 + t209) + t69 * t143 + t50 * t244 + t70 * (t208 + t211) + (t69 * (t194 + t223) * t153 + (t69 * (-rSges(4,3) - qJ(3)) + t70 * t194) * t152) * qJD(2) - (-qJD(2) * t99 + t207 - t69) * t70) * m(4) + (t7 + t11 + t6) * t202 / 0.2e1 + ((t27 + t25) * t153 + (-t165 + t251) * t152) * qJD(5) * t231; 0.2e1 * (t12 * t236 + t13 * t237) * m(6) + 0.2e1 * (t236 * t37 + t237 * t38) * m(5) + 0.2e1 * (t236 * t50 + t237 * t51) * m(4); -t156 * t239 + 0.2e1 * (m(5) * (t152 * t37 + t153 * t38) / 0.2e1 + m(6) * (t12 * t152 + t13 * t153) / 0.2e1) * t155; qJD(2) * (-t8 * t152 + t7 * t153 + (-t251 * t152 - t153 * t27) * qJD(2)) / 0.2e1 + ((t109 * t225 + t110 * t224 - t152 * t90) * qJD(2) + (-t109 * t240 + t110 * t164 + t152 * t175) * qJD(5)) * t190 + ((t107 * t225 + t108 * t224 + t153 * t90) * qJD(2) + (-t107 * t240 + t164 * t108 - t175 * t153) * qJD(5)) * t189 + ((t125 * t224 - t169 * t225) * qJD(2) + (t125 * t164 + t169 * t240) * qJD(5)) * t231 - (qJD(2) * t10 + ((t109 * t32 + t110 * t34 - t152 * t30 - t205 * t55 + t58 * t71 - t62 * t72) * t153 - t152 * (t109 * t31 + t110 * t33 - t152 * t29 - t205 * t57 + t60 * t71 + t63 * t72) + (-t16 * t152 + t153 * t183) * qJD(2)) * t246) * t152 / 0.2e1 + (t11 * qJD(2) + (-(t107 * t31 + t108 * t33 + t153 * t29 - t206 * t57 + t60 * t73 + t63 * t74) * t152 + t153 * (t107 * t32 + t108 * t34 + t153 * t30 - t206 * t55 + t58 * t73 - t62 * t74) + (-t14 * t152 - t15 * t153) * qJD(2)) * t246) * t153 / 0.2e1 - (t5 + t167) * t206 / 0.2e1 - (t6 + t166) * t205 / 0.2e1 + ((t21 * t87 - t9 * t66 + t28 * (-qJD(2) * t64 - t35) + (qJD(2) * t22 + t13) * t97) * t153 + (t22 * t87 - t9 * t64 + t28 * (qJD(2) * t66 - t36) + (-qJD(2) * t21 + t12) * t97) * t152 - (-t21 * t81 + t22 * t82) * qJD(2) - (t28 * (-t152 * t81 - t153 * t82) + (t152 * t22 + t153 * t21) * (-rSges(6,1) * t169 - rSges(6,2) * t125)) * qJD(5)) * m(6);];
tauc = t1(:);
