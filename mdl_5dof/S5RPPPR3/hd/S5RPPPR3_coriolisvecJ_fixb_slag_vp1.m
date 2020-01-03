% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:51
% EndTime: 2019-12-31 17:43:58
% DurationCPUTime: 5.26s
% Computational Cost: add. (5576->362), mult. (7744->478), div. (0->0), fcn. (7382->8), ass. (0->176)
t159 = qJ(1) + pkin(7);
t157 = cos(t159);
t161 = cos(pkin(8));
t234 = t157 * t161;
t138 = pkin(4) * t234;
t160 = sin(pkin(8));
t235 = t157 * t160;
t156 = sin(t159);
t123 = t157 * pkin(2) + t156 * qJ(3);
t164 = cos(qJ(1));
t158 = t164 * pkin(1);
t265 = t158 + t123;
t270 = pkin(3) * t234 + qJ(4) * t235 + t265;
t276 = t138 + t270;
t162 = sin(qJ(5));
t254 = cos(qJ(5));
t125 = t160 * t254 - t161 * t162;
t109 = t125 * t157;
t178 = t160 * t162 + t161 * t254;
t110 = t178 * t157;
t108 = t178 * t156;
t102 = Icges(6,4) * t108;
t107 = t125 * t156;
t60 = Icges(6,2) * t107 + Icges(6,6) * t157 + t102;
t101 = Icges(6,4) * t107;
t64 = -Icges(6,1) * t108 - Icges(6,5) * t157 - t101;
t250 = t109 * t60 - t110 * t64;
t57 = Icges(6,5) * t108 + Icges(6,6) * t107 + Icges(6,3) * t157;
t16 = -t156 * t57 + t250;
t215 = rSges(5,1) * t234 + t156 * rSges(5,2) + rSges(5,3) * t235;
t275 = t270 + t215;
t274 = t107 * t60 - t108 * t64;
t26 = -t125 * t64 - t178 * t60;
t165 = qJD(1) ^ 2;
t236 = t156 * t161;
t252 = pkin(6) * t157;
t116 = pkin(4) * t236 + t252;
t66 = rSges(6,1) * t108 + rSges(6,2) * t107 + rSges(6,3) * t157;
t272 = -t116 - t66;
t14 = t157 * t57 + t274;
t239 = Icges(6,4) * t110;
t62 = Icges(6,2) * t109 - Icges(6,6) * t156 + t239;
t103 = Icges(6,4) * t109;
t65 = Icges(6,1) * t110 - Icges(6,5) * t156 + t103;
t249 = t109 * t62 + t110 * t65;
t59 = Icges(6,5) * t110 + Icges(6,6) * t109 - Icges(6,3) * t156;
t198 = t156 * t59 - t249;
t271 = t198 - t14;
t180 = rSges(4,1) * t234 - rSges(4,2) * t235 + t156 * rSges(4,3);
t269 = t180 + t265;
t91 = Icges(6,5) * t125 - Icges(6,6) * t178;
t238 = Icges(6,4) * t125;
t93 = -Icges(6,2) * t178 + t238;
t120 = Icges(6,4) * t178;
t95 = Icges(6,1) * t125 - t120;
t172 = -t107 * t93 - t108 * t95 - t157 * t91;
t268 = t172 * qJD(1);
t267 = 2 * qJD(5);
t266 = t107 * t62 + t108 * t65;
t200 = t157 * rSges(3,1) - rSges(3,2) * t156;
t264 = t158 + t200;
t263 = qJD(5) * t125;
t261 = t156 * (-Icges(6,2) * t110 + t103 + t65) - t157 * (-Icges(6,2) * t108 + t101 - t64);
t225 = qJD(1) * t157;
t117 = t178 * qJD(5);
t226 = qJD(1) * t156;
t71 = -t117 * t157 - t125 * t226;
t72 = t157 * t263 - t178 * t226;
t246 = t72 * rSges(6,1) + t71 * rSges(6,2);
t37 = -rSges(6,3) * t225 + t246;
t73 = qJD(1) * t109 - t117 * t156;
t74 = qJD(1) * t110 + t156 * t263;
t191 = -rSges(6,1) * t74 - rSges(6,2) * t73;
t38 = -rSges(6,3) * t226 - t191;
t232 = t110 * rSges(6,1) + t109 * rSges(6,2);
t68 = -rSges(6,3) * t156 + t232;
t9 = (-t156 * t38 - t157 * t37 + (t156 * t68 - t157 * t66) * qJD(1)) * qJD(5);
t260 = m(6) * t9;
t258 = t156 / 0.2e1;
t257 = -t157 / 0.2e1;
t255 = -rSges(6,3) - pkin(6);
t163 = sin(qJ(1));
t253 = pkin(1) * t163;
t251 = -qJD(1) / 0.2e1;
t245 = -Icges(6,2) * t125 - t120 + t95;
t244 = -Icges(6,1) * t178 - t238 - t93;
t242 = rSges(4,2) * t160;
t241 = -rSges(5,3) - qJ(4);
t146 = qJD(3) * t157;
t106 = qJD(1) * t123 - t146;
t187 = pkin(3) * t161 + qJ(4) * t160;
t224 = qJD(4) * t160;
t212 = t156 * t224;
t240 = -t187 * t225 - t106 - t212;
t97 = rSges(6,1) * t125 - rSges(6,2) * t178;
t182 = qJD(5) * t97 + t224;
t22 = -t146 + t182 * t156 + (-pkin(6) * t156 + t276 + t68) * qJD(1);
t237 = qJD(1) * t22;
t132 = t156 * t242;
t231 = rSges(4,3) * t225 + qJD(1) * t132;
t130 = t157 * t224;
t145 = qJD(3) * t156;
t230 = t130 + t145;
t229 = t157 * rSges(4,3) + t132;
t228 = qJ(3) * t225 + t145;
t148 = t157 * qJ(3);
t121 = pkin(2) * t156 - t148;
t227 = -qJD(1) * t121 + t145;
t223 = qJD(5) * t156;
t222 = qJD(5) * t157;
t221 = qJD(1) * qJD(3);
t15 = t157 * t59 + t266;
t220 = t165 * t253;
t219 = t165 * t158;
t218 = rSges(4,1) * t236;
t83 = t97 * t222;
t216 = t130 + t228;
t211 = -rSges(4,1) * t161 - pkin(2);
t207 = t223 / 0.2e1;
t206 = -t222 / 0.2e1;
t204 = -t121 - t253;
t202 = t148 - t253;
t114 = t187 * t156;
t199 = -qJD(1) * t114 + t130 + t227;
t197 = t157 * t221 - t219;
t196 = -t114 + t204;
t192 = -pkin(2) + (-rSges(5,1) - pkin(3)) * t161;
t122 = rSges(3,1) * t156 + rSges(3,2) * t157;
t188 = rSges(5,1) * t161 + rSges(5,3) * t160;
t186 = t156 * (Icges(6,5) * t109 - Icges(6,6) * t110) - t157 * (Icges(6,5) * t107 - Icges(6,6) * t108);
t181 = qJD(1) * (-pkin(2) * t226 + t228) + t156 * t221 - t220;
t179 = -pkin(2) - t187;
t176 = (t14 * t157 - t15 * t156) * qJD(5);
t175 = (t156 * t198 + t157 * t16) * qJD(5);
t173 = t181 + (-t187 * t226 + 0.2e1 * t130) * qJD(1);
t171 = -(Icges(6,1) * t109 - t239 - t62) * t156 + (Icges(6,1) * t107 - t102 - t60) * t157;
t169 = -pkin(4) * t161 + t179;
t152 = t157 * rSges(5,2);
t144 = rSges(5,2) * t225;
t99 = t218 - t229;
t98 = t156 * t188 - t152;
t90 = -Icges(6,5) * t178 - Icges(6,6) * t125;
t87 = -rSges(6,1) * t117 - rSges(6,2) * t263;
t86 = -Icges(6,1) * t117 - Icges(6,4) * t263;
t85 = -Icges(6,4) * t117 - Icges(6,2) * t263;
t84 = -Icges(6,5) * t117 - Icges(6,6) * t263;
t82 = rSges(6,1) * t109 - rSges(6,2) * t110;
t81 = rSges(6,1) * t107 - rSges(6,2) * t108;
t56 = qJD(1) * t269 - t146;
t55 = t145 + (t204 - t99) * qJD(1);
t51 = (-qJD(1) * t180 - t106) * qJD(1) + t197;
t50 = qJD(1) * (-qJD(1) * t218 + t231) + t181;
t43 = qJD(1) * t275 - t146 + t212;
t42 = (t196 - t98) * qJD(1) + t230;
t36 = Icges(6,1) * t74 + Icges(6,4) * t73 - Icges(6,5) * t226;
t35 = Icges(6,1) * t72 + Icges(6,4) * t71 - Icges(6,5) * t225;
t34 = Icges(6,4) * t74 + Icges(6,2) * t73 - Icges(6,6) * t226;
t33 = Icges(6,4) * t72 + Icges(6,2) * t71 - Icges(6,6) * t225;
t32 = Icges(6,5) * t74 + Icges(6,6) * t73 - Icges(6,3) * t226;
t31 = Icges(6,5) * t72 + Icges(6,6) * t71 - Icges(6,3) * t225;
t30 = (-qJD(1) * t215 - t212 + t240) * qJD(1) + t197;
t29 = (-t188 * t226 + t144) * qJD(1) + t173;
t28 = -qJD(4) * t161 + qJD(2) + (-t156 * t66 - t157 * t68) * qJD(5);
t27 = t125 * t65 - t178 * t62;
t25 = t109 * t93 + t110 * t95 - t156 * t91;
t23 = t25 * qJD(1);
t21 = t83 + (t196 + t272) * qJD(1) + t230;
t13 = t87 * t222 + (-qJD(1) * t138 - t38 + (pkin(6) * qJD(1) - t182) * t156 + t240) * qJD(1) + t197;
t12 = t87 * t223 + (-qJD(1) * t116 + t37 + t83) * qJD(1) + t173;
t11 = t107 * t85 + t108 * t86 + t157 * t84 - t226 * t91 + t73 * t93 + t74 * t95;
t10 = t109 * t85 + t110 * t86 - t156 * t84 - t225 * t91 + t71 * t93 + t72 * t95;
t8 = -t117 * t65 + t125 * t35 - t178 * t33 - t263 * t62;
t7 = t117 * t64 + t125 * t36 - t178 * t34 - t263 * t60;
t6 = t23 + t175;
t5 = t176 - t268;
t1 = [(-t117 * t95 + t125 * t86 - t178 * t85 - t263 * t93) * qJD(1) + (t23 + (t250 * t157 + (t271 + t274) * t156) * qJD(5)) * t206 + m(3) * ((-t122 * t165 - t220) * t264 + (-t219 + (-0.2e1 * t200 - t158 + t264) * t165) * (-t122 - t253)) - (t8 + t10) * t223 / 0.2e1 + (t5 + t268 + ((t249 + t271) * t157 + t266 * t156) * qJD(5)) * t207 + (t13 * (t202 - t66 - t252) + t21 * (t146 + t191) + t12 * (t232 + t276) + t22 * (t216 + t246) + (t12 * t255 + t13 * t169 - t21 * t224) * t156 + ((-t22 * t163 - t21 * t164) * pkin(1) + (t169 * t21 + t22 * t255) * t157 + (t21 * (-qJ(3) - t255) + t22 * t169) * t156) * qJD(1) - (-t21 + t83 + t199) * t22 - (-t253 + t272) * t237) * m(6) + (t30 * (t152 + t202) + t42 * t146 + t29 * t275 + t43 * (t144 + t216) + (t30 * t192 + (-t42 * qJD(4) + t241 * t30) * t160) * t156 + ((-t43 * t163 - t42 * t164) * pkin(1) + t42 * (t160 * t241 + t192) * t157 + (t42 * (-rSges(5,2) - qJ(3)) + t43 * (t179 - t188)) * t156) * qJD(1) - (-t42 + (-t98 - t253) * qJD(1) + t199) * t43) * m(5) + (t51 * (t156 * t211 + t202 + t229) + t55 * t146 + t50 * t269 + t56 * (t228 + t231) + ((-t56 * t163 - t55 * t164) * pkin(1) + t55 * (t211 + t242) * t157 + (t55 * (-rSges(4,3) - qJ(3)) + t56 * t211) * t156) * qJD(1) - (-t55 + (-t99 - t253) * qJD(1) + t227) * t56) * m(4) + (t7 + t11 + t6) * t222 / 0.2e1 + ((t26 - t172) * t156 + (t27 + t25) * t157) * qJD(5) * t251; t260; 0.2e1 * (t12 * t257 + t13 * t258) * m(6) + 0.2e1 * (t257 * t29 + t258 * t30) * m(5) + 0.2e1 * (t257 * t50 + t258 * t51) * m(4); -t161 * t260 + 0.2e1 * (m(5) * (t156 * t29 + t157 * t30) / 0.2e1 + m(6) * (t12 * t156 + t13 * t157) / 0.2e1) * t160; qJD(1) * (-t8 * t156 + t7 * t157 + (-t26 * t156 - t157 * t27) * qJD(1)) / 0.2e1 + ((t109 * t245 + t110 * t244 - t156 * t90) * qJD(1) + (-t109 * t261 + t110 * t171 + t156 * t186) * qJD(5)) * t207 + ((t107 * t245 + t108 * t244 + t157 * t90) * qJD(1) + (-t107 * t261 + t171 * t108 - t186 * t157) * qJD(5)) * t206 + ((t125 * t244 - t178 * t245) * qJD(1) + (t125 * t171 + t178 * t261) * qJD(5)) * t251 - (t10 * qJD(1) + ((t109 * t34 + t110 * t36 - t156 * t32 - t225 * t57 + t60 * t71 - t64 * t72) * t157 - t156 * (t109 * t33 + t110 * t35 - t156 * t31 - t225 * t59 + t62 * t71 + t65 * t72) + (-t16 * t156 + t157 * t198) * qJD(1)) * t267) * t156 / 0.2e1 + (qJD(1) * t11 + (-t156 * (t107 * t33 + t108 * t35 + t157 * t31 - t226 * t59 + t62 * t73 + t65 * t74) + t157 * (t107 * t34 + t108 * t36 + t157 * t32 - t226 * t57 + t60 * t73 - t64 * t74) + (-t14 * t156 - t15 * t157) * qJD(1)) * t267) * t157 / 0.2e1 - (t5 + t176) * t226 / 0.2e1 - (t6 + t175) * t225 / 0.2e1 + ((t21 * t87 - t9 * t68 + t28 * (-qJD(1) * t66 - t37) + (t13 + t237) * t97) * t157 + (t22 * t87 - t9 * t66 + t28 * (qJD(1) * t68 - t38) + (-qJD(1) * t21 + t12) * t97) * t156 - (-t21 * t81 + t22 * t82) * qJD(1) - (t28 * (-t156 * t81 - t157 * t82) + (t22 * t156 + t21 * t157) * (-rSges(6,1) * t178 - rSges(6,2) * t125)) * qJD(5)) * m(6);];
tauc = t1(:);
