% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:46
% EndTime: 2019-12-05 17:28:56
% DurationCPUTime: 6.69s
% Computational Cost: add. (8366->419), mult. (8328->582), div. (0->0), fcn. (7746->10), ass. (0->205)
t153 = qJ(1) + pkin(7);
t150 = cos(t153);
t143 = qJD(3) * t150;
t148 = sin(t153);
t155 = sin(pkin(8));
t235 = qJD(4) * t155;
t224 = t148 * t235;
t203 = t143 - t224;
t152 = pkin(9) + qJ(5);
t147 = sin(t152);
t149 = cos(t152);
t157 = cos(pkin(8));
t248 = t148 * t157;
t100 = -t150 * t147 + t149 * t248;
t245 = t150 * t157;
t101 = t147 * t245 - t148 * t149;
t246 = t150 * t155;
t102 = t147 * t148 + t149 * t245;
t91 = Icges(6,4) * t102;
t58 = Icges(6,2) * t101 - Icges(6,6) * t246 - t91;
t90 = Icges(6,4) * t101;
t60 = Icges(6,1) * t102 + Icges(6,5) * t246 - t90;
t99 = t147 * t248 + t149 * t150;
t271 = -t100 * t60 - t58 * t99;
t195 = t101 * t58 + t102 * t60;
t54 = Icges(6,5) * t102 - Icges(6,6) * t101 + Icges(6,3) * t246;
t22 = -t157 * t54 + (t147 * t58 + t149 * t60) * t155;
t285 = t150 * qJ(3);
t156 = cos(pkin(9));
t243 = t156 * t157;
t154 = sin(pkin(9));
t244 = t154 * t157;
t250 = t148 * t154;
t284 = rSges(5,2) * (-t148 * t156 + t150 * t244) - rSges(5,1) * (t150 * t243 + t250);
t247 = t150 * t154;
t283 = -rSges(5,2) * (t148 * t244 + t150 * t156) - rSges(5,1) * (-t148 * t243 + t247);
t141 = -qJD(5) * t157 + qJD(1);
t234 = qJD(5) * t155;
t222 = t150 * t234;
t197 = -t102 * rSges(6,1) + t101 * rSges(6,2);
t228 = rSges(6,3) * t246;
t62 = -t197 + t228;
t98 = -rSges(6,3) * t157 + (rSges(6,1) * t149 - rSges(6,2) * t147) * t155;
t282 = -t141 * t62 + t98 * t222;
t159 = cos(qJ(1));
t273 = t159 * pkin(1);
t205 = -rSges(3,1) * t150 - t273;
t279 = t148 * rSges(3,2) + t205;
t249 = t148 * t155;
t53 = -Icges(6,5) * t100 + Icges(6,6) * t99 - Icges(6,3) * t249;
t254 = Icges(6,4) * t100;
t56 = Icges(6,2) * t99 - Icges(6,6) * t249 - t254;
t89 = Icges(6,4) * t99;
t59 = -Icges(6,1) * t100 - Icges(6,5) * t249 + t89;
t16 = -t101 * t56 + t102 * t59 + t53 * t246;
t237 = qJD(1) * t150;
t142 = qJD(3) * t148;
t238 = qJD(1) * t148;
t239 = -pkin(2) * t238 + t142;
t103 = qJ(3) * t237 + t239;
t133 = t150 * t235;
t236 = qJD(1) * t155;
t226 = t148 * t236;
t276 = pkin(3) * t157;
t240 = qJ(4) * t226 + t238 * t276;
t278 = -t103 - 0.2e1 * t133 + t240;
t170 = t148 * (Icges(6,2) * t100 + t59 + t89) - t150 * (-Icges(6,2) * t102 + t60 - t90);
t171 = t148 * (-Icges(6,1) * t99 - t254 + t56) - t150 * (Icges(6,1) * t101 - t58 + t91);
t160 = qJD(1) ^ 2;
t204 = -t100 * rSges(6,1) + t99 * rSges(6,2);
t168 = rSges(6,3) * t249 - t204;
t71 = qJD(1) * t99 - qJD(5) * t102;
t72 = -qJD(1) * t100 - qJD(5) * t101;
t206 = t72 * rSges(6,1) + t71 * rSges(6,2);
t38 = -rSges(6,3) * t226 + t206;
t225 = t150 * t236;
t73 = qJD(1) * t101 + qJD(5) * t100;
t74 = -qJD(1) * t102 + qJD(5) * t99;
t265 = t74 * rSges(6,1) + t73 * rSges(6,2);
t39 = -rSges(6,3) * t225 + t265;
t161 = t155 * (-t148 * t38 - t150 * t39 + (-t148 * t168 - t150 * t62) * qJD(1));
t9 = qJD(5) * t161;
t277 = m(6) * t9;
t275 = pkin(4) * t154;
t158 = sin(qJ(1));
t274 = t158 * pkin(1);
t272 = -pkin(6) - qJ(4);
t270 = t100 * t59 - t99 * t56;
t264 = t284 * qJD(1);
t262 = rSges(4,1) * t157;
t259 = rSges(6,3) * t155;
t95 = -Icges(6,3) * t157 + (Icges(6,5) * t149 - Icges(6,6) * t147) * t155;
t252 = Icges(6,4) * t149;
t96 = -Icges(6,6) * t157 + (-Icges(6,2) * t147 + t252) * t155;
t253 = Icges(6,4) * t147;
t97 = -Icges(6,5) * t157 + (Icges(6,1) * t149 - t253) * t155;
t30 = -t101 * t96 + t102 * t97 + t246 * t95;
t258 = t141 * t30;
t257 = rSges(4,3) + qJ(3);
t117 = (-Icges(6,2) * t149 - t253) * t155;
t256 = t117 + t97;
t118 = (-Icges(6,1) * t147 - t252) * t155;
t255 = t118 - t96;
t251 = t148 * qJ(3);
t145 = pkin(4) * t156 + pkin(3);
t242 = t157 * t145;
t124 = t150 * pkin(2) + t251;
t123 = qJD(1) * t124;
t241 = t123 - t143;
t233 = qJD(3) * qJD(1);
t232 = pkin(4) * t247;
t231 = t160 * t273;
t230 = qJD(1) * t274;
t229 = rSges(5,3) * t249;
t227 = -t133 - t239;
t223 = t148 * t234;
t221 = -pkin(2) - t262;
t220 = t155 * t272;
t218 = -t234 / 0.2e1;
t217 = t234 / 0.2e1;
t216 = -t124 - t273;
t214 = -pkin(2) - t242;
t213 = qJ(3) + t275;
t212 = t98 * t223;
t211 = t148 * t218;
t210 = t148 * t217;
t209 = t150 * t218;
t208 = t150 * t217;
t207 = qJD(1) * t218;
t200 = rSges(3,1) * t148 + rSges(3,2) * t150;
t199 = -rSges(4,2) * t155 + t262;
t198 = t283 * qJD(1);
t196 = t155 * qJ(4) + t276;
t15 = -t249 * t54 + t271;
t191 = t249 * t53 + t270;
t194 = -t15 * t148 + t150 * t191;
t193 = t148 * (Icges(6,5) * t99 + Icges(6,6) * t100) - t150 * (-Icges(6,5) * t101 - Icges(6,6) * t102);
t192 = -qJD(1) * t241 + t150 * t233 - t231;
t190 = t148 * t207;
t189 = t150 * t207;
t188 = -qJD(1) * (-t148 * pkin(2) + t285) - t142 + t230;
t187 = -rSges(4,1) * t245 - rSges(4,3) * t148;
t176 = qJD(1) * t196;
t186 = t192 + qJD(1) * (-t150 * t176 - t224);
t185 = -pkin(6) * t155 + (pkin(3) - t145) * t157;
t180 = (-rSges(6,1) * t147 - rSges(6,2) * t149) * t155;
t179 = (t148 * t191 + t15 * t150) * t155;
t17 = t54 * t246 + t195;
t178 = (-t148 * t16 + t150 * t17) * t155;
t177 = t220 - t242;
t116 = (-Icges(6,5) * t147 - Icges(6,6) * t149) * t155;
t120 = t196 * t150;
t175 = qJD(1) * t120 + t123 - t203;
t173 = -t273 + t284;
t172 = t148 * t176 - t133 + t188;
t169 = -t155 * rSges(5,3) - pkin(2) - t196;
t167 = -t177 + t259;
t166 = rSges(5,3) * t246 - t173;
t32 = Icges(6,5) * t72 + Icges(6,6) * t71 - Icges(6,3) * t226;
t33 = Icges(6,5) * t74 + Icges(6,6) * t73 - Icges(6,3) * t225;
t34 = Icges(6,4) * t72 + Icges(6,2) * t71 - Icges(6,6) * t226;
t35 = Icges(6,4) * t74 + Icges(6,2) * t73 - Icges(6,6) * t225;
t36 = Icges(6,1) * t72 + Icges(6,4) * t71 - Icges(6,5) * t226;
t37 = Icges(6,1) * t74 + Icges(6,4) * t73 - Icges(6,5) * t225;
t165 = (-(-t101 * t35 + t102 * t37 + t56 * t71 + t59 * t72 + (t150 * t33 - t238 * t53) * t155) * t148 + t150 * (-t101 * t34 + t102 * t36 - t58 * t71 + t60 * t72 + (t150 * t32 - t238 * t54) * t155) + (-t17 * t148 - t16 * t150) * qJD(1)) * t155;
t164 = (qJD(1) * t194 - t148 * (-t100 * t37 + t35 * t99 + t56 * t73 + t59 * t74 + (-t148 * t33 - t237 * t53) * t155) + t150 * (-t100 * t36 + t34 * t99 - t58 * t73 + t60 * t74 + (-t148 * t32 - t237 * t54) * t155)) * t155;
t21 = -t157 * t53 + (-t147 * t56 + t149 * t59) * t155;
t7 = -t157 * t33 + (-t147 * t35 + t149 * t37 + (-t147 * t59 - t149 * t56) * qJD(5)) * t155;
t8 = -t157 * t32 + (-t147 * t34 + t149 * t36 + (-t147 * t60 + t149 * t58) * qJD(5)) * t155;
t163 = (-t148 * t7 + t150 * t8 + (-t148 * t22 - t150 * t21) * qJD(1)) * t155;
t162 = t155 * (-t150 * t204 + (-t62 + t228) * t148);
t146 = t160 * t274;
t139 = rSges(3,2) * t238;
t134 = rSges(4,2) * t246;
t129 = rSges(4,2) * t225;
t126 = t220 * t237;
t109 = qJD(5) * t180;
t108 = qJD(5) * t118;
t107 = qJD(5) * t117;
t106 = qJD(5) * t116;
t88 = -t134 - t187;
t82 = -rSges(6,1) * t101 - rSges(6,2) * t102;
t81 = rSges(6,1) * t99 + rSges(6,2) * t100;
t70 = -pkin(4) * t250 + t150 * t185;
t65 = t143 + (t216 - t88) * qJD(1);
t64 = -qJD(1) * (rSges(4,3) * t150 - t199 * t148) + t188;
t49 = (qJD(1) * t187 + t129) * qJD(1) + t192;
t48 = t146 + (-rSges(4,3) * t237 - t103 + (qJD(1) * t199 - qJD(3)) * t148) * qJD(1);
t41 = (-t120 - t124 - t166) * qJD(1) + t203;
t29 = -t100 * t97 - t249 * t95 + t96 * t99;
t28 = t29 * t141;
t27 = ((-rSges(5,3) * t237 - qJD(4) * t148) * t155 + t264) * qJD(1) + t186;
t26 = t146 + ((rSges(5,3) * t236 - qJD(3)) * t148 + t198 + t278) * qJD(1);
t25 = -qJD(4) * t157 + qJD(5) * t162 + qJD(2);
t24 = -t106 * t157 + (-t107 * t147 + t108 * t149 + (-t147 * t97 - t149 * t96) * qJD(5)) * t155;
t23 = t24 * t141;
t20 = (-t120 + t216 + t70) * qJD(1) + t203 + t282;
t19 = t141 * t168 - qJD(1) * (t148 * t185 + t232) - t212 + t172;
t13 = -t100 * t108 + t107 * t99 + t73 * t96 + t74 * t97 + (-t106 * t148 - t237 * t95) * t155;
t12 = -t101 * t107 + t102 * t108 + t71 * t96 + t72 * t97 + (t106 * t150 - t238 * t95) * t155;
t11 = t109 * t223 + t141 * t39 + (t126 + (-qJD(1) * t275 - t235) * t148 + (t98 * t234 + (t196 - t242) * qJD(1)) * t150) * qJD(1) + t186;
t10 = t109 * t222 - t141 * t38 - t148 * t233 + t146 + (-t240 - t212 + (-t148 * t177 - t232) * qJD(1) + t278) * qJD(1);
t6 = qJD(5) * t178 + t258;
t5 = qJD(5) * t179 + t28;
t1 = [(t28 + (t271 * t150 + (-t17 + t191 + t195) * t148) * t234) * t209 + t23 + m(3) * ((t160 * t200 + t146) * t279 - (qJD(1) * t200 + t230) * (qJD(1) * t205 + t139) + (t231 + (rSges(3,1) * t237 + qJD(1) * t279 - t139) * qJD(1)) * (t200 + t274)) + (t6 - t258 + (-(-t16 - t271) * t148 + (-t195 - t270) * t150 + (-t54 * t148 ^ 2 + (-t148 * t53 - t150 * t54) * t150) * t155 + t194) * t234) * t210 + (t22 + t30) * t190 + (t21 + t29) * t189 + (t10 * (t197 - t273) + t20 * (-t206 + t227) + t11 * (t204 - t274) - t19 * (t126 + t143 + t265) + (t10 * (-pkin(2) - t167) + t11 * t213) * t150 + (-t10 * t213 + t11 * t214 + (t11 * (-rSges(6,3) + t272) + t19 * qJD(4)) * t155) * t148 + ((t20 * t158 + t159 * t19) * pkin(1) + (-t20 * t213 - t19 * (t214 - t259)) * t150 + (t167 * t20 + t19 * t213) * t148) * qJD(1) - (t20 + (-t70 + t273) * qJD(1) + t175 - t282) * t19) * m(6) + (t26 * t173 - t27 * t274 + (t26 * t169 + t27 * (rSges(5,1) * t154 + rSges(5,2) * t156 + qJ(3))) * t150 + (-t26 * qJ(3) + (-pkin(2) + (-rSges(5,1) * t156 + rSges(5,2) * t154 - pkin(3)) * t157 + (-rSges(5,3) - qJ(4)) * t155) * t27) * t148 + (t198 + t227 + t240 + (t229 + t274 - t285) * qJD(1)) * t41 + (-t175 - t264 - t41 + (-t169 * t150 - t166 + t251 + t273) * qJD(1) - t203) * (-qJD(1) * (-t229 - t283) + t172)) * m(5) + (t48 * (t134 - t273) - t65 * t239 - t49 * t274 - t64 * (t129 + t143) + (t221 * t48 + t257 * t49) * t150 + (-t48 * t257 + t49 * (-pkin(2) - t199)) * t148 + ((t65 * t158 + t159 * t64) * pkin(1) + (-t221 * t64 - t257 * t65) * t150 + (t199 * t65 + t257 * t64) * t148) * qJD(1) - (t65 + (t88 + t273) * qJD(1) + t241) * t64) * m(4) + (t7 + t13) * t211 + (t8 + t5 + t12) * t208; t277; m(4) * (t148 * t49 + t150 * t48) + m(5) * (t148 * t27 + t150 * t26) + m(6) * (t10 * t150 + t11 * t148); -t157 * t277 + 0.2e1 * (m(5) * (-t148 * t26 + t150 * t27) / 0.2e1 + m(6) * (-t10 * t148 + t11 * t150) / 0.2e1) * t155; -t157 * (qJD(5) * t163 + t23) / 0.2e1 + t141 * (-t157 * t24 + t163) / 0.2e1 - (qJD(5) * t164 + t13 * t141) * t249 / 0.2e1 + (-t157 * t29 + t179) * t189 + (-t13 * t157 + t164) * t211 + (qJD(5) * t165 + t12 * t141) * t246 / 0.2e1 + (-t157 * t30 + t178) * t190 + (-t12 * t157 + t165) * t208 - t141 * (-t157 * t116 * t141 + ((-t147 * t256 + t149 * t255) * t141 + ((t147 * t170 + t149 * t171) * t155 + t193 * t157) * qJD(5)) * t155) / 0.2e1 + ((-t100 * t255 - t116 * t249 + t256 * t99) * t141 + (-t100 * t171 - t170 * t99 + t193 * t249) * t234) * t210 + ((-t101 * t256 + t102 * t255 + t116 * t246) * t141 + (t101 * t170 + t102 * t171 - t193 * t246) * t234) * t209 - (t148 * t6 + t150 * t5) * t236 / 0.2e1 + (t9 * t162 + t25 * t161 + t10 * (t157 * t62 + t246 * t98) + t20 * (t157 * t38 + (t109 * t150 - t238 * t98) * t155) + t11 * (t157 * t168 + t249 * t98) - t19 * (-t157 * t39 + (t109 * t148 + t237 * t98) * t155) - (-t19 * t81 - t20 * t82) * t141 - (t25 * (-t148 * t82 - t150 * t81) + (-t19 * t148 + t20 * t150) * t180) * t234) * m(6);];
tauc = t1(:);
