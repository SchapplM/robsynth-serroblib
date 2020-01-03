% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:05
% EndTime: 2019-12-31 16:49:09
% DurationCPUTime: 2.28s
% Computational Cost: add. (13381->245), mult. (10475->344), div. (0->0), fcn. (9716->8), ass. (0->163)
t193 = qJ(3) + qJ(4);
t190 = cos(t193);
t183 = Icges(5,4) * t190;
t189 = sin(t193);
t147 = -Icges(5,2) * t189 + t183;
t148 = Icges(5,1) * t189 + t183;
t308 = t147 + t148;
t192 = qJ(1) + pkin(7);
t187 = sin(t192);
t286 = t187 / 0.2e1;
t188 = cos(t192);
t285 = -t188 / 0.2e1;
t307 = t188 / 0.2e1;
t185 = t187 ^ 2;
t186 = t188 ^ 2;
t227 = t185 + t186;
t194 = sin(qJ(3));
t196 = cos(qJ(3));
t163 = rSges(4,1) * t194 + rSges(4,2) * t196;
t141 = t163 * t187;
t142 = t163 * t188;
t261 = Icges(4,4) * t194;
t159 = Icges(4,2) * t196 + t261;
t162 = Icges(4,1) * t196 - t261;
t150 = rSges(5,1) * t189 + rSges(5,2) * t190;
t281 = pkin(3) * t194;
t204 = t150 + t281;
t302 = t204 * t188;
t303 = t204 * t187;
t253 = t187 * t190;
t254 = t187 * t189;
t105 = rSges(5,1) * t253 - rSges(5,2) * t254 - rSges(5,3) * t188;
t280 = pkin(3) * t196;
t184 = pkin(2) + t280;
t297 = -pkin(6) - pkin(5);
t229 = -t184 * t187 - t188 * t297;
t283 = sin(qJ(1)) * pkin(1);
t74 = -t105 + t229 - t283;
t172 = t187 * t297;
t248 = t188 * t189;
t222 = -rSges(5,2) * t248 + t187 * rSges(5,3);
t271 = rSges(5,1) * t190;
t282 = cos(qJ(1)) * pkin(1);
t75 = t282 - t172 + (t184 + t271) * t188 + t222;
t182 = t188 * pkin(5);
t272 = rSges(4,1) * t196;
t225 = pkin(2) + t272;
t252 = t187 * t194;
t228 = rSges(4,2) * t252 + rSges(4,3) * t188;
t78 = -t187 * t225 + t182 + t228 - t283;
t270 = rSges(4,2) * t194;
t171 = t188 * t270;
t79 = t282 - t171 + t225 * t188 + (rSges(4,3) + pkin(5)) * t187;
t305 = -(t162 / 0.2e1 - t159 / 0.2e1) * t194 - m(5) * (-t302 * t75 + t303 * t74) - m(4) * (t141 * t78 - t142 * t79);
t151 = -rSges(5,2) * t189 + t271;
t211 = Icges(5,5) * t189 + Icges(5,6) * t190;
t119 = t211 * t187;
t120 = t188 * t211;
t260 = Icges(5,4) * t189;
t149 = Icges(5,1) * t190 - t260;
t104 = Icges(5,5) * t187 + t149 * t188;
t146 = Icges(5,2) * t190 + t260;
t234 = -t146 * t188 + t104;
t102 = Icges(5,6) * t187 + t147 * t188;
t236 = -t148 * t188 - t102;
t201 = -t189 * t234 + t190 * t236;
t154 = Icges(5,4) * t254;
t103 = Icges(5,1) * t253 - Icges(5,5) * t188 - t154;
t235 = -Icges(5,2) * t253 + t103 - t154;
t101 = Icges(5,4) * t253 - Icges(5,2) * t254 - Icges(5,6) * t188;
t237 = t148 * t187 + t101;
t202 = t189 * t235 + t190 * t237;
t278 = (-t185 * t120 + (t202 * t188 + (t119 + t201) * t187) * t188) * t286 + (-t186 * t119 + (t201 * t187 + (t120 + t202) * t188) * t187) * t285;
t247 = t188 * t190;
t66 = t187 * t105 + t188 * (rSges(5,1) * t247 + t222);
t44 = -t187 * (pkin(2) * t187 - t182 + t229) + (-t187 * pkin(5) - t172 + (-pkin(2) + t184) * t188) * t188 + t66;
t125 = t150 * t187;
t126 = t150 * t188;
t73 = -t125 * t187 - t126 * t188;
t6 = t278 + m(5) * (t44 * t73 + (t187 * t303 + t188 * t302) * t151);
t304 = t6 * qJD(4);
t191 = Icges(4,4) * t196;
t160 = -Icges(4,2) * t194 + t191;
t161 = Icges(4,1) * t194 + t191;
t217 = t308 * t190 / 0.2e1 + (-t146 / 0.2e1 + t149 / 0.2e1) * t189;
t145 = Icges(5,5) * t190 - Icges(5,6) * t189;
t256 = t145 * t188;
t100 = Icges(5,3) * t187 + t256;
t99 = Icges(5,5) * t253 - Icges(5,6) * t254 - Icges(5,3) * t188;
t220 = t102 * t189 - t99;
t80 = t104 * t253;
t221 = t188 * t100 - t80;
t257 = t101 * t189;
t276 = t100 * t187 + t104 * t247;
t277 = -t103 * t247 - t187 * t99;
t53 = -t102 * t254 - t221;
t54 = -t101 * t248 - t277;
t55 = -t102 * t248 + t276;
t224 = ((t53 - t80 + (t100 + t257) * t188 + t277) * t188 + t276 * t187) * t285 + (t187 * t55 - t188 * t54) * t307 + ((t220 * t187 + t221 + t53 + t54) * t187 + ((-t103 * t190 + t257) * t187 - t276 + t55 + (t99 + t220) * t188) * t188) * t286;
t300 = 2 * qJD(3);
t299 = m(4) / 0.2e1;
t298 = m(5) / 0.2e1;
t207 = (-t187 * t75 - t188 * t74) * t151;
t292 = m(5) * (-t125 * t302 + t126 * t303 + t207);
t291 = m(5) * (t207 + (t187 * t302 - t188 * t303) * t150);
t289 = m(5) * (t125 * t74 - t126 * t75);
t288 = m(5) * t73;
t287 = -t187 / 0.2e1;
t251 = t187 * t196;
t112 = Icges(4,5) * t251 - Icges(4,6) * t252 - Icges(4,3) * t188;
t169 = Icges(4,4) * t252;
t116 = Icges(4,1) * t251 - Icges(4,5) * t188 - t169;
t240 = t196 * t116;
t274 = -t187 * t112 - t188 * t240;
t213 = Icges(4,5) * t196 - Icges(4,6) * t194;
t113 = Icges(4,3) * t187 + t188 * t213;
t117 = Icges(4,5) * t187 + t162 * t188;
t239 = t196 * t117;
t273 = t113 * t187 + t188 * t239;
t114 = Icges(4,4) * t251 - Icges(4,2) * t252 - Icges(4,6) * t188;
t242 = t194 * t114;
t115 = Icges(4,6) * t187 + t160 * t188;
t241 = t194 * t115;
t233 = t161 * t187 + t114;
t232 = -t161 * t188 - t115;
t231 = -Icges(4,2) * t251 + t116 - t169;
t230 = -t159 * t188 + t117;
t223 = -t151 - t280;
t89 = t187 * t239;
t219 = t188 * t113 - t89;
t218 = -t112 + t241;
t212 = -Icges(4,5) * t194 - Icges(4,6) * t196;
t198 = (-t146 + t149) * t190 - t308 * t189;
t205 = -t224 + (t145 * t187 + t188 * t198 + t189 * t236 + t190 * t234) * t286 + (t187 * t198 - t189 * t237 + t190 * t235 - t256) * t285;
t203 = -t217 + (t286 + t287) * (t190 * t101 + t189 * t103);
t200 = t194 * t231 + t196 * t233;
t199 = -t194 * t230 + t196 * t232;
t165 = -t270 + t272;
t136 = t212 * t188;
t135 = t212 * t187;
t109 = t223 * t188;
t107 = t223 * t187;
t77 = -t141 * t187 - t142 * t188;
t70 = qJD(4) * t288;
t63 = -t227 * t281 + t73;
t59 = -t188 * t241 + t273;
t58 = -t188 * t242 - t274;
t57 = -t187 * t241 - t219;
t38 = t187 * t59 - t188 * t58;
t37 = t187 * t57 - t188 * (-(-t240 + t242) * t187 - t188 * t112);
t34 = t217 + t289;
t29 = t291 / 0.2e1;
t24 = t292 / 0.2e1;
t19 = (t161 / 0.2e1 + t160 / 0.2e1) * t196 + t217 - t305;
t16 = (t57 - t89 + (t113 + t242) * t188 + t274) * t188 + t273 * t187;
t15 = (t188 * t218 - t273 + t59) * t188 + (t187 * t218 + t219 + t58) * t187;
t8 = m(5) * (t150 * t151 * t227 + t66 * t73) + t278;
t7 = t8 * qJD(4);
t4 = t24 - t291 / 0.2e1 + t224;
t3 = t29 - t292 / 0.2e1 + t224;
t2 = t24 + t29 + t205;
t1 = (t38 / 0.2e1 - t16 / 0.2e1) * t188 + (t15 / 0.2e1 + t37 / 0.2e1) * t187 + t224;
t5 = [qJD(3) * t19 + qJD(4) * t34, 0, t19 * qJD(1) + t2 * qJD(4) + (((-t187 * t79 - t188 * t78) * t165 + (-t141 * t188 + t142 * t187) * t163) * t299 + (t107 * t75 + t109 * t74) * t298) * t300 + (t16 * t307 + (t194 * t232 + t196 * t230) * t286 + t205 + (t185 / 0.2e1 + t186 / 0.2e1) * t213 + (t15 + t37) * t287 + (-t194 * t233 + t196 * t231 + t38) * t285) * qJD(3), t34 * qJD(1) + t2 * qJD(3) + (t205 + (t207 + (-t125 * t188 + t126 * t187) * t150) * m(5)) * qJD(4); 0, 0, (t298 * t63 + t299 * t77) * t300 + t70, qJD(3) * t288 + t70; t1 * qJD(3) + t4 * qJD(4) + (t203 - (t161 + t160) * t196 / 0.2e1 + t305) * qJD(1), 0, t1 * qJD(1) + (m(4) * (t227 * t165 * t163 + (t187 * (rSges(4,1) * t251 - t228) + t188 * (t187 * rSges(4,3) + t188 * t272 - t171)) * t77) + (t185 * t136 + (t200 * t188 + (-t135 + t199) * t187) * t188) * t286 + (t186 * t135 + (t199 * t187 + (-t136 + t200) * t188) * t187) * t285 + m(5) * (-t107 * t303 - t109 * t302 + t44 * t63) + t278) * qJD(3) + t304, t4 * qJD(1) + t6 * qJD(3) + t304; (t203 - t289) * qJD(1) + t3 * qJD(3) + t224 * qJD(4), 0, t3 * qJD(1) + ((t63 * t66 + (-t107 * t187 - t109 * t188) * t150) * m(5) + t278) * qJD(3) + t7, qJD(1) * t224 + qJD(3) * t8 + t7;];
Cq = t5;
