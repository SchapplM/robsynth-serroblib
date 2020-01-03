% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:31
% EndTime: 2019-12-31 16:32:34
% DurationCPUTime: 2.24s
% Computational Cost: add. (13307->242), mult. (10393->342), div. (0->0), fcn. (9640->6), ass. (0->161)
t191 = qJ(3) + qJ(4);
t188 = cos(t191);
t181 = Icges(5,4) * t188;
t187 = sin(t191);
t147 = -Icges(5,2) * t187 + t181;
t148 = Icges(5,1) * t187 + t181;
t303 = t147 + t148;
t190 = pkin(7) + qJ(2);
t185 = sin(t190);
t280 = t185 / 0.2e1;
t186 = cos(t190);
t279 = -t186 / 0.2e1;
t302 = t186 / 0.2e1;
t183 = t185 ^ 2;
t184 = t186 ^ 2;
t223 = t183 + t184;
t192 = sin(qJ(3));
t193 = cos(qJ(3));
t163 = rSges(4,1) * t192 + rSges(4,2) * t193;
t141 = t163 * t185;
t142 = t163 * t186;
t257 = Icges(4,4) * t192;
t159 = Icges(4,2) * t193 + t257;
t162 = Icges(4,1) * t193 - t257;
t150 = rSges(5,1) * t187 + rSges(5,2) * t188;
t277 = pkin(3) * t192;
t200 = t150 + t277;
t297 = t200 * t186;
t298 = t200 * t185;
t249 = t185 * t188;
t250 = t185 * t187;
t105 = rSges(5,1) * t249 - rSges(5,2) * t250 - rSges(5,3) * t186;
t276 = pkin(3) * t193;
t182 = pkin(2) + t276;
t291 = -pkin(6) - pkin(5);
t225 = -t182 * t185 - t186 * t291;
t76 = -t105 + t225;
t170 = t185 * t291;
t244 = t186 * t187;
t218 = -rSges(5,2) * t244 + rSges(5,3) * t185;
t267 = rSges(5,1) * t188;
t77 = -t170 + (t182 + t267) * t186 + t218;
t180 = t186 * pkin(5);
t268 = rSges(4,1) * t193;
t221 = pkin(2) + t268;
t248 = t185 * t192;
t224 = rSges(4,2) * t248 + rSges(4,3) * t186;
t81 = -t185 * t221 + t180 + t224;
t260 = t192 * rSges(4,2);
t169 = t186 * t260;
t82 = -t169 + t221 * t186 + (rSges(4,3) + pkin(5)) * t185;
t300 = -(t162 / 0.2e1 - t159 / 0.2e1) * t192 - m(5) * (-t297 * t77 + t298 * t76) - m(4) * (t141 * t81 - t142 * t82);
t151 = -rSges(5,2) * t187 + t267;
t207 = Icges(5,5) * t187 + Icges(5,6) * t188;
t119 = t207 * t185;
t120 = t186 * t207;
t256 = Icges(5,4) * t187;
t149 = Icges(5,1) * t188 - t256;
t104 = Icges(5,5) * t185 + t149 * t186;
t146 = Icges(5,2) * t188 + t256;
t230 = -t146 * t186 + t104;
t102 = Icges(5,6) * t185 + t147 * t186;
t232 = -t148 * t186 - t102;
t197 = -t187 * t230 + t188 * t232;
t154 = Icges(5,4) * t250;
t103 = Icges(5,1) * t249 - Icges(5,5) * t186 - t154;
t231 = -Icges(5,2) * t249 + t103 - t154;
t101 = Icges(5,4) * t249 - Icges(5,2) * t250 - Icges(5,6) * t186;
t233 = t148 * t185 + t101;
t198 = t187 * t231 + t188 * t233;
t274 = (-t183 * t120 + (t198 * t186 + (t119 + t197) * t185) * t186) * t280 + (-t184 * t119 + (t197 * t185 + (t120 + t198) * t186) * t185) * t279;
t243 = t186 * t188;
t66 = t185 * t105 + t186 * (rSges(5,1) * t243 + t218);
t43 = -t185 * (pkin(2) * t185 - t180 + t225) + (-t185 * pkin(5) - t170 + (-pkin(2) + t182) * t186) * t186 + t66;
t125 = t150 * t185;
t126 = t150 * t186;
t73 = -t125 * t185 - t126 * t186;
t6 = t274 + m(5) * (t43 * t73 + (t185 * t298 + t186 * t297) * t151);
t299 = t6 * qJD(4);
t189 = Icges(4,4) * t193;
t160 = -Icges(4,2) * t192 + t189;
t161 = Icges(4,1) * t192 + t189;
t213 = t303 * t188 / 0.2e1 + (-t146 / 0.2e1 + t149 / 0.2e1) * t187;
t145 = Icges(5,5) * t188 - Icges(5,6) * t187;
t252 = t145 * t186;
t100 = Icges(5,3) * t185 + t252;
t99 = Icges(5,5) * t249 - Icges(5,6) * t250 - Icges(5,3) * t186;
t216 = t102 * t187 - t99;
t78 = t104 * t249;
t217 = t186 * t100 - t78;
t253 = t101 * t187;
t272 = t100 * t185 + t104 * t243;
t273 = -t103 * t243 - t185 * t99;
t53 = -t102 * t250 - t217;
t54 = -t101 * t244 - t273;
t55 = -t102 * t244 + t272;
t220 = ((t53 - t78 + (t100 + t253) * t186 + t273) * t186 + t272 * t185) * t279 + (t185 * t55 - t186 * t54) * t302 + ((t216 * t185 + t217 + t53 + t54) * t185 + ((-t103 * t188 + t253) * t185 - t272 + t55 + (t99 + t216) * t186) * t186) * t280;
t294 = 2 * qJD(3);
t293 = m(4) / 0.2e1;
t292 = m(5) / 0.2e1;
t203 = (-t185 * t77 - t186 * t76) * t151;
t286 = m(5) * (-t125 * t297 + t126 * t298 + t203);
t285 = m(5) * (t203 + (t185 * t297 - t186 * t298) * t150);
t283 = m(5) * (t125 * t76 - t126 * t77);
t282 = m(5) * t73;
t281 = -t185 / 0.2e1;
t247 = t185 * t193;
t112 = Icges(4,5) * t247 - Icges(4,6) * t248 - Icges(4,3) * t186;
t167 = Icges(4,4) * t248;
t116 = Icges(4,1) * t247 - Icges(4,5) * t186 - t167;
t236 = t193 * t116;
t270 = -t185 * t112 - t186 * t236;
t209 = Icges(4,5) * t193 - Icges(4,6) * t192;
t113 = Icges(4,3) * t185 + t186 * t209;
t117 = Icges(4,5) * t185 + t162 * t186;
t235 = t193 * t117;
t269 = t113 * t185 + t186 * t235;
t114 = Icges(4,4) * t247 - Icges(4,2) * t248 - Icges(4,6) * t186;
t238 = t192 * t114;
t115 = Icges(4,6) * t185 + t160 * t186;
t237 = t192 * t115;
t229 = t161 * t185 + t114;
t228 = -t161 * t186 - t115;
t227 = -Icges(4,2) * t247 + t116 - t167;
t226 = -t159 * t186 + t117;
t219 = -t151 - t276;
t89 = t185 * t235;
t215 = t186 * t113 - t89;
t214 = -t112 + t237;
t208 = -Icges(4,5) * t192 - Icges(4,6) * t193;
t194 = (-t146 + t149) * t188 - t303 * t187;
t201 = -t220 + (t145 * t185 + t186 * t194 + t187 * t232 + t188 * t230) * t280 + (t185 * t194 - t187 * t233 + t188 * t231 - t252) * t279;
t199 = -t213 + (t280 + t281) * (t188 * t101 + t103 * t187);
t196 = t192 * t227 + t193 * t229;
t195 = -t192 * t226 + t193 * t228;
t164 = -t260 + t268;
t136 = t208 * t186;
t135 = t208 * t185;
t109 = t219 * t186;
t107 = t219 * t185;
t75 = -t141 * t185 - t142 * t186;
t70 = qJD(4) * t282;
t63 = -t223 * t277 + t73;
t59 = -t186 * t237 + t269;
t58 = -t186 * t238 - t270;
t57 = -t185 * t237 - t215;
t38 = t185 * t59 - t186 * t58;
t37 = t185 * t57 - t186 * (-(-t236 + t238) * t185 - t186 * t112);
t34 = t213 + t283;
t31 = t285 / 0.2e1;
t30 = t286 / 0.2e1;
t19 = (t161 / 0.2e1 + t160 / 0.2e1) * t193 + t213 - t300;
t16 = (t57 - t89 + (t113 + t238) * t186 + t270) * t186 + t269 * t185;
t15 = (t186 * t214 - t269 + t59) * t186 + (t185 * t214 + t215 + t58) * t185;
t8 = m(5) * (t150 * t151 * t223 + t66 * t73) + t274;
t7 = t8 * qJD(4);
t4 = t30 - t285 / 0.2e1 + t220;
t3 = t31 - t286 / 0.2e1 + t220;
t2 = t30 + t31 + t201;
t1 = (t38 / 0.2e1 - t16 / 0.2e1) * t186 + (t15 / 0.2e1 + t37 / 0.2e1) * t185 + t220;
t5 = [0, 0, (t292 * t63 + t293 * t75) * t294 + t70, qJD(3) * t282 + t70; 0, qJD(3) * t19 + qJD(4) * t34, t19 * qJD(2) + t2 * qJD(4) + (((-t185 * t82 - t186 * t81) * t164 + (-t141 * t186 + t142 * t185) * t163) * t293 + (t107 * t77 + t109 * t76) * t292) * t294 + (t16 * t302 + (t192 * t228 + t193 * t226) * t280 + t201 + (t184 / 0.2e1 + t183 / 0.2e1) * t209 + (t15 + t37) * t281 + (-t192 * t229 + t193 * t227 + t38) * t279) * qJD(3), t34 * qJD(2) + t2 * qJD(3) + (t201 + (t203 + (-t125 * t186 + t126 * t185) * t150) * m(5)) * qJD(4); 0, t1 * qJD(3) + t4 * qJD(4) + (t199 - (t161 + t160) * t193 / 0.2e1 + t300) * qJD(2), t1 * qJD(2) + (m(4) * (t223 * t164 * t163 + (t185 * (rSges(4,1) * t247 - t224) + t186 * (t185 * rSges(4,3) + t186 * t268 - t169)) * t75) + (t183 * t136 + (t196 * t186 + (-t135 + t195) * t185) * t186) * t280 + (t184 * t135 + (t195 * t185 + (-t136 + t196) * t186) * t185) * t279 + m(5) * (-t107 * t298 - t109 * t297 + t43 * t63) + t274) * qJD(3) + t299, t4 * qJD(2) + t6 * qJD(3) + t299; 0, (t199 - t283) * qJD(2) + t3 * qJD(3) + t220 * qJD(4), t3 * qJD(2) + ((t63 * t66 + (-t107 * t185 - t109 * t186) * t150) * m(5) + t274) * qJD(3) + t7, qJD(2) * t220 + qJD(3) * t8 + t7;];
Cq = t5;
