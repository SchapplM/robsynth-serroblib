% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:08
% EndTime: 2019-12-31 17:22:12
% DurationCPUTime: 2.19s
% Computational Cost: add. (17524->214), mult. (11510->288), div. (0->0), fcn. (9848->8), ass. (0->145)
t200 = qJ(1) + qJ(2);
t199 = qJ(3) + t200;
t194 = sin(t199);
t195 = cos(t199);
t161 = -t194 * rSges(4,1) - t195 * rSges(4,2);
t196 = sin(t200);
t266 = pkin(2) * t196;
t149 = t161 - t266;
t162 = t195 * rSges(4,1) - t194 * rSges(4,2);
t197 = cos(t200);
t265 = pkin(2) * t197;
t150 = t162 + t265;
t232 = t162 * t149 - t150 * t161;
t203 = cos(qJ(4));
t258 = t203 * rSges(5,1);
t224 = pkin(3) + t258;
t201 = sin(qJ(4));
t244 = t194 * t201;
t225 = rSges(5,2) * t244 + t195 * rSges(5,3);
t121 = t195 * pkin(7) - t224 * t194 + t225;
t115 = t121 - t266;
t259 = t201 * rSges(5,2);
t177 = t195 * t259;
t122 = -t177 + t224 * t195 + (rSges(5,3) + pkin(7)) * t194;
t116 = t122 + t265;
t260 = t122 * t115 - t116 * t121;
t295 = m(5) / 0.2e1;
t296 = m(4) / 0.2e1;
t264 = sin(qJ(1)) * pkin(1);
t113 = t115 - t264;
t263 = cos(qJ(1)) * pkin(1);
t114 = t116 + t263;
t52 = -t122 * t113 + t114 * t121;
t147 = t149 - t264;
t148 = t150 + t263;
t93 = -t162 * t147 + t148 * t161;
t261 = (t260 + t52) * t295 + (t232 + t93) * t296;
t262 = (-t260 + t52) * t295 + (-t232 + t93) * t296;
t5 = t261 - t262;
t310 = t5 * qJD(1);
t198 = Icges(5,4) * t203;
t180 = -Icges(5,2) * t201 + t198;
t181 = Icges(5,1) * t201 + t198;
t309 = t180 + t181;
t269 = m(3) * (t263 * (-t196 * rSges(3,1) - t197 * rSges(3,2)) + (t197 * rSges(3,1) - t196 * rSges(3,2)) * t264);
t256 = Icges(5,4) * t201;
t179 = Icges(5,2) * t203 + t256;
t182 = Icges(5,1) * t203 - t256;
t220 = t309 * t203 / 0.2e1 + (-t179 / 0.2e1 + t182 / 0.2e1) * t201;
t183 = t201 * rSges(5,1) + t203 * rSges(5,2);
t157 = t183 * t194;
t104 = t113 * t157;
t158 = t183 * t195;
t63 = -t114 * t158 + t104;
t67 = t121 * t157 - t122 * t158;
t285 = m(5) * (t67 + t63);
t306 = t220 + t285 / 0.2e1;
t64 = t115 * t157 - t116 * t158;
t284 = m(5) * (t67 + t64);
t305 = t220 + t284 / 0.2e1;
t83 = -t150 * t147 + t148 * t149;
t303 = t194 ^ 2;
t302 = t195 ^ 2;
t301 = 0.4e1 * qJD(1);
t298 = 2 * qJD(3);
t292 = m(4) * t83;
t290 = m(4) * t93;
t289 = m(4) * t232;
t286 = m(5) * (t64 + t63);
t283 = m(5) * (t104 + (-t115 * t194 + (-t114 + t116) * t195) * t183);
t282 = m(5) * (t104 + (-t121 * t194 + (-t114 + t122) * t195) * t183);
t281 = m(5) * ((-t116 + t122) * t195 + (t115 - t121) * t194) * t183;
t49 = -t116 * t113 + t114 * t115;
t280 = m(5) * t49;
t278 = m(5) * t52;
t277 = m(5) * t260;
t275 = m(5) * t63;
t274 = m(5) * t64;
t273 = m(5) * t67;
t272 = -t194 / 0.2e1;
t271 = t194 / 0.2e1;
t270 = -t195 / 0.2e1;
t243 = t194 * t203;
t178 = Icges(5,5) * t203 - Icges(5,6) * t201;
t241 = t195 * t178;
t140 = Icges(5,4) * t243 - Icges(5,2) * t244 - Icges(5,6) * t195;
t240 = t201 * t140;
t141 = Icges(5,6) * t194 + t180 * t195;
t239 = t201 * t141;
t175 = Icges(5,4) * t244;
t142 = Icges(5,1) * t243 - Icges(5,5) * t195 - t175;
t236 = t203 * t142;
t143 = Icges(5,5) * t194 + t182 * t195;
t235 = t203 * t143;
t138 = Icges(5,5) * t243 - Icges(5,6) * t244 - Icges(5,3) * t195;
t231 = -t194 * t138 - t195 * t236;
t139 = Icges(5,3) * t194 + t241;
t230 = t194 * t139 + t195 * t235;
t229 = t181 * t194 + t140;
t228 = -t181 * t195 - t141;
t227 = -Icges(5,2) * t243 + t142 - t175;
t226 = -t179 * t195 + t143;
t129 = t194 * t235;
t222 = t195 * t139 - t129;
t221 = -t138 + t239;
t219 = t286 / 0.2e1 + t220;
t215 = Icges(5,5) * t201 + Icges(5,6) * t203;
t211 = (-t157 * t195 + t158 * t194) * t183;
t75 = -t195 * t240 - t231;
t76 = -t195 * t239 + t230;
t17 = (t221 * t195 - t230 + t76) * t195 + (t221 * t194 + t222 + t75) * t194;
t74 = -t194 * t239 - t222;
t18 = (t74 - t129 + (t139 + t240) * t195 + t231) * t195 + t230 * t194;
t206 = (-t179 + t182) * t203 - t309 * t201;
t46 = t74 * t194 - (-(-t236 + t240) * t194 - t195 * t138) * t195;
t47 = t76 * t194 - t75 * t195;
t210 = t195 * t18 / 0.2e1 + (t17 + t46) * t272 + (t194 * t178 + t206 * t195 + t228 * t201 + t226 * t203) * t271 + (t206 * t194 - t229 * t201 + t227 * t203 - t241 + t47) * t270;
t209 = -t220 + (t271 + t272) * (t203 * t140 + t201 * t142);
t208 = t227 * t201 + t229 * t203;
t207 = -t226 * t201 + t228 * t203;
t185 = t258 - t259;
t152 = t195 * t215;
t151 = t215 * t194;
t53 = t220 + t273;
t50 = t220 + t274;
t48 = t220 + t275;
t40 = t281 / 0.2e1;
t38 = t282 / 0.2e1;
t35 = -t277 - t289;
t34 = t283 / 0.2e1;
t30 = t278 + t290;
t21 = t269 + t280 + t292;
t14 = -t281 / 0.2e1 + t305;
t13 = t40 + t305;
t12 = -t282 / 0.2e1 + t306;
t11 = t38 + t306;
t10 = -t283 / 0.2e1 + t219;
t9 = t34 + t219;
t8 = t40 - t284 / 0.2e1 + t209;
t7 = t38 - t285 / 0.2e1 + t209;
t6 = t34 - t286 / 0.2e1 + t209;
t3 = t261 + t262;
t2 = (t47 / 0.2e1 - t18 / 0.2e1) * t195 + (t17 / 0.2e1 + t46 / 0.2e1) * t194;
t1 = t2 * qJD(4);
t4 = [t21 * qJD(2) + t30 * qJD(3) + t48 * qJD(4), t21 * qJD(1) + t3 * qJD(3) + t9 * qJD(4) + 0.2e1 * (t269 / 0.2e1 + t49 * t295 + t83 * t296) * qJD(2), t30 * qJD(1) + t3 * qJD(2) + t11 * qJD(4) + (t52 * t295 + t93 * t296) * t298, t48 * qJD(1) + t9 * qJD(2) + t11 * qJD(3) + (((-t113 * t195 - t114 * t194) * t185 + t211) * m(5) + t210) * qJD(4); -t5 * qJD(3) + t10 * qJD(4) + (-t269 / 0.4e1 - t292 / 0.4e1 - t280 / 0.4e1) * t301, t35 * qJD(3) + t50 * qJD(4), -t310 + t35 * qJD(2) + t13 * qJD(4) + (-t232 * t296 - t260 * t295) * t298, t10 * qJD(1) + t50 * qJD(2) + t13 * qJD(3) + (((-t115 * t195 - t116 * t194) * t185 + t211) * m(5) + t210) * qJD(4); t5 * qJD(2) + t12 * qJD(4) + (-t290 / 0.4e1 - t278 / 0.4e1) * t301, t310 + t14 * qJD(4) + 0.4e1 * (t289 / 0.4e1 + t277 / 0.4e1) * qJD(2), t53 * qJD(4), t12 * qJD(1) + t14 * qJD(2) + t53 * qJD(3) + (((-t121 * t195 - t122 * t194) * t185 + t211) * m(5) + t210) * qJD(4); (t209 - t275) * qJD(1) + t6 * qJD(2) + t7 * qJD(3) + t1, t6 * qJD(1) + (t209 - t274) * qJD(2) + t8 * qJD(3) + t1, t7 * qJD(1) + t8 * qJD(2) + (t209 - t273) * qJD(3) + t1, (m(5) * ((t194 * (rSges(5,1) * t243 - t225) + t195 * (t194 * rSges(5,3) + t195 * t258 - t177)) * (-t194 * t157 - t195 * t158) + (t302 + t303) * t185 * t183) + (-t303 * t152 + (t208 * t195 + (t151 + t207) * t194) * t195) * t271 + (-t302 * t151 + (t207 * t194 + (t152 + t208) * t195) * t194) * t270) * qJD(4) + (qJD(1) + qJD(2) + qJD(3)) * t2;];
Cq = t4;
