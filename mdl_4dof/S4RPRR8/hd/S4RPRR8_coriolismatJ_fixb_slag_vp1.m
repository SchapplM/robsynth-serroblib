% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:06
% EndTime: 2019-12-31 16:55:10
% DurationCPUTime: 2.40s
% Computational Cost: add. (7805->260), mult. (10814->375), div. (0->0), fcn. (9964->6), ass. (0->168)
t203 = cos(qJ(1));
t201 = sin(qJ(1));
t199 = qJ(3) + qJ(4);
t186 = cos(t199);
t185 = sin(t199);
t254 = Icges(5,4) * t185;
t218 = Icges(5,2) * t186 + t254;
t119 = Icges(5,6) * t203 + t218 * t201;
t244 = t186 * t201;
t180 = Icges(5,4) * t244;
t247 = t185 * t201;
t121 = Icges(5,1) * t247 + Icges(5,5) * t203 + t180;
t212 = -t119 * t186 - t121 * t185;
t313 = t212 * t203;
t222 = t185 * rSges(5,1) + t186 * rSges(5,2);
t200 = sin(qJ(3));
t266 = pkin(3) * t200;
t306 = t222 + t266;
t115 = t306 * t203;
t253 = Icges(5,4) * t186;
t166 = -Icges(5,2) * t185 + t253;
t220 = Icges(5,1) * t185 + t253;
t312 = t166 + t220;
t197 = t201 ^ 2;
t290 = m(4) / 0.2e1;
t289 = m(5) / 0.2e1;
t310 = -t201 / 0.2e1;
t269 = t201 / 0.2e1;
t267 = t203 / 0.2e1;
t214 = Icges(5,5) * t185 + Icges(5,6) * t186;
t309 = t214 * t201;
t308 = t214 * t203;
t251 = t222 * t201;
t250 = t222 * t203;
t198 = t203 ^ 2;
t231 = t197 + t198;
t93 = t231 * t222;
t202 = cos(qJ(3));
t255 = Icges(4,4) * t202;
t173 = -Icges(4,2) * t200 + t255;
t221 = Icges(4,1) * t200 + t255;
t307 = (t221 / 0.2e1 + t173 / 0.2e1) * t202;
t122 = -Icges(5,5) * t201 + t220 * t203;
t238 = t166 * t203 + t122;
t239 = -Icges(5,2) * t247 + t121 + t180;
t120 = -Icges(5,6) * t201 + t218 * t203;
t168 = Icges(5,1) * t186 - t254;
t240 = -t168 * t203 + t120;
t241 = -t168 * t201 + t119;
t305 = -(t240 * t201 - t241 * t203) * t185 + (t238 * t201 - t239 * t203) * t186;
t143 = -Icges(4,5) * t201 + t221 * t203;
t234 = t173 * t203 + t143;
t242 = t201 * t202;
t182 = Icges(4,4) * t242;
t243 = t200 * t201;
t142 = Icges(4,1) * t243 + Icges(4,5) * t203 + t182;
t235 = -Icges(4,2) * t243 + t142 + t182;
t256 = Icges(4,4) * t200;
t219 = Icges(4,2) * t202 + t256;
t141 = -Icges(4,6) * t201 + t219 * t203;
t175 = Icges(4,1) * t202 - t256;
t236 = -t175 * t203 + t141;
t140 = Icges(4,6) * t203 + t219 * t201;
t237 = -t175 * t201 + t140;
t304 = -(t236 * t201 - t237 * t203) * t200 + (t234 * t201 - t235 * t203) * t202;
t170 = rSges(5,1) * t186 - rSges(5,2) * t185;
t183 = pkin(3) * t242;
t114 = t170 * t201 + t183;
t265 = pkin(3) * t202;
t224 = (t170 + t265) * t203;
t215 = Icges(5,5) * t186 - Icges(5,6) * t185;
t145 = t201 * t215;
t146 = t215 * t203;
t264 = (-t197 * t146 + (t201 * t145 + t305) * t203) * t269 + (t198 * t145 + (-t203 * t146 - t305) * t201) * t267;
t194 = t201 * rSges(5,3);
t106 = t203 * (-t194 + t250);
t123 = rSges(5,3) * t203 + t251;
t286 = -pkin(6) - pkin(5);
t184 = t203 * t286;
t53 = -t203 * (-t201 * pkin(6) + t203 * t266) - t106 + (-pkin(3) * t243 + pkin(5) * t203 - t123 + t184) * t201;
t151 = rSges(5,1) * t244 - rSges(5,2) * t247;
t152 = t170 * t203;
t84 = -t151 * t201 - t203 * t152;
t6 = t264 + m(5) * (-t114 * t251 - t224 * t250 + t53 * t84);
t302 = t6 * qJD(4);
t301 = t201 * t224;
t300 = (t141 * t202 + t143 * t200) * t203;
t299 = (t120 * t186 + t122 * t185) * t203;
t111 = t151 + t183;
t73 = -t111 * t203 + t301;
t178 = rSges(4,1) * t202 - rSges(4,2) * t200;
t162 = t178 * t201;
t163 = t178 * t203;
t91 = -t162 * t203 + t201 * t163;
t263 = t73 * t289 + t91 * t290;
t297 = t312 * t185 - (-t218 + t168) * t186;
t227 = -t312 * t186 / 0.2e1 + (t218 / 0.2e1 - t168 / 0.2e1) * t185;
t118 = -Icges(5,3) * t201 + t308;
t54 = t203 * (Icges(5,3) * t203 + t309) + t119 * t244 + t121 * t247;
t55 = -t203 * t118 - t120 * t244 - t122 * t247;
t57 = -t201 * t118 + t299;
t5 = (t55 * t201 + t203 * t54) * t310 + ((t55 - t313) * t201 + (t57 - t299 + (t118 + t212) * t201 + t54) * t203) * t269 + (t197 * t118 + t57 * t201 + (t55 + (t118 - t212) * t203 + t313) * t203) * t267;
t292 = 4 * qJD(1);
t291 = 2 * qJD(3);
t83 = -t151 * t203 + t201 * t152;
t288 = t83 / 0.2e1;
t287 = pkin(1) + pkin(5);
t187 = t203 * qJ(2);
t285 = m(3) * ((rSges(3,3) * t203 + t187) * t203 + (rSges(3,3) + qJ(2)) * t197);
t223 = rSges(4,1) * t200 + rSges(4,2) * t202;
t204 = -t201 * rSges(4,3) + t223 * t203;
t94 = -t287 * t201 + t187 + t204;
t95 = (rSges(4,3) + t287) * t203 + (qJ(2) + t223) * t201;
t284 = m(4) * (t162 * t95 + t163 * t94);
t283 = m(4) * (t95 * t201 + t203 * t94);
t85 = t187 - t194 + (-pkin(1) + t286) * t201 + t115;
t86 = -t184 + (rSges(5,3) + pkin(1)) * t203 + (qJ(2) + t306) * t201;
t225 = t86 * t250 - t85 * t251;
t278 = m(5) * (t114 * t152 - t151 * t224 + t225);
t277 = m(5) * (t73 * t170 + t225);
t276 = m(5) * (t111 * t86 + t224 * t85);
t275 = m(5) * (t151 * t86 + t152 * t85);
t274 = m(5) * (t86 * t201 + t203 * t85);
t272 = m(5) * (t114 * t203 - t301);
t268 = -t203 / 0.2e1;
t262 = m(5) * qJD(4);
t230 = -m(5) * t83 / 0.2e1;
t216 = Icges(4,5) * t200 + Icges(4,6) * t202;
t138 = Icges(4,3) * t203 + t216 * t201;
t58 = t203 * t138 + t140 * t242 + t142 * t243;
t139 = -Icges(4,3) * t201 + t216 * t203;
t59 = -t203 * t139 - t141 * t242 - t143 * t243;
t228 = t231 * t223;
t217 = Icges(4,5) * t202 - Icges(4,6) * t200;
t113 = t306 * t201;
t213 = -t113 * t201 - t115 * t203;
t209 = -t140 * t202 - t142 * t200;
t207 = -t5 + (t238 * t185 + t240 * t186 + t203 * t297 - t309) * t269 + (-t239 * t185 - t241 * t186 - t201 * t297 - t308) * t267;
t205 = -t227 + (t267 + t268) * (t120 * t185 - t122 * t186);
t157 = t217 * t203;
t156 = t201 * t217;
t125 = t201 * t138;
t90 = t93 * t262;
t78 = t262 * t288;
t77 = -t201 * t123 - t106;
t72 = -t231 * t265 + t84;
t70 = t272 / 0.2e1;
t61 = -t201 * t139 + t300;
t60 = t209 * t203 + t125;
t41 = t61 * t201 + t203 * t60;
t40 = t59 * t201 + t203 * t58;
t30 = t227 + t275;
t29 = t277 / 0.2e1;
t28 = t278 / 0.2e1;
t26 = -t272 / 0.2e1 + t263;
t25 = t70 + t263;
t24 = t70 - t263;
t23 = t274 + t283 + t285;
t18 = -t307 + (-t175 / 0.2e1 + t219 / 0.2e1) * t200 + t284 + t276 + t227;
t16 = t197 * t139 + (-t125 + t59 + (t139 - t209) * t203) * t203;
t15 = (-t60 + t125 + t59) * t201 + (t61 - t300 + (t139 + t209) * t201 + t58) * t203;
t8 = m(5) * (-t170 * t93 + t77 * t84) + t264;
t7 = t8 * qJD(4);
t4 = t28 - t277 / 0.2e1 + t5;
t3 = t29 - t278 / 0.2e1 + t5;
t2 = t28 + t29 + t207;
t1 = (t41 / 0.2e1 + t16 / 0.2e1) * t203 + (t15 / 0.2e1 - t40 / 0.2e1) * t201 + t5;
t9 = [t23 * qJD(2) + t18 * qJD(3) + t30 * qJD(4), qJD(1) * t23 + qJD(3) * t25 + t78, t18 * qJD(1) + t25 * qJD(2) + t2 * qJD(4) + ((t91 * t178 - (t201 * t94 - t203 * t95) * t223) * t290 + (-t113 * t85 + t115 * t86 + (-t111 + t114) * t224) * t289) * t291 + (t15 * t310 + (-t235 * t200 - t237 * t202) * t267 + t207 - (t198 / 0.2e1 + t197 / 0.2e1) * t216 + (t234 * t200 + t236 * t202 + t40) * t269 + (t41 + t16) * t268) * qJD(3), t30 * qJD(1) + t2 * qJD(3) + t207 * qJD(4) + (qJD(2) * t288 + (t83 * t170 + t225) * qJD(4)) * m(5); t26 * qJD(3) + t78 + (-t285 / 0.4e1 - t283 / 0.4e1 - t274 / 0.4e1) * t292, 0, t26 * qJD(1) + (t213 * t289 - t228 * t290) * t291 - t90, -t90 + 0.2e1 * (t83 * qJD(1) / 0.4e1 - t93 * qJD(3) / 0.2e1) * m(5); (t205 + (-t219 + t175) * t200 / 0.2e1 + t307) * qJD(1) + t24 * qJD(2) + t1 * qJD(3) + t4 * qJD(4) + (-t284 / 0.4e1 - t276 / 0.4e1) * t292, t24 * qJD(1), t1 * qJD(1) + (m(4) * ((-t203 * t204 + (-t203 * rSges(4,3) - t223 * t201) * t201) * (-t201 * t162 - t163 * t203) - t178 * t228) + (t198 * t156 + (-t203 * t157 - t304) * t201) * t267 + (-t197 * t157 + (t201 * t156 + t304) * t203) * t269 + m(5) * (-t113 * t114 - t115 * t224 + t53 * t72) + t264) * qJD(3) + t302, t4 * qJD(1) + t6 * qJD(3) + t302; (t205 - t275) * qJD(1) + qJD(2) * t230 + t3 * qJD(3) + t5 * qJD(4), qJD(1) * t230, t3 * qJD(1) + ((t213 * t170 + t77 * t72) * m(5) + t264) * qJD(3) + t7, qJD(1) * t5 + qJD(3) * t8 + t7;];
Cq = t9;
