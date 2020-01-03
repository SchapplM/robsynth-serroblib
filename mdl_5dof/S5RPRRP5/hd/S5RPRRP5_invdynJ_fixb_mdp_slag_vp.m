% Calculate vector of inverse dynamics joint torques for
% S5RPRRP5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:41:01
% EndTime: 2019-12-31 18:41:04
% DurationCPUTime: 1.27s
% Computational Cost: add. (1330->224), mult. (2219->252), div. (0->0), fcn. (1205->12), ass. (0->119)
t250 = qJDD(1) + qJDD(3);
t258 = sin(qJ(3));
t261 = cos(qJ(3));
t256 = cos(pkin(8));
t241 = pkin(1) * t256 + pkin(2);
t227 = t241 * qJD(1);
t255 = sin(pkin(8));
t339 = pkin(1) * t255;
t346 = qJD(3) * t227 + qJDD(1) * t339;
t300 = qJD(1) * t339;
t347 = -qJD(3) * t300 + t241 * qJDD(1);
t341 = -t347 * t258 - t346 * t261;
t348 = pkin(7) * t250 + qJD(2) * qJD(4) - t341;
t260 = cos(qJ(4));
t202 = t227 * t258 + t261 * t300;
t251 = qJD(1) + qJD(3);
t198 = pkin(7) * t251 + t202;
t257 = sin(qJ(4));
t329 = t198 * t257;
t190 = qJD(2) * t260 - t329;
t345 = qJD(5) - t190;
t297 = t257 * qJDD(2) + t348 * t260;
t302 = qJDD(4) * qJ(5);
t176 = t302 + (qJD(5) - t329) * qJD(4) + t297;
t305 = qJD(4) * t260;
t286 = -t260 * qJDD(2) + t198 * t305 + t348 * t257;
t332 = (qJDD(4) * pkin(4));
t342 = qJDD(5) - t332;
t177 = t286 + t342;
t344 = t176 * t260 + t177 * t257;
t187 = -qJD(4) * pkin(4) + t345;
t191 = qJD(2) * t257 + t198 * t260;
t188 = qJD(4) * qJ(5) + t191;
t336 = pkin(4) * t260;
t279 = qJ(5) * t257 + t336;
t215 = -pkin(3) - t279;
t252 = qJ(1) + pkin(8);
t248 = qJ(3) + t252;
t239 = sin(t248);
t240 = cos(t248);
t343 = -t240 * pkin(7) - t215 * t239;
t313 = g(1) * t240 + g(2) * t239;
t288 = t241 * t261 - t258 * t339;
t312 = t258 * t241 + t261 * t339;
t276 = t187 * t257 + t188 * t260;
t306 = qJD(4) * t257;
t340 = t187 * t305 - t188 * t306 + t344;
t338 = pkin(3) * t250;
t337 = pkin(3) * t251;
t263 = qJD(4) ^ 2;
t335 = pkin(7) * t263;
t235 = g(1) * t239;
t334 = g(2) * t240;
t333 = pkin(7) * qJDD(4);
t201 = t227 * t261 - t258 * t300;
t328 = t201 * t251;
t327 = t202 * t251;
t204 = t288 * qJD(3);
t326 = t204 * t251;
t205 = t312 * qJD(3);
t325 = t205 * t251;
t207 = pkin(7) + t312;
t324 = t207 * t263;
t323 = t215 * t250;
t322 = t215 * t251;
t321 = t239 * t257;
t320 = t240 * t257;
t318 = t250 * t257;
t317 = t251 * t257;
t316 = t257 * t260;
t304 = qJD(5) * t257;
t208 = pkin(4) * t306 - qJ(5) * t305 - t304;
t315 = t208 - t202;
t314 = g(1) * t321 - g(2) * t320;
t253 = t257 ^ 2;
t254 = t260 ^ 2;
t311 = t253 - t254;
t310 = t253 + t254;
t307 = qJD(4) * t251;
t301 = qJDD(4) * t207;
t249 = t251 ^ 2;
t298 = t249 * t316;
t222 = t260 * t235;
t296 = t201 * t306 + t260 * t327 + t222;
t284 = -t346 * t258 + t347 * t261;
t183 = -t284 - t338;
t294 = -t183 - t334;
t293 = t310 * t250;
t291 = t190 + t329;
t199 = t215 - t288;
t290 = t199 * t251 - t204;
t197 = -t201 - t337;
t285 = t183 * t257 + t197 * t305 - t314;
t283 = t239 * pkin(7) + qJ(5) * t320 + (pkin(3) + t336) * t240;
t282 = t335 - t338;
t259 = sin(qJ(1));
t262 = cos(qJ(1));
t280 = g(1) * t259 - g(2) * t262;
t278 = pkin(4) * t257 - qJ(5) * t260;
t216 = qJDD(4) * t257 + t260 * t263;
t217 = qJDD(4) * t260 - t257 * t263;
t277 = 0.2e1 * (t250 * t316 - t311 * t307) * MDP(9) + (t250 * t253 + 0.2e1 * t305 * t317) * MDP(8) + t216 * MDP(10) + t217 * MDP(11) + t250 * MDP(5);
t178 = (t278 * qJD(4) - t304) * t251 + t323 - t284;
t275 = -t178 - t323 - t335;
t272 = t235 + t284 - t334;
t206 = -pkin(3) - t288;
t271 = t206 * t250 + t324 + t325;
t270 = g(1) * t320 + g(2) * t321 - g(3) * t260 - t286;
t269 = -t301 + (t206 * t251 - t204) * qJD(4);
t189 = t205 + t208;
t268 = -t189 * t251 - t199 * t250 - t178 - t324;
t267 = -t313 + t340;
t266 = qJD(4) * t191 + t270;
t265 = t313 + t341;
t209 = t278 * t251;
t192 = t197 * t306;
t186 = -t201 + t322;
t184 = t186 * t306;
t1 = [qJDD(1) * MDP(1) + t280 * MDP(2) + (g(1) * t262 + g(2) * t259) * MDP(3) + (t280 + (t255 ^ 2 + t256 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t288 * t250 + t272 - t325) * MDP(6) + (-t312 * t250 + t265 - t326) * MDP(7) + (t192 + t222 + t269 * t257 + (-t271 + t294) * t260) * MDP(13) + (t271 * t257 + t269 * t260 + t285) * MDP(14) + (t184 + t222 + (t290 * qJD(4) - t301) * t257 + (t268 - t334) * t260) * MDP(15) + (t207 * t293 + t310 * t326 + t267) * MDP(16) + ((t301 + (-t186 - t290) * qJD(4)) * t260 + t268 * t257 + t314) * MDP(17) + (t178 * t199 + t186 * t189 - g(1) * (-pkin(2) * sin(t252) - t259 * pkin(1) - t343) - g(2) * (pkin(2) * cos(t252) + t262 * pkin(1) + t283) + t276 * t204 + t340 * t207) * MDP(18) + t277; (qJDD(2) - g(3)) * MDP(4) + (t276 * qJD(4) + t176 * t257 - t177 * t260 - g(3)) * MDP(18) + (MDP(13) + MDP(15)) * t217 + (-MDP(14) + MDP(17)) * t216; (t272 + t327) * MDP(6) + (t265 + t328) * MDP(7) + (t192 + (-pkin(3) * t307 - t333) * t257 + (-t282 + t294) * t260 + t296) * MDP(13) + ((-t333 + (t201 - t337) * qJD(4)) * t260 + (t282 - t327) * t257 + t285) * MDP(14) + (t184 + (t215 * t307 - t333) * t257 + (-t208 * t251 + t275 - t334) * t260 + t296) * MDP(15) + (pkin(7) * t293 - t310 * t328 + t267) * MDP(16) + ((t333 + (-t186 - t201 - t322) * qJD(4)) * t260 + (-t315 * t251 + t275) * t257 + t314) * MDP(17) + (t178 * t215 - g(2) * t283 - t276 * t201 + t315 * t186 + ((t187 * t260 - t188 * t257) * qJD(4) + t344) * pkin(7) + t343 * g(1)) * MDP(18) + t277; -MDP(8) * t298 + t311 * MDP(9) * t249 + MDP(10) * t318 + qJDD(4) * MDP(12) + (-t197 * t317 + t266) * MDP(13) + (g(3) * t257 + t291 * qJD(4) + (-t197 * t251 + t313) * t260 - t297) * MDP(14) + ((2 * t332) - qJDD(5) + (-t186 * t257 + t209 * t260) * t251 + t266) * MDP(15) + (0.2e1 * t302 + (t209 * t251 - g(3)) * t257 + (t186 * t251 - t313) * t260 + (0.2e1 * qJD(5) - t291) * qJD(4) + t297) * MDP(17) + (-t177 * pkin(4) - g(3) * t279 + t176 * qJ(5) - t186 * t209 - t187 * t191 + t345 * t188 + t313 * t278) * MDP(18) + (t260 * MDP(11) - t278 * MDP(16)) * t250; (-qJDD(4) - t298) * MDP(15) + MDP(16) * t318 + (-t249 * t253 - t263) * MDP(17) + (-qJD(4) * t188 + t186 * t317 - t270 + t342) * MDP(18);];
tau = t1;
