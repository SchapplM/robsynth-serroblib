% Calculate Coriolis joint torque vector for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:15:01
% EndTime: 2021-01-15 19:15:17
% DurationCPUTime: 3.95s
% Computational Cost: add. (2561->297), mult. (6656->382), div. (0->0), fcn. (4830->6), ass. (0->115)
t274 = sin(qJ(4));
t306 = qJD(4) * t274;
t272 = sin(pkin(8));
t273 = cos(pkin(8));
t275 = sin(qJ(3));
t277 = cos(qJ(3));
t345 = -t272 * t275 + t277 * t273;
t349 = t345 * qJD(1);
t337 = pkin(6) + qJ(2);
t258 = t337 * t272;
t255 = qJD(1) * t258;
t259 = t337 * t273;
t256 = qJD(1) * t259;
t227 = -t275 * t255 + t277 * t256;
t348 = qJD(3) * t227;
t254 = t272 * t277 + t273 * t275;
t249 = t254 * qJD(1);
t276 = cos(qJ(4));
t305 = qJD(4) * t276;
t307 = qJD(3) * t274;
t201 = t249 * t305 + (qJD(4) + t349) * t307;
t347 = (MDP(6) * qJ(2) + MDP(5)) * (t272 ^ 2 + t273 ^ 2);
t235 = t249 * t276 + t307;
t250 = t345 * qJD(3);
t282 = qJD(1) * t250;
t304 = t276 * qJD(3);
t200 = -qJD(4) * t304 + t249 * t306 - t276 * t282;
t266 = -pkin(2) * t273 - pkin(1);
t257 = qJD(1) * t266 + qJD(2);
t204 = -pkin(3) * t349 - pkin(7) * t249 + t257;
t222 = qJD(3) * pkin(7) + t227;
t190 = t204 * t274 + t222 * t276;
t284 = t345 * qJD(2);
t342 = -t255 * t277 - t275 * t256;
t195 = qJD(1) * t284 + qJD(3) * t342;
t251 = t254 * qJD(3);
t239 = qJD(1) * t251;
t210 = t239 * pkin(3) - pkin(7) * t282;
t209 = t276 * t210;
t281 = -qJD(4) * t190 - t195 * t274 + t209;
t279 = qJ(5) * t200 + t281;
t338 = pkin(4) * t239;
t176 = -qJD(5) * t235 + t279 + t338;
t233 = t249 * t274 - t304;
t184 = -qJ(5) * t233 + t190;
t242 = qJD(4) - t349;
t334 = t184 * t242;
t346 = t176 + t334;
t230 = t258 * t277 + t275 * t259;
t340 = t235 ^ 2;
t339 = pkin(4) * t233;
t336 = -qJ(5) - pkin(7);
t333 = t200 * t274;
t332 = t233 * t242;
t331 = t233 * t349;
t330 = t233 * t249;
t329 = t235 * t242;
t328 = t235 * t249;
t327 = t235 * t274;
t326 = t349 * t274;
t325 = t349 * t276;
t324 = t254 * t274;
t323 = t254 * t276;
t317 = t274 * t239;
t231 = -t258 * t275 + t259 * t277;
t228 = t276 * t231;
t238 = t276 * t239;
t189 = t276 * t204 - t222 * t274;
t183 = -qJ(5) * t235 + t189;
t180 = pkin(4) * t242 + t183;
t315 = t180 - t183;
t223 = pkin(3) * t249 - pkin(7) * t349;
t215 = t276 * t223;
t298 = qJD(4) * t336;
t314 = pkin(4) * t249 - qJ(5) * t325 - t276 * t298 + t215 + (qJD(5) - t342) * t274;
t311 = t274 * t223 + t276 * t342;
t313 = -qJ(5) * t326 - qJD(5) * t276 - t274 * t298 + t311;
t312 = -t274 * t201 - t233 * t305;
t225 = -pkin(3) * t345 - pkin(7) * t254 + t266;
t310 = t274 * t225 + t228;
t302 = qJD(1) * qJD(2);
t205 = -qJD(3) * t230 + t284;
t224 = pkin(3) * t251 - pkin(7) * t250;
t301 = t276 * t205 + t274 * t224 + t225 * t305;
t300 = t254 * t305;
t297 = -qJD(5) - t339;
t296 = t242 * t276;
t295 = t276 * t195 + t204 * t305 + t274 * t210 - t222 * t306;
t196 = t254 * t302 + t348;
t286 = qJ(5) * t201 - t295;
t177 = -qJD(5) * t233 - t286;
t293 = -t180 * t242 + t177;
t292 = -qJ(5) * t250 - qJD(5) * t254;
t290 = t238 + (-t306 + t326) * t242;
t221 = -qJD(3) * pkin(3) - t342;
t289 = t274 * t250 + t300;
t288 = t276 * t250 - t254 * t306;
t186 = pkin(4) * t201 + t196;
t287 = -pkin(7) * t239 + t221 * t242;
t206 = qJD(2) * t254 + qJD(3) * t231;
t268 = -pkin(4) * t276 - pkin(3);
t261 = t336 * t276;
t260 = t336 * t274;
t232 = t233 ^ 2;
t218 = t276 * t225;
t216 = t276 * t224;
t207 = pkin(4) * t324 + t230;
t197 = pkin(4) * t326 + t227;
t193 = t221 - t297;
t192 = pkin(4) * t289 + t206;
t191 = -qJ(5) * t324 + t310;
t187 = -pkin(4) * t345 - qJ(5) * t323 - t231 * t274 + t218;
t179 = -qJ(5) * t300 + (-qJD(4) * t231 + t292) * t274 + t301;
t178 = pkin(4) * t251 - t205 * t274 + t216 + t292 * t276 + (-t228 + (qJ(5) * t254 - t225) * t274) * qJD(4);
t1 = [(t249 * t250 + t254 * t282) * MDP(7) + (-t254 * t239 - t249 * t251 + t250 * t349 + t282 * t345) * MDP(8) + (t239 * t266 + t251 * t257) * MDP(12) + t257 * t250 * MDP(13) + (-t200 * t323 + t235 * t288) * MDP(14) + ((-t233 * t276 - t327) * t250 + (t333 - t201 * t276 + (t233 * t274 - t235 * t276) * qJD(4)) * t254) * MDP(15) + (t200 * t345 + t235 * t251 + t238 * t254 + t242 * t288) * MDP(16) + (t201 * t345 - t233 * t251 - t242 * t289 - t254 * t317) * MDP(17) + (-t239 * t345 + t242 * t251) * MDP(18) + ((-t231 * t305 + t216) * t242 + t218 * t239 - (-t222 * t305 + t209) * t345 + t189 * t251 + t206 * t233 + t230 * t201 + t221 * t300 + ((-qJD(4) * t225 - t205) * t242 - t231 * t239 - (-qJD(4) * t204 - t195) * t345 + t196 * t254 + t221 * t250) * t274) * MDP(19) + (-(-t231 * t306 + t301) * t242 - t310 * t239 + t295 * t345 - t190 * t251 + t206 * t235 - t230 * t200 + t196 * t323 + t288 * t221) * MDP(20) + (-t176 * t345 + t178 * t242 + t180 * t251 + t186 * t324 + t187 * t239 + t192 * t233 + t193 * t289 + t201 * t207) * MDP(21) + (t177 * t345 - t179 * t242 - t184 * t251 + t186 * t323 - t191 * t239 + t192 * t235 + t193 * t288 - t200 * t207) * MDP(22) + (-t178 * t235 - t179 * t233 + t187 * t200 - t191 * t201 + (-t180 * t276 - t184 * t274) * t250 + (-t176 * t276 - t177 * t274 + (t180 * t274 - t184 * t276) * qJD(4)) * t254) * MDP(23) + (t176 * t187 + t177 * t191 + t178 * t180 + t179 * t184 + t186 * t207 + t192 * t193) * MDP(24) + 0.2e1 * t302 * t347 + (t250 * MDP(9) - t251 * MDP(10) - t206 * MDP(12) + (t266 * t349 - t205) * MDP(13)) * qJD(3); t312 * MDP(23) - t193 * t249 * MDP(24) + (MDP(20) + MDP(22)) * (-t242 ^ 2 * t276 - t317 - t328) + (MDP(19) + MDP(21)) * (t290 - t330) + ((t200 + t331) * MDP(23) + t346 * MDP(24)) * t276 + (MDP(23) * t329 + MDP(24) * t293) * t274 + (t249 * MDP(12) + t349 * MDP(13) + (MDP(12) * t254 + MDP(13) * t345) * qJD(1)) * qJD(3) - qJD(1) ^ 2 * t347; (-t196 + t348) * MDP(12) + (t235 * t296 - t333) * MDP(14) + ((-t200 + t331) * t276 - t242 * t327 + t312) * MDP(15) + (t242 * t296 + t317 - t328) * MDP(16) + (t290 + t330) * MDP(17) + (-pkin(3) * t201 - t196 * t276 - t227 * t233 + (-pkin(7) * t305 - t215) * t242 + (t242 * t342 + t287) * t274) * MDP(19) + (pkin(3) * t200 + t196 * t274 - t227 * t235 + (pkin(7) * t306 + t311) * t242 + t287 * t276) * MDP(20) + (-t186 * t276 - t197 * t233 + t201 * t268 + t239 * t260 - t314 * t242 + (t193 + t339) * t306) * MDP(21) + (-t193 * t325 + t186 * t274 - t197 * t235 - t200 * t268 + t239 * t261 + t313 * t242 + (pkin(4) * t327 + t193 * t276) * qJD(4)) * MDP(22) + (t200 * t260 + t201 * t261 + t313 * t233 + t314 * t235 - t346 * t274 + t293 * t276) * MDP(23) + (t176 * t260 - t177 * t261 + t186 * t268 + (pkin(4) * t306 - t197) * t193 - t313 * t184 - t314 * t180) * MDP(24) + (-t257 * MDP(12) - t242 * MDP(18) - t189 * MDP(19) + t190 * MDP(20) - t180 * MDP(21) + t184 * MDP(22) + MDP(8) * t249) * t249 + ((-qJD(2) - t257) * MDP(13) - t193 * t274 * MDP(21) - t249 * MDP(7) - MDP(8) * t349) * t349; t235 * t233 * MDP(14) + (-t232 + t340) * MDP(15) + (-t200 + t332) * MDP(16) + (-t201 + t329) * MDP(17) + t239 * MDP(18) + (t190 * t242 - t221 * t235 + t281) * MDP(19) + (t189 * t242 + t221 * t233 - t295) * MDP(20) + (0.2e1 * t338 + t334 + (-t193 + t297) * t235 + t279) * MDP(21) + (-pkin(4) * t340 + t183 * t242 + (qJD(5) + t193) * t233 + t286) * MDP(22) + (pkin(4) * t200 - t233 * t315) * MDP(23) + (t315 * t184 + (-t193 * t235 + t176) * pkin(4)) * MDP(24); (-t200 - t332) * MDP(22) + (-t232 - t340) * MDP(23) + (t180 * t235 + t184 * t233 + t186) * MDP(24) + (t329 + t201) * MDP(21);];
tauc = t1;
