% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:12
% EndTime: 2019-12-31 19:13:17
% DurationCPUTime: 2.23s
% Computational Cost: add. (1679->249), mult. (3602->351), div. (0->0), fcn. (2410->6), ass. (0->116)
t249 = qJD(3) + qJD(4);
t257 = cos(qJ(5));
t255 = sin(qJ(4));
t258 = cos(qJ(4));
t259 = cos(qJ(3));
t300 = qJD(1) * t259;
t286 = t258 * t300;
t256 = sin(qJ(3));
t301 = qJD(1) * t256;
t215 = t255 * t301 - t286;
t254 = sin(qJ(5));
t311 = t215 * t254;
t204 = -t249 * t257 - t311;
t221 = t255 * t259 + t256 * t258;
t216 = t221 * qJD(1);
t329 = qJD(5) + t216;
t332 = t204 * t329;
t269 = t215 * t257 - t249 * t254;
t331 = t269 * t329;
t297 = qJD(4) * t255;
t299 = qJD(3) * t256;
t330 = -t255 * t299 - t256 * t297;
t328 = MDP(8) * (t256 ^ 2 - t259 ^ 2);
t327 = t216 * t249;
t260 = -pkin(1) - pkin(6);
t230 = qJD(1) * t260 + qJD(2);
t213 = -pkin(7) * t300 + t230 * t259;
t209 = qJD(3) * pkin(3) + t213;
t283 = pkin(7) * qJD(1) - t230;
t298 = qJD(3) * t259;
t211 = t283 * t298;
t326 = t258 * (qJD(4) * t209 - t211);
t325 = qJ(2) * MDP(6) + MDP(5);
t210 = t283 * t299;
t212 = -pkin(7) * t301 + t230 * t256;
t281 = t210 * t255 - t212 * t297;
t173 = t281 + t326;
t312 = t212 * t258;
t190 = t209 * t255 + t312;
t282 = -t210 * t258 - t211 * t255;
t174 = qJD(4) * t190 + t282;
t323 = pkin(7) - t260;
t224 = t323 * t256;
t225 = t323 * t259;
t202 = -t224 * t255 + t225 * t258;
t218 = t323 * t299;
t219 = qJD(3) * t225;
t180 = -qJD(4) * t202 + t218 * t255 - t219 * t258;
t313 = t212 * t255;
t189 = t209 * t258 - t313;
t186 = -pkin(4) * t249 - t189;
t226 = pkin(3) * t301 + qJD(1) * qJ(2);
t191 = pkin(4) * t216 + pkin(8) * t215 + t226;
t308 = t258 * t259;
t272 = t249 * t308;
t304 = t330 * qJD(1);
t197 = qJD(1) * t272 + t304;
t220 = t255 * t256 - t308;
t241 = pkin(3) * t256 + qJ(2);
t199 = pkin(4) * t221 + pkin(8) * t220 + t241;
t296 = qJD(4) * t258;
t200 = -t255 * t298 - t256 * t296 - t258 * t299 - t259 * t297;
t203 = -t224 * t258 - t225 * t255;
t324 = -(qJD(5) * t199 + t180) * t329 - t203 * t197 - (qJD(5) * t191 + t173) * t221 - t174 * t220 + t186 * t200;
t295 = qJD(5) * t254;
t294 = qJD(5) * t257;
t305 = t249 * t294 - t257 * t327;
t178 = t215 * t295 + t305;
t321 = t178 * t220;
t320 = t178 * t254;
t319 = t186 * t216;
t318 = t186 * t220;
t317 = t327 * t254;
t316 = t199 * t197;
t315 = t200 * t249;
t201 = t272 + t330;
t314 = t201 * t249;
t310 = t254 * t197;
t309 = t257 * t197;
t262 = qJD(1) ^ 2;
t307 = t259 * t262;
t261 = qJD(3) ^ 2;
t306 = t260 * t261;
t250 = qJD(1) * qJD(2);
t292 = qJD(1) * qJD(3);
t285 = t259 * t292;
t223 = pkin(3) * t285 + t250;
t231 = pkin(3) * t298 + qJD(2);
t291 = 0.2e1 * qJD(1);
t290 = pkin(3) * t300;
t284 = -pkin(3) * t249 - t209;
t280 = t257 * t329;
t198 = -pkin(4) * t215 + pkin(8) * t216;
t243 = pkin(3) * t255 + pkin(8);
t276 = qJD(5) * t243 + t198 + t290;
t192 = t213 * t255 + t312;
t275 = pkin(3) * t297 - t192;
t193 = t213 * t258 - t313;
t274 = -pkin(3) * t296 + t193;
t187 = pkin(8) * t249 + t190;
t171 = t187 * t257 + t191 * t254;
t273 = -t171 * t215 + t174 * t254 + t186 * t294;
t271 = -t197 * t243 + t319;
t270 = t187 * t254 - t191 * t257;
t268 = -t174 * t257 + t186 * t295 - t215 * t270;
t267 = t215 * t226 - t282;
t266 = t216 * t226 - t281;
t265 = t200 * t257 + t220 * t295;
t179 = -qJD(5) * t269 - t317;
t263 = ((t178 - t332) * t257 + (-t179 + t331) * t254) * MDP(22) + (-t269 * t280 + t320) * MDP(21) + (-t254 * t329 ^ 2 - t204 * t215 + t309) * MDP(24) + (-t215 * t269 + t280 * t329 + t310) * MDP(23) + (-t304 + (-t215 - t286) * t249) * MDP(17) + (t215 ^ 2 - t216 ^ 2) * MDP(15) + (-MDP(14) * t216 + MDP(25) * t329) * t215;
t244 = -pkin(3) * t258 - pkin(4);
t181 = qJD(4) * t203 - t218 * t258 - t219 * t255;
t177 = pkin(4) * t201 - pkin(8) * t200 + t231;
t176 = pkin(4) * t197 + pkin(8) * t327 + t223;
t175 = t257 * t176;
t1 = [0.2e1 * t292 * t328 - t261 * t259 * MDP(10) + qJ(2) * t298 * t291 * MDP(12) + (-t259 * t306 + (-qJ(2) * t299 + qJD(2) * t259) * t291) * MDP(13) + (-t200 * t215 + t220 * t327) * MDP(14) + (t197 * t220 - t200 * t216 + t201 * t215 + t221 * t327) * MDP(15) + MDP(16) * t315 - MDP(17) * t314 + (-t181 * t249 + t197 * t241 + t201 * t226 + t216 * t231 + t221 * t223) * MDP(19) + (-t180 * t249 + t200 * t226 - t215 * t231 - t220 * t223 - t241 * t327) * MDP(20) + (-t257 * t321 - t265 * t269) * MDP(21) + ((-t204 * t257 + t254 * t269) * t200 + (t320 + t179 * t257 + (-t204 * t254 - t257 * t269) * qJD(5)) * t220) * MDP(22) + (t178 * t221 - t201 * t269 - t220 * t309 + t265 * t329) * MDP(23) + (t220 * t310 - t179 * t221 - t201 * t204 + (-t200 * t254 + t220 * t294) * t329) * MDP(24) + (t197 * t221 + t201 * t329) * MDP(25) + (-t270 * t201 + t175 * t221 + t202 * t179 + t181 * t204 + (t177 * t329 + t316 + (-t187 * t221 - t203 * t329 - t318) * qJD(5)) * t257 + t324 * t254) * MDP(26) + (-t171 * t201 + t202 * t178 - t181 * t269 + (-(-qJD(5) * t203 + t177) * t329 - t316 - (-qJD(5) * t187 + t176) * t221 + qJD(5) * t318) * t254 + t324 * t257) * MDP(27) + 0.2e1 * t325 * t250 + (-0.2e1 * MDP(7) * t285 - t261 * MDP(9) + (qJD(2) * t291 - t306) * MDP(12)) * t256; (-qJD(1) * t216 + t315) * MDP(19) + (qJD(1) * t215 - t314) * MDP(20) + (t179 * t220 - t200 * t204 - t221 * t310) * MDP(26) + (t200 * t269 - t221 * t309 + t321) * MDP(27) - t325 * t262 + ((-qJD(1) * t257 - t201 * t254 - t221 * t294) * MDP(26) + (qJD(1) * t254 - t201 * t257 + t221 * t295) * MDP(27)) * t329 + (MDP(12) * t256 + MDP(13) * t259) * (-t261 - t262); (-t216 * t290 + t192 * t249 + (t255 * t284 - t312) * qJD(4) + t267) * MDP(19) + (t215 * t290 + t193 * t249 + (qJD(4) * t284 + t211) * t258 + t266) * MDP(20) + (t244 * t179 + t271 * t254 + t275 * t204 + (t254 * t274 - t257 * t276) * t329 + t268) * MDP(26) + (t244 * t178 + t271 * t257 - t275 * t269 + (t254 * t276 + t257 * t274) * t329 + t273) * MDP(27) + t263 - t262 * t328 + t256 * MDP(7) * t307 + (MDP(13) * t256 * t262 - MDP(12) * t307) * qJ(2); (t267 + (-qJD(4) + t249) * t190) * MDP(19) + (t189 * t249 + t266 - t326) * MDP(20) + (-pkin(4) * t179 - (-t189 * t254 + t198 * t257) * t329 - t190 * t204 + t254 * t319 + (-t294 * t329 - t310) * pkin(8) + t268) * MDP(26) + (-pkin(4) * t178 + (t189 * t257 + t198 * t254) * t329 + t190 * t269 + t257 * t319 + (t295 * t329 - t309) * pkin(8) + t273) * MDP(27) + t263; -t269 * t204 * MDP(21) + (-t204 ^ 2 + t269 ^ 2) * MDP(22) + (t305 + t332) * MDP(23) + (t317 - t331) * MDP(24) + t197 * MDP(25) + (t171 * t329 - t173 * t254 + t186 * t269 + t175) * MDP(26) + (-t173 * t257 - t176 * t254 + t186 * t204 - t270 * t329) * MDP(27) + (MDP(23) * t311 + MDP(24) * t269 - MDP(26) * t171 + MDP(27) * t270) * qJD(5);];
tauc = t1;
