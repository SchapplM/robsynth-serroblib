% Calculate Coriolis joint torque vector for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:17:02
% EndTime: 2022-01-23 09:17:05
% DurationCPUTime: 1.78s
% Computational Cost: add. (1362->214), mult. (3944->320), div. (0->0), fcn. (3018->8), ass. (0->116)
t253 = sin(pkin(9));
t255 = cos(pkin(9));
t258 = sin(qJ(4));
t260 = cos(qJ(4));
t234 = t253 * t260 + t255 * t258;
t230 = t234 * qJD(4);
t254 = sin(pkin(8));
t214 = t254 * t230;
t204 = qJD(1) * t214;
t298 = qJD(1) * t254;
t285 = t253 * t298;
t275 = t258 * t285;
t236 = qJD(4) * t275;
t305 = t255 * t260;
t287 = t254 * t305;
t276 = qJD(4) * t287;
t205 = qJD(1) * t276 - t236;
t212 = qJD(1) * t287 - t275;
t257 = sin(qJ(5));
t259 = cos(qJ(5));
t292 = qJD(5) * t257;
t265 = qJD(1) * t234;
t209 = t254 * t265;
t303 = t259 * t209;
t171 = -qJD(5) * t303 - t259 * t204 - t257 * t205 - t212 * t292;
t269 = -t209 * t257 + t259 * t212;
t172 = t269 * qJD(5) - t204 * t257 + t259 * t205;
t178 = t212 * t257 + t303;
t256 = cos(pkin(8));
t297 = qJD(1) * t256;
t318 = qJD(4) - t297;
t242 = -qJD(5) - t318;
t309 = t178 * t242;
t310 = t269 * t242;
t319 = t178 * t269 * MDP(17) + (-t172 - t310) * MDP(20) + (-t178 ^ 2 + t269 ^ 2) * MDP(18) + (t171 - t309) * MDP(19);
t237 = -pkin(2) * t256 - qJ(3) * t254 - pkin(1);
t225 = t237 * qJD(1) + qJD(2);
t217 = t255 * t225;
t264 = -pkin(6) * t254 * t255 + (-qJ(2) * t253 - pkin(3)) * t256;
t184 = t264 * qJD(1) + t217;
t286 = qJ(2) * t297;
t194 = t253 * t225 + t255 * t286;
t189 = -pkin(6) * t285 + t194;
t271 = -t184 * t258 - t189 * t260;
t295 = qJD(3) * t254;
t296 = qJD(2) * t256;
t227 = -t253 * t296 - t255 * t295;
t221 = t227 * qJD(1);
t228 = -t253 * t295 + t255 * t296;
t222 = t228 * qJD(1);
t278 = t260 * t221 - t258 * t222;
t263 = t271 * qJD(4) + t278;
t164 = pkin(7) * t204 + t263;
t170 = -pkin(7) * t209 - t271;
t168 = t170 * t292;
t243 = qJ(2) * t298 + qJD(3);
t226 = pkin(3) * t285 + t243;
t187 = pkin(4) * t209 + t226;
t317 = t187 * t178 + t168 + (t170 * t242 - t164) * t257;
t314 = qJD(5) + t242;
t293 = qJD(4) * t260;
t294 = qJD(4) * t258;
t267 = -t184 * t293 + t189 * t294 - t258 * t221 - t260 * t222;
t163 = -pkin(7) * t205 - t267;
t281 = -t257 * t163 + t259 * t164;
t313 = -t187 * t269 + t281;
t312 = pkin(4) * t212;
t311 = qJ(2) * t256;
t308 = t209 * t318;
t307 = t212 * t318;
t306 = t253 * t254;
t304 = t259 * t170;
t302 = -t256 * t265 + t230;
t233 = -t253 * t258 + t305;
t301 = t318 * t233;
t300 = t253 * t237 + t255 * t311;
t235 = pkin(3) * t306 + t254 * qJ(2);
t251 = t254 ^ 2;
t252 = t256 ^ 2;
t299 = t251 + t252;
t290 = t254 * qJD(2);
t289 = qJD(1) * qJD(2);
t288 = 0.2e1 * qJD(2) * t251;
t284 = qJ(2) * t289;
t280 = t260 * t184 - t189 * t258;
t169 = -pkin(7) * t212 + t280;
t167 = pkin(4) * t318 + t169;
t283 = pkin(4) * t242 - t167;
t277 = t260 * t227 - t228 * t258;
t274 = qJD(5) * t234 + t302;
t273 = qJD(5) * t233 + t301;
t272 = -t167 * t257 - t304;
t232 = t255 * t237;
t190 = t232 + t264;
t195 = -pkin(6) * t306 + t300;
t270 = -t190 * t258 - t195 * t260;
t223 = t234 * t254;
t224 = t233 * t254;
t185 = t223 * t259 + t224 * t257;
t186 = -t223 * t257 + t224 * t259;
t266 = t190 * t293 - t195 * t294 + t258 * t227 + t260 * t228;
t261 = qJD(1) ^ 2;
t248 = t254 * t289;
t245 = t251 * t284;
t215 = -t294 * t306 + t276;
t196 = pkin(4) * t215 + t290;
t193 = -t253 * t286 + t217;
t192 = pkin(4) * t205 + t248;
t191 = pkin(4) * t223 + t235;
t176 = -pkin(7) * t223 - t270;
t175 = -pkin(4) * t256 - pkin(7) * t224 + t190 * t260 - t195 * t258;
t174 = t186 * qJD(5) - t214 * t257 + t259 * t215;
t173 = -t185 * qJD(5) - t214 * t259 - t215 * t257;
t166 = pkin(7) * t214 + t270 * qJD(4) + t277;
t165 = -pkin(7) * t215 + t266;
t1 = [0.2e1 * t299 * MDP(5) * t289 + 0.2e1 * (t252 * t284 + t245) * MDP(6) + (-t221 * t256 + (-t227 * t256 + t253 * t288) * qJD(1)) * MDP(7) + (t222 * t256 + (t228 * t256 + t255 * t288) * qJD(1)) * MDP(8) + (t222 * t300 + t194 * t228 + t221 * (-t253 * t311 + t232) + t193 * t227 + t245 + t243 * t290) * MDP(9) + (-t204 * t224 - t212 * t214) * MDP(10) + (t204 * t223 - t205 * t224 + t209 * t214 - t212 * t215) * MDP(11) + (t204 * t256 - t214 * t318) * MDP(12) + (t205 * t256 - t215 * t318) * MDP(13) + (t277 * t318 - t278 * t256 + t235 * t205 + t226 * t215 + (-t271 * t256 + t270 * t318) * qJD(4) + (qJD(1) * t223 + t209) * t290) * MDP(15) + (-t266 * t318 - t267 * t256 - t235 * t204 - t226 * t214 + (qJD(1) * t224 + t212) * t290) * MDP(16) + (t171 * t186 + t173 * t269) * MDP(17) + (-t171 * t185 - t172 * t186 - t173 * t178 - t174 * t269) * MDP(18) + (-t171 * t256 - t173 * t242) * MDP(19) + (t172 * t256 + t174 * t242) * MDP(20) + (-(-t165 * t257 + t166 * t259) * t242 - t281 * t256 + t196 * t178 + t191 * t172 + t192 * t185 + t187 * t174 + (-(-t175 * t257 - t176 * t259) * t242 - t272 * t256) * qJD(5)) * MDP(22) + (-t168 * t256 + t191 * t171 + t187 * t173 + t196 * t269 + t192 * t186 + ((-qJD(5) * t176 + t166) * t242 + t164 * t256) * t257 + ((qJD(5) * t175 + t165) * t242 + (qJD(5) * t167 + t163) * t256) * t259) * MDP(23); (t221 * t255 + t222 * t253 + (-t243 * t254 + (t193 * t253 - t194 * t255) * t256) * qJD(1)) * MDP(9) + (-t209 * t298 - t302 * t318) * MDP(15) + (-t212 * t298 - t301 * t318) * MDP(16) + (-t178 * t298 + (t273 * t257 + t274 * t259) * t242) * MDP(22) + (-t269 * t298 + (-t274 * t257 + t273 * t259) * t242) * MDP(23) - (qJ(2) * MDP(6) + t253 * MDP(7) + t255 * MDP(8) + MDP(5)) * t299 * t261; t248 * MDP(9) + (-t236 + t307) * MDP(15) - MDP(16) * t308 + (t172 - t310) * MDP(22) + (t171 + t309) * MDP(23) + ((-MDP(7) * t255 + MDP(8) * t253) * t261 * t256 + ((t193 * t255 + t194 * t253) * MDP(9) + (MDP(15) * t305 - t234 * MDP(16)) * qJD(4)) * qJD(1)) * t254; t212 * t209 * MDP(10) + (-t209 ^ 2 + t212 ^ 2) * MDP(11) + (-t204 + t308) * MDP(12) + (-t205 + t307) * MDP(13) + (-t226 * t212 - t271 * t318 + t263) * MDP(15) + (t226 * t209 + t280 * t318 + t267) * MDP(16) + ((-t169 * t257 - t304) * t242 - t178 * t312 + (t283 * t257 - t304) * qJD(5) + t313) * MDP(22) + (-t269 * t312 + (t283 * qJD(5) - t169 * t242 - t163) * t259 + t317) * MDP(23) + t319; (t314 * t272 + t313) * MDP(22) + ((-t314 * t167 - t163) * t259 + t317) * MDP(23) + t319;];
tauc = t1;
