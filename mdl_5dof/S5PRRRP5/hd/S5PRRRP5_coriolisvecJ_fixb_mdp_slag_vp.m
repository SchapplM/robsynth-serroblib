% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:22
% EndTime: 2019-12-05 16:49:25
% DurationCPUTime: 1.16s
% Computational Cost: add. (990->178), mult. (2485->250), div. (0->0), fcn. (1672->6), ass. (0->102)
t214 = cos(qJ(3));
t212 = sin(qJ(2));
t247 = t212 * qJD(1);
t197 = qJD(2) * pkin(6) + t247;
t237 = pkin(7) * qJD(2) + t197;
t174 = t237 * t214;
t210 = sin(qJ(4));
t167 = t210 * t174;
t211 = sin(qJ(3));
t173 = t237 * t211;
t170 = qJD(3) * pkin(3) - t173;
t213 = cos(qJ(4));
t234 = t213 * t170 - t167;
t190 = t210 * t214 + t211 * t213;
t181 = t190 * qJD(2);
t267 = t181 * qJ(5);
t281 = t267 - t234;
t207 = qJD(3) + qJD(4);
t279 = t207 * t190;
t158 = t279 * qJD(2);
t245 = qJD(2) * qJD(3);
t280 = -0.2e1 * t245;
t278 = MDP(5) * t211;
t263 = t213 * t214;
t265 = t210 * t211;
t189 = -t263 + t265;
t176 = t189 * t212;
t252 = qJD(3) * t211;
t215 = cos(qJ(2));
t254 = qJD(1) * t215;
t163 = -t197 * t252 + (-pkin(7) * t252 + t214 * t254) * qJD(2);
t277 = (qJD(4) * t170 + t163) * t213;
t276 = (t211 ^ 2 - t214 ^ 2) * MDP(6);
t271 = pkin(6) + pkin(7);
t242 = qJD(3) * t271;
t191 = t211 * t242;
t192 = t214 * t242;
t193 = t271 * t211;
t194 = t271 * t214;
t228 = t193 * t210 - t194 * t213;
t275 = t228 * qJD(4) + t190 * t254 + t210 * t191 - t213 * t192;
t224 = t189 * t215;
t249 = qJD(4) * t213;
t250 = qJD(4) * t210;
t274 = -qJD(1) * t224 + t213 * t191 + t210 * t192 + t193 * t249 + t194 * t250;
t273 = t211 * MDP(10) + t214 * MDP(11);
t272 = t181 ^ 2;
t270 = qJD(2) * pkin(2);
t240 = qJD(2) * t263;
t253 = qJD(2) * t211;
t241 = t210 * t253;
t179 = -t240 + t241;
t269 = qJ(5) * t179;
t206 = -pkin(3) * t214 - pkin(2);
t186 = t206 * qJD(2) - t254;
t160 = pkin(4) * t179 + qJD(5) + t186;
t268 = t160 * t181;
t266 = t186 * t181;
t216 = qJD(3) ^ 2;
t264 = t211 * t216;
t169 = t213 * t174;
t262 = t214 * t216;
t261 = -qJ(5) * t279 - qJD(5) * t189 - t274;
t230 = t207 * t265;
t251 = qJD(3) * t214;
t161 = -t213 * t251 - t214 * t249 + t230;
t260 = qJ(5) * t161 - qJD(5) * t190 + t275;
t144 = pkin(4) * t207 - t281;
t259 = t144 + t281;
t258 = -t213 * t173 - t167;
t239 = t214 * t245;
t257 = -qJD(4) * t240 - t213 * t239;
t243 = pkin(3) * t252;
t187 = (t243 + t247) * qJD(2);
t244 = pkin(3) * t253;
t238 = -pkin(3) * t207 - t170;
t164 = -t197 * t251 + (-pkin(7) * t251 - t211 * t254) * qJD(2);
t236 = -t210 * t163 + t213 * t164;
t235 = t210 * t164 - t174 * t250;
t233 = t173 * t210 - t169;
t153 = pkin(4) * t158 + t187;
t229 = -t170 * t210 - t169;
t226 = t243 - t247;
t225 = t186 * t179 - t235;
t178 = t179 ^ 2;
t223 = t181 * t179 * MDP(12) + (-t257 + (t179 - t241) * t207) * MDP(14) + (t181 * t207 - t158) * MDP(15) + (-t178 + t272) * MDP(13);
t221 = -0.2e1 * qJD(3) * t270;
t220 = t229 * qJD(4) + t236;
t217 = qJD(2) ^ 2;
t205 = pkin(3) * t213 + pkin(4);
t175 = t190 * t212;
t157 = qJD(2) * t230 + t257;
t156 = -qJ(5) * t189 - t228;
t155 = -qJ(5) * t190 - t193 * t213 - t194 * t210;
t150 = t207 * t176 - t215 * t181;
t149 = -qJD(2) * t224 - t212 * t279;
t148 = -t267 + t258;
t147 = t233 + t269;
t146 = -t229 - t269;
t141 = qJ(5) * t157 - qJD(5) * t181 + t220;
t140 = -qJ(5) * t158 - qJD(5) * t179 + t235 + t277;
t1 = [(-t149 * t179 - t150 * t181 - t157 * t175 + t158 * t176) * MDP(19) + (-t140 * t176 - t141 * t175 + t144 * t150 + t146 * t149) * MDP(20) + (MDP(17) * t150 - MDP(18) * t149) * t207 + (-t158 * MDP(17) + t157 * MDP(18) - t153 * MDP(20) - t217 * MDP(4) + t273 * t280) * t215 + (-t217 * MDP(3) + (t179 * MDP(17) + t181 * MDP(18) + t160 * MDP(20)) * qJD(2) + (-MDP(10) * t214 + MDP(11) * t211) * (t216 + t217)) * t212; 0.2e1 * t239 * t278 + t276 * t280 + MDP(7) * t262 - MDP(8) * t264 + (-pkin(6) * t262 + t211 * t221) * MDP(10) + (pkin(6) * t264 + t214 * t221) * MDP(11) + (-t157 * t190 - t161 * t181) * MDP(12) + (t157 * t189 - t158 * t190 + t161 * t179 - t181 * t279) * MDP(13) + (t206 * t158 + t226 * t179 + t186 * t279 + t187 * t189) * MDP(17) + (-t206 * t157 - t186 * t161 + t226 * t181 + t187 * t190) * MDP(18) + (-t140 * t189 - t141 * t190 + t144 * t161 - t146 * t279 + t155 * t157 - t156 * t158 - t261 * t179 - t260 * t181) * MDP(19) + (t140 * t156 + t141 * t155 + t153 * (pkin(4) * t189 + t206) + (pkin(4) * t279 + t226) * t160 + t261 * t146 + t260 * t144) * MDP(20) + (-t161 * MDP(14) - MDP(15) * t279 + t275 * MDP(17) + t274 * MDP(18)) * t207; (-t233 * t207 - t179 * t244 - t266 + (t238 * t210 - t169) * qJD(4) + t236) * MDP(17) + (t258 * t207 - t181 * t244 + (t238 * qJD(4) - t163) * t213 + t225) * MDP(18) + (t157 * t205 + (t146 + t147) * t181 + (-t144 + t148) * t179 + (-t158 * t210 + (-t179 * t213 + t181 * t210) * qJD(4)) * pkin(3)) * MDP(19) + (-pkin(4) * t268 + t141 * t205 - t144 * t147 - t146 * t148 + (-t160 * t253 + t140 * t210 + (-t144 * t210 + t146 * t213) * qJD(4)) * pkin(3)) * MDP(20) + t223 + t273 * qJD(2) * t270 + (-t214 * t278 + t276) * t217; (-t229 * t207 + t220 - t266) * MDP(17) + (t234 * t207 + t225 - t277) * MDP(18) + (pkin(4) * t157 - t259 * t179) * MDP(19) + (t259 * t146 + (t141 - t268) * pkin(4)) * MDP(20) + t223; (-t178 - t272) * MDP(19) + (t144 * t181 + t146 * t179 + t153) * MDP(20);];
tauc = t1;
