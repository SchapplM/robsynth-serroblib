% Calculate Coriolis joint torque vector for
% S5RPRRP3
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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:46:10
% EndTime: 2021-01-15 12:46:15
% DurationCPUTime: 1.28s
% Computational Cost: add. (1337->187), mult. (3159->251), div. (0->0), fcn. (2024->6), ass. (0->104)
t211 = sin(qJ(3));
t212 = cos(qJ(3));
t200 = sin(pkin(8)) * pkin(1) + pkin(6);
t267 = pkin(7) + t200;
t236 = t267 * qJD(1);
t171 = qJD(2) * t211 + t236 * t212;
t210 = sin(qJ(4));
t165 = t210 * t171;
t170 = t212 * qJD(2) - t236 * t211;
t266 = qJD(3) * pkin(3);
t168 = t170 + t266;
t269 = cos(qJ(4));
t234 = t269 * t168 - t165;
t189 = t210 * t212 + t269 * t211;
t250 = qJD(1) * t189;
t262 = t250 * qJ(5);
t142 = -t262 + t234;
t205 = qJD(3) + qJD(4);
t274 = t205 * MDP(19);
t273 = MDP(5) * t212;
t272 = (t211 ^ 2 - t212 ^ 2) * MDP(6);
t240 = t269 * qJD(4);
t271 = t269 * qJD(3) + t240;
t270 = t250 ^ 2;
t268 = pkin(3) * t205;
t259 = t210 * t211;
t230 = t205 * t259;
t242 = t269 * t212;
t231 = qJD(1) * t242;
t253 = t205 * t231;
t155 = qJD(1) * t230 - t253;
t265 = t155 * qJ(5);
t249 = qJD(1) * t211;
t181 = t210 * t249 - t231;
t264 = t181 * qJ(5);
t263 = t181 * t205;
t201 = -cos(pkin(8)) * pkin(1) - pkin(2);
t190 = -pkin(3) * t212 + t201;
t185 = t190 * qJD(1);
t260 = t185 * t250;
t213 = qJD(3) ^ 2;
t258 = t211 * t213;
t257 = t212 * t213;
t141 = pkin(4) * t205 + t142;
t256 = t141 - t142;
t161 = t205 * t189;
t156 = t161 * qJD(1);
t160 = -t212 * t271 + t230;
t255 = -t189 * t156 + t160 * t181;
t254 = t269 * t170 - t165;
t251 = MDP(11) * t212;
t192 = qJD(1) * t201;
t247 = qJD(4) * t210;
t237 = pkin(4) * t181 + qJD(5);
t157 = t185 + t237;
t246 = qJD(5) + t157;
t245 = qJD(1) * qJD(3);
t244 = pkin(3) * t249;
t243 = t211 * t266;
t167 = t269 * t171;
t239 = t211 * t245;
t199 = pkin(3) * t239;
t153 = pkin(4) * t156 + t199;
t238 = qJD(3) * t267;
t163 = t170 * qJD(3);
t164 = t171 * qJD(3);
t235 = -t210 * t163 - t269 * t164;
t233 = -t170 * t210 - t167;
t232 = t205 * t211;
t188 = -t242 + t259;
t228 = -t155 * t188 + t161 * t250;
t227 = 0.2e1 * qJD(3) * t192;
t226 = -t210 * t168 - t167;
t186 = t267 * t211;
t187 = t267 * t212;
t225 = t210 * t186 - t269 * t187;
t180 = t181 ^ 2;
t224 = t250 * t181 * MDP(12) + (-t210 * qJD(1) * t232 + t253 + t263) * MDP(14) + (-t180 + t270) * MDP(13);
t178 = t211 * t238;
t179 = t212 * t238;
t223 = -t269 * t178 - t210 * t179 - t186 * t240 - t187 * t247;
t222 = t226 * qJD(4) + t235;
t221 = t225 * qJD(4) + t210 * t178 - t269 * t179;
t220 = t269 * t163 - t210 * t164 + t168 * t240 - t171 * t247;
t219 = t222 + t265;
t218 = t185 * t181 - t220;
t217 = -t156 * qJ(5) + t220;
t216 = (-t167 + (-t168 - t268) * t210) * qJD(4) + t235;
t215 = t246 * t181 - t217;
t203 = t269 * pkin(3) + pkin(4);
t193 = t240 * t268;
t172 = pkin(4) * t250 + t244;
t169 = pkin(4) * t188 + t190;
t154 = pkin(4) * t161 + t243;
t149 = -t188 * qJ(5) - t225;
t148 = -t189 * qJ(5) - t269 * t186 - t210 * t187;
t145 = -t262 + t254;
t144 = t233 + t264;
t143 = -t226 - t264;
t140 = t160 * qJ(5) - t189 * qJD(5) + t221;
t139 = -qJ(5) * t161 - qJD(5) * t188 + t223;
t138 = -qJD(5) * t250 + t219;
t137 = -t181 * qJD(5) + t217;
t1 = [0.2e1 * t239 * t273 - 0.2e1 * t245 * t272 + MDP(7) * t257 - MDP(8) * t258 + (-t200 * t257 + t211 * t227) * MDP(10) + (t200 * t258 + t212 * t227) * MDP(11) + (-t155 * t189 - t160 * t250) * MDP(12) + (-t228 + t255) * MDP(13) + (t190 * t156 + t185 * t161 + t181 * t243 + t188 * t199) * MDP(17) + (-t190 * t155 - t185 * t160 + 0.2e1 * t250 * t243) * MDP(18) + (t153 * t188 + t154 * t181 + t156 * t169 + t157 * t161) * MDP(19) + (t153 * t189 + t154 * t250 - t155 * t169 - t157 * t160) * MDP(20) + (-t137 * t188 - t138 * t189 - t139 * t181 - t140 * t250 + t141 * t160 - t143 * t161 + t148 * t155 - t149 * t156) * MDP(21) + (t137 * t149 + t138 * t148 + t139 * t143 + t140 * t141 + t153 * t169 + t154 * t157) * MDP(22) + (-t160 * MDP(14) - t161 * MDP(15) + t221 * MDP(17) - t223 * MDP(18) + t140 * MDP(19) - t139 * MDP(20)) * t205; (t228 + t255) * MDP(21) + (t137 * t189 - t138 * t188 - t141 * t161 - t143 * t160) * MDP(22) + (-MDP(10) * t211 - t251) * t213 + ((-MDP(17) - MDP(19)) * t161 + (MDP(18) + MDP(20)) * t160) * t205; (-t181 * t244 - t233 * t205 + t216 - t260) * MDP(17) + (t254 * t205 - t244 * t250 - t193 + t218) * MDP(18) + (-t144 * t205 - t172 * t181 - t246 * t250 + t216 + t265) * MDP(19) + (t145 * t205 - t172 * t250 - t193 + t215) * MDP(20) + (t203 * t155 + (t143 + t144) * t250 + (-t141 + t145) * t181 + (-t156 * t210 + (-t269 * t181 + t210 * t250) * qJD(4)) * pkin(3)) * MDP(21) + (t138 * t203 - t141 * t144 - t143 * t145 - t157 * t172 + (t137 * t210 + (-t141 * t210 + t269 * t143) * qJD(4)) * pkin(3)) * MDP(22) + t224 + (-t211 * t273 + t272) * qJD(1) ^ 2 + (-MDP(10) * t249 - qJD(1) * t251) * t192; (-t226 * t205 + t222 - t260) * MDP(17) + (t234 * t205 + t218) * MDP(18) + (t143 * t205 + (-t157 - t237) * t250 + t219) * MDP(19) + (-t270 * pkin(4) + t142 * t205 + t215) * MDP(20) + (pkin(4) * t155 - t256 * t181) * MDP(21) + (t256 * t143 + (-t157 * t250 + t138) * pkin(4)) * MDP(22) + t224; t250 * t274 + (t253 - t263) * MDP(20) + (-t180 - t270) * MDP(21) + (t141 * t250 + t143 * t181 + t153) * MDP(22) + (t271 * t211 * MDP(19) + (-MDP(20) * t232 + t212 * t274) * t210) * qJD(1);];
tauc = t1;
