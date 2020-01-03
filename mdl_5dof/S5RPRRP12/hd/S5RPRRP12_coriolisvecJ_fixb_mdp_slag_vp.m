% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:23
% EndTime: 2019-12-31 18:57:28
% DurationCPUTime: 1.71s
% Computational Cost: add. (1103->233), mult. (2335->345), div. (0->0), fcn. (1276->4), ass. (0->112)
t223 = 2 * qJD(1);
t193 = cos(qJ(3));
t189 = t193 ^ 2;
t191 = sin(qJ(3));
t267 = MDP(8) * (t191 ^ 2 - t189);
t266 = qJ(2) * MDP(6) + MDP(5);
t178 = pkin(3) * t191 - pkin(7) * t193 + qJ(2);
t190 = sin(qJ(4));
t192 = cos(qJ(4));
t194 = -pkin(1) - pkin(6);
t247 = t191 * t194;
t239 = t190 * t178 + t192 * t247;
t235 = qJD(1) * t193;
t217 = t192 * t235;
t234 = qJD(3) * t190;
t173 = t217 + t234;
t265 = t173 ^ 2;
t264 = -qJ(5) - pkin(7);
t224 = qJD(3) * qJD(4);
t186 = t192 * t224;
t227 = t192 * qJD(3);
t213 = t191 * t227;
t229 = qJD(4) * t193;
t215 = t190 * t229;
t199 = t213 + t215;
t148 = qJD(1) * t199 - t186;
t262 = t148 * t190;
t226 = qJD(1) * qJD(3);
t212 = t191 * t226;
t182 = t190 * t212;
t149 = qJD(4) * t173 - t182;
t261 = t149 * t192;
t184 = qJD(1) * t194 + qJD(2);
t254 = t184 * t193;
t164 = -qJD(3) * pkin(3) - t254;
t260 = t164 * t190;
t218 = t190 * t235;
t171 = t218 - t227;
t259 = t171 * t190;
t258 = t171 * t192;
t236 = qJD(1) * t191;
t185 = qJD(4) + t236;
t257 = t173 * t185;
t256 = t173 * t190;
t255 = t173 * t192;
t253 = t185 * t190;
t252 = t185 * t191;
t251 = t185 * t192;
t250 = t190 * t193;
t249 = t190 * t194;
t177 = t191 * t184;
t248 = t191 * t192;
t246 = t192 * t193;
t196 = qJD(1) ^ 2;
t245 = t193 * t196;
t195 = qJD(3) ^ 2;
t244 = t194 * t195;
t160 = t178 * qJD(1);
t163 = qJD(3) * pkin(7) + t177;
t143 = t192 * t160 - t163 * t190;
t138 = -qJ(5) * t173 + t143;
t137 = pkin(4) * t185 + t138;
t243 = t137 - t138;
t209 = qJD(4) * t264;
t228 = qJD(5) * t192;
t206 = pkin(3) * t193 + pkin(7) * t191;
t176 = t206 * qJD(1);
t240 = t190 * t176 + t184 * t246;
t242 = t228 - t240 + (-qJ(5) * t236 + t209) * t190;
t162 = t192 * t176;
t221 = t184 * t250;
t241 = -qJD(5) * t190 + t192 * t209 + t221 - t162 - (pkin(4) * t193 + qJ(5) * t248) * qJD(1);
t237 = -t195 - t196;
t233 = qJD(3) * t191;
t232 = qJD(3) * t193;
t231 = qJD(4) * t190;
t230 = qJD(4) * t192;
t225 = qJD(3) * MDP(18);
t220 = t190 * t247;
t170 = qJD(3) * t206 + qJD(2);
t216 = t193 * t227;
t219 = t190 * t170 + t178 * t230 + t194 * t216;
t214 = t192 * t229;
t141 = pkin(4) * t149 + t184 * t233;
t211 = pkin(4) - t249;
t210 = pkin(4) * t190 - t194;
t208 = t171 + t227;
t207 = -t173 + t234;
t144 = t160 * t190 + t163 * t192;
t139 = -qJ(5) * t171 + t144;
t205 = t137 * t192 + t139 * t190;
t204 = t137 * t190 - t139 * t192;
t203 = qJD(1) * t189 - t252;
t202 = -pkin(7) * t232 + t164 * t191;
t154 = t170 * qJD(1);
t200 = -t190 * t154 - t160 * t230 + t163 * t231 - t184 * t216;
t151 = t192 * t154;
t198 = -t144 * qJD(4) + t151;
t197 = (t255 + t259) * MDP(21) - t205 * MDP(22) + (-t192 * MDP(19) + t190 * MDP(20)) * t185;
t181 = t264 * t192;
t180 = t264 * t190;
t169 = t171 ^ 2;
t168 = t192 * t178;
t158 = t192 * t170;
t147 = -qJ(5) * t250 + t239;
t146 = pkin(4) * t171 + qJD(5) + t164;
t145 = -qJ(5) * t246 + t191 * t211 + t168;
t136 = -qJ(5) * t214 + (-qJD(5) * t193 + (qJ(5) * qJD(3) - qJD(4) * t194) * t191) * t190 + t219;
t135 = qJ(5) * t213 + t158 - t239 * qJD(4) + (qJ(5) * t231 + qJD(3) * t211 - t228) * t193;
t134 = -qJ(5) * t149 - qJD(5) * t171 - t200;
t133 = qJ(5) * t148 - qJD(5) * t173 + (pkin(4) * qJD(1) - t184 * t190) * t232 + t198;
t1 = [-0.2e1 * t193 * MDP(7) * t212 + 0.2e1 * t226 * t267 + (-t191 * t244 + (qJ(2) * t232 + qJD(2) * t191) * t223) * MDP(12) + (-t193 * t244 + (-qJ(2) * t233 + qJD(2) * t193) * t223) * MDP(13) + (-t148 * t246 - t173 * t199) * MDP(14) + ((t256 + t258) * t233 + (t262 - t261 + (-t255 + t259) * qJD(4)) * t193) * MDP(15) + (-t185 * t215 - t148 * t191 + (t173 * t193 + t192 * t203) * qJD(3)) * MDP(16) + (-t185 * t214 - t149 * t191 + (-t171 * t193 - t190 * t203) * qJD(3)) * MDP(17) + (t185 + t236) * t193 * t225 + (-t193 * t194 * t149 + t151 * t191 + t158 * t185 + (-t144 * t191 + t164 * t246 - t185 * t239) * qJD(4) + ((t171 * t194 - t260) * t191 + (-t185 * t249 + (t168 - t220) * qJD(1) + t143) * t193) * qJD(3)) * MDP(19) + (-(-qJD(4) * t220 + t219) * t185 + t200 * t191 + (t194 * t148 - t164 * t231) * t193 + ((-qJD(1) * t239 - t144) * t193 + (t194 * t173 + (-t164 + t254) * t192) * t191) * qJD(3)) * MDP(20) + (-t135 * t173 - t136 * t171 + t145 * t148 - t147 * t149 + t205 * t233 + (qJD(4) * t204 - t133 * t192 - t134 * t190) * t193) * MDP(21) + (t133 * t145 + t134 * t147 + t137 * t135 + t139 * t136 - t146 * t210 * t233 + (pkin(4) * t146 * t230 + t141 * t210) * t193) * MDP(22) + t266 * qJD(2) * t223 + (-MDP(10) * t193 - MDP(9) * t191) * t195; -t266 * t196 + t197 * qJD(1) + (t237 * MDP(13) - t149 * MDP(19) + t148 * MDP(20) - t141 * MDP(22) + ((t256 - t258) * MDP(21) - t204 * MDP(22) + (-t190 * MDP(19) - t192 * MDP(20)) * t185) * qJD(3)) * t193 + (t237 * MDP(12) + (-t261 - t262) * MDP(21) + (-t133 * t190 + t134 * t192) * MDP(22) + t197 * qJD(4) + ((t171 - t218) * MDP(19) + (t173 - t217) * MDP(20) + t146 * MDP(22)) * qJD(3)) * t191; t191 * MDP(7) * t245 - t196 * t267 + (t173 * t251 - t262) * MDP(14) + ((-t171 * t185 - t148) * t192 + (-t149 - t257) * t190) * MDP(15) + (t185 * t230 + (t185 * t248 + t193 * t207) * qJD(1)) * MDP(16) + (-t185 * t231 + (-t190 * t252 + t193 * t208) * qJD(1)) * MDP(17) - t185 * MDP(18) * t235 + (-pkin(3) * t149 - t162 * t185 + (t185 * t250 - t191 * t208) * t184 + (-pkin(7) * t251 + t260) * qJD(4) + (-t143 * t193 + t190 * t202) * qJD(1)) * MDP(19) + (pkin(3) * t148 + t240 * t185 + t207 * t177 + (pkin(7) * t253 + t164 * t192) * qJD(4) + (t144 * t193 + t192 * t202) * qJD(1)) * MDP(20) + (t148 * t180 + t149 * t181 - t241 * t173 - t242 * t171 + (-t137 * t185 + t134) * t192 + (-t139 * t185 - t133) * t190) * MDP(21) + (-t134 * t181 + t133 * t180 + t141 * (-pkin(4) * t192 - pkin(3)) + (pkin(4) * t253 - t177) * t146 + t242 * t139 + t241 * t137) * MDP(22) + (MDP(13) * t191 * t196 - MDP(12) * t245) * qJ(2); (-t169 + t265) * MDP(15) + t186 * MDP(16) + (-t190 * t224 + t182 + t257) * MDP(17) + (-qJD(3) * t221 + t144 * t185 - t164 * t173 + t198) * MDP(19) + (t143 * t185 + t200) * MDP(20) + t243 * MDP(22) * t139 + (t148 * MDP(21) + (-t146 * t173 + t133) * MDP(22)) * pkin(4) + (t173 * MDP(14) + t185 * MDP(16) + t164 * MDP(20) - MDP(21) * t243) * t171 + (-MDP(16) * t213 + (t225 + (-t190 * MDP(16) - t192 * MDP(17)) * qJD(4)) * t193) * qJD(1); (-t169 - t265) * MDP(21) + (t137 * t173 + t139 * t171 + t141) * MDP(22);];
tauc = t1;
