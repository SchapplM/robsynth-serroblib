% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:43
% EndTime: 2019-12-31 20:15:45
% DurationCPUTime: 0.88s
% Computational Cost: add. (834->170), mult. (1312->228), div. (0->0), fcn. (753->6), ass. (0->96)
t199 = sin(qJ(5));
t200 = sin(qJ(4));
t231 = qJD(5) * t199;
t233 = qJD(4) * t200;
t260 = -t199 * t233 - t200 * t231;
t203 = cos(qJ(4));
t259 = -t203 * t200 * MDP(10) + MDP(11) * (t200 ^ 2 - t203 ^ 2);
t258 = -MDP(5) + MDP(7);
t193 = t200 * pkin(4);
t190 = qJ(3) + t193;
t196 = qJD(1) + qJD(2);
t205 = -pkin(2) - pkin(7);
t204 = cos(qJ(2));
t251 = pkin(1) * qJD(1);
t228 = t204 * t251;
t216 = qJD(3) - t228;
t158 = t196 * t205 + t216;
t257 = -pkin(8) * t196 + t158;
t195 = qJD(4) + qJD(5);
t255 = qJD(5) - t195;
t194 = t196 ^ 2;
t223 = -pkin(1) * t204 - pkin(2);
t188 = -pkin(7) + t223;
t253 = -pkin(8) + t188;
t252 = -pkin(8) + t205;
t250 = pkin(1) * qJD(2);
t249 = qJ(3) * t196;
t202 = cos(qJ(5));
t172 = t199 * t203 + t200 * t202;
t145 = t195 * t172;
t248 = t145 * t195;
t242 = t202 * t203;
t215 = t195 * t242;
t146 = t215 + t260;
t247 = t146 * t195;
t201 = sin(qJ(2));
t246 = t195 * t201;
t245 = t196 * t203;
t244 = t199 * t200;
t147 = t257 * t200;
t243 = t202 * t147;
t206 = qJD(4) ^ 2;
t241 = t205 * t206;
t232 = qJD(4) * t203;
t192 = pkin(4) * t232;
t183 = qJD(3) + t192;
t225 = qJD(1) * t250;
t187 = t204 * t225;
t151 = t183 * t196 + t187;
t229 = t201 * t251;
t157 = t190 * t196 + t229;
t173 = t242 - t244;
t240 = -t157 * t145 + t151 * t173;
t239 = t157 * t146 + t151 * t172;
t169 = qJD(3) * t196 + t187;
t174 = t229 + t249;
t238 = t169 * t200 + t174 * t232;
t237 = t260 * t196;
t230 = pkin(4) * t245;
t227 = t201 * t250;
t226 = t204 * t250;
t224 = t196 * t242;
t165 = t253 * t203;
t178 = t252 * t203;
t189 = pkin(1) * t201 + qJ(3);
t148 = -pkin(8) * t245 + t203 * t158;
t144 = qJD(4) * pkin(4) + t148;
t219 = -pkin(4) * t195 - t144;
t217 = t201 * t225;
t182 = qJD(3) + t226;
t213 = t182 * t196 - t188 * t206;
t212 = t189 * t196 + t227;
t180 = t203 * t217;
t137 = -t233 * t257 + t180;
t138 = t200 * t217 + t232 * t257;
t154 = t196 * t244 - t224;
t211 = t202 * t137 - t199 * t138 + t157 * t154;
t135 = t145 * t196;
t155 = t172 * t196;
t210 = -t154 * t155 * MDP(17) + (t155 * t195 - t135) * MDP(19) + (-t237 + (-t154 - t224) * t195) * MDP(20) + (t154 ^ 2 - t155 ^ 2) * MDP(18);
t209 = -t174 * t196 + t217;
t208 = t147 * t231 + (-t147 * t195 - t137) * t199 + t157 * t155;
t136 = t196 * t215 + t237;
t207 = (t135 * t172 - t136 * t173 + t145 * t155 + t146 * t154) * MDP(18) + (-t135 * t173 + t145 * t154) * MDP(17) - MDP(19) * t248 - MDP(20) * t247 + (-MDP(12) * t200 - MDP(13) * t203) * t206 + 0.2e1 * t259 * qJD(4) * t196;
t191 = pkin(8) * t233;
t179 = t189 + t193;
t177 = t252 * t200;
t171 = -pkin(2) * t196 + t216;
t168 = t182 + t192;
t167 = qJD(4) * t178;
t166 = -t205 * t233 + t191;
t164 = t253 * t200;
t162 = t169 * t203;
t150 = qJD(4) * t165 + t200 * t227;
t149 = -t188 * t233 + t203 * t227 + t191;
t1 = [t207 + (t168 * t155 + t179 * t136 + (t149 * t202 - t150 * t199 + (-t164 * t202 - t165 * t199) * qJD(5)) * t195 + t239) * MDP(22) + (t162 + t213 * t203 + (-t174 - t212) * t233) * MDP(16) + (t200 * t213 + t212 * t232 + t238) * MDP(15) + (-t168 * t154 - t179 * t135 - (t149 * t199 + t150 * t202 + (-t164 * t199 + t165 * t202) * qJD(5)) * t195 + t240) * MDP(23) + (t169 * t189 + t174 * t182 + (qJD(1) * t223 + t171) * t227) * MDP(9) + (-t196 * t226 - t187) * MDP(6) + (t187 + (qJD(3) + t182) * t196) * MDP(8) + t258 * (qJD(1) + t196) * t227; t207 + (t162 + (t196 * t216 - t241) * t203 + (-t174 + t229 - t249) * t233) * MDP(16) + (t183 * t155 + t190 * t136 + (t166 * t202 - t167 * t199 + (-t177 * t202 - t178 * t199) * qJD(5)) * t195 + (-t204 * t155 - t173 * t246) * t251 + t239) * MDP(22) + (-t229 * t232 - t200 * t241 + (qJ(3) * t232 + t200 * t216) * t196 + t238) * MDP(15) + (qJ(3) * t169 + qJD(3) * t174 + (-t174 * t204 + (-pkin(2) * qJD(2) - t171) * t201) * t251) * MDP(9) + (-t183 * t154 - t190 * t135 - (t166 * t199 + t167 * t202 + (-t177 * t199 + t178 * t202) * qJD(5)) * t195 + (t204 * t154 + t172 * t246) * t251 + t240) * MDP(23) + (t196 * t228 - t187) * MDP(6) + (t187 + (0.2e1 * qJD(3) - t228) * t196) * MDP(8) + t258 * (qJD(2) - t196) * t229; -t194 * MDP(8) + t209 * MDP(9) + (-t155 * t196 - t248) * MDP(22) + (t154 * t196 - t247) * MDP(23) + (t200 * MDP(15) + t203 * MDP(16)) * (-t194 - t206); (-t174 * t245 + t180) * MDP(15) - t209 * MDP(16) * t200 + (-t155 * t230 - (-t148 * t199 - t243) * t195 + (t199 * t219 - t243) * qJD(5) + t211) * MDP(22) + (t154 * t230 + (qJD(5) * t219 + t148 * t195 - t138) * t202 + t208) * MDP(23) + t210 - t259 * t194; (t211 + t255 * (-t144 * t199 - t243)) * MDP(22) + ((-t255 * t144 - t138) * t202 + t208) * MDP(23) + t210;];
tauc = t1;
