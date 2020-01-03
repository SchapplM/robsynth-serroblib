% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRR2
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
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:21
% EndTime: 2019-12-31 17:23:23
% DurationCPUTime: 0.84s
% Computational Cost: add. (638->129), mult. (1145->203), div. (0->0), fcn. (698->6), ass. (0->85)
t239 = pkin(6) + pkin(7);
t180 = sin(qJ(4));
t181 = sin(qJ(3));
t183 = cos(qJ(4));
t184 = cos(qJ(3));
t155 = t180 * t184 + t183 * t181;
t177 = qJD(1) + qJD(2);
t146 = t155 * t177;
t176 = qJD(3) + qJD(4);
t235 = qJD(4) - t176;
t238 = MDP(7) * t181;
t237 = (t181 ^ 2 - t184 ^ 2) * MDP(8);
t236 = t181 * MDP(12) + t184 * MDP(13);
t185 = cos(qJ(2));
t233 = t185 * pkin(1);
t182 = sin(qJ(2));
t171 = pkin(1) * t182 + pkin(6);
t232 = -pkin(7) - t171;
t231 = pkin(1) * qJD(1);
t230 = pkin(1) * qJD(2);
t229 = t176 * t185;
t228 = t177 * t181;
t227 = t180 * t181;
t186 = qJD(3) ^ 2;
t225 = t181 * t186;
t198 = t177 * t239 + t182 * t231;
t142 = t198 * t184;
t224 = t183 * t142;
t222 = t183 * t184;
t221 = t184 * t186;
t154 = -t222 + t227;
t133 = t176 * t154;
t173 = -pkin(3) * t184 - pkin(2);
t206 = t185 * t231;
t147 = t173 * t177 - t206;
t205 = qJD(1) * t230;
t196 = t182 * t205;
t216 = qJD(3) * t181;
t201 = t177 * t216;
t148 = pkin(3) * t201 + t196;
t220 = -t133 * t147 + t148 * t155;
t134 = t176 * t155;
t219 = t134 * t147 + t148 * t154;
t160 = -pkin(2) * t177 - t206;
t215 = qJD(3) * t184;
t218 = t160 * t215 + t181 * t196;
t214 = qJD(3) * t185;
t211 = -qJD(1) - t177;
t210 = -qJD(2) + t177;
t209 = pkin(3) * t228;
t208 = t185 * t230;
t207 = pkin(3) * t216;
t204 = t177 * t227;
t203 = t177 * t222;
t202 = qJD(3) * t239;
t200 = t177 * t215;
t141 = t198 * t181;
t138 = qJD(3) * pkin(3) - t141;
t199 = -pkin(3) * t176 - t138;
t197 = qJD(3) * t232;
t195 = t185 * t205;
t193 = qJD(3) * t198;
t131 = -t181 * t193 + t184 * t195;
t132 = -t181 * t195 - t184 * t193;
t191 = -t180 * t131 + t132 * t183 - t147 * t146;
t126 = qJD(4) * t203 - t176 * t204 + t183 * t200;
t144 = -t203 + t204;
t190 = t146 * t144 * MDP(14) + (-t144 ^ 2 + t146 ^ 2) * MDP(15) + (t144 * t176 + t126) * MDP(16);
t188 = t147 * t144 + (t142 * t235 - t132) * t180;
t127 = t134 * t177;
t187 = -MDP(10) * t225 + (-t126 * t154 - t127 * t155 + t133 * t144 - t134 * t146) * MDP(15) + (t126 * t155 - t133 * t146) * MDP(14) - 0.2e1 * t177 * qJD(3) * t237 + 0.2e1 * t200 * t238 + MDP(9) * t221 + (-MDP(16) * t133 - MDP(17) * t134) * t176;
t174 = t184 * pkin(7);
t172 = -pkin(2) - t233;
t169 = pkin(6) * t184 + t174;
t168 = t239 * t181;
t163 = t173 - t233;
t158 = t182 * t230 + t207;
t157 = t184 * t202;
t156 = t181 * t202;
t152 = t171 * t184 + t174;
t151 = t232 * t181;
t149 = t160 * t216;
t140 = -t181 * t208 + t184 * t197;
t139 = t181 * t197 + t184 * t208;
t1 = [(-t171 * t221 + t172 * t201 + t149) * MDP(12) + (t171 * t225 + t172 * t200 + t218) * MDP(13) + (t158 * t144 + t163 * t127 + (-t180 * t139 + t183 * t140 + (-t151 * t180 - t152 * t183) * qJD(4)) * t176 + t219) * MDP(19) + (t158 * t146 + t163 * t126 - (t183 * t139 + t180 * t140 + (t151 * t183 - t152 * t180) * qJD(4)) * t176 + t220) * MDP(20) + ((t211 * MDP(6) - qJD(3) * t236) * t185 + (MDP(13) * t228 + (MDP(12) * t184 + MDP(5)) * t211) * t182) * t230 + t187; (-pkin(2) * t201 - pkin(6) * t221 + t149 + (t182 * t184 * t210 + t181 * t214) * t231) * MDP(12) + (-pkin(2) * t200 + pkin(6) * t225 + (-t182 * t228 + t184 * t214) * t231 + t218) * MDP(13) + (t144 * t207 + t173 * t127 + (t180 * t156 - t183 * t157 + (t168 * t180 - t169 * t183) * qJD(4)) * t176 + (-t182 * t144 + t155 * t229) * t231 + t219) * MDP(19) + (t146 * t207 + t173 * t126 - (-t183 * t156 - t180 * t157 + (-t168 * t183 - t169 * t180) * qJD(4)) * t176 + (-t182 * t146 - t154 * t229) * t231 + t220) * MDP(20) + t187 + (MDP(5) * t182 + MDP(6) * t185) * t210 * t231; (-t144 * t209 - (t141 * t180 - t224) * t176 + (t180 * t199 - t224) * qJD(4) + t191) * MDP(19) + (-t146 * t209 + (qJD(4) * t199 - t141 * t176 - t131) * t183 + t188) * MDP(20) + t190 + t236 * (-t160 * t177 - t195) + (-t184 * t238 + t237) * t177 ^ 2; (t191 + t235 * (-t138 * t180 - t224)) * MDP(19) + ((-t138 * t235 - t131) * t183 + t188) * MDP(20) + t190;];
tauc = t1;
