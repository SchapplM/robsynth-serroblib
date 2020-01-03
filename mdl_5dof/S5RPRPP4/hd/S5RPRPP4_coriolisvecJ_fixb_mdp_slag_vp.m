% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:52
% EndTime: 2019-12-31 18:14:54
% DurationCPUTime: 0.76s
% Computational Cost: add. (944->166), mult. (1925->217), div. (0->0), fcn. (1099->4), ass. (0->78)
t192 = sin(qJ(3));
t193 = cos(qJ(3));
t240 = -t193 * t192 * MDP(7) + (t192 ^ 2 - t193 ^ 2) * MDP(8) + (t193 * MDP(12) - t192 * MDP(13)) * qJ(2);
t239 = MDP(14) + MDP(17);
t236 = t192 * MDP(12) + t193 * MDP(13);
t194 = -pkin(1) - pkin(6);
t227 = qJ(4) - t194;
t207 = t227 * t193;
t217 = t192 * qJD(4);
t154 = -qJD(3) * t207 - t217;
t190 = sin(pkin(7));
t191 = cos(pkin(7));
t215 = t193 * qJD(4);
t221 = qJD(3) * t192;
t200 = t227 * t221 - t215;
t133 = t190 * t154 - t191 * t200;
t134 = t191 * t154 + t190 * t200;
t170 = t227 * t192;
t144 = -t190 * t170 + t191 * t207;
t145 = -t191 * t170 - t190 * t207;
t212 = qJD(1) * qJD(3);
t208 = t193 * t212;
t171 = t191 * t208;
t223 = qJD(1) * t192;
t209 = t190 * t223;
t152 = qJD(3) * t209 - t171;
t205 = t190 * t193 + t191 * t192;
t153 = t205 * t212;
t159 = t205 * qJD(1);
t222 = qJD(1) * t193;
t163 = t191 * t222 - t209;
t234 = t133 * t163 - t134 * t159 - t144 * t153 + t145 * t152;
t158 = t163 ^ 2;
t172 = t194 * qJD(1) + qJD(2);
t220 = qJD(3) * t193;
t146 = t172 * t220 + (-qJ(4) * t220 - t217) * qJD(1);
t197 = -t172 * t221 + (qJ(4) * t221 - t215) * qJD(1);
t123 = t190 * t146 - t191 * t197;
t233 = t123 * t144;
t166 = -t190 * t192 + t191 * t193;
t232 = t123 * t166;
t155 = (-qJ(4) * qJD(1) + t172) * t192;
t229 = t190 * t155;
t149 = t191 * t155;
t228 = t192 * pkin(3) + qJ(2);
t124 = t191 * t146 + t190 * t197;
t156 = -qJ(4) * t222 + t193 * t172;
t151 = qJD(3) * pkin(3) + t156;
t131 = t190 * t151 + t149;
t226 = pkin(3) * t208 + qJD(1) * qJD(2);
t169 = pkin(3) * t223 + qJD(1) * qJ(2) + qJD(4);
t219 = t169 * qJD(1);
t214 = pkin(3) * t220 + qJD(2);
t136 = t191 * t156 - t229;
t213 = qJD(5) - t136;
t211 = qJD(3) * qJD(5);
t206 = -qJ(2) * MDP(6) - MDP(5);
t130 = t191 * t151 - t229;
t203 = -t152 * pkin(4) + t153 * qJ(5) + t226;
t201 = -t190 * t222 - t191 * t223;
t122 = t211 + t124;
t127 = -qJD(3) * pkin(4) + qJD(5) - t130;
t128 = qJD(3) * qJ(5) + t131;
t161 = t190 * t221 - t191 * t220;
t162 = -t190 * t220 - t191 * t221;
t199 = t122 * t205 - t127 * t162 - t128 * t161 - t232;
t198 = t124 * t205 + t130 * t162 - t131 * t161 - t232;
t196 = qJD(1) ^ 2;
t195 = qJD(3) ^ 2;
t182 = -t191 * pkin(3) - pkin(4);
t178 = t190 * pkin(3) + qJ(5);
t138 = pkin(4) * t205 - t166 * qJ(5) + t228;
t137 = pkin(3) * t222 + t163 * pkin(4) + t159 * qJ(5);
t135 = t190 * t156 + t149;
t132 = t159 * pkin(4) - t163 * qJ(5) + t169;
t126 = -t161 * pkin(4) - t162 * qJ(5) - t166 * qJD(5) + t214;
t121 = -t163 * qJD(5) + t203;
t1 = [(-t198 + t234) * MDP(14) + (t124 * t145 - t130 * t133 + t131 * t134 + t169 * t214 + t226 * t228 + t233) * MDP(15) + (-t133 * qJD(3) + t121 * t205 + t126 * t159 - t132 * t161 - t138 * t152) * MDP(16) + (-t199 + t234) * MDP(17) + (t134 * qJD(3) - t121 * t166 - t126 * t163 - t132 * t162 + t138 * t153) * MDP(18) + (t121 * t138 + t122 * t145 + t132 * t126 + t127 * t133 + t128 * t134 + t233) * MDP(19) + ((-MDP(13) * t194 - MDP(10)) * t193 + (-MDP(12) * t194 - MDP(9)) * t192) * t195 + (0.2e1 * (-t206 + t236) * qJD(2) + 0.2e1 * t240 * qJD(3)) * qJD(1); (t198 - t219) * MDP(15) + (-qJD(1) * t159 + t162 * qJD(3)) * MDP(16) + (qJD(1) * t163 - t161 * qJD(3)) * MDP(18) + (-t132 * qJD(1) + t199) * MDP(19) + t206 * t196 + t236 * (-t195 - t196) + t239 * (t152 * t205 + t166 * t153 + t161 * t159 - t162 * t163); (t130 * t135 - t131 * t136) * MDP(15) + (t135 * qJD(3) - t123) * MDP(16) + (t178 * t152 - t182 * t153) * MDP(17) + (-t136 * qJD(3) + t124 + 0.2e1 * t211) * MDP(18) + (t122 * t178 + t123 * t182 - t127 * t135 + t213 * t128 - t132 * t137) * MDP(19) + ((t131 - t135) * MDP(14) - t132 * MDP(16) + (t128 - t135) * MDP(17) + t137 * MDP(18)) * t163 + ((-t130 + t136) * MDP(14) - t137 * MDP(16) + (t127 - t213) * MDP(17) - t132 * MDP(18)) * t159 + ((t152 * t190 + t153 * t191) * MDP(14) + (-t123 * t191 + t124 * t190 - t193 * t219) * MDP(15)) * pkin(3) - t240 * t196; (t130 * t163 + t131 * t159 + t226) * MDP(15) + t171 * MDP(16) + (t128 * t159 + (-qJD(5) - t127) * t163 + t203) * MDP(19) + ((t163 - t209) * MDP(16) + (t159 - t201) * MDP(18)) * qJD(3) + t239 * (-t159 ^ 2 - t158); t163 * t159 * MDP(16) + (-t158 - t195) * MDP(18) + (t132 * t163 + t123) * MDP(19) + ((t159 + t201) * MDP(17) - t128 * MDP(19)) * qJD(3);];
tauc = t1;
