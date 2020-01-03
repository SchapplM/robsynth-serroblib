% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:24
% EndTime: 2020-01-03 11:36:28
% DurationCPUTime: 0.66s
% Computational Cost: add. (631->111), mult. (1331->177), div. (0->0), fcn. (719->8), ass. (0->70)
t143 = qJD(1) + qJD(3);
t138 = cos(pkin(8)) * pkin(1) + pkin(2);
t134 = t138 * qJD(1);
t151 = sin(qJ(3));
t201 = pkin(1) * sin(pkin(8));
t173 = qJD(3) * t201;
t167 = qJD(1) * t173;
t153 = cos(qJ(3));
t179 = qJD(3) * t153;
t162 = t134 * t179 - t151 * t167;
t114 = t143 * qJD(4) + t162;
t146 = sin(pkin(9));
t141 = t146 ^ 2;
t148 = cos(pkin(9));
t182 = t148 ^ 2 + t141;
t208 = t182 * t114;
t207 = t148 * MDP(8) + MDP(6);
t150 = sin(qJ(5));
t152 = cos(qJ(5));
t193 = t141 * t152;
t206 = -t150 * MDP(12) * t193 + (t150 ^ 2 - t152 ^ 2) * MDP(13) * t141;
t202 = -t138 * t153 + t151 * t201;
t174 = qJD(1) * t201;
t122 = t153 * t134 - t151 * t174;
t161 = qJD(4) - t122;
t186 = t151 * t134;
t123 = t153 * t174 + t186;
t118 = t143 * qJ(4) + t123;
t106 = -t148 * qJD(2) + t118 * t146;
t200 = t106 * t146;
t157 = t138 * t179 - t151 * t173;
t124 = qJD(4) + t157;
t198 = t124 * t143;
t195 = t141 * t143;
t194 = t141 * t150;
t192 = t143 * t146;
t191 = t143 * t148;
t189 = t148 * t150;
t188 = t148 * t152;
t135 = -qJD(5) + t191;
t187 = t150 * t135;
t185 = t152 * t135;
t131 = -t148 * pkin(4) - t146 * pkin(7) - pkin(3);
t104 = t131 * t143 + t161;
t107 = qJD(2) * t146 + t118 * t148;
t121 = qJD(3) * t186 + t153 * t167;
t184 = (t114 * t188 + t150 * t121 + (t104 * t152 - t107 * t150) * qJD(5)) * t148 + t114 * t193;
t177 = qJD(5) * t150;
t176 = qJD(5) * t152;
t175 = qJD(5) + t135;
t172 = t106 * t192;
t171 = t143 * t193;
t170 = t146 * t176;
t169 = t146 * t177;
t129 = t169 * t191;
t166 = (t135 * t169 + t129) * MDP(14) + (t135 + t191) * MDP(15) * t170 + 0.2e1 * t206 * qJD(5) * t143;
t160 = -t104 * t150 - t107 * t152;
t103 = t160 * qJD(5) - t114 * t189 + t152 * t121;
t164 = -t103 * t148 + t106 * t170 + t114 * t194;
t159 = t107 * t148 + t200;
t158 = t135 * t148 + t195;
t156 = t151 * t138 + t153 * t201;
t155 = qJ(4) * t176 + t161 * t150;
t140 = t143 ^ 2;
t126 = qJ(4) + t156;
t125 = t156 * qJD(3);
t119 = t131 + t202;
t116 = t121 * t146;
t115 = -t143 * pkin(3) + t161;
t1 = [(-t157 * t143 - t162) * MDP(7) + (t125 * t192 + t116) * MDP(9) + (t182 * t198 + t208) * MDP(10) + (t121 * (-pkin(3) + t202) + t115 * t125 + t159 * t124 + t126 * t208) * MDP(11) + (-(-t124 * t189 + t125 * t152) * t135 + t194 * t198 + (-(-t119 * t150 - t126 * t188) * t135 + t126 * t171) * qJD(5) + t164) * MDP(17) + ((t124 * t188 + t125 * t150) * t135 + t124 * t171 + (t119 * t185 + (-t158 * t126 - t200) * t150) * qJD(5) + t184) * MDP(18) + t166 + t207 * (-t125 * t143 - t121); t129 * MDP(18) + (-MDP(18) * t187 + (t135 - t191) * MDP(17) * t152) * t146 * qJD(5); (t122 * t143 - t162) * MDP(7) + (-t123 * t192 + t116) * MDP(9) + (t143 * t161 * t182 + t208) * MDP(10) + (-pkin(3) * t121 + qJ(4) * t208 - t115 * t123 + t159 * t161) * MDP(11) + ((t152 * t123 + t131 * t177 + t155 * t148) * t135 + t155 * t195 + t164) * MDP(17) + (-t123 * t187 + t161 * t158 * t152 + (t131 * t185 + (-t158 * qJ(4) - t200) * t150) * qJD(5) + t184) * MDP(18) + t166 + t207 * (t123 * t143 - t121); (-t159 * t143 + t121) * MDP(11) - t182 * MDP(10) * t140 + (t150 * MDP(17) + t152 * MDP(18)) * (-t135 ^ 2 - t141 * t140); (t160 * t135 - t152 * t172 + t103) * MDP(17) + ((-t175 * t104 - t148 * t114) * t152 + (t175 * t107 - t121 + t172) * t150) * MDP(18) - (t150 * MDP(14) + t152 * MDP(15)) * t175 * t192 - t206 * t140;];
tauc = t1;
