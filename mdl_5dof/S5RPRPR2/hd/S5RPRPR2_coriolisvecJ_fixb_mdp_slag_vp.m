% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR2
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
%   see S5RPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:05
% EndTime: 2020-01-03 11:34:09
% DurationCPUTime: 0.70s
% Computational Cost: add. (601->102), mult. (1237->151), div. (0->0), fcn. (750->8), ass. (0->61)
t155 = sin(pkin(9));
t157 = cos(pkin(9));
t178 = t155 ^ 2 + t157 ^ 2;
t154 = qJD(1) + qJD(3);
t146 = cos(pkin(8)) * pkin(1) + pkin(2);
t143 = t146 * qJD(1);
t160 = sin(qJ(3));
t192 = pkin(1) * sin(pkin(8));
t174 = qJD(3) * t192;
t169 = qJD(1) * t174;
t162 = cos(qJ(3));
t177 = qJD(3) * t162;
t167 = t143 * t177 - t160 * t169;
t111 = t154 * qJD(4) + t167;
t198 = t178 * t111;
t197 = t157 * MDP(8) + MDP(6);
t159 = sin(qJ(5));
t161 = cos(qJ(5));
t135 = t155 * t161 + t157 * t159;
t126 = t135 * t154;
t196 = t154 * t178;
t175 = qJD(1) * t192;
t121 = t143 * t162 - t160 * t175;
t193 = t121 - qJD(4);
t191 = t157 * pkin(4);
t147 = -pkin(3) - t191;
t107 = t147 * t154 - t193;
t181 = t160 * t143;
t117 = qJD(3) * t181 + t162 * t169;
t182 = t157 * t161;
t186 = t155 * t159;
t134 = -t182 + t186;
t131 = t134 * qJD(5);
t190 = -t107 * t131 + t117 * t135;
t187 = t154 * t155;
t132 = t135 * qJD(5);
t180 = t107 * t132 + t117 * t134;
t173 = t154 * t186;
t172 = t154 * t182;
t136 = qJD(5) * t172;
t118 = -qJD(5) * t173 + t136;
t119 = t154 * t132;
t124 = -t172 + t173;
t170 = (-t118 * t134 - t119 * t135 + t124 * t131 - t126 * t132) * MDP(13) + (t118 * t135 - t126 * t131) * MDP(12) + (-t131 * MDP(14) - t132 * MDP(15)) * qJD(5);
t168 = -t146 * t162 + t160 * t192 - pkin(3);
t122 = t162 * t175 + t181;
t165 = t178 * (t154 * qJ(4) + t122);
t164 = t146 * t177 - t160 * t174;
t163 = t160 * t146 + t162 * t192;
t150 = t157 * pkin(7);
t139 = qJ(4) * t157 + t150;
t138 = (-pkin(7) - qJ(4)) * t155;
t130 = qJ(4) + t163;
t127 = t163 * qJD(3);
t123 = qJD(4) + t164;
t120 = t168 - t191;
t116 = t130 * t157 + t150;
t115 = (-pkin(7) - t130) * t155;
t113 = t117 * t155;
t112 = -t154 * pkin(3) - t193;
t1 = [(-t164 * t154 - t167) * MDP(7) + (t127 * t187 + t113) * MDP(9) + (t123 * t196 + t198) * MDP(10) + (t112 * t127 + t117 * t168 + t165 * t123 + t130 * t198) * MDP(11) + (t120 * t119 + t127 * t124 + ((-t115 * t159 - t116 * t161) * qJD(5) - t135 * t123) * qJD(5) + t180) * MDP(17) + (t120 * t118 + t127 * t126 + ((-t115 * t161 + t116 * t159) * qJD(5) + t134 * t123) * qJD(5) + t190) * MDP(18) + t170 + t197 * (-t127 * t154 - t117); (-MDP(17) * t132 + MDP(18) * t131) * qJD(5); (t121 * t154 - t167) * MDP(7) + (-t122 * t187 + t113) * MDP(9) + (-t193 * t196 + t198) * MDP(10) + (-pkin(3) * t117 + qJ(4) * t198 - t112 * t122 - t193 * t165) * MDP(11) + (t147 * t119 - t122 * t124 + ((-t138 * t159 - t139 * t161) * qJD(5) + t193 * t135) * qJD(5) + t180) * MDP(17) + (t147 * t118 - t122 * t126 + ((-t138 * t161 + t139 * t159) * qJD(5) - t193 * t134) * qJD(5) + t190) * MDP(18) + t170 + t197 * (t122 * t154 - t117); (-t165 * t154 + t117) * MDP(11) + t136 * MDP(18) - t178 * MDP(10) * t154 ^ 2 + (0.2e1 * t126 * MDP(17) + (-t124 - t173) * MDP(18)) * qJD(5); t126 * t124 * MDP(12) + (-t124 ^ 2 + t126 ^ 2) * MDP(13) + (t136 + (t124 - t173) * qJD(5)) * MDP(14) + (-t107 * t126 - t135 * t111) * MDP(17) + (t107 * t124 + t134 * t111) * MDP(18);];
tauc = t1;
