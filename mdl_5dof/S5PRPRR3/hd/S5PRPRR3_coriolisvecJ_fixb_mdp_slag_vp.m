% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:36
% EndTime: 2019-12-05 15:47:41
% DurationCPUTime: 0.94s
% Computational Cost: add. (529->128), mult. (1303->199), div. (0->0), fcn. (934->8), ass. (0->80)
t165 = cos(qJ(2));
t195 = qJD(1) * t165;
t146 = qJD(2) * pkin(2) + t195;
t158 = sin(pkin(9));
t162 = sin(qJ(2));
t196 = qJD(1) * t162;
t147 = t158 * t196;
t159 = cos(pkin(9));
t123 = t146 * t159 - t147;
t164 = cos(qJ(4));
t186 = -pkin(4) * t164 - pkin(3);
t119 = t186 * qJD(2) - t123;
t132 = t159 * t195 - t147;
t211 = t119 + t132;
t160 = sin(qJ(5));
t161 = sin(qJ(4));
t163 = cos(qJ(5));
t141 = t160 * t164 + t161 * t163;
t155 = qJD(4) + qJD(5);
t209 = t155 * t141;
t116 = t209 * qJD(2);
t206 = qJD(5) - t155;
t210 = MDP(6) * t161;
t140 = t160 * t161 - t163 * t164;
t171 = t155 * t140;
t208 = (t161 ^ 2 - t164 ^ 2) * MDP(7);
t148 = t159 * t196;
t124 = t158 * t146 + t148;
t180 = t124 + (pkin(6) + pkin(7)) * qJD(2);
t113 = t164 * qJD(3) - t180 * t161;
t207 = MDP(11) * t161 + MDP(12) * t164;
t114 = qJD(3) * t161 + t180 * t164;
t205 = pkin(2) * t159;
t152 = pkin(2) * t158 + pkin(6);
t204 = pkin(7) + t152;
t139 = t158 * t165 + t159 * t162;
t203 = t139 * t155;
t166 = qJD(4) ^ 2;
t202 = t161 * t166;
t201 = t163 * t114;
t200 = t164 * t166;
t194 = qJD(2) * t161;
t193 = qJD(2) * t164;
t191 = qJD(4) * t161;
t190 = qJD(4) * t164;
t189 = qJD(2) * qJD(4);
t188 = pkin(4) * t194;
t187 = pkin(4) * t191;
t185 = t160 * t194;
t184 = t163 * t193;
t183 = t164 * t189;
t112 = qJD(4) * pkin(4) + t113;
t182 = -pkin(4) * t155 - t112;
t181 = qJD(4) * t204;
t130 = t158 * t195 + t148;
t178 = -t130 + t187;
t138 = t158 * t162 - t159 * t165;
t131 = t138 * qJD(2);
t126 = qJD(1) * t131;
t106 = qJD(4) * t113 - t164 * t126;
t107 = -qJD(4) * t114 + t161 * t126;
t135 = -t160 * t193 - t163 * t194;
t175 = -t160 * t106 + t163 * t107 + t119 * t135;
t129 = t139 * qJD(2);
t125 = qJD(1) * t129;
t174 = qJD(2) * t130 - t152 * t166 - t125;
t121 = -qJD(2) * pkin(3) - t123;
t173 = qJD(4) * (qJD(2) * (-pkin(3) - t205) + t121 + t132);
t115 = qJD(5) * t184 - t155 * t185 + t163 * t183;
t133 = -t184 + t185;
t170 = -t135 * t133 * MDP(13) + (t133 * t155 + t115) * MDP(15) + (-t135 * t155 - t116) * MDP(16) + (-t133 ^ 2 + t135 ^ 2) * MDP(14);
t169 = t119 * t133 + (t206 * t114 - t107) * t160;
t167 = qJD(2) ^ 2;
t143 = t186 - t205;
t137 = t204 * t164;
t136 = t204 * t161;
t128 = t164 * t181;
t127 = t161 * t181;
t120 = (t139 * qJD(1) + t187) * qJD(2);
t1 = [(-t123 * t129 - t124 * t131 + t125 * t138 - t126 * t139) * MDP(5) + (t131 * t191 - t139 * t200 + (-t129 * t164 + t138 * t191) * qJD(2)) * MDP(11) + (t131 * t190 + t139 * t202 + (t129 * t161 + t138 * t190) * qJD(2)) * MDP(12) + (t138 * t116 + t129 * t133 + t131 * t209 + t171 * t203) * MDP(18) + (t138 * t115 - t129 * t135 - t131 * t171 + t203 * t209) * MDP(19) + (-t162 * MDP(3) - t165 * MDP(4)) * t167; (t123 * t130 - t124 * t132 + (-t125 * t159 - t126 * t158) * pkin(2)) * MDP(5) + 0.2e1 * t183 * t210 - 0.2e1 * t189 * t208 + MDP(8) * t200 - MDP(9) * t202 + (t161 * t173 + t174 * t164) * MDP(11) + (-t174 * t161 + t164 * t173) * MDP(12) + (t115 * t141 + t135 * t171) * MDP(13) + (-t115 * t140 - t116 * t141 + t133 * t171 + t135 * t209) * MDP(14) + (t143 * t116 + t120 * t140 + t178 * t133 + t211 * t209) * MDP(18) + (t143 * t115 + t120 * t141 - t178 * t135 - t211 * t171) * MDP(19) + (-t171 * MDP(15) - t209 * MDP(16) + (t127 * t160 - t128 * t163 + (t136 * t160 - t137 * t163) * qJD(5)) * MDP(18) + (t127 * t163 + t128 * t160 - (-t136 * t163 - t137 * t160) * qJD(5)) * MDP(19)) * t155; -t207 * t166 + (-MDP(18) * t209 + t171 * MDP(19)) * t155; (-(-t113 * t160 - t201) * t155 - t133 * t188 + (t182 * t160 - t201) * qJD(5) + t175) * MDP(18) + (t135 * t188 + (t182 * qJD(5) + t113 * t155 - t106) * t163 + t169) * MDP(19) + t170 + t207 * (-qJD(2) * t121 + t126) + (-t164 * t210 + t208) * t167; (t175 + t206 * (-t112 * t160 - t201)) * MDP(18) + ((-t112 * t206 - t106) * t163 + t169) * MDP(19) + t170;];
tauc = t1;
