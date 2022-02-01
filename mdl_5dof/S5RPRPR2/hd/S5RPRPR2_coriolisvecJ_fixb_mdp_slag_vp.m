% Calculate Coriolis joint torque vector for
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
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:08
% EndTime: 2022-01-23 09:19:10
% DurationCPUTime: 0.53s
% Computational Cost: add. (587->98), mult. (1198->145), div. (0->0), fcn. (730->8), ass. (0->59)
t152 = sin(pkin(9));
t154 = cos(pkin(9));
t175 = t152 ^ 2 + t154 ^ 2;
t151 = qJD(1) + qJD(3);
t143 = cos(pkin(8)) * pkin(1) + pkin(2);
t140 = t143 * qJD(1);
t157 = sin(qJ(3));
t188 = pkin(1) * sin(pkin(8));
t171 = qJD(3) * t188;
t166 = qJD(1) * t171;
t159 = cos(qJ(3));
t174 = qJD(3) * t159;
t164 = t140 * t174 - t157 * t166;
t109 = t151 * qJD(4) + t164;
t194 = t175 * t109;
t193 = t154 * MDP(8) + MDP(6);
t156 = sin(qJ(5));
t158 = cos(qJ(5));
t132 = t158 * t152 + t156 * t154;
t123 = t132 * t151;
t192 = t151 * t175;
t172 = qJD(1) * t188;
t118 = t159 * t140 - t157 * t172;
t189 = t118 - qJD(4);
t187 = t154 * pkin(4);
t144 = -pkin(3) - t187;
t105 = t144 * t151 - t189;
t179 = t157 * t140;
t114 = qJD(3) * t179 + t159 * t166;
t129 = t132 * qJD(5);
t177 = t158 * t154;
t181 = t156 * t152;
t131 = -t177 + t181;
t186 = t105 * t129 + t114 * t131;
t128 = t131 * qJD(5);
t185 = -t105 * t128 + t114 * t132;
t170 = t151 * t181;
t169 = t151 * t177;
t133 = qJD(5) * t169;
t115 = -qJD(5) * t170 + t133;
t116 = t151 * t129;
t121 = -t169 + t170;
t167 = (-t115 * t131 - t132 * t116 + t128 * t121 - t123 * t129) * MDP(12) + (t115 * t132 - t123 * t128) * MDP(11) + (-t128 * MDP(13) - t129 * MDP(14)) * qJD(5);
t165 = -t159 * t143 + t157 * t188 - pkin(3);
t119 = t159 * t172 + t179;
t162 = t175 * (t151 * qJ(4) + t119);
t161 = t143 * t174 - t157 * t171;
t160 = t157 * t143 + t159 * t188;
t147 = t154 * pkin(7);
t136 = t154 * qJ(4) + t147;
t135 = (-pkin(7) - qJ(4)) * t152;
t127 = qJ(4) + t160;
t124 = t160 * qJD(3);
t120 = qJD(4) + t161;
t117 = t165 - t187;
t113 = t154 * t127 + t147;
t112 = (-pkin(7) - t127) * t152;
t110 = -t151 * pkin(3) - t189;
t1 = [(-t161 * t151 - t164) * MDP(7) + (t120 * t192 + t194) * MDP(9) + (t110 * t124 + t114 * t165 + t162 * t120 + t127 * t194) * MDP(10) + (t117 * t116 + t124 * t121 + ((-t112 * t156 - t113 * t158) * qJD(5) - t132 * t120) * qJD(5) + t186) * MDP(16) + (t117 * t115 + t124 * t123 + ((-t112 * t158 + t113 * t156) * qJD(5) + t131 * t120) * qJD(5) + t185) * MDP(17) + t167 + t193 * (-t124 * t151 - t114); (-MDP(16) * t129 + MDP(17) * t128) * qJD(5); (t118 * t151 - t164) * MDP(7) + (-t189 * t192 + t194) * MDP(9) + (-t114 * pkin(3) + qJ(4) * t194 - t110 * t119 - t189 * t162) * MDP(10) + (t144 * t116 - t119 * t121 + ((-t135 * t156 - t136 * t158) * qJD(5) + t189 * t132) * qJD(5) + t186) * MDP(16) + (t144 * t115 - t119 * t123 + ((-t135 * t158 + t136 * t156) * qJD(5) - t189 * t131) * qJD(5) + t185) * MDP(17) + t167 + t193 * (t119 * t151 - t114); (-t162 * t151 + t114) * MDP(10) + t133 * MDP(17) - t175 * MDP(9) * t151 ^ 2 + (0.2e1 * t123 * MDP(16) + (-t121 - t170) * MDP(17)) * qJD(5); t123 * t121 * MDP(11) + (-t121 ^ 2 + t123 ^ 2) * MDP(12) + (t133 + (t121 - t170) * qJD(5)) * MDP(13) + (-t105 * t123 - t132 * t109) * MDP(16) + (t105 * t121 + t131 * t109) * MDP(17);];
tauc = t1;
