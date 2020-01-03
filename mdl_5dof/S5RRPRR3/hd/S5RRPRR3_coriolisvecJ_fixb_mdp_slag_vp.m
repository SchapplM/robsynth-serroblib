% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:36
% EndTime: 2020-01-03 12:00:37
% DurationCPUTime: 0.49s
% Computational Cost: add. (707->105), mult. (1505->158), div. (0->0), fcn. (862->8), ass. (0->73)
t140 = sin(pkin(9));
t147 = cos(qJ(2));
t141 = cos(pkin(9));
t144 = sin(qJ(2));
t185 = t141 * t144;
t153 = pkin(1) * (-t140 * t147 - t185);
t123 = qJD(1) * t153;
t186 = t140 * t144;
t152 = pkin(1) * (t141 * t147 - t186);
t125 = qJD(1) * t152;
t143 = sin(qJ(4));
t146 = cos(qJ(4));
t132 = t141 * pkin(2) + pkin(3);
t192 = pkin(2) * t140;
t155 = t143 * t132 + t146 * t192;
t196 = -t155 * qJD(4) - t146 * t123 + t143 * t125;
t142 = sin(qJ(5));
t145 = cos(qJ(5));
t174 = t145 * MDP(17);
t195 = t142 * MDP(16) + t174;
t194 = t144 * MDP(5) + t147 * MDP(6);
t193 = (t142 ^ 2 - t145 ^ 2) * MDP(12);
t137 = qJD(1) + qJD(2);
t190 = pkin(1) * qJD(1);
t130 = t137 * pkin(2) + t147 * t190;
t168 = t144 * t190;
t112 = t141 * t130 - t140 * t168;
t110 = t137 * pkin(3) + t112;
t126 = qJD(2) * t152;
t117 = qJD(1) * t126;
t113 = t140 * t130 + t141 * t168;
t124 = qJD(2) * t153;
t116 = qJD(1) * t124;
t178 = qJD(4) * t143;
t165 = -t113 * t178 + t143 * t116;
t151 = -(qJD(4) * t110 + t117) * t146 - t165;
t105 = t143 * t110 + t146 * t113;
t164 = -t146 * t116 + t143 * t117;
t97 = t105 * qJD(4) + t164;
t136 = qJD(4) + t137;
t191 = t136 * pkin(4);
t134 = t147 * pkin(1) + pkin(2);
t162 = -pkin(1) * t186 + t141 * t134;
t120 = pkin(3) + t162;
t127 = pkin(1) * t185 + t140 * t134;
t156 = t143 * t120 + t146 * t127;
t189 = (t156 * qJD(4) - t146 * t124 + t143 * t126) * t136;
t104 = t146 * t110 - t143 * t113;
t102 = -t104 - t191;
t173 = t145 * qJD(5);
t188 = t102 * t173 + t97 * t142;
t187 = t105 * t136;
t148 = qJD(5) ^ 2;
t183 = t145 * t148;
t179 = MDP(11) * t145;
t177 = qJD(4) * t146;
t175 = t142 * qJD(5);
t172 = t148 * MDP(14);
t167 = MDP(13) * t183 + 0.2e1 * (t179 * t142 - t193) * qJD(5) * t136;
t166 = -t102 * t136 + t151;
t161 = pkin(8) * t148 - t187;
t160 = qJD(5) * (t104 - t191);
t159 = (pkin(8) + t156) * t148 + t189;
t157 = t146 * t120 - t143 * t127;
t98 = t157 * qJD(4) + t143 * t124 + t146 * t126;
t158 = qJD(5) * ((-pkin(4) - t157) * t136 - t98);
t154 = t146 * t132 - t143 * t192;
t100 = t102 * t175;
t150 = t100 * MDP(16) + t188 * MDP(17) + t167;
t135 = t136 ^ 2;
t122 = pkin(8) + t155;
t121 = -pkin(4) - t154;
t1 = [(t112 * t124 + t113 * t126 + t116 * t162 + t117 * t127) * MDP(7) + (-t97 - t189) * MDP(9) + (-t98 * t136 + t151) * MDP(10) + ((-t159 - t97) * MDP(16) + MDP(17) * t158) * t145 + (MDP(16) * t158 + t159 * MDP(17) - t172) * t142 + t150 + t194 * pkin(1) * qJD(2) * (-qJD(1) - t137); (-t112 * t123 - t113 * t125 + (t116 * t141 + t117 * t140) * pkin(2)) * MDP(7) + (-t110 * t178 - t113 * t177 - t164) * MDP(9) + (-t110 * t177 - t146 * t117 - t165) * MDP(10) - t142 * t172 + (-t122 * t183 - t97 * t145 + t100) * MDP(16) + (t148 * t142 * t122 + t188) * MDP(17) + (t196 * MDP(9) + (t121 * t175 + t196 * t145) * MDP(16) + (t121 * t173 - t196 * t142) * MDP(17)) * t136 + t167 + (MDP(10) * t136 + t195 * qJD(5)) * (-t154 * qJD(4) + t143 * t123 + t146 * t125) + t194 * t190 * (-qJD(2) + t137); -t195 * t148; (-t97 + t187) * MDP(9) + (t104 * t136 + t151) * MDP(10) + ((-t161 - t97) * MDP(16) + MDP(17) * t160) * t145 + (MDP(16) * t160 + t161 * MDP(17) - t172) * t142 + t150; t166 * t174 + t135 * t193 + (t166 * MDP(16) - t135 * t179) * t142;];
tauc = t1;
