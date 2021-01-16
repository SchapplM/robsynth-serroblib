% Calculate Coriolis joint torque vector for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:14:43
% EndTime: 2021-01-15 16:14:46
% DurationCPUTime: 0.73s
% Computational Cost: add. (581->153), mult. (1032->212), div. (0->0), fcn. (474->4), ass. (0->83)
t141 = sin(qJ(4));
t143 = cos(qJ(4));
t138 = qJD(2) + qJD(3);
t142 = sin(qJ(3));
t185 = pkin(2) * qJD(2);
t163 = t142 * t185;
t120 = pkin(7) * t138 + t163;
t157 = qJ(5) * t138 + t120;
t150 = t157 * t143;
t103 = qJD(1) * t141 + t150;
t191 = qJD(4) * t103;
t139 = t141 ^ 2;
t140 = t143 ^ 2;
t190 = (t139 - t140) * MDP(9);
t144 = cos(qJ(3));
t189 = pkin(2) * t144;
t188 = pkin(4) * t139;
t187 = pkin(4) * t143;
t186 = -qJ(5) - pkin(7);
t184 = pkin(2) * qJD(3);
t183 = qJD(4) * pkin(4);
t135 = t143 * qJD(1);
t102 = -t157 * t141 + t135;
t99 = t102 + t183;
t182 = -t102 + t99;
t162 = t144 * t185;
t121 = -pkin(3) * t138 - t162;
t181 = t121 * t138;
t137 = t138 ^ 2;
t180 = t137 * t143;
t179 = t138 * t141;
t178 = t138 * t143;
t145 = qJD(4) ^ 2;
t177 = t143 * t145;
t176 = t144 * MDP(7);
t131 = pkin(2) * t142 + pkin(7);
t175 = -qJ(5) - t131;
t133 = -pkin(3) - t187;
t107 = t133 * t138 + qJD(5) - t162;
t160 = qJD(2) * t184;
t129 = t142 * t160;
t168 = qJD(4) * t141;
t159 = t138 * t168;
t111 = pkin(4) * t159 + t129;
t167 = qJD(4) * t143;
t174 = t107 * t167 + t111 * t141;
t173 = t121 * t167 + t141 * t129;
t153 = qJD(4) * t162;
t172 = t141 * t153 + t163 * t178;
t170 = MDP(16) * t138;
t166 = qJD(5) * t138;
t165 = t143 * MDP(13);
t164 = -MDP(14) - MDP(16);
t161 = t144 * t184;
t158 = qJD(4) * t186;
t156 = t107 * t168 - t111 * t143;
t155 = qJD(4) * t175;
t154 = (-t139 - t140) * MDP(17);
t152 = t144 * t160;
t151 = t103 * t143 - t141 * t99;
t149 = 0.2e1 * t143 * MDP(8) * t159 - 0.2e1 * t138 * qJD(4) * t190 + MDP(10) * t177 - t129 * MDP(6);
t148 = (-pkin(3) - t189) * t138 - t161;
t147 = qJD(4) * qJD(1) + t152;
t146 = (-qJD(5) - t107) * t138 - t147;
t136 = t143 * qJ(5);
t134 = t143 * qJD(5);
t127 = pkin(7) * t143 + t136;
t126 = t186 * t141;
t124 = t143 * t153;
t122 = t133 - t189;
t118 = pkin(4) * t168 + t142 * t184;
t116 = t131 * t143 + t136;
t115 = t175 * t141;
t114 = t120 * t168;
t112 = t121 * t168;
t110 = -qJD(5) * t141 + t143 * t158;
t109 = t141 * t158 + t134;
t101 = (-qJD(5) - t161) * t141 + t143 * t155;
t100 = t141 * t155 + t143 * t161 + t134;
t98 = (-t152 - t166) * t141 - t191;
t97 = -qJ(5) * t159 - t114 + (t147 + t166) * t143;
t96 = t97 * t143;
t1 = [(t151 * qJD(4) + t97 * t141 + t98 * t143) * MDP(18) + (t164 * t143 + (-MDP(13) - MDP(15)) * t141) * t145; (-t131 * t177 + t112) * MDP(13) + t173 * MDP(14) + (-t118 * t178 + t156) * MDP(15) + t174 * MDP(16) + (t100 * t178 + t96) * MDP(17) + (t100 * t103 + t101 * t99 + t107 * t118 + t111 * t122 + t115 * t98 + t116 * t97) * MDP(18) + (t118 * t170 + (-t101 * t138 - t98) * MDP(17) + (t131 * MDP(14) - MDP(11)) * t145) * t141 + ((-qJD(2) - t138) * t176 + (-qJD(2) * t165 + (t141 * MDP(14) - MDP(6) - t165) * t138) * t142) * t184 + (t101 * MDP(15) - t100 * MDP(16) + (t148 * MDP(14) + t122 * t170 + (-t115 * t138 - t99) * MDP(17)) * t143 + (t148 * MDP(13) + t122 * t138 * MDP(15) + (-t116 * t138 - t103) * MDP(17)) * t141) * qJD(4) + t149; (-pkin(7) * t177 + t112 + t172) * MDP(13) + (t124 + t173) * MDP(14) + (qJD(4) * t110 + t156 + t172) * MDP(15) + (-qJD(4) * t109 + t124 + t174) * MDP(16) + (-t99 * t167 + t96) * MDP(17) + (t103 * t109 + t110 * t99 + t111 * t133 + t126 * t98 + t127 * t97) * MDP(18) + ((-t98 - t191) * MDP(17) + t107 * MDP(18) * t183 + (pkin(7) * MDP(14) - MDP(11)) * t145) * t141 + ((-t107 * t142 - t151 * t144) * MDP(18) + (-t142 * t165 - t176) * qJD(3)) * t185 + ((t109 * t143 - t110 * t141) * MDP(17) + ((MDP(7) + t154) * t144 + (t164 * t141 + MDP(6)) * t142) * t185 + (MDP(16) * t188 + (-pkin(3) * MDP(14) + t133 * MDP(16) - t126 * MDP(17)) * t143 + (-pkin(3) * MDP(13) + (t133 - t187) * MDP(15) - t127 * MDP(17)) * t141) * qJD(4)) * t138 + t149; -t141 * MDP(8) * t180 + t137 * t190 + (-t152 - t181) * t141 * MDP(13) + (t114 + (-t120 * t141 + t135) * qJD(4) + (-t147 - t181) * t143) * MDP(14) + ((t103 - t150) * qJD(4) + (pkin(4) * t180 + t146) * t141) * MDP(15) + (-t137 * t188 + t114 + (qJ(5) * t179 + t102) * qJD(4) + t146 * t143) * MDP(16) + (t182 - t183) * MDP(17) * t178 + (t182 * t103 + (-t107 * t179 + t98) * pkin(4)) * MDP(18); t111 * MDP(18) + t137 * t154 + (-t151 * MDP(18) + 0.2e1 * (t141 * MDP(15) + t143 * MDP(16)) * qJD(4)) * t138;];
tauc = t1;
