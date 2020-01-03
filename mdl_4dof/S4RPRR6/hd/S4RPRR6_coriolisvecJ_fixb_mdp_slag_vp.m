% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:43
% EndTime: 2019-12-31 16:52:45
% DurationCPUTime: 0.94s
% Computational Cost: add. (645->131), mult. (1785->191), div. (0->0), fcn. (1354->6), ass. (0->68)
t152 = cos(pkin(7));
t156 = cos(qJ(3));
t179 = qJD(3) * t156;
t145 = t152 * qJD(1) * t179;
t151 = sin(pkin(7));
t154 = sin(qJ(3));
t180 = qJD(1) * t154;
t174 = t151 * t180;
t130 = -qJD(3) * t174 + t145;
t184 = t151 * t156;
t139 = t152 * t154 + t184;
t136 = t139 * qJD(3);
t131 = qJD(1) * t136;
t153 = sin(qJ(4));
t155 = cos(qJ(4));
t183 = t152 * t156;
t165 = t151 * t154 - t183;
t133 = t165 * qJD(1);
t134 = t139 * qJD(1);
t169 = -t133 * t153 + t155 * t134;
t100 = t169 * qJD(4) + t130 * t153 + t155 * t131;
t127 = t155 * t133;
t115 = t134 * t153 + t127;
t150 = qJD(3) + qJD(4);
t186 = t115 * t150;
t187 = t169 * t150;
t178 = qJD(4) * t153;
t99 = -qJD(4) * t127 + t155 * t130 - t153 * t131 - t134 * t178;
t200 = (t99 + t186) * MDP(17) + t115 * t169 * MDP(15) + (-t100 + t187) * MDP(18) + (-t115 ^ 2 + t169 ^ 2) * MDP(16);
t198 = (qJ(2) * MDP(7) + MDP(6)) * (t151 ^ 2 + t152 ^ 2);
t162 = t139 * qJD(2);
t161 = qJD(1) * t162;
t189 = pkin(5) + qJ(2);
t143 = t189 * t151;
t140 = qJD(1) * t143;
t144 = t189 * t152;
t141 = qJD(1) * t144;
t167 = t140 * t154 - t141 * t156;
t105 = -pkin(6) * t130 + t167 * qJD(3) - t161;
t111 = -pkin(6) * t133 - t167;
t147 = -pkin(2) * t152 - pkin(1);
t142 = t147 * qJD(1) + qJD(2);
t121 = pkin(3) * t133 + t142;
t197 = t121 * t115 + t111 * t178 + (-t111 * t150 - t105) * t153;
t193 = -t156 * t140 - t141 * t154;
t192 = qJD(4) - t150;
t104 = -pkin(6) * t131 - qJD(2) * t133 + t193 * qJD(3);
t191 = -t153 * t104 + t155 * t105 - t121 * t169;
t190 = pkin(3) * t134;
t182 = t155 * t111;
t176 = qJD(1) * qJD(2);
t110 = -pkin(6) * t134 + t193;
t109 = qJD(3) * pkin(3) + t110;
t173 = -pkin(3) * t150 - t109;
t168 = -t139 * t153 - t155 * t165;
t120 = t139 * t155 - t153 * t165;
t166 = t143 * t154 - t144 * t156;
t160 = -t143 * t179 + qJD(2) * t183 + (-qJD(2) * t151 - qJD(3) * t144) * t154;
t158 = t166 * qJD(3) - t162;
t135 = t165 * qJD(3);
t125 = pkin(3) * t165 + t147;
t113 = -pkin(6) * t165 - t166;
t112 = -pkin(6) * t139 - t143 * t156 - t154 * t144;
t107 = pkin(6) * t135 + t158;
t106 = -pkin(6) * t136 + t160;
t103 = t120 * qJD(4) - t135 * t153 + t155 * t136;
t102 = t168 * qJD(4) - t135 * t155 - t136 * t153;
t1 = [(t130 * t139 - t134 * t135) * MDP(8) + (-t130 * t165 - t131 * t139 + t133 * t135 - t134 * t136) * MDP(9) + (t147 * t131 + t142 * t136) * MDP(13) + (t147 * t130 - t142 * t135) * MDP(14) + (t102 * t169 + t120 * t99) * MDP(15) + (-t100 * t120 - t102 * t115 - t103 * t169 + t168 * t99) * MDP(16) + (t125 * t100 + t121 * t103 + (t115 * t136 - t131 * t168) * pkin(3)) * MDP(20) + (t125 * t99 + t121 * t102 + (t120 * t131 + t136 * t169) * pkin(3)) * MDP(21) + 0.2e1 * t176 * t198 + (t102 * MDP(17) - t103 * MDP(18) + (-t106 * t153 + t107 * t155 + (-t112 * t153 - t113 * t155) * qJD(4)) * MDP(20) + (-t106 * t155 - t107 * t153 - (t112 * t155 - t113 * t153) * qJD(4)) * MDP(21)) * t150 + (-t135 * MDP(10) - t136 * MDP(11) + t158 * MDP(13) - t160 * MDP(14)) * qJD(3); t145 * MDP(14) + (t100 + t187) * MDP(20) + (t99 - t186) * MDP(21) + ((qJD(1) * t184 + t152 * t180 + t134) * MDP(13) + (-t133 - t174) * MDP(14)) * qJD(3) - qJD(1) ^ 2 * t198; t134 * t133 * MDP(8) + (-t133 ^ 2 + t134 ^ 2) * MDP(9) + (t145 + (t133 - t174) * qJD(3)) * MDP(10) + (-t142 * t134 - t161) * MDP(13) + (t142 * t133 + t165 * t176) * MDP(14) + (-t115 * t190 - (-t110 * t153 - t182) * t150 + (t173 * t153 - t182) * qJD(4) + t191) * MDP(20) + (-t169 * t190 + (t173 * qJD(4) + t110 * t150 - t104) * t155 + t197) * MDP(21) + t200; (t192 * (-t109 * t153 - t182) + t191) * MDP(20) + ((-t192 * t109 - t104) * t155 + t197) * MDP(21) + t200;];
tauc = t1;
