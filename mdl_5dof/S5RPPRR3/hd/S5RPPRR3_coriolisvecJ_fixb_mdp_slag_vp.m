% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:28:21
% EndTime: 2020-01-03 11:28:29
% DurationCPUTime: 1.03s
% Computational Cost: add. (824->144), mult. (2069->211), div. (0->0), fcn. (1564->8), ass. (0->75)
t172 = sin(pkin(9));
t177 = sin(qJ(4));
t174 = cos(pkin(9));
t179 = cos(qJ(4));
t207 = t174 * t179;
t187 = t172 * t177 - t207;
t153 = t187 * qJD(1);
t178 = cos(qJ(5));
t147 = t178 * t153;
t200 = t179 * qJD(4);
t204 = qJD(1) * t174;
t164 = t200 * t204;
t203 = qJD(1) * t177;
t198 = t172 * t203;
t150 = -qJD(4) * t198 + t164;
t208 = t172 * t179;
t161 = t174 * t177 + t208;
t156 = t161 * qJD(4);
t151 = qJD(1) * t156;
t154 = t161 * qJD(1);
t176 = sin(qJ(5));
t202 = qJD(5) * t176;
t116 = -qJD(5) * t147 + t178 * t150 - t176 * t151 - t154 * t202;
t190 = -t153 * t176 + t178 * t154;
t117 = qJD(5) * t190 + t150 * t176 + t178 * t151;
t129 = t154 * t176 + t147;
t171 = qJD(4) + qJD(5);
t210 = t129 * t171;
t211 = t190 * t171;
t221 = t129 * t190 * MDP(16) + (-t117 + t211) * MDP(19) + (-t129 ^ 2 + t190 ^ 2) * MDP(17) + (t116 + t210) * MDP(18);
t185 = t161 * qJD(3);
t184 = qJD(1) * t185;
t166 = sin(pkin(8)) * pkin(1) + qJ(3);
t163 = qJD(1) * t166;
t168 = t174 * qJD(2);
t138 = t168 + (-pkin(6) * qJD(1) - t163) * t172;
t146 = t172 * qJD(2) + t174 * t163;
t139 = pkin(6) * t204 + t146;
t192 = -t138 * t177 - t139 * t179;
t115 = -pkin(7) * t150 + qJD(4) * t192 - t184;
t121 = -pkin(7) * t153 - t192;
t162 = -cos(pkin(8)) * pkin(1) - pkin(3) * t174 - pkin(2);
t152 = qJD(1) * t162 + qJD(3);
t133 = pkin(4) * t153 + t152;
t220 = t133 * t129 + t121 * t202 + (-t121 * t171 - t115) * t176;
t217 = t179 * t138 - t139 * t177;
t216 = qJD(3) * t153;
t215 = qJD(5) - t171;
t114 = -pkin(7) * t151 + t217 * qJD(4) - t216;
t214 = -t176 * t114 + t178 * t115 - t133 * t190;
t213 = pkin(4) * t154;
t212 = pkin(6) + t166;
t206 = t178 * t121;
t205 = t172 ^ 2 + t174 ^ 2;
t120 = -pkin(7) * t154 + t217;
t119 = qJD(4) * pkin(4) + t120;
t197 = -pkin(4) * t171 - t119;
t194 = qJD(1) * t205;
t191 = (-t163 * t172 + t168) * t172 - t146 * t174;
t157 = t212 * t172;
t158 = t212 * t174;
t189 = t157 * t177 - t158 * t179;
t188 = -t161 * t176 - t178 * t187;
t135 = t161 * t178 - t176 * t187;
t183 = -t157 * t200 + qJD(3) * t207 + (-qJD(3) * t172 - qJD(4) * t158) * t177;
t181 = qJD(4) * t189 - t185;
t155 = t187 * qJD(4);
t137 = pkin(4) * t187 + t162;
t127 = -pkin(7) * t187 - t189;
t126 = -pkin(7) * t161 - t157 * t179 - t177 * t158;
t125 = pkin(7) * t155 + t181;
t124 = -pkin(7) * t156 + t183;
t123 = qJD(5) * t135 - t155 * t176 + t178 * t156;
t122 = qJD(5) * t188 - t155 * t178 - t156 * t176;
t1 = [(t150 * t161 - t154 * t155) * MDP(9) + (-t150 * t187 - t151 * t161 + t153 * t155 - t154 * t156) * MDP(10) + (t162 * t151 + t152 * t156) * MDP(14) + (t162 * t150 - t152 * t155) * MDP(15) + (t116 * t135 + t122 * t190) * MDP(16) + (t116 * t188 - t117 * t135 - t122 * t129 - t123 * t190) * MDP(17) + (t137 * t117 + t133 * t123 + (t129 * t156 - t151 * t188) * pkin(4)) * MDP(21) + (t137 * t116 + t133 * t122 + (t135 * t151 + t156 * t190) * pkin(4)) * MDP(22) + (t122 * MDP(18) - t123 * MDP(19) + (-t124 * t176 + t125 * t178 + (-t126 * t176 - t127 * t178) * qJD(5)) * MDP(21) + (-t124 * t178 - t125 * t176 - (t126 * t178 - t127 * t176) * qJD(5)) * MDP(22)) * t171 + (0.2e1 * MDP(7) * t194 + (t166 * t194 - t191) * MDP(8)) * qJD(3) + (-t155 * MDP(11) - t156 * MDP(12) + MDP(14) * t181 - MDP(15) * t183) * qJD(4); (-t123 * MDP(21) - t122 * MDP(22)) * t171 + (-t156 * MDP(14) + t155 * MDP(15)) * qJD(4); t164 * MDP(15) + (t117 + t211) * MDP(21) + (t116 - t210) * MDP(22) - t205 * MDP(7) * qJD(1) ^ 2 + t191 * MDP(8) * qJD(1) + ((qJD(1) * t208 + t174 * t203 + t154) * MDP(14) + (-t153 - t198) * MDP(15)) * qJD(4); t154 * t153 * MDP(9) + (-t153 ^ 2 + t154 ^ 2) * MDP(10) + (t164 + (t153 - t198) * qJD(4)) * MDP(11) + (-t152 * t154 - t184) * MDP(14) + (t152 * t153 + t216) * MDP(15) + (-t129 * t213 - (-t120 * t176 - t206) * t171 + (t176 * t197 - t206) * qJD(5) + t214) * MDP(21) + (-t190 * t213 + (qJD(5) * t197 + t120 * t171 - t114) * t178 + t220) * MDP(22) + t221; (t215 * (-t119 * t176 - t206) + t214) * MDP(21) + ((-t215 * t119 - t114) * t178 + t220) * MDP(22) + t221;];
tauc = t1;
