% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:07
% EndTime: 2019-12-31 17:17:10
% DurationCPUTime: 0.88s
% Computational Cost: add. (846->159), mult. (2171->213), div. (0->0), fcn. (1332->4), ass. (0->77)
t214 = -MDP(16) - MDP(18);
t172 = qJD(2) + qJD(3);
t176 = sin(qJ(2));
t178 = cos(qJ(2));
t216 = (t176 ^ 2 - t178 ^ 2) * MDP(5);
t212 = -pkin(6) - pkin(5);
t159 = t212 * t176;
t154 = qJD(1) * t159;
t177 = cos(qJ(3));
t198 = qJD(3) * t177;
t160 = t212 * t178;
t156 = qJD(1) * t160;
t175 = sin(qJ(3));
t209 = t175 * t156;
t215 = -pkin(2) * t198 + t177 * t154 + t209;
t153 = t175 * t178 + t177 * t176;
t201 = qJD(1) * t153;
t213 = t201 ^ 2;
t211 = qJD(2) * pkin(2);
t206 = t177 * t178;
t189 = qJD(1) * t206;
t200 = qJD(1) * t176;
t190 = t175 * t200;
t145 = -t189 + t190;
t170 = -t178 * pkin(2) - pkin(1);
t158 = t170 * qJD(1);
t121 = t145 * pkin(3) - qJ(4) * t201 + t158;
t210 = t121 * t201;
t208 = t175 * t176;
t207 = t177 * t156;
t205 = t178 * MDP(4);
t204 = qJD(4) - t215;
t203 = t172 * t189;
t199 = qJD(3) * t175;
t197 = t178 * MDP(10);
t149 = t154 + t211;
t129 = t177 * t149 + t209;
t196 = qJD(4) - t129;
t195 = 0.2e1 * qJD(1);
t193 = pkin(2) * t200;
t119 = t203 + (t145 - t190) * t172;
t192 = t119 * MDP(13) + (-t145 ^ 2 + t213) * MDP(12);
t191 = qJD(2) * t212;
t187 = qJD(1) * t191;
t150 = t176 * t187;
t151 = t178 * t187;
t188 = t149 * t198 + t177 * t150 + t175 * t151 + t156 * t199;
t114 = t149 * t199 + t175 * t150 - t177 * t151 - t156 * t198;
t186 = pkin(2) * t199 - t175 * t154 + t207;
t185 = t172 * t208;
t127 = pkin(3) * t201 + t145 * qJ(4);
t171 = t172 * qJD(4);
t113 = t171 + t188;
t130 = t175 * t149 - t207;
t184 = t177 * t159 + t175 * t160;
t136 = t175 * t159 - t177 * t160;
t183 = t129 * t172 - t188;
t182 = t130 * t172 - t114;
t181 = MDP(11) * t201 + t158 * MDP(17) - t121 * MDP(20);
t134 = t172 * t153;
t169 = -t177 * pkin(2) - pkin(3);
t167 = t175 * pkin(2) + qJ(4);
t157 = t178 * t191;
t155 = t176 * t191;
t152 = -t206 + t208;
t133 = -qJD(2) * t206 - t178 * t198 + t185;
t128 = t152 * pkin(3) - t153 * qJ(4) + t170;
t126 = t134 * qJD(1);
t125 = qJD(1) * t185 - t203;
t124 = t172 * qJ(4) + t130;
t123 = t127 + t193;
t122 = -t172 * pkin(3) + t196;
t116 = t136 * qJD(3) + t175 * t155 - t177 * t157;
t115 = t184 * qJD(3) + t177 * t155 + t175 * t157;
t112 = t134 * pkin(3) + t133 * qJ(4) - t153 * qJD(4) + t176 * t211;
t111 = t126 * pkin(3) + t125 * qJ(4) + qJD(2) * t193 - qJD(4) * t201;
t1 = [(-t125 * t153 - t133 * t201) * MDP(11) + (t125 * t152 - t153 * t126 + t133 * t145 - t134 * t201) * MDP(12) + (t170 * t126 + t158 * t134) * MDP(16) + (-t170 * t125 - t158 * t133) * MDP(17) + (t111 * t152 + t112 * t145 + t121 * t134 + t128 * t126) * MDP(18) + (-t113 * t152 + t114 * t153 - t115 * t145 + t116 * t201 - t122 * t133 - t124 * t134 + t125 * t184 - t136 * t126) * MDP(19) + (-t111 * t153 - t112 * t201 + t121 * t133 + t128 * t125) * MDP(20) + (t111 * t128 + t121 * t112 + t113 * t136 - t114 * t184 + t124 * t115 + t122 * t116) * MDP(21) + (-t133 * MDP(13) - t134 * MDP(14) + t214 * t116 + (-MDP(17) + MDP(20)) * t115) * t172 + ((-pkin(1) * t197 - t216) * t195 + ((-pkin(1) * MDP(9) + t205) * t195 + ((qJD(1) * t152 + t145) * MDP(16) + 0.2e1 * t201 * MDP(17)) * pkin(2)) * t176 + (t178 * MDP(6) - t176 * MDP(7) + (MDP(10) * t176 - MDP(9) * t178) * pkin(5)) * qJD(2)) * qJD(2); -t188 * MDP(17) + (-t169 * t125 - t167 * t126) * MDP(19) + t113 * MDP(20) + (t113 * t167 - t121 * t123 + t186 * t122 + t204 * t124) * MDP(21) + (-t158 * MDP(16) - MDP(17) * t193 - t121 * MDP(18) + (t124 + t186) * MDP(19) + t123 * MDP(20)) * t201 + (-MDP(16) * t193 - t123 * MDP(18) + (t122 - t204) * MDP(19) + t181) * t145 + (MDP(17) * t215 + t204 * MDP(20) + t186 * t214) * t172 + (-t176 * t205 + t216 + (t176 * MDP(9) + t197) * pkin(1)) * qJD(1) ^ 2 + t192 + (t169 * MDP(21) + t214) * t114; (-t158 * t201 + t182) * MDP(16) + t183 * MDP(17) + (t182 - t210) * MDP(18) + (pkin(3) * t125 - t126 * qJ(4) - (-t124 + t130) * t201) * MDP(19) + (t127 * t201 + 0.2e1 * t171 - t183) * MDP(20) + (-t114 * pkin(3) + t113 * qJ(4) - t121 * t127 - t122 * t130 + t196 * t124) * MDP(21) + (-t127 * MDP(18) + (t122 - t196) * MDP(19) + t181) * t145 + t192; t201 * t145 * MDP(18) + t119 * MDP(19) + (-t172 ^ 2 - t213) * MDP(20) + (-t124 * t172 + t114 + t210) * MDP(21);];
tauc = t1;
