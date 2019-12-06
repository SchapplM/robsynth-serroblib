% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:31
% EndTime: 2019-12-05 16:04:38
% DurationCPUTime: 1.56s
% Computational Cost: add. (578->199), mult. (1380->296), div. (0->0), fcn. (905->8), ass. (0->96)
t148 = cos(qJ(4));
t209 = MDP(8) * t148;
t145 = sin(qJ(4));
t190 = qJD(2) * t145;
t137 = qJD(5) + t190;
t144 = sin(qJ(5));
t147 = cos(qJ(5));
t162 = t147 * MDP(20) - t144 * MDP(21);
t157 = t162 * t137;
t141 = t148 ^ 2;
t208 = (t145 ^ 2 - t141) * MDP(9);
t150 = -pkin(2) - pkin(7);
t207 = qJD(2) * pkin(2);
t149 = cos(qJ(2));
t142 = sin(pkin(5));
t194 = qJD(1) * t142;
t177 = t149 * t194;
t165 = qJD(3) - t177;
t126 = t150 * qJD(2) + t165;
t143 = cos(pkin(5));
t146 = sin(qJ(2));
t176 = qJD(2) * t142 * t146;
t158 = -qJD(4) * t143 + t176;
t187 = qJD(4) * t145;
t192 = qJD(1) * t148;
t112 = t126 * t187 - t158 * t192;
t206 = t112 * t144;
t205 = t112 * t147;
t181 = t147 * qJD(4);
t138 = qJD(5) * t181;
t184 = qJD(5) * t144;
t175 = t148 * t184;
t156 = -t145 * t181 - t175;
t121 = t156 * qJD(2) + t138;
t204 = t121 * t144;
t189 = qJD(2) * t148;
t173 = t144 * t189;
t130 = t173 - t181;
t203 = t130 * t137;
t172 = t147 * t189;
t188 = qJD(4) * t144;
t132 = t172 + t188;
t202 = t132 * t137;
t201 = t142 * t149;
t200 = t144 * t137;
t199 = t144 * t150;
t198 = t147 * t137;
t197 = t147 * t150;
t151 = qJD(4) ^ 2;
t152 = qJD(2) ^ 2;
t195 = -t151 - t152;
t193 = qJD(1) * t145;
t191 = qJD(2) * qJ(3);
t186 = qJD(4) * t148;
t185 = qJD(4) * t150;
t183 = qJD(5) * t147;
t159 = -t126 * t148 + t143 * t193;
t114 = -qJD(4) * pkin(4) + t159;
t182 = t114 * qJD(5);
t180 = t148 * MDP(19);
t179 = qJD(2) * qJD(4);
t178 = t146 * t194;
t174 = t137 * t183;
t171 = t145 * t179;
t170 = qJD(4) * t180;
t119 = t126 * t145 + t143 * t192;
t115 = qJD(4) * pkin(8) + t119;
t169 = t137 * t150 + t115;
t168 = pkin(4) * t148 + pkin(8) * t145;
t127 = t168 * qJD(4) + qJD(3);
t167 = -t127 + t177;
t134 = t178 + t191;
t166 = -t134 + t178;
t135 = pkin(4) * t145 - pkin(8) * t148 + qJ(3);
t123 = t135 * qJD(2) + t178;
t110 = t115 * t147 + t123 * t144;
t164 = t115 * t144 - t123 * t147;
t163 = qJD(2) * t141 - t137 * t145;
t161 = -t144 * MDP(20) - t147 * MDP(21);
t160 = -pkin(8) * t186 + t114 * t145;
t124 = t143 * t145 + t148 * t201;
t125 = t143 * t148 - t145 * t201;
t155 = t166 - t191;
t128 = (qJD(3) + t177) * qJD(2);
t154 = t165 * qJD(2) - t150 * t151 + t128;
t111 = t126 * t186 + t158 * t193;
t153 = -qJD(4) * t114 - qJD(5) * t123 + t137 * t178 - t111;
t136 = t144 * t171;
t133 = t168 * qJD(2);
t129 = t165 - t207;
t122 = t132 * qJD(5) - t136;
t120 = (t127 + t177) * qJD(2);
t117 = t125 * qJD(4) - t148 * t176;
t116 = -t124 * qJD(4) + t145 * t176;
t113 = t147 * t120;
t1 = [((-t116 * t144 - t125 * t183) * t137 + t117 * t130 + t124 * t122) * MDP(20) + (-(t116 * t147 - t125 * t184) * t137 + t117 * t132 + t124 * t121) * MDP(21) + (t161 * t125 * t189 - t117 * MDP(13) - t116 * MDP(14)) * qJD(4) + (((MDP(7) * t134 + t157) * qJD(2) + (t145 * MDP(13) + t148 * MDP(14) - MDP(4) + MDP(6)) * t152) * t149 + (t128 * MDP(7) + (-MDP(3) + MDP(5)) * t152 + t161 * t137 * qJD(5) + ((t129 - t177) * MDP(7) + (-MDP(14) * t145 + (MDP(13) + t162) * t148) * qJD(4)) * qJD(2)) * t146) * t142; 0.2e1 * qJD(2) * qJD(3) * MDP(6) + (qJ(3) * t128 + qJD(3) * t134 + (-t134 * t149 + (-t129 - t207) * t146) * t194) * MDP(7) - 0.2e1 * t171 * t209 + 0.2e1 * t179 * t208 + (t154 * t145 - t155 * t186) * MDP(13) + (t154 * t148 + t155 * t187) * MDP(14) + (t121 * t147 * t148 + t156 * t132) * MDP(15) + ((t130 * t147 + t132 * t144) * t187 + (-t204 - t122 * t147 + (t130 * t144 - t132 * t147) * qJD(5)) * t148) * MDP(16) + (-t137 * t175 + t121 * t145 + (t132 * t148 + t163 * t147) * qJD(4)) * MDP(17) + (-t148 * t174 - t122 * t145 + (-t130 * t148 - t163 * t144) * qJD(4)) * MDP(18) + (t137 + t190) * t170 + ((-t135 * t184 - t167 * t147) * t137 + (t130 * t185 + t153 * t144 - t169 * t183 + t113) * t145 + (t130 * t178 + t147 * t182 + t206 - t150 * t122 + (-t137 * t199 + (t135 * t147 - t145 * t199) * qJD(2) - t164) * qJD(4)) * t148) * MDP(20) + ((-t135 * t183 + t167 * t144) * t137 + (t132 * t185 + (t169 * qJD(5) - t120) * t144 + t153 * t147) * t145 + (t132 * t178 - t144 * t182 + t205 - t150 * t121 + (-t137 * t197 - (t135 * t144 + t145 * t197) * qJD(2) - t110) * qJD(4)) * t148) * MDP(21) + (-t145 * MDP(10) - t148 * MDP(11)) * t151; -t152 * MDP(6) + (t166 * MDP(7) - t157) * qJD(2) + (t195 * MDP(14) + (-t137 * t188 - t122) * MDP(20) + (-t137 * t181 - t121) * MDP(21)) * t148 + (t195 * MDP(13) - qJD(5) * t157 + ((t130 - t173) * MDP(20) + (t132 - t172) * MDP(21)) * qJD(4)) * t145; (t132 * t198 + t204) * MDP(15) + ((t121 - t203) * t147 + (-t122 - t202) * t144) * MDP(16) + (t174 + (t145 * t198 + (-t132 + t188) * t148) * qJD(2)) * MDP(17) + (-t137 * t184 + (-t145 * t200 + (t130 + t181) * t148) * qJD(2)) * MDP(18) - t137 * qJD(2) * t180 + (-pkin(4) * t122 - t205 - (t133 * t147 + t144 * t159) * t137 - t119 * t130 + (-pkin(8) * t198 + t114 * t144) * qJD(5) + (t160 * t144 + t148 * t164) * qJD(2)) * MDP(20) + (-pkin(4) * t121 + t206 + (t144 * t133 - t147 * t159) * t137 - t119 * t132 + (pkin(8) * t200 + t114 * t147) * qJD(5) + (t110 * t148 + t160 * t147) * qJD(2)) * MDP(21) + (t145 * t209 - t208) * t152 + (MDP(13) * t189 - t190 * MDP(14)) * t166; t132 * t130 * MDP(15) + (-t130 ^ 2 + t132 ^ 2) * MDP(16) + (-t147 * t171 + t138 + t203) * MDP(17) + (t136 + t202) * MDP(18) + qJD(2) * t170 + (t110 * t137 - t144 * t111 - t114 * t132 + t113) * MDP(20) + (-t147 * t111 + t114 * t130 - t144 * t120 - t137 * t164) * MDP(21) + (-MDP(17) * t173 - t132 * MDP(18) - t110 * MDP(20) + t164 * MDP(21)) * qJD(5);];
tauc = t1;
