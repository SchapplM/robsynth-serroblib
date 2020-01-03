% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:53
% EndTime: 2019-12-31 17:47:56
% DurationCPUTime: 1.00s
% Computational Cost: add. (472->151), mult. (1255->252), div. (0->0), fcn. (821->6), ass. (0->85)
t215 = pkin(3) + qJ(2);
t213 = MDP(6) + MDP(8);
t164 = cos(pkin(8));
t165 = cos(pkin(7));
t200 = qJD(1) * t165;
t143 = t164 * t200 + qJD(5);
t211 = t143 * t164;
t161 = t165 ^ 2;
t168 = qJD(1) ^ 2;
t210 = t161 * t168;
t163 = sin(pkin(7));
t209 = t163 * t164;
t166 = sin(qJ(5));
t208 = t163 * t166;
t167 = cos(qJ(5));
t207 = t163 * t167;
t206 = t164 * t165;
t205 = t165 * t166;
t204 = t165 * t167;
t183 = -qJ(3) * t163 - pkin(1);
t134 = (-pkin(2) - qJ(4)) * t165 + t183;
t120 = t134 * qJD(1) + qJD(2);
t195 = t163 * qJD(1);
t149 = qJ(2) * t195 + qJD(3);
t135 = pkin(3) * t195 + t149;
t162 = sin(pkin(8));
t112 = t164 * t120 + t162 * t135;
t144 = t215 * t163;
t203 = t164 * t134 + t162 * t144;
t190 = t162 * t200;
t181 = qJD(5) * t190;
t184 = qJD(5) * t195;
t123 = t166 * t181 + t167 * t184;
t202 = t215 * t165;
t159 = t163 ^ 2;
t201 = t159 + t161;
t199 = qJD(2) * t161;
t198 = qJD(2) * t163;
t197 = qJD(2) * t165;
t196 = qJD(5) * t143;
t194 = t163 * qJD(3);
t193 = qJD(1) * qJD(2);
t192 = qJD(1) * qJD(3);
t191 = 0.2e1 * t199;
t189 = t166 * t196;
t188 = t167 * t196;
t136 = t215 * t200 + qJD(4);
t187 = qJ(2) * t193;
t186 = t163 * t193;
t185 = t165 * t193;
t182 = (-t162 ^ 2 - t164 ^ 2) * MDP(14);
t110 = pkin(6) * t195 + t112;
t173 = (pkin(4) * t164 + pkin(6) * t162) * t165;
t115 = qJD(1) * t173 + t136;
t180 = -t110 * t167 - t115 * t166;
t179 = t110 * t166 - t115 * t167;
t111 = -t120 * t162 + t135 * t164;
t178 = -t111 * t162 + t112 * t164;
t142 = -qJD(4) * t165 - t194;
t139 = t142 * qJD(1);
t116 = t139 * t162 - t164 * t186;
t117 = t139 * t164 + t162 * t186;
t177 = -t116 * t162 - t117 * t164;
t176 = -t134 * t162 + t144 * t164;
t175 = -MDP(12) * t162 - MDP(13) * t164;
t174 = -pkin(2) * t165 + t183;
t172 = t162 * t204 - t208;
t132 = t162 * t205 + t207;
t124 = -t166 * t184 + t167 * t181;
t171 = t167 * t117 + t166 * t185;
t170 = t180 * qJD(5) - t166 * t117 + t167 * t185;
t126 = t132 * qJD(1);
t129 = -t166 * t195 + t167 * t190;
t169 = (-t126 * t162 - t166 * t211) * MDP(21) + (-t129 * t162 - t167 * t211) * MDP(22);
t151 = t159 * t187;
t131 = t174 * qJD(1) + qJD(2);
t128 = t172 * qJD(5);
t127 = t132 * qJD(5);
t122 = t142 * t164 + t162 * t198;
t121 = t142 * t162 - t164 * t198;
t119 = t173 + t202;
t114 = pkin(6) * t163 + t203;
t113 = -pkin(4) * t163 - t176;
t109 = -pkin(4) * t195 - t111;
t1 = [0.2e1 * (t161 * t187 + t151) * MDP(7) + 0.2e1 * t159 * MDP(10) * t192 + (t151 + (t149 * qJD(2) - t131 * qJD(3)) * t163 + (qJ(2) * t191 - t174 * t194) * qJD(1)) * MDP(11) + (-t116 * t163 + (-t121 * t163 + t164 * t191) * qJD(1)) * MDP(12) + (-t117 * t163 + (-t122 * t163 - 0.2e1 * t162 * t199) * qJD(1)) * MDP(13) + (t117 * t203 + t112 * t122 - t116 * t176 - t111 * t121 + (qJD(1) * t202 + t136) * t197) * MDP(15) + (-t123 * t172 - t127 * t129) * MDP(16) + (t123 * t132 - t124 * t172 + t126 * t127 - t128 * t129) * MDP(17) + (t123 * t206 + t127 * t143) * MDP(18) + (t124 * t206 + t128 * t143) * MDP(19) + ((t167 * t197 - t122 * t166 + (-t114 * t167 - t119 * t166) * qJD(5)) * t143 + t170 * t206 - t121 * t126 - t113 * t124 - t116 * t132 - t109 * t128) * MDP(21) + (-(t166 * t197 + t122 * t167 + (-t114 * t166 + t119 * t167) * qJD(5)) * t143 - (-t179 * qJD(5) + t171) * t206 - t121 * t129 + t113 * t123 - t116 * t172 + t109 * t127) * MDP(22) + 0.2e1 * t213 * t201 * t193 + (-0.2e1 * t163 * MDP(9) * t192 + ((-t121 * t162 - t122 * t164) * qJD(1) + t177) * MDP(14)) * t165; (-qJ(2) * t210 + (-qJD(3) - t149) * t195) * MDP(11) + ((-t136 * t165 + (-t111 * t164 - t112 * t162) * t163) * qJD(1) - t177) * MDP(15) + (-t164 * t188 - t162 * t124 + (-(-t162 * t208 + t204) * t143 - t126 * t209) * qJD(1)) * MDP(21) + (t164 * t189 + t162 * t123 + ((t162 * t207 + t205) * t143 - t129 * t209) * qJD(1)) * MDP(22) + (-t164 * MDP(12) + MDP(13) * t162 - qJ(2) * MDP(7) - t213) * t201 * t168; (-t116 * t164 + t117 * t162) * MDP(15) + (t124 * t164 - t162 * t188) * MDP(21) + (-t123 * t164 + t162 * t189) * MDP(22) + (-MDP(10) + t175) * t168 * t159 + ((MDP(9) + t182) * t168 * t165 + ((qJD(2) + t131) * MDP(11) + t178 * MDP(15) + t169) * qJD(1)) * t163; t182 * t210 + (-t166 * MDP(21) - t167 * MDP(22)) * t196 + (t175 * t168 * t163 + ((qJD(2) + t178) * MDP(15) + t169) * qJD(1)) * t165; t129 * t126 * MDP(16) + (-t126 ^ 2 + t129 ^ 2) * MDP(17) + (-t126 * t143 + t123) * MDP(18) + (-t129 * t143 + t124) * MDP(19) + (t109 * t129 - t180 * t143 + t170) * MDP(21) + (-t109 * t126 - t171 + (qJD(5) - t143) * t179) * MDP(22);];
tauc = t1;
