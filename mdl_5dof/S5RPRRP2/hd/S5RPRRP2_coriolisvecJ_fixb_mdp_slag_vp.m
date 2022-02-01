% Calculate Coriolis joint torque vector for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:14
% EndTime: 2022-01-23 09:28:16
% DurationCPUTime: 0.74s
% Computational Cost: add. (849->162), mult. (1683->218), div. (0->0), fcn. (876->6), ass. (0->89)
t168 = sin(qJ(4));
t170 = cos(qJ(4));
t156 = cos(pkin(8)) * pkin(1) + pkin(2);
t150 = t156 * qJD(1);
t169 = sin(qJ(3));
t171 = cos(qJ(3));
t213 = pkin(1) * sin(pkin(8));
t187 = qJD(1) * t213;
t138 = t150 * t169 + t171 * t187;
t163 = qJD(1) + qJD(3);
t128 = pkin(7) * t163 + t138;
t183 = qJ(5) * t163 + t128;
t177 = t183 * t170;
t116 = qJD(2) * t168 + t177;
t216 = qJD(4) * t116;
t164 = t168 ^ 2;
t165 = t170 ^ 2;
t215 = (t164 - t165) * MDP(9);
t214 = t156 * t171 - t169 * t213;
t212 = pkin(4) * t164;
t211 = pkin(4) * t170;
t210 = -qJ(5) - pkin(7);
t209 = qJD(4) * pkin(4);
t208 = t116 * t170;
t137 = t150 * t171 - t169 * t187;
t127 = -pkin(3) * t163 - t137;
t207 = t127 * t163;
t174 = t156 * t169 + t171 * t213;
t141 = t174 * qJD(3);
t206 = t141 * t163;
t162 = t163 ^ 2;
t204 = t162 * t170;
t203 = t163 * t168;
t172 = qJD(4) ^ 2;
t202 = t170 * t172;
t143 = pkin(7) + t174;
t201 = -qJ(5) - t143;
t160 = t170 * qJD(2);
t115 = -t168 * t183 + t160;
t114 = t115 + t209;
t200 = t114 - t115;
t158 = -pkin(3) - t211;
t120 = t158 * t163 + qJD(5) - t137;
t178 = qJD(3) * t187;
t192 = qJD(3) * t150;
t135 = t169 * t192 + t171 * t178;
t191 = qJD(4) * t168;
t185 = t163 * t191;
t121 = pkin(4) * t185 + t135;
t190 = qJD(4) * t170;
t199 = t120 * t190 + t121 * t168;
t198 = t127 * t190 + t135 * t168;
t197 = t138 * t170 * t163 + t137 * t191;
t195 = MDP(16) * t163;
t194 = MDP(17) * t163;
t189 = qJD(5) * t163;
t188 = -MDP(14) - MDP(16);
t186 = 0.2e1 * t170 * MDP(8) * t185 - 0.2e1 * t163 * qJD(4) * t215 + MDP(10) * t202;
t184 = qJD(4) * t210;
t140 = t214 * qJD(3);
t142 = -pkin(3) - t214;
t182 = t142 * t163 - t140;
t181 = qJD(4) * t201;
t180 = (-t164 - t165) * MDP(17);
t134 = -t169 * t178 + t171 * t192;
t179 = qJD(4) * qJD(2) + t134;
t176 = t114 * t168 - t208;
t175 = t143 * t172 + t206;
t173 = (-qJD(5) - t120) * t163 - t179;
t161 = t170 * qJ(5);
t159 = t170 * qJD(5);
t152 = pkin(7) * t170 + t161;
t151 = t210 * t168;
t145 = -qJD(5) * t168 + t170 * t184;
t144 = t168 * t184 + t159;
t136 = t142 - t211;
t133 = pkin(4) * t191 + t141;
t132 = t143 * t170 + t161;
t131 = t201 * t168;
t130 = t137 * t190;
t124 = t128 * t191;
t122 = t127 * t191;
t117 = t120 * t191;
t113 = (-qJD(5) - t140) * t168 + t170 * t181;
t112 = t140 * t170 + t168 * t181 + t159;
t111 = (-t134 - t189) * t168 - t216;
t110 = -qJ(5) * t185 - t124 + (t179 + t189) * t170;
t109 = t110 * t170;
t1 = [(-t135 - t206) * MDP(6) + (-t140 * t163 - t134) * MDP(7) + t122 * MDP(13) + t198 * MDP(14) + t117 * MDP(15) + t199 * MDP(16) + t109 * MDP(17) + (t110 * t132 + t111 * t131 + t112 * t116 + t113 * t114 + t120 * t133 + t121 * t136) * MDP(18) + ((-t135 - t175) * MDP(13) + (-t133 * t163 - t121) * MDP(15) + t112 * t194) * t170 + (-t172 * MDP(11) + t175 * MDP(14) + t133 * t195 + (-t113 * t163 - t111) * MDP(17)) * t168 + (t113 * MDP(15) - t112 * MDP(16) + (t182 * MDP(14) + t136 * t195 + (-t131 * t163 - t114) * MDP(17)) * t170 + (t182 * MDP(13) + t136 * t163 * MDP(15) + (-t132 * t163 - t116) * MDP(17)) * t168) * qJD(4) + t186; (-qJD(4) * t176 + t110 * t168 + t111 * t170) * MDP(18) + (t188 * t170 + (-MDP(13) - MDP(15)) * t168) * t172; -t135 * MDP(6) - t134 * MDP(7) + (-pkin(7) * t202 - t135 * t170 + t122 + t197) * MDP(13) + (t130 + t198) * MDP(14) + (qJD(4) * t145 - t121 * t170 + t117 + t197) * MDP(15) + (-qJD(4) * t144 + t130 + t199) * MDP(16) + (-t114 * t190 + t109) * MDP(17) + (t110 * t152 + t111 * t151 + t114 * t145 + t116 * t144 - t120 * t138 + t121 * t158 - t137 * t208) * MDP(18) + ((-t111 - t216) * MDP(17) + (t114 * t137 + t120 * t209) * MDP(18) + (pkin(7) * MDP(14) - MDP(11)) * t172) * t168 + ((t144 * t170 - t145 * t168) * MDP(17) + (MDP(7) + t180) * t137 + (t168 * t188 + MDP(6)) * t138 + (MDP(16) * t212 + (-pkin(3) * MDP(14) + t158 * MDP(16) - t151 * MDP(17)) * t170 + (-pkin(3) * MDP(13) + (t158 - t211) * MDP(15) - t152 * MDP(17)) * t168) * qJD(4)) * t163 + t186; -t168 * MDP(8) * t204 + t162 * t215 + (-t134 - t207) * t168 * MDP(13) + (t124 + (-t128 * t168 + t160) * qJD(4) + (-t179 - t207) * t170) * MDP(14) + ((t116 - t177) * qJD(4) + (pkin(4) * t204 + t173) * t168) * MDP(15) + (-t162 * t212 + t124 + (qJ(5) * t203 + t115) * qJD(4) + t173 * t170) * MDP(16) + (t200 - t209) * t170 * t194 + (t200 * t116 + (-t120 * t203 + t111) * pkin(4)) * MDP(18); t121 * MDP(18) + t162 * t180 + (t176 * MDP(18) + 0.2e1 * (t168 * MDP(15) + t170 * MDP(16)) * qJD(4)) * t163;];
tauc = t1;
