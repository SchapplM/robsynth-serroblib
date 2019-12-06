% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:41:00
% EndTime: 2019-12-05 18:41:01
% DurationCPUTime: 0.53s
% Computational Cost: add. (795->125), mult. (1724->184), div. (0->0), fcn. (960->8), ass. (0->84)
t167 = sin(qJ(3));
t171 = cos(qJ(2));
t168 = sin(qJ(2));
t170 = cos(qJ(3));
t213 = t168 * t170;
t182 = -t167 * t171 - t213;
t201 = qJD(3) * t170;
t192 = t168 * t201;
t174 = (t182 * qJD(2) - t192) * pkin(1);
t173 = qJD(1) * t174;
t161 = qJD(1) + qJD(2);
t222 = pkin(1) * qJD(1);
t195 = t171 * t222;
t149 = pkin(2) * t161 + t195;
t202 = qJD(3) * t167;
t193 = t149 * t202;
t226 = t173 - t193;
t225 = t168 * MDP(5) + t171 * MDP(6);
t166 = sin(qJ(5));
t169 = cos(qJ(5));
t224 = (t166 ^ 2 - t169 ^ 2) * MDP(12);
t160 = qJD(3) + t161;
t172 = qJD(5) ^ 2;
t157 = pkin(2) * t170 + pkin(3);
t164 = sin(pkin(9));
t165 = cos(pkin(9));
t216 = t165 * t167;
t205 = pkin(2) * t216 + t164 * t157;
t144 = t182 * t222;
t212 = t170 * t171;
t215 = t167 * t168;
t181 = t212 - t215;
t145 = t181 * t222;
t221 = pkin(2) * qJD(3);
t208 = -t144 * t165 + t145 * t164 - (t164 * t170 + t216) * t221;
t223 = -t208 * t160 + (pkin(8) + t205) * t172;
t196 = t168 * t222;
t136 = t149 * t167 + t170 * t196;
t220 = t136 * t164;
t218 = t164 * t167;
t217 = t165 * t136;
t189 = qJD(2) * t195;
t190 = t167 * t196;
t206 = (qJD(2) + qJD(3)) * t190;
t124 = (qJD(3) * t149 + t189) * t170 - t206;
t111 = t124 * t164 - t165 * t226;
t135 = t170 * t149 - t190;
t132 = pkin(3) * t160 + t135;
t119 = t132 * t165 - t220;
t117 = -pkin(4) * t160 - t119;
t200 = qJD(5) * t117;
t210 = t111 * t166 + t169 * t200;
t158 = pkin(1) * t171 + pkin(2);
t143 = -pkin(1) * t215 + t158 * t170 + pkin(3);
t146 = pkin(1) * t213 + t158 * t167;
t209 = t164 * t143 + t165 * t146;
t207 = -t144 * t164 - t145 * t165 + (t165 * t170 - t218) * t221;
t203 = MDP(11) * t169;
t199 = t169 * MDP(17);
t198 = t172 * MDP(14);
t194 = t172 * t169 * MDP(13) + 0.2e1 * (t203 * t166 - t224) * qJD(5) * t160;
t112 = t165 * t124 + t226 * t164;
t191 = -t117 * t160 - t112;
t188 = MDP(8) * t192;
t186 = (-pkin(2) * t160 - t149) * qJD(3);
t129 = t158 * t201 + (t181 * qJD(2) - t168 * t202) * pkin(1);
t130 = -t158 * t202 + t174;
t113 = t129 * t164 - t130 * t165;
t185 = t113 * t160 + (pkin(8) + t209) * t172;
t121 = t135 * t164 + t217;
t184 = -t121 * t160 + (pkin(3) * t164 + pkin(8)) * t172;
t183 = t143 * t165 - t146 * t164;
t114 = t129 * t165 + t130 * t164;
t180 = qJD(5) * ((-pkin(4) - t183) * t160 - t114);
t122 = t135 * t165 - t220;
t179 = qJD(5) * ((-pkin(3) * t165 - pkin(4)) * t160 + t122);
t178 = -t149 * t201 + t206;
t177 = -pkin(2) * t218 + t157 * t165;
t176 = qJD(5) * ((-pkin(4) - t177) * t160 - t207);
t115 = t166 * t200;
t175 = t115 * MDP(16) + t210 * MDP(17) + t194;
t159 = t160 ^ 2;
t120 = t164 * t132 + t217;
t1 = [(t130 * t160 - t193) * MDP(8) + (-t129 * t160 + t178) * MDP(9) + (-t111 * t183 + t112 * t209 - t119 * t113 + t120 * t114) * MDP(10) + ((-t111 - t185) * MDP(16) + MDP(17) * t180) * t169 + (MDP(16) * t180 + t185 * MDP(17) - t198) * t166 + (-qJD(1) * t188 + (-t225 * t161 + ((-t170 * MDP(8) - MDP(5)) * t168 + (-t167 * MDP(8) - t170 * MDP(9) - MDP(6)) * t171) * qJD(1)) * qJD(2)) * pkin(1) + t175; (-t144 * t160 + t167 * t186 + t173) * MDP(8) + (t145 * t160 + (t186 - t189) * t170 + t206) * MDP(9) + (-t111 * t177 + t112 * t205 + t208 * t119 + t207 * t120) * MDP(10) - t166 * t198 + (t115 + t166 * t176 + (-t111 - t223) * t169) * MDP(16) + (t166 * t223 + t169 * t176 + t210) * MDP(17) + t194 + t225 * (-qJD(2) + t161) * t222; (t136 * t160 - t193) * MDP(8) + (t135 * t160 + t178) * MDP(9) + (t119 * t121 - t120 * t122 + (-t111 * t165 + t112 * t164) * pkin(3)) * MDP(10) + (-t188 + (t182 * MDP(8) - MDP(9) * t212) * qJD(2)) * t222 + ((-t111 - t184) * MDP(16) + MDP(17) * t179) * t169 + (MDP(16) * t179 + t184 * MDP(17) - t198) * t166 + t175; (-t166 * MDP(16) - t199) * t172; t191 * t199 + t159 * t224 + (t191 * MDP(16) - t159 * t203) * t166;];
tauc = t1;
