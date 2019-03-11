% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:16
% EndTime: 2019-03-09 01:30:19
% DurationCPUTime: 1.25s
% Computational Cost: add. (611->189), mult. (1236->270), div. (0->0), fcn. (673->6), ass. (0->95)
t140 = sin(pkin(9)) * pkin(1) + qJ(3);
t134 = qJD(1) * t140;
t218 = t134 * MDP(7);
t217 = MDP(6) + MDP(8);
t152 = cos(qJ(5));
t146 = t152 ^ 2;
t150 = sin(qJ(5));
t216 = (t150 ^ 2 - t146) * MDP(12);
t133 = qJD(4) + t134;
t125 = -qJD(1) * pkin(7) + t133;
t120 = qJD(2) * t152 + t125 * t150;
t183 = qJD(3) * qJD(1);
t112 = qJD(5) * t120 - t152 * t183;
t149 = sin(qJ(6));
t214 = t112 * t149;
t151 = cos(qJ(6));
t213 = t112 * t151;
t187 = t151 * qJD(5);
t141 = qJD(6) * t187;
t191 = qJD(6) * t149;
t175 = t152 * t191;
t156 = -t150 * t187 - t175;
t117 = t156 * qJD(1) + t141;
t212 = t117 * t149;
t196 = qJD(1) * t152;
t172 = t151 * t196;
t194 = qJD(5) * t149;
t131 = t172 + t194;
t184 = qJD(1) * qJD(5);
t171 = t150 * t184;
t136 = t149 * t171;
t118 = t131 * qJD(6) - t136;
t211 = t118 * t150;
t173 = t149 * t196;
t129 = t173 - t187;
t197 = qJD(1) * t150;
t139 = qJD(6) + t197;
t210 = t129 * t139;
t209 = t129 * t152;
t208 = t131 * t139;
t207 = t131 * t152;
t206 = t139 * t149;
t205 = t139 * t151;
t204 = t149 * t150;
t203 = t150 * t151;
t192 = qJD(5) * t152;
t202 = t117 * t150 + t131 * t192;
t153 = qJD(5) ^ 2;
t154 = qJD(1) ^ 2;
t200 = -t153 - t154;
t138 = cos(pkin(9)) * pkin(1) + pkin(2) + qJ(4);
t199 = qJD(1) * t138;
t198 = qJD(1) * t146;
t137 = -pkin(7) + t140;
t195 = qJD(5) * t137;
t193 = qJD(5) * t150;
t190 = qJD(6) * t151;
t189 = t150 * MDP(17);
t188 = t151 * MDP(23);
t186 = t152 * MDP(16);
t126 = -qJD(3) + t199;
t185 = qJD(3) - t126;
t181 = 0.2e1 * qJD(4) * qJD(1);
t180 = t139 * t204;
t179 = t139 * t203;
t178 = t149 * t198;
t177 = t139 * t194;
t176 = t139 * t191;
t174 = t139 * t190;
t170 = t152 * t184;
t169 = MDP(22) * t192;
t115 = qJD(5) * pkin(8) + t120;
t168 = t137 * t139 + t115;
t167 = t139 * t175;
t166 = pkin(5) * t152 + pkin(8) * t150;
t165 = t149 * t170 + t174;
t164 = qJD(3) + t126 + t199;
t124 = pkin(5) * t150 - pkin(8) * t152 + t138;
t116 = t124 * qJD(1) - qJD(3);
t110 = t115 * t151 + t116 * t149;
t163 = t115 * t149 - t116 * t151;
t162 = -t139 * t150 + t198;
t161 = qJD(2) * t150 - t125 * t152;
t160 = -t137 * t153 + t181;
t159 = t150 * t177 - t152 * t174;
t114 = -qJD(5) * pkin(5) + t161;
t158 = -pkin(8) * t192 + t114 * t150;
t157 = (t149 * MDP(24) - t188) * t139;
t128 = t166 * qJD(5) + qJD(4);
t111 = -t161 * qJD(5) + t150 * t183;
t155 = -qJD(3) * t139 - t114 * qJD(5) - qJD(6) * t116 - t111;
t132 = t166 * qJD(1);
t122 = t128 * qJD(1);
t121 = t151 * t122;
t1 = [0.2e1 * qJD(3) * t218 + MDP(9) * t181 + (qJD(3) * t133 + qJD(4) * t126 + (qJD(3) * t140 + qJD(4) * t138) * qJD(1)) * MDP(10) + 0.2e1 * t184 * t216 - t153 * t152 * MDP(14) + t164 * t192 * MDP(16) + (t160 * t152 - t164 * t193) * MDP(17) + (t117 * t151 * t152 + t156 * t131) * MDP(18) + ((t129 * t151 + t131 * t149) * t193 + (-t212 - t118 * t151 + (t129 * t149 - t131 * t151) * qJD(6)) * t152) * MDP(19) + (t162 * t187 - t167 + t202) * MDP(20) + (-t211 + (-t178 - t209) * qJD(5) + t159) * MDP(21) + (t139 + t197) * t169 + ((-t124 * t191 + t128 * t151) * t139 + (t114 * t190 - qJD(3) * t129 + t214 - t137 * t118 + (-t137 * t206 + (t124 * t151 - t137 * t204) * qJD(1) - t163) * qJD(5)) * t152) * MDP(23) + (-(t124 * t190 + t128 * t149) * t139 + (-t114 * t191 - qJD(3) * t131 + t213 - t137 * t117 + (-t137 * t205 - (t124 * t149 + t137 * t203) * qJD(1) - t110) * qJD(5)) * t152) * MDP(24) + 0.2e1 * t217 * t183 + (-0.2e1 * MDP(11) * t170 - t153 * MDP(13) + t160 * MDP(16) + (t129 * t195 + t155 * t149 - t168 * t190 + t121) * MDP(23) + (t131 * t195 + (t168 * qJD(6) - t122) * t149 + t155 * t151) * MDP(24)) * t150; (t159 + t211) * MDP(23) + (t167 + t202) * MDP(24) + (-t186 + t189) * t153 + ((-t178 + t209) * MDP(23) - t162 * MDP(24) * t151) * qJD(5); MDP(23) * t176 + t165 * MDP(24) - t217 * t154 + (-t218 + (-qJD(4) - t133) * MDP(10) + (t180 + t209) * MDP(23) + (t179 + t207) * MDP(24) + (0.2e1 * t189 + (-0.2e1 * MDP(16) - t188) * t152) * qJD(5)) * qJD(1); -t154 * MDP(9) + (t185 * MDP(10) + t157) * qJD(1) + (t200 * MDP(17) + (-t118 - t177) * MDP(23) + (-t139 * t187 - t117) * MDP(24)) * t152 + (t200 * MDP(16) + qJD(6) * t157 + ((t129 - t173) * MDP(23) + (t131 - t172) * MDP(24)) * qJD(5)) * t150; (t131 * t205 + t212) * MDP(18) + ((t117 - t210) * t151 + (-t118 - t208) * t149) * MDP(19) + ((t179 - t207) * qJD(1) + t165) * MDP(20) + (-t176 + (-t180 + (t129 + t187) * t152) * qJD(1)) * MDP(21) - t139 * MDP(22) * t196 + (-pkin(5) * t118 - t213 - (t132 * t151 + t149 * t161) * t139 - t120 * t129 + (-pkin(8) * t205 + t114 * t149) * qJD(6) + (t158 * t149 + t152 * t163) * qJD(1)) * MDP(23) + (-pkin(5) * t117 + t214 + (t132 * t149 - t151 * t161) * t139 - t120 * t131 + (pkin(8) * t206 + t114 * t151) * qJD(6) + (t110 * t152 + t158 * t151) * qJD(1)) * MDP(24) + (MDP(11) * t150 * t152 - t216) * t154 + (-t197 * MDP(17) + qJD(1) * t186) * t185; t131 * t129 * MDP(18) + (-t129 ^ 2 + t131 ^ 2) * MDP(19) + (-t151 * t171 + t141 + t210) * MDP(20) + (t136 + t208) * MDP(21) + qJD(1) * t169 + (t110 * t139 - t149 * t111 - t114 * t131 + t121) * MDP(23) + (-t151 * t111 + t114 * t129 - t149 * t122 - t139 * t163) * MDP(24) + (-MDP(20) * t173 - t131 * MDP(21) - t110 * MDP(23) + t163 * MDP(24)) * qJD(6);];
tauc  = t1;
