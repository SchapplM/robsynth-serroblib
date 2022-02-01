% Calculate Coriolis joint torque vector for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:53
% EndTime: 2022-01-20 10:05:56
% DurationCPUTime: 0.72s
% Computational Cost: add. (621->122), mult. (1322->189), div. (0->0), fcn. (742->8), ass. (0->81)
t156 = sin(pkin(8));
t160 = sin(qJ(2));
t209 = pkin(1) * qJD(1);
t184 = t160 * t209;
t141 = t156 * t184;
t158 = cos(pkin(8));
t162 = cos(qJ(2));
t208 = pkin(1) * qJD(2);
t143 = t158 * t162 * t208;
t127 = qJD(1) * t143 - qJD(2) * t141;
t152 = qJD(1) + qJD(2);
t122 = qJD(4) * t152 + t127;
t155 = sin(pkin(9));
t150 = t155 ^ 2;
t157 = cos(pkin(9));
t191 = t157 ^ 2 + t150;
t217 = t122 * t191;
t159 = sin(qJ(5));
t161 = cos(qJ(5));
t203 = t150 * t161;
t216 = -t159 * MDP(11) * t203 + (t159 ^ 2 - t161 ^ 2) * MDP(12) * t150;
t215 = t160 * MDP(5) + t162 * MDP(6);
t183 = t162 * t209;
t132 = t158 * t183 - t141;
t211 = qJD(4) - t132;
t210 = pkin(2) * t158;
t137 = t152 * pkin(2) + t183;
t142 = t158 * t184;
t124 = t156 * t137 + t142;
t119 = qJ(4) * t152 + t124;
t110 = -t157 * qJD(3) + t119 * t155;
t207 = t110 * t155;
t144 = t156 * t160 * pkin(1);
t165 = -qJD(2) * t144 + t143;
t128 = qJD(4) + t165;
t206 = t128 * t152;
t205 = t150 * t152;
t204 = t150 * t159;
t202 = t152 * t155;
t201 = t152 * t157;
t200 = t157 * MDP(8);
t199 = t157 * t159;
t198 = t157 * t161;
t197 = t158 * t160;
t140 = -qJD(5) + t201;
t196 = t159 * t140;
t195 = t161 * t140;
t166 = -pkin(4) * t157 - pkin(7) * t155 - pkin(3);
t123 = t137 * t158 - t141;
t170 = qJD(4) - t123;
t108 = t166 * t152 + t170;
t111 = qJD(3) * t155 + t119 * t157;
t131 = (t156 * t162 + t197) * t208;
t126 = qJD(1) * t131;
t194 = (t122 * t198 + t159 * t126 + (t108 * t161 - t111 * t159) * qJD(5)) * t157 + t122 * t203;
t147 = pkin(1) * t162 + pkin(2);
t192 = pkin(1) * t197 + t156 * t147;
t187 = qJD(5) * t159;
t186 = qJD(5) * t161;
t185 = qJD(5) + t140;
t182 = t110 * t202;
t181 = t152 * t203;
t180 = t155 * t186;
t179 = t155 * t187;
t177 = t147 * t158 - t144;
t136 = t179 * t201;
t176 = (t140 * t179 + t136) * MDP(13) + (t140 + t201) * MDP(14) * t180 + 0.2e1 * t216 * qJD(5) * t152;
t169 = -t108 * t159 - t111 * t161;
t107 = t169 * qJD(5) - t122 * t199 + t161 * t126;
t172 = -t107 * t157 + t110 * t180 + t122 * t204;
t168 = t111 * t157 + t207;
t167 = t140 * t157 + t205;
t146 = pkin(2) * t156 + qJ(4);
t164 = t146 * t186 + t211 * t159;
t149 = t152 ^ 2;
t134 = t166 - t210;
t130 = t156 * t183 + t142;
t129 = qJ(4) + t192;
t120 = t166 - t177;
t116 = -pkin(3) * t152 + t170;
t1 = [(-t123 * t131 + t124 * t165 - t126 * t177 + t127 * t192) * MDP(7) + (-t131 * t152 - t126) * t200 + (t191 * t206 + t217) * MDP(9) + (t126 * (-pkin(3) - t177) + t116 * t131 + t168 * t128 + t129 * t217) * MDP(10) + (-(-t128 * t199 + t131 * t161) * t140 + t204 * t206 + (-(-t120 * t159 - t129 * t198) * t140 + t129 * t181) * qJD(5) + t172) * MDP(16) + ((t128 * t198 + t131 * t159) * t140 + t128 * t181 + (t120 * t195 + (-t167 * t129 - t207) * t159) * qJD(5) + t194) * MDP(17) + t176 + t215 * (-qJD(1) - t152) * t208; (t123 * t130 - t124 * t132 + (-t126 * t158 + t127 * t156) * pkin(2)) * MDP(7) + (t130 * t152 - t126) * t200 + (t211 * t152 * t191 + t217) * MDP(9) + (t126 * (-pkin(3) - t210) - t116 * t130 + t146 * t217 + t211 * t168) * MDP(10) + ((t130 * t161 + t134 * t187 + t164 * t157) * t140 + t164 * t205 + t172) * MDP(16) + (-t130 * t196 + t211 * t167 * t161 + (t134 * t195 + (-t167 * t146 - t207) * t159) * qJD(5) + t194) * MDP(17) + t176 + t215 * (-qJD(2) + t152) * t209; t136 * MDP(17) + (-MDP(17) * t196 + (t140 - t201) * MDP(16) * t161) * t155 * qJD(5); (-t168 * t152 + t126) * MDP(10) - t191 * MDP(9) * t149 + (MDP(16) * t159 + MDP(17) * t161) * (-t140 ^ 2 - t150 * t149); (t169 * t140 - t161 * t182 + t107) * MDP(16) + ((-t185 * t108 - t157 * t122) * t161 + (t185 * t111 - t126 + t182) * t159) * MDP(17) - (t159 * MDP(13) + t161 * MDP(14)) * t185 * t202 - t216 * t149;];
tauc = t1;
