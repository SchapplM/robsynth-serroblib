% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:51
% EndTime: 2019-12-05 15:14:55
% DurationCPUTime: 0.96s
% Computational Cost: add. (461->112), mult. (1205->178), div. (0->0), fcn. (884->8), ass. (0->71)
t153 = sin(pkin(9));
t154 = cos(pkin(9));
t157 = sin(qJ(3));
t160 = cos(qJ(3));
t211 = -t153 * t157 + t154 * t160;
t125 = t211 * qJD(1);
t159 = cos(qJ(4));
t148 = -pkin(4) * t159 - pkin(3);
t120 = t148 * qJD(3) - t125;
t212 = t120 + t125;
t133 = t153 * t160 + t154 * t157;
t128 = t133 * qJD(3);
t155 = sin(qJ(5));
t156 = sin(qJ(4));
t158 = cos(qJ(5));
t135 = t155 * t159 + t156 * t158;
t150 = qJD(4) + qJD(5);
t209 = t150 * t135;
t114 = t209 * qJD(3);
t204 = pkin(6) + pkin(7);
t205 = qJD(5) - t150;
t210 = MDP(6) * t156;
t134 = t155 * t156 - t158 * t159;
t166 = t150 * t134;
t208 = (t156 ^ 2 - t159 ^ 2) * MDP(7);
t126 = t133 * qJD(1);
t175 = t204 * qJD(3) + t126;
t115 = t159 * qJD(2) - t175 * t156;
t207 = t156 * MDP(11) + t159 * MDP(12);
t124 = qJD(1) * t128;
t206 = qJD(3) * t126 - t124;
t116 = qJD(2) * t156 + t159 * t175;
t203 = qJD(3) * pkin(3);
t202 = t133 * t150;
t161 = qJD(4) ^ 2;
t197 = t156 * t161;
t196 = t158 * t116;
t195 = t159 * t161;
t190 = qJD(3) * t156;
t189 = qJD(3) * t159;
t188 = qJD(4) * t156;
t187 = qJD(4) * t159;
t183 = qJD(3) * qJD(4);
t182 = pkin(4) * t190;
t181 = pkin(4) * t188;
t180 = qJD(4) * t204;
t179 = t155 * t190;
t178 = t158 * t189;
t177 = t159 * t183;
t112 = qJD(4) * pkin(4) + t115;
t176 = -pkin(4) * t150 - t112;
t173 = -t126 + t181;
t170 = pkin(6) * t161 - t206;
t121 = -t125 - t203;
t169 = qJD(4) * (t121 + t125 - t203);
t127 = t211 * qJD(3);
t123 = qJD(1) * t127;
t106 = t115 * qJD(4) + t159 * t123;
t107 = -t116 * qJD(4) - t156 * t123;
t131 = -t155 * t189 - t158 * t190;
t168 = -t155 * t106 + t158 * t107 + t120 * t131;
t113 = qJD(5) * t178 - t150 * t179 + t158 * t177;
t129 = -t178 + t179;
t165 = -t131 * t129 * MDP(13) + (t129 * t150 + t113) * MDP(15) + (-t131 * t150 - t114) * MDP(16) + (-t129 ^ 2 + t131 ^ 2) * MDP(14);
t164 = t120 * t129 + (t205 * t116 - t107) * t155;
t139 = t204 * t159;
t138 = t204 * t156;
t137 = t159 * t180;
t136 = t156 * t180;
t119 = qJD(3) * t181 + t124;
t1 = [(-t127 * t188 - t133 * t195) * MDP(11) + (-t127 * t187 + t133 * t197) * MDP(12) + (-t114 * t211 - t127 * t209 + t128 * t129 + t202 * t166) * MDP(18) + (-t113 * t211 + t127 * t166 - t128 * t131 + t202 * t209) * MDP(19) + (-t128 * MDP(4) - t127 * MDP(5) + (-t128 * t159 - t188 * t211) * MDP(11) + (t128 * t156 - t187 * t211) * MDP(12)) * qJD(3); -t207 * t161 + (-MDP(18) * t209 + t166 * MDP(19)) * t150; t206 * MDP(4) + 0.2e1 * t177 * t210 - 0.2e1 * t183 * t208 + MDP(8) * t195 - MDP(9) * t197 + (t156 * t169 - t159 * t170) * MDP(11) + (t156 * t170 + t159 * t169) * MDP(12) + (t113 * t135 + t131 * t166) * MDP(13) + (-t113 * t134 - t114 * t135 + t129 * t166 + t131 * t209) * MDP(14) + (t148 * t114 + t119 * t134 + t173 * t129 + t212 * t209) * MDP(18) + (t148 * t113 + t119 * t135 - t173 * t131 - t212 * t166) * MDP(19) + (-t166 * MDP(15) - t209 * MDP(16) + (t136 * t155 - t137 * t158 + (t138 * t155 - t139 * t158) * qJD(5)) * MDP(18) + (t136 * t158 + t137 * t155 - (-t138 * t158 - t139 * t155) * qJD(5)) * MDP(19)) * t150; (-(-t115 * t155 - t196) * t150 - t129 * t182 + (t155 * t176 - t196) * qJD(5) + t168) * MDP(18) + (t131 * t182 + (qJD(5) * t176 + t115 * t150 - t106) * t158 + t164) * MDP(19) + t165 + t207 * (-qJD(3) * t121 - t123) + (-t159 * t210 + t208) * qJD(3) ^ 2; (t168 + t205 * (-t112 * t155 - t196)) * MDP(18) + ((-t205 * t112 - t106) * t158 + t164) * MDP(19) + t165;];
tauc = t1;
