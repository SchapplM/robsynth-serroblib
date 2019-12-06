% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRRR3
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
%   see S5PPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:17:02
% EndTime: 2019-12-05 15:17:06
% DurationCPUTime: 1.04s
% Computational Cost: add. (467->131), mult. (1219->215), div. (0->0), fcn. (870->8), ass. (0->81)
t170 = cos(qJ(4));
t160 = -t170 * pkin(4) - pkin(3);
t168 = sin(qJ(3));
t171 = cos(qJ(3));
t164 = sin(pkin(9));
t210 = qJD(1) * t164;
t223 = t171 * qJD(2) - t168 * t210;
t132 = t160 * qJD(3) - t223;
t228 = t132 + t223;
t218 = pkin(6) + pkin(7);
t143 = t168 * qJD(2) + t171 * t210;
t166 = sin(qJ(5));
t167 = sin(qJ(4));
t169 = cos(qJ(5));
t145 = t166 * t170 + t169 * t167;
t161 = qJD(4) + qJD(5);
t226 = t145 * t161;
t124 = qJD(3) * t226;
t219 = qJD(5) - t161;
t198 = qJD(3) * qJD(4);
t227 = -0.2e1 * t198;
t144 = t166 * t167 - t169 * t170;
t176 = t144 * t161;
t172 = qJD(4) ^ 2;
t173 = qJD(3) ^ 2;
t225 = (t172 + t173) * t168;
t224 = (t167 ^ 2 - t170 ^ 2) * MDP(7);
t204 = qJD(4) * t170;
t221 = -qJD(5) * t170 - t204;
t205 = qJD(4) * t167;
t220 = qJD(5) * t167 + t205;
t217 = qJD(3) * pkin(3);
t216 = t164 * t168;
t215 = t164 * t171;
t165 = cos(pkin(9));
t209 = qJD(1) * t165;
t154 = t167 * t209;
t188 = t218 * qJD(3) + t143;
t126 = t188 * t170 - t154;
t214 = t169 * t126;
t213 = t173 * MDP(5);
t208 = qJD(3) * t167;
t207 = qJD(3) * t170;
t206 = qJD(3) * t171;
t197 = pkin(4) * t208;
t196 = t173 * t215;
t195 = qJD(4) * t218;
t194 = qJD(3) * t216;
t193 = t166 * t208;
t192 = t169 * t207;
t190 = t170 * t198;
t125 = -t188 * t167 - t170 * t209;
t122 = qJD(4) * pkin(4) + t125;
t189 = -pkin(4) * t161 - t122;
t133 = -t223 - t217;
t136 = t223 * qJD(3);
t187 = -t133 * qJD(3) - t136;
t186 = t171 * t227;
t185 = t167 * t194;
t184 = t170 * t194;
t183 = pkin(4) * t205 - t143;
t181 = pkin(6) * t172;
t180 = qJD(4) * (t133 + t223 - t217);
t138 = -t165 * t167 + t170 * t215;
t137 = -t165 * t170 - t167 * t215;
t118 = qJD(4) * t125 + t170 * t136;
t119 = t154 * qJD(4) - t167 * t136 - t188 * t204;
t141 = -t166 * t207 - t169 * t208;
t179 = -t166 * t118 + t169 * t119 + t132 * t141;
t123 = qJD(5) * t192 - t161 * t193 + t169 * t190;
t139 = -t192 + t193;
t175 = -t141 * t139 * MDP(13) + (t139 * t161 + t123) * MDP(15) + (-t141 * t161 - t124) * MDP(16) + (-t139 ^ 2 + t141 ^ 2) * MDP(14);
t174 = t132 * t139 + (t219 * t126 - t119) * t166;
t149 = t218 * t170;
t148 = t218 * t167;
t147 = t170 * t195;
t146 = t167 * t195;
t131 = t143 * qJD(3) + qJD(4) * t197;
t130 = t137 * qJD(4) - t184;
t129 = -t138 * qJD(4) + t185;
t1 = [-MDP(4) * t196 + t213 * t216 + (-t170 * t196 + (t129 + t185) * qJD(4)) * MDP(11) + (t167 * t196 + (-t130 + t184) * qJD(4)) * MDP(12) + ((t169 * t129 - t166 * t130 + (-t137 * t166 - t138 * t169) * qJD(5)) * t161 + (t168 * t124 + t139 * t206) * t164) * MDP(18) + (-(t166 * t129 + t169 * t130 + (t137 * t169 - t138 * t166) * qJD(5)) * t161 + (t168 * t123 - t141 * t206) * t164) * MDP(19); -t173 * t168 * MDP(4) - t171 * t213 + (t167 * t186 - t170 * t225) * MDP(11) + (t167 * t225 + t170 * t186) * MDP(12) + (-0.2e1 * t171 * t124 + ((t220 * t166 + t221 * t169) * t161 + qJD(3) * t139) * t168) * MDP(18) + ((qJD(3) * t176 - t123) * t171 + (-(t221 * t166 - t220 * t169) * t161 - qJD(3) * t141) * t168) * MDP(19); t224 * t227 + (t123 * t145 + t141 * t176) * MDP(13) + (-t123 * t144 - t145 * t124 + t139 * t176 + t141 * t226) * MDP(14) + (t160 * t124 + t131 * t144 + t183 * t139 + t226 * t228) * MDP(18) + (t160 * t123 + t131 * t145 - t183 * t141 - t176 * t228) * MDP(19) + (-t181 * MDP(11) + t180 * MDP(12) + t172 * MDP(8)) * t170 + (t180 * MDP(11) + t181 * MDP(12) + 0.2e1 * MDP(6) * t190 - t172 * MDP(9)) * t167 + (-t176 * MDP(15) - t226 * MDP(16) + (t166 * t146 - t169 * t147 + (t148 * t166 - t149 * t169) * qJD(5)) * MDP(18) + (t169 * t146 + t166 * t147 - (-t148 * t169 - t149 * t166) * qJD(5)) * MDP(19)) * t161; t173 * t224 + t187 * t170 * MDP(12) + (-(-t166 * t125 - t214) * t161 - t139 * t197 + (t189 * t166 - t214) * qJD(5) + t179) * MDP(18) + (t141 * t197 + (t189 * qJD(5) + t125 * t161 - t118) * t169 + t174) * MDP(19) + t175 + (-t173 * t170 * MDP(6) + t187 * MDP(11)) * t167; (t179 + t219 * (-t166 * t122 - t214)) * MDP(18) + ((-t219 * t122 - t118) * t169 + t174) * MDP(19) + t175;];
tauc = t1;
