% Calculate joint inertia matrix for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR7_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:59:48
% EndTime: 2019-03-09 13:59:50
% DurationCPUTime: 0.67s
% Computational Cost: add. (811->169), mult. (1451->232), div. (0->0), fcn. (1555->8), ass. (0->94)
t182 = sin(qJ(4));
t186 = cos(qJ(4));
t188 = -pkin(2) - pkin(3);
t162 = qJ(3) * t186 + t182 * t188;
t160 = -pkin(9) + t162;
t181 = sin(qJ(5));
t185 = cos(qJ(5));
t198 = MDP(27) * t181 + MDP(28) * t185;
t200 = t181 * MDP(24) + t185 * MDP(25);
t233 = pkin(10) - t160;
t139 = t233 * t181;
t140 = t233 * t185;
t180 = sin(qJ(6));
t184 = cos(qJ(6));
t154 = t180 * t181 - t184 * t185;
t157 = t180 * t185 + t181 * t184;
t245 = -t157 * MDP(31) + t154 * MDP(32);
t205 = (t139 * t184 + t140 * t180) * MDP(34) - (t139 * t180 - t140 * t184) * MDP(35) + t245;
t247 = -t160 * t198 - t200 + t205;
t237 = pkin(9) + pkin(10);
t164 = t237 * t181;
t166 = t237 * t185;
t206 = (-t164 * t184 - t166 * t180) * MDP(34) - (-t164 * t180 + t166 * t184) * MDP(35) - t245;
t246 = -pkin(9) * t198 + t200 + t206;
t221 = MDP(27) * t185;
t199 = -MDP(28) * t181 + t221;
t196 = MDP(20) + t199;
t218 = MDP(34) * t154;
t244 = t186 * (-MDP(35) * t157 + t196 - t218) - MDP(21) * t182;
t183 = sin(qJ(2));
t187 = cos(qJ(2));
t156 = -t182 * t187 + t183 * t186;
t130 = t157 * t156;
t131 = t154 * t156;
t238 = pkin(7) - pkin(8);
t165 = t238 * t183;
t167 = t238 * t187;
t136 = -t165 * t186 + t167 * t182;
t138 = t165 * t182 + t167 * t186;
t176 = t181 ^ 2;
t178 = t185 ^ 2;
t219 = MDP(29) * t157;
t243 = -(t130 * t157 - t131 * t154) * MDP(30) - t196 * t136 - t138 * MDP(21) - (t176 - t178) * t156 * MDP(23) - t131 * t219;
t161 = t182 * qJ(3) - t186 * t188;
t159 = pkin(4) + t161;
t235 = pkin(5) * t185;
t145 = t159 + t235;
t242 = -0.2e1 * t145;
t163 = -t187 * pkin(2) - t183 * qJ(3) - pkin(1);
t151 = t187 * pkin(3) - t163;
t241 = 0.2e1 * t151;
t240 = -0.2e1 * MDP(30);
t239 = 0.2e1 * MDP(35);
t155 = t182 * t183 + t186 * t187;
t236 = pkin(5) * t155;
t234 = pkin(4) + t159;
t126 = pkin(4) * t155 - pkin(9) * t156 + t151;
t229 = t138 * t185;
t118 = t229 + (-pkin(10) * t156 + t126) * t181;
t232 = t118 * t184;
t228 = t156 * t181;
t127 = pkin(5) * t228 + t136;
t231 = t127 * t154;
t230 = t127 * t157;
t227 = t156 * t185;
t143 = t157 * t182;
t144 = t154 * t182;
t226 = -t143 * MDP(34) + t144 * MDP(35);
t170 = -pkin(4) - t235;
t225 = t145 - t170;
t177 = t183 ^ 2;
t224 = t187 ^ 2 + t177;
t128 = t130 * MDP(32);
t129 = t131 * MDP(31);
t215 = t161 * MDP(20);
t214 = t162 * MDP(21);
t213 = t185 * MDP(23);
t212 = t176 * MDP(22) + MDP(19);
t211 = MDP(26) + MDP(33);
t210 = t155 * MDP(33) - t128 - t129;
t209 = t181 * t213;
t208 = 0.2e1 * t209 + t212;
t207 = -pkin(2) * MDP(14) - MDP(11);
t119 = t185 * t126 - t138 * t181;
t117 = -pkin(10) * t227 + t119 + t236;
t114 = t184 * t117 - t118 * t180;
t204 = -t185 * t181 * MDP(22) - MDP(17);
t201 = MDP(24) * t185 - MDP(25) * t181;
t197 = t154 * t240 + t219;
t195 = (MDP(34) * t184 - MDP(35) * t180) * pkin(5);
t194 = 0.2e1 * t199;
t120 = t126 * t181 + t229;
t115 = t117 * t180 + t232;
t1 = [t177 * MDP(4) + (pkin(7) ^ 2 * t224 + t163 ^ 2) * MDP(14) + MDP(1) + t211 * t155 ^ 2 - (-MDP(29) * t131 + t130 * t240) * t131 + (MDP(20) * t241 - 0.2e1 * t128 - 0.2e1 * t129) * t155 + 0.2e1 * (t114 * t155 + t127 * t130) * MDP(34) + (-t115 * t155 - t127 * t131) * t239 + 0.2e1 * (t119 * t155 + t136 * t228) * MDP(27) + 0.2e1 * (-t120 * t155 + t136 * t227) * MDP(28) + 0.2e1 * t224 * MDP(12) * pkin(7) + 0.2e1 * (-t163 * MDP(11) + pkin(1) * MDP(9)) * t187 + 0.2e1 * (-pkin(1) * MDP(10) - t163 * MDP(13) + t187 * MDP(5)) * t183 + (MDP(21) * t241 + 0.2e1 * (-MDP(16) + t201) * t155 + (t178 * MDP(22) + MDP(15) - 0.2e1 * t209) * t156) * t156; t183 * MDP(6) + t187 * MDP(7) + (-pkin(2) * t183 + qJ(3) * t187) * MDP(12) + (t130 * t145 - t231) * MDP(34) + (-t131 * t145 - t230) * MDP(35) + (t159 * t198 + t204) * t156 + (MDP(18) + t247) * t155 + ((qJ(3) * MDP(14) - MDP(10) + MDP(13)) * t187 + (-MDP(9) + t207) * t183) * pkin(7) - t243; MDP(8) + 0.2e1 * pkin(2) * MDP(11) + 0.2e1 * qJ(3) * MDP(13) + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(14) + t218 * t242 + t159 * t194 + 0.2e1 * t215 + 0.2e1 * t214 + (MDP(35) * t242 + t197) * t157 + t208; (-t130 * t186 - t143 * t155) * MDP(34) + (t131 * t186 + t144 * t155) * MDP(35) + (MDP(14) * pkin(7) + MDP(12)) * t183 + t198 * (-t155 * t182 - t156 * t186); t207 - t244; MDP(14); (t130 * t170 + t231) * MDP(34) + (-t131 * t170 + t230) * MDP(35) + (-pkin(4) * t198 - t204) * t156 + (-MDP(18) + t246) * t155 + t243; -t215 - t214 - t234 * t221 + t225 * t218 + (MDP(28) * t234 - 0.2e1 * t213) * t181 + (MDP(35) * t225 - t197) * t157 - t212; t244; 0.2e1 * t170 * t218 + pkin(4) * t194 + (t170 * t239 + t197) * t157 + t208; t155 * MDP(26) + t119 * MDP(27) - t120 * MDP(28) + (t184 * t236 + t114) * MDP(34) + (-t232 + (-t117 - t236) * t180) * MDP(35) + t201 * t156 + t210; t247; -t182 * t198 + t226; t246; 0.2e1 * t195 + t211; t114 * MDP(34) - t115 * MDP(35) + t210; t205; t226; t206; MDP(33) + t195; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
