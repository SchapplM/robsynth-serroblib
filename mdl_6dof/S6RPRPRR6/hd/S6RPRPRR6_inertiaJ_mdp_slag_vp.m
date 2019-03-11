% Calculate joint inertia matrix for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR6_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:52:55
% EndTime: 2019-03-09 03:52:56
% DurationCPUTime: 0.60s
% Computational Cost: add. (1222->150), mult. (2376->212), div. (0->0), fcn. (2790->10), ass. (0->85)
t154 = sin(pkin(11));
t156 = cos(pkin(11));
t159 = sin(qJ(5));
t162 = cos(qJ(5));
t140 = t154 * t162 + t156 * t159;
t155 = sin(pkin(10));
t157 = cos(pkin(10));
t160 = sin(qJ(3));
t199 = cos(qJ(3));
t141 = t155 * t199 + t160 * t157;
t115 = t140 * t141;
t138 = t154 * t159 - t162 * t156;
t116 = t138 * t141;
t139 = t155 * t160 - t157 * t199;
t149 = -pkin(2) * t157 - pkin(1);
t119 = pkin(3) * t139 - qJ(4) * t141 + t149;
t196 = pkin(7) + qJ(2);
t143 = t196 * t155;
t145 = t196 * t157;
t130 = -t160 * t143 + t145 * t199;
t110 = t154 * t119 + t156 * t130;
t192 = t141 * t154;
t103 = -pkin(8) * t192 + t110;
t109 = t156 * t119 - t130 * t154;
t191 = t141 * t156;
t98 = pkin(4) * t139 - pkin(8) * t191 + t109;
t95 = -t103 * t159 + t162 * t98;
t96 = t103 * t162 + t159 * t98;
t209 = -t116 * MDP(21) - t115 * MDP(22) + t95 * MDP(24) - t96 * MDP(25);
t158 = sin(qJ(6));
t161 = cos(qJ(6));
t106 = t161 * t115 - t116 * t158;
t107 = -t115 * t158 - t116 * t161;
t207 = t107 * MDP(28) - t106 * MDP(29);
t195 = pkin(8) + qJ(4);
t142 = t195 * t154;
t144 = t195 * t156;
t127 = -t162 * t142 - t144 * t159;
t129 = -t142 * t159 + t144 * t162;
t113 = -pkin(9) * t140 + t127;
t114 = -pkin(9) * t138 + t129;
t124 = t161 * t138 + t140 * t158;
t125 = -t138 * t158 + t140 * t161;
t177 = t125 * MDP(28) - t124 * MDP(29) + (t113 * t161 - t114 * t158) * MDP(31) - (t113 * t158 + t114 * t161) * MDP(32);
t206 = t140 * MDP(21) - t138 * MDP(22) + t127 * MDP(24) - t129 * MDP(25) + t177;
t183 = t138 * MDP(24);
t120 = t124 * MDP(31);
t188 = -t125 * MDP(32) - t120;
t205 = -t140 * MDP(25) - t183 + t188;
t187 = t154 ^ 2 + t156 ^ 2;
t204 = MDP(17) * t187;
t148 = -pkin(4) * t156 - pkin(3);
t131 = pkin(5) * t138 + t148;
t203 = 0.2e1 * t131;
t202 = 0.2e1 * t148;
t201 = -2 * MDP(20);
t200 = -2 * MDP(27);
t198 = pkin(1) * MDP(7);
t197 = pkin(5) * t139;
t194 = pkin(3) * MDP(18);
t94 = -pkin(9) * t115 + t96;
t193 = t161 * t94;
t190 = t155 * MDP(5);
t189 = t157 * MDP(4);
t185 = t107 * MDP(26);
t184 = t116 * MDP(19);
t182 = t154 * MDP(16);
t181 = t156 * MDP(15);
t180 = MDP(23) + MDP(30);
t179 = t139 * MDP(30) + t207;
t93 = pkin(9) * t116 + t197 + t95;
t90 = -t158 * t94 + t161 * t93;
t178 = t187 * MDP(18);
t128 = t143 * t199 + t160 * t145;
t175 = t90 * MDP(31) - (t158 * t93 + t193) * MDP(32);
t174 = t109 * t156 + t110 * t154;
t173 = -t109 * t154 + t110 * t156;
t172 = MDP(15) * t154 + MDP(16) * t156;
t170 = t115 * MDP(24) - t116 * MDP(25);
t168 = t106 * MDP(31) + t107 * MDP(32);
t111 = pkin(4) * t192 + t128;
t167 = (MDP(31) * t161 - MDP(32) * t158) * pkin(5);
t165 = t181 - t182 + t205;
t108 = t115 * pkin(5) + t111;
t1 = [(t109 ^ 2 + t110 ^ 2 + t128 ^ 2) * MDP(18) + MDP(1) + (0.2e1 * t149 * MDP(14) + MDP(8) * t141) * t141 + t180 * t139 ^ 2 - (t115 * t201 - t184) * t116 + (t106 * t200 + t185) * t107 + (0.2e1 * t189 - 0.2e1 * t190 + t198) * pkin(1) + 0.2e1 * t170 * t111 + 0.2e1 * t168 * t108 + 0.2e1 * (-MDP(17) * t174 + t128 * t172) * t141 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t155 ^ 2 + t157 ^ 2) * qJ(2) + 0.2e1 * (t149 * MDP(13) + t109 * MDP(15) - t110 * MDP(16) - t141 * MDP(9) + t175 + t207 + t209) * t139; -t189 + t190 - t198 + t174 * MDP(18) + (MDP(14) - t204) * t141 + (MDP(13) + t165) * t139; MDP(7) + t178; t141 * MDP(10) - t128 * MDP(13) - t130 * MDP(14) + (-pkin(3) * t192 - t128 * t156) * MDP(15) + (-pkin(3) * t191 + t128 * t154) * MDP(16) + t173 * MDP(17) + (-pkin(3) * t128 + qJ(4) * t173) * MDP(18) - t140 * t184 + (-t115 * t140 + t116 * t138) * MDP(20) + (t111 * t138 + t115 * t148) * MDP(24) + (t111 * t140 - t116 * t148) * MDP(25) + t125 * t185 + (-t106 * t125 - t107 * t124) * MDP(27) + (t106 * t131 + t108 * t124) * MDP(31) + (t107 * t131 + t108 * t125) * MDP(32) + (-qJ(4) * t172 - MDP(11) + t206) * t139; 0; t183 * t202 + t120 * t203 + MDP(12) + (0.2e1 * t181 - 0.2e1 * t182 + t194) * pkin(3) + (MDP(19) * t140 + MDP(25) * t202 + t138 * t201) * t140 + (MDP(26) * t125 + MDP(32) * t203 + t124 * t200) * t125 + (qJ(4) * t178 + 0.2e1 * t204) * qJ(4); t128 * MDP(18) + t141 * t172 + t168 + t170; 0; -t165 - t194; MDP(18); t139 * MDP(23) + (t161 * t197 + t90) * MDP(31) + (-t193 + (-t93 - t197) * t158) * MDP(32) + t179 + t209; t205; t206; 0; 0.2e1 * t167 + t180; t175 + t179; t188; t177; 0; MDP(30) + t167; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
