% Calculate joint inertia matrix for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR3_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:13:25
% EndTime: 2019-03-09 18:13:27
% DurationCPUTime: 0.61s
% Computational Cost: add. (863->168), mult. (1473->216), div. (0->0), fcn. (1582->8), ass. (0->83)
t153 = sin(qJ(6));
t182 = t153 * MDP(35);
t205 = -MDP(27) + t182;
t204 = MDP(16) + MDP(18);
t157 = cos(qJ(6));
t180 = t157 * MDP(34);
t203 = t180 - t182;
t202 = MDP(34) * t153 + MDP(35) * t157;
t155 = sin(qJ(3));
t156 = sin(qJ(2));
t159 = cos(qJ(3));
t160 = cos(qJ(2));
t122 = t155 * t156 - t159 * t160;
t123 = t155 * t160 + t159 * t156;
t154 = sin(qJ(5));
t158 = cos(qJ(5));
t106 = t154 * t122 + t158 * t123;
t151 = t153 ^ 2;
t152 = t157 ^ 2;
t195 = -pkin(8) - pkin(7);
t132 = t195 * t156;
t133 = t195 * t160;
t112 = t155 * t132 - t159 * t133;
t100 = t122 * pkin(9) + t112;
t111 = -t159 * t132 - t155 * t133;
t99 = -t123 * pkin(9) + t111;
t94 = t154 * t100 - t158 * t99;
t91 = t94 * t157;
t95 = t158 * t100 + t154 * t99;
t201 = -t95 * MDP(28) - t91 * MDP(34) - (t151 - t152) * t106 * MDP(30);
t200 = -t154 * MDP(28) + (t180 - t205) * t158;
t142 = -t160 * pkin(2) - pkin(1);
t103 = t122 * pkin(3) - t123 * qJ(4) + t142;
t98 = -t122 * pkin(4) - t103;
t199 = 0.2e1 * t98;
t198 = 0.2e1 * t160;
t197 = 2 * MDP(18);
t196 = 2 * MDP(20);
t140 = t159 * pkin(2) + pkin(3);
t137 = -pkin(4) - t140;
t148 = t155 * pkin(2);
t138 = t148 + qJ(4);
t116 = t154 * t137 + t158 * t138;
t161 = -pkin(3) - pkin(4);
t128 = t158 * qJ(4) + t154 * t161;
t191 = MDP(29) * t157;
t105 = -t158 * t122 + t154 * t123;
t188 = t105 * MDP(32);
t187 = t105 * MDP(33);
t130 = t154 * t138;
t115 = -t158 * t137 + t130;
t186 = t115 * MDP(27);
t185 = t116 * MDP(28);
t146 = t154 * qJ(4);
t127 = -t158 * t161 + t146;
t184 = t127 * MDP(27);
t183 = t128 * MDP(28);
t181 = t155 * MDP(17);
t179 = t151 * MDP(29) + MDP(26);
t178 = t153 * t157 * MDP(30);
t177 = 0.2e1 * t178 + t179;
t135 = -0.2e1 * t178;
t176 = t135 - t179;
t174 = MDP(15) + t177;
t113 = pkin(5) + t115;
t114 = -pkin(10) + t116;
t173 = -t105 * t114 + t106 * t113;
t125 = pkin(5) + t127;
t126 = -pkin(10) + t128;
t172 = -t105 * t126 + t106 * t125;
t170 = MDP(31) * t157 - MDP(32) * t153;
t169 = -t153 * MDP(31) - t157 * MDP(32);
t167 = pkin(3) * t197 + t174;
t166 = 0.2e1 * t203;
t165 = -t105 * MDP(31) - t94 * MDP(35) - t106 * t191;
t164 = -MDP(18) - t200;
t163 = -pkin(10) * t202 - t169;
t162 = t123 * MDP(13) - t122 * MDP(14) - t106 * MDP(24) + t105 * MDP(25) + t94 * MDP(27) - t201 + (-MDP(17) + MDP(20)) * t112 - t204 * t111;
t150 = pkin(5) * t153;
t90 = t105 * pkin(5) - t106 * pkin(10) + t98;
t89 = t153 * t90 + t157 * t95;
t88 = -t153 * t95 + t157 * t90;
t1 = [(t103 ^ 2 + t111 ^ 2 + t112 ^ 2) * MDP(21) + t106 * MDP(28) * t199 + pkin(1) * MDP(9) * t198 + MDP(1) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t156 + MDP(5) * t198) * t156 + (t152 * MDP(29) + MDP(22) + t135) * t106 ^ 2 + (MDP(11) * t123 - 0.2e1 * t122 * MDP(12) + 0.2e1 * t142 * MDP(17) - 0.2e1 * t103 * MDP(20)) * t123 + (MDP(27) * t199 + t187) * t105 + 0.2e1 * (t111 * t123 - t112 * t122) * MDP(19) + 0.2e1 * (t94 * t153 * t106 + t88 * t105) * MDP(34) + 0.2e1 * (-t89 * t105 + t106 * t91) * MDP(35) + 0.2e1 * (t142 * MDP(16) + t103 * MDP(18)) * t122 + 0.2e1 * (-MDP(23) + t170) * t106 * t105; t156 * MDP(6) + t160 * MDP(7) + (-t138 * t122 - t140 * t123) * MDP(19) + (-t111 * t140 + t112 * t138) * MDP(21) + t162 + (t173 * MDP(34) + t165) * t153 + (t173 * MDP(35) - t188) * t157 + (-t160 * MDP(10) - t156 * MDP(9)) * pkin(7); MDP(8) + (t138 ^ 2 + t140 ^ 2) * MDP(21) + t113 * t166 + 0.2e1 * (t159 * MDP(16) - t181) * pkin(2) + t140 * t197 + t138 * t196 + 0.2e1 * t186 + 0.2e1 * t185 + t174; (t172 * MDP(34) + t165) * t153 + (-pkin(3) * t123 - t122 * qJ(4)) * MDP(19) + (-t111 * pkin(3) + t112 * qJ(4)) * MDP(21) + (t172 * MDP(35) - t188) * t157 + t162; (0.2e1 * qJ(4) + t148) * MDP(20) + (t140 * pkin(3) + t138 * qJ(4)) * MDP(21) + (t130 + t146 + (-t137 - t161) * t158) * MDP(27) + (t128 + t116) * MDP(28) + (t204 * t159 - t181) * pkin(2) + t167 + t203 * (t113 + t125); qJ(4) * t196 + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(21) + t125 * t166 + 0.2e1 * t184 + 0.2e1 * t183 + t167; t123 * MDP(19) + t111 * MDP(21) + t202 * (-t105 * t154 - t106 * t158); -t140 * MDP(21) + t164; -pkin(3) * MDP(21) + t164; MDP(21); t205 * t94 + (-pkin(5) * t202 + t153 * t191 + MDP(24)) * t106 + (-MDP(25) + t163) * t105 + t201; -t186 - t185 + (t113 * t153 + t150) * MDP(35) + (-pkin(5) - t113) * t180 + t176; -t184 - t183 + (t125 * t153 + t150) * MDP(35) + (-pkin(5) - t125) * t180 + t176; t200; pkin(5) * t166 + t177; t88 * MDP(34) - t89 * MDP(35) + t106 * t170 + t187; -t114 * t202 + t169; -t126 * t202 + t169; -t202 * t154; t163; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
