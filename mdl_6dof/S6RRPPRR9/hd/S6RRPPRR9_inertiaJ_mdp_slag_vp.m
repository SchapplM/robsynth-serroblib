% Calculate joint inertia matrix for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR9_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:32:09
% EndTime: 2019-03-09 09:32:11
% DurationCPUTime: 0.78s
% Computational Cost: add. (655->200), mult. (1415->267), div. (0->0), fcn. (1365->8), ass. (0->89)
t142 = cos(qJ(5));
t194 = 0.2e1 * t142;
t133 = sin(pkin(6));
t143 = cos(qJ(2));
t185 = t133 * t143;
t138 = sin(qJ(6));
t141 = cos(qJ(6));
t193 = t138 * MDP(31) + t141 * MDP(32);
t192 = 0.2e1 * MDP(31);
t191 = 0.2e1 * MDP(32);
t134 = cos(pkin(6));
t190 = pkin(1) * t134;
t136 = pkin(2) + qJ(4);
t135 = qJ(3) - pkin(9);
t140 = sin(qJ(2));
t156 = -t136 * t143 - pkin(1);
t100 = (-t135 * t140 + t156) * t133;
t139 = sin(qJ(5));
t114 = pkin(8) * t185 + t140 * t190;
t126 = t134 * qJ(3);
t108 = -t126 - t114;
t120 = pkin(3) * t185;
t105 = t120 - t108;
t99 = pkin(4) * t185 - t134 * pkin(9) + t105;
t95 = -t139 * t100 + t142 * t99;
t93 = -pkin(5) * t185 - t95;
t189 = t93 * t138;
t188 = t93 * t141;
t187 = qJ(3) * t140;
t186 = t133 * t140;
t184 = t134 * MDP(8);
t183 = t135 * t138;
t182 = t135 * t141;
t181 = t138 * t139;
t180 = t139 * t141;
t179 = pkin(8) * t186 - t143 * t190;
t177 = MDP(23) * t143;
t112 = t134 * t142 + t139 * t186;
t103 = t112 * t141 + t138 * t185;
t176 = MDP(28) * t103;
t102 = t112 * t138 - t141 * t185;
t175 = MDP(29) * t102;
t111 = t134 * t139 - t142 * t186;
t174 = MDP(30) * t111;
t173 = t103 * MDP(26);
t172 = t112 * MDP(20);
t171 = t112 * MDP(21);
t170 = t112 * MDP(25);
t169 = t135 * MDP(25);
t168 = t136 * MDP(24);
t167 = t136 * MDP(25);
t165 = t139 * MDP(25);
t164 = t141 * MDP(26);
t162 = MDP(12) - MDP(17);
t127 = t134 * pkin(2);
t110 = -t127 + t179;
t161 = 0.2e1 * t126 + t114;
t160 = t138 * t141 * MDP(27);
t159 = t133 * (MDP(11) + MDP(15));
t158 = t135 * MDP(24) + MDP(21);
t157 = t134 * qJ(4) - t110;
t96 = t142 * t100 + t139 * t99;
t155 = t95 * MDP(24) - t96 * MDP(25);
t154 = -t114 * MDP(10) - t179 * MDP(9);
t153 = MDP(28) * t141 - MDP(29) * t138;
t115 = t139 * pkin(5) - t142 * pkin(10) + t136;
t106 = t141 * t115 - t135 * t181;
t107 = t138 * t115 + t135 * t180;
t152 = t106 * MDP(31) - t107 * MDP(32);
t151 = t141 * MDP(31) - t138 * MDP(32);
t149 = -MDP(20) + t153;
t148 = MDP(24) + t151;
t98 = (-pkin(3) - pkin(4)) * t186 + t157;
t147 = t138 * MDP(28) + t141 * MDP(29) - pkin(10) * t193;
t94 = pkin(10) * t185 + t96;
t97 = t111 * pkin(5) - t112 * pkin(10) + t98;
t91 = -t138 * t94 + t141 * t97;
t92 = t138 * t97 + t141 * t94;
t146 = MDP(31) * t91 - MDP(32) * t92 + t174 - t175 + t176;
t145 = -MDP(22) + t147;
t144 = qJ(3) ^ 2;
t132 = t142 ^ 2;
t131 = t141 ^ 2;
t130 = t139 ^ 2;
t129 = t138 ^ 2;
t109 = (-pkin(2) * t143 - pkin(1) - t187) * t133;
t104 = (t156 - t187) * t133;
t101 = pkin(3) * t186 - t157;
t1 = [(t101 ^ 2 + t104 ^ 2 + t105 ^ 2) * MDP(18) + (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) * MDP(14) + t112 ^ 2 * MDP(19) + MDP(1) + (0.2e1 * MDP(6) * t186 + t184) * t134 + (-0.2e1 * t102 * MDP(27) + t173) * t103 + 0.2e1 * (MDP(7) * t134 + t171) * t185 + (-0.2e1 * MDP(22) * t185 - 0.2e1 * t172 + t174 - 0.2e1 * t175 + 0.2e1 * t176) * t111 + (t102 * t93 + t111 * t91) * t192 + (t103 * t93 - t111 * t92) * t191 + 0.2e1 * (t111 * MDP(24) + t170) * t98 + 0.2e1 * (t110 * MDP(12) - t108 * MDP(13) + t105 * MDP(16) - t101 * MDP(17) + t154) * t134 + 0.2e1 * ((t110 * MDP(11) - t109 * MDP(13) + t101 * MDP(15) - t104 * MDP(16)) * t140 + (-t108 * MDP(11) + t109 * MDP(12) + t105 * MDP(15) - t104 * MDP(17) + t155) * t143) * t133 + (t140 ^ 2 * MDP(4) + (0.2e1 * MDP(5) * t140 + t177) * t143 + 0.2e1 * (-MDP(10) * t140 + MDP(9) * t143) * pkin(1)) * t133 ^ 2; (-pkin(2) * t110 - qJ(3) * t108) * MDP(14) + (t120 + t161) * MDP(16) + t161 * MDP(13) + (-0.2e1 * t127 + t179) * MDP(12) + t184 + (t105 * qJ(3) - t101 * t136) * MDP(18) + (t136 * t134 + t157) * MDP(17) + t112 * t167 + t143 * qJ(3) * t159 + (t143 * MDP(7) + (-MDP(11) * pkin(2) - MDP(15) * t136 - MDP(17) * pkin(3) + MDP(6)) * t140) * t133 + (t152 + t168) * t111 + (-t172 + t98 * MDP(24) + (-MDP(22) - t169) * t185 + t146) * t139 + (t103 * t164 + (-t102 * t135 + t189) * MDP(31) + (-t102 * t141 - t103 * t138) * MDP(27) + (-t103 * t135 + t188) * MDP(32) + t98 * MDP(25) + t112 * MDP(19) + t158 * t185 + t149 * t111) * t142 + t154; MDP(8) - 0.2e1 * pkin(2) * MDP(12) + (pkin(2) ^ 2 + t144) * MDP(14) + (t136 ^ 2 + t144) * MDP(18) + t130 * MDP(30) + 0.2e1 * t136 * MDP(17) + t167 * t194 + 0.2e1 * (MDP(16) + MDP(13)) * qJ(3) + (t131 * MDP(26) - t182 * t191 - t183 * t192 + MDP(19) - 0.2e1 * t160) * t132 + (t106 * t192 - t107 * t191 + t149 * t194 + 0.2e1 * t168) * t139; t110 * MDP(14) + t101 * MDP(18) - t148 * t111 + t162 * t134 + t140 * t159 - t170; -pkin(2) * MDP(14) - t136 * MDP(18) - t142 * MDP(25) - t148 * t139 + t162; MDP(14) + MDP(18); t134 * MDP(16) + t105 * MDP(18) + (-t142 * t102 - t111 * t181) * MDP(31) + (-t142 * t103 - t111 * t180) * MDP(32) + (t142 * MDP(24) + MDP(15) - t165) * t185; qJ(3) * MDP(18) + MDP(16) + t193 * (-t130 - t132); 0; MDP(18); t171 + t133 * t177 + t138 * t173 + (-t138 * t102 + t103 * t141) * MDP(27) + (-pkin(5) * t102 - t188) * MDP(31) + (-pkin(5) * t103 + t189) * MDP(32) + t145 * t111 + t155; (t145 - t169) * t139 + (t138 * t164 + (-t129 + t131) * MDP(27) + (-pkin(5) * t138 + t182) * MDP(31) + (-pkin(5) * t141 - t183) * MDP(32) + t158) * t142; 0; t148 * t142 - t165; t129 * MDP(26) + 0.2e1 * pkin(5) * t151 + MDP(23) + 0.2e1 * t160; t146; t139 * MDP(30) + t153 * t142 + t152; -t151; -t193 * t139; t147; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
