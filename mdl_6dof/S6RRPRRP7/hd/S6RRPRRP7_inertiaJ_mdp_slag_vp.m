% Calculate joint inertia matrix for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP7_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:20:09
% EndTime: 2019-03-09 12:20:11
% DurationCPUTime: 0.64s
% Computational Cost: add. (883->182), mult. (1460->237), div. (0->0), fcn. (1413->6), ass. (0->82)
t138 = sin(qJ(5));
t133 = t138 ^ 2;
t141 = cos(qJ(5));
t136 = t141 ^ 2;
t186 = t133 + t136;
t164 = t186 * MDP(32);
t202 = pkin(9) * t164;
t140 = sin(qJ(2));
t143 = cos(qJ(2));
t123 = -t143 * pkin(2) - t140 * qJ(3) - pkin(1);
t112 = t143 * pkin(3) - t123;
t139 = sin(qJ(4));
t142 = cos(qJ(4));
t114 = t140 * t139 + t143 * t142;
t115 = -t143 * t139 + t140 * t142;
t106 = t114 * pkin(4) - t115 * pkin(9) + t112;
t195 = pkin(7) - pkin(8);
t125 = t195 * t140;
t126 = t195 * t143;
t110 = t139 * t125 + t142 * t126;
t102 = t138 * t106 + t141 * t110;
t190 = t114 * qJ(6);
t99 = t102 + t190;
t201 = t99 * t141;
t144 = -pkin(2) - pkin(3);
t121 = t142 * qJ(3) + t139 * t144;
t119 = -pkin(9) + t121;
t162 = t119 * t164;
t200 = t186 * MDP(30);
t173 = MDP(27) + MDP(29);
t172 = MDP(28) - MDP(31);
t199 = -t172 * t138 + t173 * t141 + MDP(20);
t197 = 2 * MDP(29);
t196 = 2 * MDP(31);
t194 = pkin(9) * t114;
t193 = t114 * pkin(5);
t120 = t139 * qJ(3) - t142 * t144;
t118 = pkin(4) + t120;
t192 = pkin(4) + t118;
t191 = MDP(32) * pkin(9);
t189 = t114 * t119;
t188 = t115 * t141;
t161 = -t141 * pkin(5) - t138 * qJ(6);
t122 = -pkin(4) + t161;
t111 = t120 - t122;
t187 = -t111 + t122;
t135 = t140 ^ 2;
t185 = t143 ^ 2 + t135;
t184 = MDP(32) * t139;
t167 = -t141 * t106 + t138 * t110;
t100 = t167 - t193;
t183 = t100 * MDP(32);
t182 = t111 * MDP(32);
t181 = t114 * MDP(26);
t180 = t139 * t200;
t179 = t119 * MDP(32);
t178 = t120 * MDP(20);
t177 = t121 * MDP(21);
t176 = t122 * MDP(32);
t175 = t141 * MDP(23);
t174 = t133 * MDP(22) + MDP(19);
t171 = t138 * t175;
t170 = 0.2e1 * t171 + t174;
t169 = -pkin(2) * MDP(14) - MDP(11);
t168 = -MDP(32) * pkin(5) - MDP(29);
t163 = -pkin(4) * t115 - t194;
t160 = pkin(5) * t138 - t141 * qJ(6);
t159 = -t115 * t122 + t194;
t158 = t111 * t115 - t189;
t157 = t115 * t118 - t189;
t109 = -t142 * t125 + t139 * t126;
t155 = t141 * MDP(24) - t138 * MDP(25);
t154 = -t167 * MDP(27) - t102 * MDP(28);
t153 = t138 * t196 + t141 * t197;
t152 = 0.2e1 * t141 * MDP(27) - 0.2e1 * t138 * MDP(28);
t151 = -t138 * MDP(24) - t141 * MDP(25) + t160 * MDP(30);
t103 = t160 * t115 + t109;
t150 = -t114 * MDP(25) + t109 * MDP(27) + t103 * MDP(29);
t149 = -MDP(22) * t188 - t114 * MDP(24) - t109 * MDP(28) + t103 * MDP(31);
t148 = (MDP(32) * qJ(6) - t172) * t141 + (-MDP(27) + t168) * t138;
t147 = -t114 * MDP(18) - t109 * MDP(20) - t110 * MDP(21) + (t100 * t138 + t201) * MDP(30) + (MDP(17) - (t133 - t136) * MDP(23)) * t115;
t1 = [(t100 ^ 2 + t103 ^ 2 + t99 ^ 2) * MDP(32) + t135 * MDP(4) + (t185 * pkin(7) ^ 2 + t123 ^ 2) * MDP(14) + MDP(1) + (t136 * MDP(22) + MDP(15) - 0.2e1 * t171) * t115 ^ 2 + (0.2e1 * t112 * MDP(20) + t181) * t114 + 0.2e1 * (-t100 * MDP(29) + t99 * MDP(31) + t154) * t114 + 0.2e1 * t185 * MDP(12) * pkin(7) + 0.2e1 * (-t123 * MDP(11) + pkin(1) * MDP(9)) * t143 + 0.2e1 * (-pkin(1) * MDP(10) - t123 * MDP(13) + t143 * MDP(5)) * t140 + 0.2e1 * (t112 * MDP(21) + (t100 * t141 - t138 * t99) * MDP(30) + (t138 * MDP(27) + t141 * MDP(28)) * t109 + (t138 * MDP(29) - t141 * MDP(31)) * t103 + (-MDP(16) + t155) * t114) * t115; t140 * MDP(6) + t143 * MDP(7) + (-t140 * pkin(2) + t143 * qJ(3)) * MDP(12) + t103 * t182 + (t157 * MDP(28) - t158 * MDP(31) + t99 * t179 + t150) * t141 + (t157 * MDP(27) + t158 * MDP(29) + t100 * t179 + t149) * t138 + ((MDP(14) * qJ(3) - MDP(10) + MDP(13)) * t143 + (-MDP(9) + t169) * t140) * pkin(7) - t147; MDP(8) + 0.2e1 * pkin(2) * MDP(11) + 0.2e1 * qJ(3) * MDP(13) + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(14) + t118 * t152 + (t153 + t182) * t111 + 0.2e1 * t178 + 0.2e1 * t177 + t170 + (-0.2e1 * t200 + t162) * t119; -t103 * t142 * MDP(32) + (pkin(7) * MDP(14) + MDP(12)) * t140 + t184 * t201 + t139 * t183 * t138 + (-t138 * t173 - t141 * t172) * (t114 * t139 + t115 * t142); -t180 + (MDP(21) + t162) * t139 + (-t182 - t199) * t142 + t169; MDP(14) + (t186 * t139 ^ 2 + t142 ^ 2) * MDP(32); t103 * t176 + (t163 * MDP(28) + t159 * MDP(31) + t99 * t191 - t150) * t141 + (t163 * MDP(27) - t159 * MDP(29) + pkin(9) * t183 - t149) * t138 + t147; t111 * t176 - t178 - t177 + t119 * t200 + (-t200 + t162) * pkin(9) + (-t192 * MDP(27) + t187 * MDP(29)) * t141 + (t192 * MDP(28) + t187 * MDP(31) - 0.2e1 * t175) * t138 - t174; t180 + (-MDP(21) + t202) * t139 + (-t176 + t199) * t142; (-t153 + t176) * t122 + pkin(4) * t152 + t170 + (0.2e1 * t200 + t202) * pkin(9); t181 + (-t167 + 0.2e1 * t193) * MDP(29) + (t102 + 0.2e1 * t190) * MDP(31) + (-t100 * pkin(5) + t99 * qJ(6)) * MDP(32) + (t161 * MDP(30) + t155) * t115 + t154; t119 * t148 + t151; t148 * t139; pkin(9) * t148 - t151; MDP(26) + pkin(5) * t197 + qJ(6) * t196 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32); -t114 * MDP(29) + MDP(30) * t188 + t183; (-MDP(30) + t179) * t138; t138 * t184; (MDP(30) + t191) * t138; t168; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
