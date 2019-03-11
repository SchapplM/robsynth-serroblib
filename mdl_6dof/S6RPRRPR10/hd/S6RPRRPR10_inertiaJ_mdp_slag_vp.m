% Calculate joint inertia matrix for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRPR10_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:37:12
% EndTime: 2019-03-09 05:37:14
% DurationCPUTime: 0.73s
% Computational Cost: add. (506->161), mult. (890->218), div. (0->0), fcn. (784->6), ass. (0->74)
t122 = sin(qJ(6));
t125 = cos(qJ(6));
t124 = sin(qJ(3));
t127 = cos(qJ(3));
t109 = t124 * pkin(3) - t127 * pkin(8) + qJ(2);
t126 = cos(qJ(4));
t123 = sin(qJ(4));
t128 = -pkin(1) - pkin(7);
t169 = t123 * t128;
t161 = -t126 * t109 + t124 * t169;
t166 = t126 * t127;
t176 = pkin(4) + pkin(5);
t85 = -pkin(9) * t166 - t176 * t124 + t161;
t117 = t124 * qJ(5);
t165 = t126 * t128;
t93 = t123 * t109 + t124 * t165;
t88 = t117 + t93;
t86 = t123 * t127 * pkin(9) + t88;
t168 = t125 * t123;
t97 = t122 * t166 - t127 * t168;
t102 = t122 * t123 + t125 * t126;
t98 = t102 * t127;
t185 = t98 * MDP(27) - t97 * MDP(28) - (t122 * t86 - t125 * t85) * MDP(30) - (t122 * t85 + t125 * t86) * MDP(31);
t184 = -t161 * MDP(19) - t93 * MDP(20) - t185;
t104 = -t122 * t126 + t168;
t172 = pkin(8) * MDP(24);
t181 = MDP(22) + t172;
t150 = -MDP(20) + MDP(23);
t151 = MDP(19) + MDP(21);
t180 = -t151 * t123 + t150 * t126;
t170 = t123 * qJ(5);
t99 = t176 * t126 + pkin(3) + t170;
t179 = 0.2e1 * t99;
t178 = 0.2e1 * t127;
t177 = -2 * MDP(26);
t175 = pkin(8) - pkin(9);
t174 = (pkin(1) * MDP(6));
t173 = pkin(4) * MDP(24);
t167 = t126 * qJ(5);
t118 = t123 ^ 2;
t120 = t126 ^ 2;
t160 = t118 + t120;
t158 = MDP(25) * t104;
t157 = t102 * MDP(30);
t156 = (t122 * qJ(5) + t125 * t176) * MDP(30);
t155 = (t125 * qJ(5) - t122 * t176) * MDP(31);
t145 = -t126 * pkin(4) - t170;
t110 = -pkin(3) + t145;
t154 = t110 * MDP(24);
t153 = t126 * MDP(15);
t152 = MDP(18) + MDP(29);
t149 = 0.2e1 * pkin(4) * MDP(21);
t148 = t175 * t123;
t147 = t160 * MDP(22);
t144 = -pkin(4) * t123 + t167;
t89 = -t124 * pkin(4) + t161;
t142 = t89 * t123 + t88 * t126;
t95 = t104 * t124;
t96 = t102 * t124;
t138 = t95 * MDP(30) - t96 * MDP(31);
t137 = t126 * MDP(16) - t123 * MDP(17);
t136 = t123 * MDP(21) - t126 * MDP(23);
t135 = t125 * MDP(30) - t122 * MDP(31);
t134 = -MDP(21) - t135;
t133 = -MDP(29) - t155 - t156;
t111 = t175 * t126;
t132 = t104 * MDP(27) - t102 * MDP(28) - (t122 * t111 - t125 * t148) * MDP(30) - (t125 * t111 + t122 * t148) * MDP(31);
t131 = (MDP(24) * qJ(5) + t150) * t126 + (-t151 - t173) * t123;
t130 = t123 * MDP(16) + t126 * MDP(17) - t132;
t121 = t127 ^ 2;
t119 = t124 ^ 2;
t94 = (-t128 - t144) * t127;
t87 = (-t176 * t123 + t128 + t167) * t127;
t1 = [(t88 ^ 2 + t89 ^ 2 + t94 ^ 2) * MDP(24) + MDP(1) + (t98 * MDP(25) + t97 * t177) * t98 + t152 * t119 + ((-2 * MDP(4) + t174) * pkin(1)) + (MDP(13) * t178 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (t97 * MDP(30) + t98 * MDP(31)) * t87 + ((-t123 * t88 + t126 * t89) * MDP(22) + t136 * t94) * t178 + (t120 * MDP(14) - 0.2e1 * t123 * t153 + MDP(7) + 0.2e1 * (-t123 * MDP(19) - t126 * MDP(20)) * t128) * t121 + 0.2e1 * (qJ(2) * MDP(12) + (-MDP(8) + t137) * t127 - t89 * MDP(21) + t88 * MDP(23) + t184) * t124; MDP(4) - t174 + (t124 * t142 - t94 * t127) * MDP(24) + (-t95 * t124 + t127 * t97) * MDP(30) + (t96 * t124 + t127 * t98) * MDP(31) + t180 * (t119 + t121); MDP(6) + (t160 * t119 + t121) * MDP(24); t98 * t158 + (-t98 * t102 - t104 * t97) * MDP(26) + (t87 * t102 + t99 * t97) * MDP(30) + (t87 * t104 + t99 * t98) * MDP(31) + (-t126 * MDP(21) - t123 * MDP(23) + t154) * t94 + (-t128 * MDP(13) + t180 * pkin(8) - MDP(10) + t130) * t124 + (MDP(9) + t128 * MDP(12) + t126 * t123 * MDP(14) + (-t118 + t120) * MDP(15) + (-pkin(3) * t123 + t165) * MDP(19) + (-pkin(3) * t126 - t169) * MDP(20) + t136 * t110) * t127 + t181 * t142; (t160 * t172 - MDP(13) + t147) * t124 + (t104 * MDP(31) + t150 * t123 + t151 * t126 + MDP(12) - t154 + t157) * t127; MDP(11) + t118 * MDP(14) + (t160 * pkin(8) ^ 2 + t110 ^ 2) * MDP(24) + t157 * t179 + 0.2e1 * pkin(8) * t147 + (MDP(31) * t179 + t102 * t177 + t158) * t104 + 0.2e1 * (pkin(3) * MDP(19) - t110 * MDP(21)) * t126 + 0.2e1 * (-pkin(3) * MDP(20) - t110 * MDP(23) + t153) * t123; -t161 * MDP(21) + (0.2e1 * t117 + t93) * MDP(23) + (-t89 * pkin(4) + t88 * qJ(5)) * MDP(24) + (MDP(18) - t133 + t149) * t124 + (MDP(22) * t145 + t137) * t127 + t184; t124 * t131 - t138; MDP(22) * t144 + pkin(8) * t131 + t130; t149 + 0.2e1 * qJ(5) * MDP(23) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(24) + 0.2e1 * t156 + 0.2e1 * t155 + t152; MDP(22) * t166 + t89 * MDP(24) + t124 * t134; t123 * t124 * MDP(24); t181 * t123; t134 - t173; MDP(24); -t124 * MDP(29) + t185; t138; t132; t133; t135; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
