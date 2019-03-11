% Calculate joint inertia matrix for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP4_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:41
% EndTime: 2019-03-09 11:55:45
% DurationCPUTime: 1.05s
% Computational Cost: add. (1492->198), mult. (2694->290), div. (0->0), fcn. (3006->8), ass. (0->80)
t153 = sin(pkin(10));
t143 = pkin(2) * t153 + pkin(8);
t156 = sin(qJ(4));
t159 = cos(qJ(4));
t186 = pkin(9) + t143;
t130 = t186 * t159;
t155 = sin(qJ(5));
t158 = cos(qJ(5));
t170 = t186 * t156;
t112 = t130 * t155 + t158 * t170;
t113 = t158 * t130 - t155 * t170;
t138 = t155 * t159 + t156 * t158;
t173 = -MDP(26) + MDP(29);
t174 = MDP(25) + MDP(27);
t193 = -t155 * t156 + t158 * t159;
t163 = t138 * MDP(22) + MDP(23) * t193 - t112 * t174 + t113 * t173;
t164 = MDP(18) * t156 + MDP(19) * t159;
t200 = t156 * MDP(15) + t159 * MDP(16) - t143 * t164 + t163;
t197 = -qJ(3) - pkin(7);
t154 = cos(pkin(10));
t157 = sin(qJ(2));
t188 = cos(qJ(2));
t134 = t153 * t188 + t154 * t157;
t104 = t138 * t134;
t105 = t193 * t134;
t196 = MDP(22) * t105 - MDP(23) * t104;
t194 = -MDP(25) * t193 + MDP(26) * t138;
t192 = -2 * MDP(21);
t191 = 2 * MDP(27);
t190 = 2 * MDP(28);
t189 = 2 * MDP(29);
t187 = pkin(4) * t158;
t133 = t153 * t157 - t154 * t188;
t123 = t133 * pkin(5);
t150 = t155 * pkin(4);
t182 = t134 * t159;
t149 = -pkin(2) * t188 - pkin(1);
t114 = pkin(3) * t133 - pkin(8) * t134 + t149;
t140 = t197 * t157;
t141 = t197 * t188;
t117 = t140 * t153 - t141 * t154;
t98 = t114 * t159 - t117 * t156;
t95 = pkin(4) * t133 - pkin(9) * t182 + t98;
t184 = t117 * t159;
t97 = t184 + (-pkin(9) * t134 + t114) * t156;
t91 = t155 * t95 + t158 * t97;
t185 = t105 * t193;
t183 = t134 * t156;
t180 = t156 * t159;
t178 = t105 * MDP(20);
t124 = t193 * MDP(27);
t127 = t138 * MDP(29);
t177 = t155 * MDP(26);
t176 = 0.2e1 * t188;
t175 = MDP(24) + MDP(17);
t122 = t133 * qJ(6);
t88 = t122 + t91;
t172 = MDP(24) * t133 + t196;
t145 = -pkin(2) * t154 - pkin(3);
t171 = MDP(14) * t180;
t169 = t155 * t97 - t158 * t95;
t168 = pkin(5) * t191 + MDP(24);
t115 = -t140 * t154 - t141 * t153;
t167 = t127 + t124 - t194;
t89 = -t123 + t169;
t101 = pkin(4) * t183 + t115;
t166 = pkin(5) * t193 + qJ(6) * t138;
t139 = -pkin(4) * t159 + t145;
t165 = t159 * MDP(18) - t156 * MDP(19);
t162 = (MDP(15) * t159 - MDP(16) * t156) * t134;
t152 = t159 ^ 2;
t151 = t156 ^ 2;
t147 = pkin(5) + t187;
t144 = t150 + qJ(6);
t135 = t138 ^ 2;
t106 = t139 - t166;
t100 = t138 * t104;
t99 = t114 * t156 + t184;
t92 = pkin(5) * t104 - qJ(6) * t105 + t101;
t1 = [(t88 ^ 2 + t89 ^ 2 + t92 ^ 2) * MDP(30) + (t115 ^ 2 + t117 ^ 2 + t149 ^ 2) * MDP(12) + pkin(1) * MDP(9) * t176 + MDP(1) + (MDP(13) * t152 - 0.2e1 * t171) * t134 ^ 2 + t175 * t133 ^ 2 + (t104 * t192 + t178) * t105 + (-0.2e1 * MDP(10) * pkin(1) + MDP(4) * t157 + MDP(5) * t176) * t157 + 0.2e1 * (t162 + t196) * t133 + (-t104 * t88 + t105 * t89) * t190 + 0.2e1 * (t101 * t104 - t133 * t169) * MDP(25) + (-t105 * t92 + t133 * t88) * t189 + (t104 * t92 - t133 * t89) * t191 + 0.2e1 * (t101 * t105 - t133 * t91) * MDP(26) + 0.2e1 * (t115 * t134 - t117 * t133) * MDP(11) + 0.2e1 * (t115 * t183 + t133 * t98) * MDP(18) + 0.2e1 * (t115 * t182 - t133 * t99) * MDP(19); t157 * MDP(6) + t188 * MDP(7) + t138 * t178 + (-t100 + t185) * MDP(21) + (-t101 * t193 + t104 * t139) * MDP(25) + (t101 * t138 + t105 * t139) * MDP(26) + (t104 * t106 - t193 * t92) * MDP(27) + (-t104 * t113 + t105 * t112 + t138 * t89 + t193 * t88) * MDP(28) + (-t105 * t106 - t138 * t92) * MDP(29) + (t106 * t92 + t112 * t89 + t113 * t88) * MDP(30) - t165 * t115 + (-MDP(10) * t188 - MDP(9) * t157) * pkin(7) + (MDP(13) * t180 + (-t151 + t152) * MDP(14) + t164 * t145) * t134 + t200 * t133 + ((-t133 * t153 - t134 * t154) * MDP(11) + (-t115 * t154 + t117 * t153) * MDP(12)) * pkin(2); MDP(8) + t151 * MDP(13) + 0.2e1 * t171 + t135 * MDP(20) - t138 * t193 * t192 + (t112 * t138 + t113 * t193) * t190 + (t112 ^ 2 + t113 ^ 2) * MDP(30) + (t153 ^ 2 + t154 ^ 2) * MDP(12) * pkin(2) ^ 2 + (MDP(30) * t106 - 0.2e1 * t124 - 0.2e1 * t127) * t106 - 0.2e1 * t165 * t145 + 0.2e1 * t194 * t139; t149 * MDP(12) + (-t100 - t185) * MDP(28) + (t138 * t88 - t193 * t89) * MDP(30) + (t138 * t173 + t174 * t193 + t165) * t133; (-t112 * t193 + t113 * t138) * MDP(30); MDP(12) + (t193 ^ 2 + t135) * MDP(30); t133 * MDP(17) + t98 * MDP(18) - t99 * MDP(19) + (t133 * t187 - t169) * MDP(25) + (-t133 * t150 - t91) * MDP(26) + (t133 * t147 - t89) * MDP(27) + (-t104 * t144 - t105 * t147) * MDP(28) + (t133 * t144 + t88) * MDP(29) + (t144 * t88 - t147 * t89) * MDP(30) + t162 + t172; (-t138 * t147 + t144 * t193) * MDP(28) + (-t112 * t147 + t113 * t144) * MDP(30) + t200; (t138 * t144 + t147 * t193) * MDP(30) + t165 + t167; (t144 ^ 2 + t147 ^ 2) * MDP(30) + 0.2e1 * (MDP(25) * t158 - t177) * pkin(4) + t147 * t191 + t144 * t189 + t175; -t169 * MDP(25) - t91 * MDP(26) + (-t89 + t123) * MDP(27) + (-pkin(5) * t105 - qJ(6) * t104) * MDP(28) + (0.2e1 * t122 + t91) * MDP(29) + (-pkin(5) * t89 + qJ(6) * t88) * MDP(30) + t172; (-pkin(5) * t138 + qJ(6) * t193) * MDP(28) + (-pkin(5) * t112 + qJ(6) * t113) * MDP(30) + t163; MDP(30) * t166 + t167; (0.2e1 * qJ(6) + t150) * MDP(29) + (pkin(5) * t147 + qJ(6) * t144) * MDP(30) + (t158 * t174 - t177) * pkin(4) + t168; qJ(6) * t189 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(30) + t168; -MDP(27) * t133 + MDP(28) * t105 + MDP(30) * t89; MDP(28) * t138 + MDP(30) * t112; -t193 * MDP(30); -MDP(30) * t147 - MDP(27); -MDP(30) * pkin(5) - MDP(27); MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
