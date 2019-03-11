% Calculate joint inertia matrix for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRRP1_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:41:16
% EndTime: 2019-03-09 11:41:18
% DurationCPUTime: 0.54s
% Computational Cost: add. (1170->154), mult. (2133->218), div. (0->0), fcn. (2441->8), ass. (0->75)
t142 = sin(qJ(5));
t145 = cos(qJ(5));
t168 = t142 * MDP(22) + t145 * MDP(23);
t147 = cos(qJ(2));
t133 = -t147 * pkin(2) - pkin(1);
t140 = sin(pkin(10));
t141 = cos(pkin(10));
t144 = sin(qJ(2));
t152 = t140 * t144 - t141 * t147;
t114 = t152 * pkin(3) + t133;
t178 = 0.2e1 * t114;
t177 = 0.2e1 * t147;
t176 = 2 * MDP(27);
t175 = pkin(5) * t145;
t174 = t140 * pkin(2);
t131 = pkin(2) * t141 + pkin(3);
t143 = sin(qJ(4));
t146 = cos(qJ(4));
t117 = t131 * t146 - t143 * t174;
t115 = -pkin(4) - t117;
t173 = pkin(4) - t115;
t172 = -qJ(3) - pkin(7);
t171 = MDP(28) * pkin(5);
t126 = t172 * t144;
t128 = t172 * t147;
t109 = t141 * t126 + t128 * t140;
t120 = t140 * t147 + t141 * t144;
t98 = -pkin(8) * t120 + t109;
t110 = t140 * t126 - t141 * t128;
t99 = -t152 * pkin(8) + t110;
t96 = t143 * t98 + t146 * t99;
t170 = t145 * t96;
t169 = t142 * t145;
t137 = t145 * qJ(6);
t138 = t142 ^ 2;
t139 = t145 ^ 2;
t167 = t138 + t139;
t166 = MDP(23) * t142;
t106 = t120 * t143 + t146 * t152;
t165 = MDP(24) * t106;
t164 = MDP(25) * t145;
t163 = MDP(26) * t142;
t162 = t117 * MDP(18);
t118 = -t131 * t143 - t146 * t174;
t161 = t118 * MDP(19);
t160 = t142 * MDP(27);
t159 = MDP(21) * t169;
t158 = t138 * MDP(20) + MDP(17) + 0.2e1 * t159;
t107 = t146 * t120 - t143 * t152;
t94 = t106 * pkin(4) - t107 * pkin(9) + t114;
t87 = -t142 * t96 + t145 * t94;
t157 = -pkin(4) * t107 - pkin(9) * t106;
t84 = pkin(5) * t106 - t107 * t137 + t87;
t86 = t170 + (-qJ(6) * t107 + t94) * t142;
t156 = t142 * t86 + t145 * t84;
t95 = t143 * t99 - t146 * t98;
t155 = t87 * MDP(25) - (t142 * t94 + t170) * MDP(26);
t116 = pkin(9) - t118;
t154 = -t106 * t116 + t107 * t115;
t111 = (-qJ(6) - t116) * t142;
t112 = t116 * t145 + t137;
t153 = t111 * t145 + t112 * t142;
t151 = -t163 + t164;
t150 = t142 * MDP(25) + t145 * MDP(26);
t149 = -t95 * MDP(18) - t96 * MDP(19) + (MDP(20) * t169 + MDP(15) + (-t138 + t139) * MDP(21)) * t107 + (-MDP(16) + t168) * t106;
t132 = -pkin(4) - t175;
t127 = pkin(9) * t145 + t137;
t125 = (-qJ(6) - pkin(9)) * t142;
t124 = t127 * t145;
t113 = t115 - t175;
t108 = t112 * t145;
t91 = t95 * t142;
t89 = pkin(5) * t107 * t142 + t95;
t85 = t86 * t145;
t1 = [MDP(1) + pkin(1) * MDP(9) * t177 + (t109 ^ 2 + t110 ^ 2 + t133 ^ 2) * MDP(12) + (t84 ^ 2 + t86 ^ 2 + t89 ^ 2) * MDP(28) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t144 + MDP(5) * t177) * t144 + 0.2e1 * (-t109 * t120 - t110 * t152) * MDP(11) + (MDP(19) * t178 - 0.2e1 * t156 * MDP(27) + 0.2e1 * t150 * t95 + (t139 * MDP(20) + MDP(13) - 0.2e1 * t159) * t107) * t107 + (MDP(18) * t178 + t165 + 0.2e1 * (MDP(22) * t145 - MDP(14) - t166) * t107 + 0.2e1 * t155) * t106; (t111 * t84 + t112 * t86 + t113 * t89) * MDP(28) + t144 * MDP(6) + t147 * MDP(7) + t149 + (t154 * t145 + t91) * MDP(26) + (t154 * t142 - t145 * t95) * MDP(25) + (-t153 * t107 - t142 * t84 + t85) * MDP(27) + (-t147 * MDP(10) - t144 * MDP(9)) * pkin(7) + ((t109 * t141 + t110 * t140) * MDP(12) + (-t141 * t120 - t140 * t152) * MDP(11)) * pkin(2); MDP(8) + (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) * MDP(28) + (t140 ^ 2 + t141 ^ 2) * MDP(12) * pkin(2) ^ 2 - 0.2e1 * t151 * t115 + 0.2e1 * t162 + 0.2e1 * t161 + (-t111 * t142 + t108) * t176 + t158; t133 * MDP(12) + t156 * MDP(28) + (-t167 * MDP(27) + MDP(19)) * t107 + (MDP(18) + t151) * t106; t153 * MDP(28); t167 * MDP(28) + MDP(12); t91 * MDP(26) + t85 * MDP(27) + (t125 * t84 + t127 * t86 + t132 * t89) * MDP(28) + (-t107 * t125 * MDP(27) - t95 * MDP(25) + t157 * MDP(26)) * t145 + (t157 * MDP(25) + (-t107 * t127 - t84) * MDP(27)) * t142 + t149; t162 + t161 + (t108 + t124) * MDP(27) + (t111 * t125 + t112 * t127 + t113 * t132) * MDP(28) + t173 * t164 + (-t173 * MDP(26) + (-t111 - t125) * MDP(27)) * t142 + t158; (t125 * t145 + t127 * t142) * MDP(28); (-t125 * t142 + t124) * t176 + (t125 ^ 2 + t127 ^ 2 + t132 ^ 2) * MDP(28) + 0.2e1 * t151 * pkin(4) + t158; t84 * t171 + t165 + (-t166 + (-MDP(27) * pkin(5) + MDP(22)) * t145) * t107 + t155; -t150 * t116 + (t111 * MDP(28) - t160) * pkin(5) + t168; -t163 + (MDP(25) + t171) * t145; -t150 * pkin(9) + (MDP(28) * t125 - t160) * pkin(5) + t168; MDP(28) * pkin(5) ^ 2 + MDP(24); t89 * MDP(28); t113 * MDP(28); 0; t132 * MDP(28); 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
