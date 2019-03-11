% Calculate joint inertia matrix for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRRPRP1_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:33:03
% EndTime: 2019-03-09 16:33:05
% DurationCPUTime: 0.58s
% Computational Cost: add. (1241->160), mult. (2240->233), div. (0->0), fcn. (2538->8), ass. (0->76)
t145 = sin(qJ(5));
t148 = cos(qJ(5));
t174 = t145 * MDP(22) + t148 * MDP(23);
t146 = sin(qJ(3));
t147 = sin(qJ(2));
t149 = cos(qJ(3));
t150 = cos(qJ(2));
t125 = t146 * t150 + t147 * t149;
t143 = sin(pkin(10));
t144 = cos(pkin(10));
t158 = t146 * t147 - t149 * t150;
t107 = t144 * t125 - t143 * t158;
t185 = 0.2e1 * t107;
t168 = t145 * MDP(26);
t136 = -t150 * pkin(2) - pkin(1);
t184 = 0.2e1 * t136;
t183 = 2 * MDP(27);
t182 = pkin(7) + pkin(8);
t181 = pkin(2) * t146;
t180 = pkin(5) * t148;
t179 = MDP(28) * pkin(5);
t128 = t182 * t147;
t129 = t182 * t150;
t159 = t146 * t128 - t149 * t129;
t102 = -t158 * qJ(4) - t159;
t160 = -t128 * t149 - t129 * t146;
t153 = -qJ(4) * t125 + t160;
t98 = t144 * t102 + t143 * t153;
t178 = t148 * t98;
t177 = t145 * t148;
t140 = t148 * qJ(6);
t96 = t102 * t143 - t144 * t153;
t176 = t96 * MDP(25);
t135 = pkin(2) * t149 + pkin(3);
t118 = t135 * t144 - t143 * t181;
t115 = -pkin(4) - t118;
t134 = -pkin(3) * t144 - pkin(4);
t175 = t115 + t134;
t119 = t143 * t135 + t144 * t181;
t141 = t145 ^ 2;
t142 = t148 ^ 2;
t173 = t141 + t142;
t172 = MDP(23) * t145;
t106 = t125 * t143 + t144 * t158;
t171 = MDP(24) * t106;
t170 = MDP(27) * t107;
t169 = MDP(27) * t145;
t167 = t148 * MDP(25);
t166 = MDP(21) * t177;
t165 = t141 * MDP(20) + MDP(15) + 0.2e1 * t166;
t114 = t158 * pkin(3) + t136;
t95 = t106 * pkin(4) - t107 * pkin(9) + t114;
t90 = -t145 * t98 + t148 * t95;
t87 = pkin(5) * t106 - t107 * t140 + t90;
t89 = t178 + (-qJ(6) * t107 + t95) * t145;
t164 = t145 * t89 + t148 * t87;
t163 = t90 * MDP(25) - (t145 * t95 + t178) * MDP(26);
t116 = pkin(9) + t119;
t162 = -t106 * t116 + t107 * t115;
t133 = pkin(3) * t143 + pkin(9);
t161 = -t106 * t133 + t107 * t134;
t157 = t167 - t168;
t156 = -t145 * MDP(25) - t148 * MDP(26);
t155 = (MDP(16) * t149 - MDP(17) * t146) * pkin(2);
t154 = 0.2e1 * t157;
t152 = t89 * t148 * MDP(27) + t125 * MDP(13) - t158 * MDP(14) + t160 * MDP(16) + t159 * MDP(17) + t96 * t168 + (MDP(20) * t177 + (-t141 + t142) * MDP(21)) * t107 + t174 * t106;
t126 = t134 - t180;
t122 = t133 * t148 + t140;
t121 = (-qJ(6) - t133) * t145;
t117 = t122 * t148;
t113 = t115 - t180;
t112 = t116 * t148 + t140;
t111 = (-qJ(6) - t116) * t145;
t108 = t112 * t148;
t92 = pkin(5) * t107 * t145 + t96;
t1 = [MDP(1) + t158 * MDP(16) * t184 + (t114 ^ 2 + t96 ^ 2 + t98 ^ 2) * MDP(19) + (t87 ^ 2 + t89 ^ 2 + t92 ^ 2) * MDP(28) + (t142 * MDP(20) - 0.2e1 * t166) * t107 ^ 2 + (MDP(11) * t125 - 0.2e1 * t158 * MDP(12) + MDP(17) * t184) * t125 + (-t164 * MDP(27) + (MDP(18) - t156) * t96) * t185 + (MDP(4) * t147 + 0.2e1 * t150 * MDP(5)) * t147 + 0.2e1 * (-t147 * MDP(10) + t150 * MDP(9)) * pkin(1) + (t171 - 0.2e1 * t98 * MDP(18) + 0.2e1 * t163 + (MDP(22) * t148 - t172) * t185) * t106; t147 * MDP(6) + t150 * MDP(7) + (t111 * t87 + t112 * t89 + t113 * t92) * MDP(28) + (-t106 * t119 - t107 * t118) * MDP(18) + (-t118 * t96 + t119 * t98) * MDP(19) + t152 + (-MDP(10) * t150 - MDP(9) * t147) * pkin(7) + (t162 * MDP(25) + (-t107 * t112 - t87) * MDP(27)) * t145 + (t162 * MDP(26) - t111 * t170 - t176) * t148; MDP(8) + (t118 ^ 2 + t119 ^ 2) * MDP(19) + (-t111 * t145 + t108) * t183 + (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) * MDP(28) - t115 * t154 + 0.2e1 * t155 + t165; ((-t106 * t143 - t107 * t144) * MDP(18) + (t143 * t98 - t144 * t96) * MDP(19)) * pkin(3) + (t161 * MDP(25) + (-t107 * t122 - t87) * MDP(27)) * t145 + (t121 * t87 + t122 * t89 + t126 * t92) * MDP(28) + (t161 * MDP(26) - t121 * t170 - t176) * t148 + t152; (t108 + t117) * MDP(27) + (t111 * t121 + t112 * t122 + t113 * t126) * MDP(28) - t175 * t167 + (t118 * t144 + t119 * t143) * MDP(19) * pkin(3) + t155 + (t175 * MDP(26) + (-t111 - t121) * MDP(27)) * t145 + t165; (-t121 * t145 + t117) * t183 + (t121 ^ 2 + t122 ^ 2 + t126 ^ 2) * MDP(28) + (t143 ^ 2 + t144 ^ 2) * MDP(19) * pkin(3) ^ 2 - t134 * t154 + t165; t114 * MDP(19) + t164 * MDP(28) + t157 * t106 - t173 * t170; (t111 * t148 + t112 * t145) * MDP(28); (t121 * t148 + t122 * t145) * MDP(28); t173 * MDP(28) + MDP(19); t87 * t179 + t171 + (-t172 + (-MDP(27) * pkin(5) + MDP(22)) * t148) * t107 + t163; t156 * t116 + (MDP(28) * t111 - t169) * pkin(5) + t174; t156 * t133 + (MDP(28) * t121 - t169) * pkin(5) + t174; -t168 + (MDP(25) + t179) * t148; MDP(28) * pkin(5) ^ 2 + MDP(24); t92 * MDP(28); t113 * MDP(28); t126 * MDP(28); 0; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
