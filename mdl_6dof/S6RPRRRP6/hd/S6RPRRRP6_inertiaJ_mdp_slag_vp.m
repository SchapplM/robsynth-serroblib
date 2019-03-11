% Calculate joint inertia matrix for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRRP6_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:16:46
% EndTime: 2019-03-09 06:16:48
% DurationCPUTime: 0.59s
% Computational Cost: add. (1064->156), mult. (2012->228), div. (0->0), fcn. (2271->8), ass. (0->77)
t142 = sin(qJ(4));
t144 = cos(qJ(4));
t149 = -t142 * MDP(20) - t144 * MDP(21);
t180 = pkin(8) + pkin(9);
t128 = t180 * t142;
t129 = t180 * t144;
t141 = sin(qJ(5));
t178 = cos(qJ(5));
t106 = -t178 * t128 - t129 * t141;
t107 = -t141 * t128 + t129 * t178;
t124 = t141 * t144 + t142 * t178;
t185 = -t141 * t142 + t178 * t144;
t151 = t124 * MDP(24) + MDP(25) * t185 + t106 * MDP(27) - t107 * MDP(28);
t186 = t142 * MDP(17) + t144 * MDP(18) + pkin(8) * t149 + t151;
t184 = t144 * MDP(17) - t142 * MDP(18);
t140 = cos(pkin(10));
t131 = -pkin(2) * t140 - pkin(1);
t183 = 0.2e1 * t131;
t182 = -2 * MDP(23);
t181 = 2 * MDP(29);
t179 = cos(qJ(3));
t177 = pkin(1) * MDP(7);
t139 = sin(pkin(10));
t143 = sin(qJ(3));
t119 = t139 * t143 - t140 * t179;
t176 = pkin(4) * t119;
t175 = pkin(4) * t141;
t174 = pkin(7) + qJ(2);
t173 = MDP(30) * pkin(5);
t120 = t139 * t179 + t143 * t140;
t99 = t185 * t120;
t172 = t185 * t99;
t126 = t174 * t139;
t127 = t174 * t140;
t103 = -t143 * t126 + t127 * t179;
t171 = t103 * t144;
t170 = t120 * t142;
t169 = t120 * t144;
t168 = t139 * MDP(5);
t167 = t140 * MDP(4);
t165 = t142 * t144;
t98 = t124 * t120;
t94 = t98 * MDP(25);
t164 = t99 * MDP(22);
t95 = t99 * MDP(24);
t163 = MDP(27) * t185 - t124 * MDP(28);
t161 = t120 * MDP(14);
t158 = MDP(19) + MDP(26);
t157 = t119 * MDP(26) - t94 + t95;
t156 = t178 * pkin(4);
t101 = pkin(3) * t119 - pkin(8) * t120 + t131;
t88 = t171 + (-pkin(9) * t120 + t101) * t142;
t155 = t178 * t88;
t134 = -pkin(4) * t144 - pkin(3);
t153 = MDP(16) * t165;
t152 = MDP(27) * t178;
t90 = t144 * t101 - t103 * t142;
t87 = -pkin(9) * t169 + t176 + t90;
t84 = -t141 * t88 + t178 * t87;
t102 = t126 * t179 + t143 * t127;
t150 = t144 * MDP(20) - t142 * MDP(21);
t93 = pkin(4) * t170 + t102;
t85 = t141 * t87 + t155;
t147 = -MDP(13) - t150;
t138 = t144 ^ 2;
t137 = t142 ^ 2;
t133 = t156 + pkin(5);
t121 = t124 ^ 2;
t108 = -pkin(5) * t185 + t134;
t97 = qJ(6) * t185 + t107;
t96 = -qJ(6) * t124 + t106;
t92 = t124 * t98;
t91 = t101 * t142 + t171;
t89 = t98 * pkin(5) + t93;
t83 = -t98 * qJ(6) + t85;
t82 = pkin(5) * t119 - qJ(6) * t99 + t84;
t1 = [(t82 ^ 2 + t83 ^ 2 + t89 ^ 2) * MDP(30) + t161 * t183 + MDP(1) + (t182 * t98 + t164) * t99 + t158 * t119 ^ 2 + (0.2e1 * t167 - 0.2e1 * t168 + t177) * pkin(1) + (t138 * MDP(15) + MDP(8) - 0.2e1 * t153) * t120 ^ 2 + (MDP(13) * t183 + 0.2e1 * t95 - 0.2e1 * t94 + 0.2e1 * (-MDP(9) + t184) * t120) * t119 + (-t82 * t99 - t83 * t98) * t181 + 0.2e1 * (t119 * t84 + t93 * t98) * MDP(27) + 0.2e1 * (-t119 * t85 + t93 * t99) * MDP(28) + 0.2e1 * (t102 * t170 + t119 * t90) * MDP(20) + 0.2e1 * (t102 * t169 - t119 * t91) * MDP(21) + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t139 ^ 2 + t140 ^ 2) * qJ(2); -t167 + t168 - t177 + t161 + (-t92 - t172) * MDP(29) + (t124 * t83 + t185 * t82) * MDP(30) + (-t147 + t163) * t119; MDP(7) + (t185 ^ 2 + t121) * MDP(30); -t103 * MDP(14) + t124 * t164 + (-t92 + t172) * MDP(23) + (t134 * t98 - t185 * t93) * MDP(27) + (t124 * t93 + t134 * t99) * MDP(28) + (-t124 * t82 + t185 * t83 - t96 * t99 - t97 * t98) * MDP(29) + (t108 * t89 + t82 * t96 + t83 * t97) * MDP(30) + t147 * t102 + (MDP(10) + MDP(15) * t165 + (-t137 + t138) * MDP(16) + t149 * pkin(3)) * t120 + (-MDP(11) + t186) * t119; (t124 * t97 + t185 * t96) * MDP(30); MDP(12) + t137 * MDP(15) + 0.2e1 * t153 + t121 * MDP(22) - t124 * t185 * t182 + (-t124 * t96 + t185 * t97) * t181 + (t108 ^ 2 + t96 ^ 2 + t97 ^ 2) * MDP(30) + 0.2e1 * pkin(3) * t150 - 0.2e1 * t134 * t163; t119 * MDP(19) + t90 * MDP(20) - t91 * MDP(21) + (t119 * t156 + t84) * MDP(27) + (-t155 + (-t87 - t176) * t141) * MDP(28) + (-t133 * t99 - t175 * t98) * MDP(29) + (t133 * t82 + t175 * t83) * MDP(30) + t157 + t184 * t120; (t124 * t175 + t133 * t185) * MDP(30) + t150 + t163; (-t124 * t133 + t175 * t185) * MDP(29) + (t133 * t96 + t175 * t97) * MDP(30) + t186; t133 ^ 2 * MDP(30) + (0.2e1 * t152 + (MDP(30) * t175 - 0.2e1 * MDP(28)) * t141) * pkin(4) + t158; t84 * MDP(27) - t85 * MDP(28) + (-t99 * MDP(29) + t82 * MDP(30)) * pkin(5) + t157; t173 * t185 + t163; (-t124 * MDP(29) + t96 * MDP(30)) * pkin(5) + t151; t133 * t173 + MDP(26) + (-MDP(28) * t141 + t152) * pkin(4); MDP(30) * pkin(5) ^ 2 + MDP(26); t89 * MDP(30); 0; t108 * MDP(30); 0; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
