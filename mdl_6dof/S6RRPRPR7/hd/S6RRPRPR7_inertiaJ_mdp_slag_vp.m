% Calculate joint inertia matrix for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR7_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:47:45
% EndTime: 2019-03-09 10:47:46
% DurationCPUTime: 0.46s
% Computational Cost: add. (766->145), mult. (1317->219), div. (0->0), fcn. (1429->8), ass. (0->67)
t131 = sin(qJ(6));
t134 = cos(qJ(6));
t144 = -MDP(29) * t131 - MDP(30) * t134;
t132 = sin(qJ(4));
t135 = cos(qJ(4));
t136 = cos(qJ(2));
t170 = pkin(7) - pkin(8);
t154 = t170 * t136;
t133 = sin(qJ(2));
t155 = t170 * t133;
t100 = t132 * t154 - t135 * t155;
t101 = t132 * t155 + t135 * t154;
t109 = -t132 * t136 + t133 * t135;
t125 = t131 ^ 2;
t127 = t134 ^ 2;
t164 = t133 * t132 + t136 * t135;
t129 = sin(pkin(10));
t130 = cos(pkin(10));
t94 = t130 * t109 - t129 * t164;
t172 = -t164 * MDP(18) + t109 * MDP(17) - t100 * MDP(20) - t101 * MDP(21) - (t125 - t127) * t94 * MDP(25);
t113 = -t136 * pkin(2) - t133 * qJ(3) - pkin(1);
t105 = t136 * pkin(3) - t113;
t171 = 0.2e1 * t105;
t140 = -t109 * qJ(5) - t100;
t91 = -t164 * qJ(5) + t101;
t86 = t129 * t91 - t130 * t140;
t169 = t86 * t94;
t168 = MDP(23) * pkin(4);
t167 = t134 * t94;
t117 = -pkin(4) * t130 - pkin(5);
t137 = -pkin(2) - pkin(3);
t111 = qJ(3) * t132 - t135 * t137;
t110 = -pkin(4) - t111;
t112 = qJ(3) * t135 + t132 * t137;
t97 = t110 * t130 - t112 * t129;
t95 = pkin(5) - t97;
t166 = t117 - t95;
t93 = t109 * t129 + t130 * t164;
t165 = t93 * MDP(28);
t98 = t129 * t110 + t130 * t112;
t126 = t133 ^ 2;
t163 = t136 ^ 2 + t126;
t161 = MDP(29) * t134;
t159 = t111 * MDP(20);
t158 = t112 * MDP(21);
t157 = t134 * MDP(25);
t156 = t125 * MDP(24) + MDP(19);
t153 = t131 * t157;
t152 = 0.2e1 * t153 + t156;
t151 = -pkin(2) * MDP(14) - MDP(11);
t96 = -pkin(9) + t98;
t150 = -t93 * t96 + t94 * t95;
t116 = pkin(4) * t129 + pkin(9);
t148 = -t116 * t93 + t117 * t94;
t147 = -t93 * MDP(27) + t86 * MDP(29);
t146 = t135 * MDP(20) - t132 * MDP(21);
t145 = -MDP(30) * t131 + t161;
t143 = (MDP(26) * t134 - MDP(27) * t131) * t94;
t99 = t164 * pkin(4) + t105;
t142 = -MDP(24) * t167 - t93 * MDP(26) - t86 * MDP(30);
t108 = t129 * t135 + t130 * t132;
t106 = t129 * t132 - t130 * t135;
t88 = t129 * t140 + t130 * t91;
t85 = t93 * pkin(5) - t94 * pkin(9) + t99;
t84 = t131 * t85 + t134 * t88;
t83 = -t131 * t88 + t134 * t85;
t1 = [MDP(1) + t126 * MDP(4) + (t163 * pkin(7) ^ 2 + t113 ^ 2) * MDP(14) + t164 * MDP(20) * t171 + (t86 ^ 2 + t88 ^ 2 + t99 ^ 2) * MDP(23) + (MDP(24) * t127 - 0.2e1 * t153) * t94 ^ 2 + t165 * t93 + (MDP(15) * t109 - 0.2e1 * t164 * MDP(16) + MDP(21) * t171) * t109 + 0.2e1 * (-t88 * t93 + t169) * MDP(22) + 0.2e1 * (t131 * t169 + t83 * t93) * MDP(29) + 0.2e1 * (t86 * t167 - t84 * t93) * MDP(30) + 0.2e1 * t163 * MDP(12) * pkin(7) + 0.2e1 * (-MDP(11) * t113 + MDP(9) * pkin(1)) * t136 + 0.2e1 * t143 * t93 + 0.2e1 * (-MDP(10) * pkin(1) - MDP(13) * t113 + MDP(5) * t136) * t133; t133 * MDP(6) + t136 * MDP(7) + (-pkin(2) * t133 + qJ(3) * t136) * MDP(12) + (-t93 * t98 - t94 * t97) * MDP(22) + (-t86 * t97 + t88 * t98) * MDP(23) + (t150 * MDP(30) + t147) * t134 + (t150 * MDP(29) + t142) * t131 + ((MDP(14) * qJ(3) - MDP(10) + MDP(13)) * t136 + (-MDP(9) + t151) * t133) * pkin(7) - t172; MDP(8) + 0.2e1 * pkin(2) * MDP(11) + 0.2e1 * qJ(3) * MDP(13) + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(14) + (t97 ^ 2 + t98 ^ 2) * MDP(23) + 0.2e1 * t145 * t95 + 0.2e1 * t159 + 0.2e1 * t158 + t152; (t106 * t86 + t108 * t88) * MDP(23) + (MDP(14) * pkin(7) + MDP(12)) * t133 + (MDP(22) - t144) * (t106 * t94 - t108 * t93); t98 * t108 * MDP(23) + (-MDP(23) * t97 + t145) * t106 - t146 + t151; MDP(14) + (t106 ^ 2 + t108 ^ 2) * MDP(23); (t148 * MDP(30) - t147) * t134 + (t148 * MDP(29) - t142) * t131 + ((-t129 * t93 - t130 * t94) * MDP(22) + (t129 * t88 - t130 * t86) * MDP(23)) * pkin(4) + t172; -t159 - t158 + t166 * t161 + (t129 * t98 + t130 * t97) * t168 + (-t166 * MDP(30) - 0.2e1 * t157) * t131 - t156; -t145 * t106 + (-t106 * t130 + t108 * t129) * t168 + t146; (t129 ^ 2 + t130 ^ 2) * MDP(23) * pkin(4) ^ 2 - 0.2e1 * t145 * t117 + t152; MDP(23) * t99 + t145 * t93; 0; 0; 0; MDP(23); MDP(29) * t83 - MDP(30) * t84 + t143 + t165; (-MDP(30) * t96 - MDP(27)) * t134 + (-MDP(29) * t96 - MDP(26)) * t131; t144 * t108; MDP(26) * t131 + MDP(27) * t134 + t144 * t116; t145; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
