% Calculate joint inertia matrix for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP6_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:29
% EndTime: 2019-03-09 04:48:31
% DurationCPUTime: 0.64s
% Computational Cost: add. (845->180), mult. (1496->255), div. (0->0), fcn. (1468->6), ass. (0->63)
t131 = sin(qJ(4));
t133 = cos(qJ(4));
t171 = MDP(19) * t131 + MDP(20) * t133;
t170 = MDP(21) + MDP(24);
t152 = MDP(22) + MDP(26);
t169 = 2 * MDP(21);
t168 = 2 * MDP(23);
t167 = 2 * MDP(24);
t166 = (pkin(1) * MDP(6));
t165 = -qJ(5) - pkin(8);
t129 = sin(pkin(9));
t130 = cos(pkin(9));
t132 = sin(qJ(3));
t134 = cos(qJ(3));
t115 = pkin(3) * t132 - pkin(8) * t134 + qJ(2);
t110 = t133 * t115;
t161 = t133 * t134;
t135 = -pkin(1) - pkin(7);
t162 = t131 * t135;
t94 = -qJ(5) * t161 + t110 + (pkin(4) - t162) * t132;
t160 = t133 * t135;
t150 = t132 * t160;
t97 = t150 + (-qJ(5) * t134 + t115) * t131;
t85 = t129 * t94 + t130 * t97;
t164 = t131 * t133;
t163 = t131 * t134;
t111 = t129 * t131 - t130 * t133;
t112 = t129 * t133 + t130 * t131;
t123 = -pkin(4) * t133 - pkin(3);
t93 = pkin(5) * t111 - qJ(6) * t112 + t123;
t159 = t93 * MDP(26);
t155 = t111 * MDP(23);
t154 = t112 * MDP(25);
t119 = pkin(4) * t129 + qJ(6);
t153 = t119 * MDP(25);
t116 = t165 * t133;
t148 = t165 * t131;
t101 = -t130 * t116 + t129 * t148;
t99 = -t116 * t129 - t130 * t148;
t151 = t101 ^ 2 + t99 ^ 2;
t82 = t132 * qJ(6) + t85;
t149 = MDP(15) * t164;
t107 = t112 * t134;
t109 = -t129 * t163 + t130 * t161;
t146 = -t101 * t107 + t109 * t99;
t84 = -t129 * t97 + t130 * t94;
t113 = pkin(4) * t163 - t134 * t135;
t141 = MDP(16) * t133 - MDP(17) * t131;
t140 = t133 * MDP(19) - t131 * MDP(20);
t138 = MDP(22) * t123 - t154 + t155 + t159;
t137 = t131 * MDP(16) + t133 * MDP(17) - t99 * MDP(23) + t101 * MDP(25) - pkin(8) * t171;
t128 = t134 ^ 2;
t127 = t133 ^ 2;
t126 = t132 ^ 2;
t125 = t131 ^ 2;
t121 = pkin(4) * t130 + pkin(5);
t108 = t111 * t132;
t105 = t112 * t132;
t103 = t115 * t131 + t150;
t102 = -t132 * t162 + t110;
t86 = pkin(5) * t107 - qJ(6) * t109 + t113;
t83 = -pkin(5) * t132 - t84;
t1 = [MDP(1) + t126 * MDP(18) + (t113 ^ 2 + t84 ^ 2 + t85 ^ 2) * MDP(22) + (t82 ^ 2 + t83 ^ 2 + t86 ^ 2) * MDP(26) + ((-2 * MDP(4) + t166) * pkin(1)) + (0.2e1 * t134 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (t127 * MDP(14) + MDP(7) - 0.2e1 * t149) * t128 + 0.2e1 * (qJ(2) * MDP(12) + (-MDP(8) + t141) * t134) * t132 + 0.2e1 * (t102 * t132 - t128 * t162) * MDP(19) + 0.2e1 * (-t103 * t132 - t128 * t160) * MDP(20) + (-t107 * t85 - t109 * t84) * t169 + (t107 * t86 - t132 * t83) * t168 + (-t107 * t82 + t109 * t83) * t167 + 0.2e1 * (-t109 * t86 + t132 * t82) * MDP(25); MDP(4) - t166 + (-t105 * t84 - t108 * t85 - t113 * t134) * MDP(22) + (-t105 * t132 - t107 * t134) * MDP(23) + (-t108 * t132 + t109 * t134) * MDP(25) + (t105 * t83 - t108 * t82 - t134 * t86) * MDP(26) + t171 * (-t126 - t128) + t170 * (t105 * t109 + t108 * t107); MDP(6) + t152 * (t105 ^ 2 + t108 ^ 2 + t128); (-t111 * t85 - t84 * t112 + t146) * MDP(21) + (t101 * t85 + t113 * t123 - t84 * t99) * MDP(22) + (t107 * t93 + t111 * t86) * MDP(23) + (-t111 * t82 + t112 * t83 + t146) * MDP(24) + (-t109 * t93 - t112 * t86) * MDP(25) + (t101 * t82 + t83 * t99 + t86 * t93) * MDP(26) + (-MDP(13) * t135 - MDP(10) + t137) * t132 + (MDP(9) + t135 * MDP(12) + MDP(14) * t164 + (-t125 + t127) * MDP(15) + (-pkin(3) * t131 + t160) * MDP(19) + (-pkin(3) * t133 - t162) * MDP(20)) * t134; -t132 * MDP(13) + (MDP(12) - t138 + t140) * t134 + t152 * (-t108 * t101 + t105 * t99) + t170 * (t105 * t112 + t108 * t111); MDP(11) + t125 * MDP(14) + 0.2e1 * t149 + (t123 ^ 2 + t151) * MDP(22) + t151 * MDP(26) + (-0.2e1 * t154 + 0.2e1 * t155 + t159) * t93 + 0.2e1 * t140 * pkin(3) + (t169 + t167) * (-t101 * t111 + t112 * t99); t102 * MDP(19) - t103 * MDP(20) + t84 * MDP(23) + (-t107 * t119 - t109 * t121) * MDP(24) + t82 * MDP(25) + (t119 * t82 - t121 * t83) * MDP(26) + t141 * t134 + (MDP(18) + (pkin(5) + t121) * MDP(23) + t153) * t132 + ((-t107 * t129 - t109 * t130) * MDP(21) + (t129 * t85 + t130 * t84) * MDP(22)) * pkin(4); -t105 * MDP(23) - t108 * MDP(25) + (-t105 * t121 - t108 * t119) * MDP(26) - t171 * t132 + (-t105 * t130 - t108 * t129) * MDP(22) * pkin(4); (-t111 * t119 - t112 * t121) * MDP(24) + (t101 * t119 - t121 * t99) * MDP(26) + ((-t111 * t129 - t112 * t130) * MDP(21) + (t101 * t129 - t130 * t99) * MDP(22)) * pkin(4) + t137; MDP(18) + (t119 ^ 2 + t121 ^ 2) * MDP(26) + (t129 ^ 2 + t130 ^ 2) * MDP(22) * pkin(4) ^ 2 + t121 * t168 + 0.2e1 * t153; MDP(22) * t113 + t107 * MDP(23) - t109 * MDP(25) + t86 * MDP(26); -t152 * t134; t138; 0; t152; -t132 * MDP(23) + t109 * MDP(24) + t83 * MDP(26); t105 * MDP(26); MDP(24) * t112 + MDP(26) * t99; -MDP(26) * t121 - MDP(23); 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
