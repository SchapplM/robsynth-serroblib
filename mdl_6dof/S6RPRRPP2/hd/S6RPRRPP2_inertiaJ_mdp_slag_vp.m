% Calculate joint inertia matrix for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP2_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:22
% EndTime: 2019-03-09 04:33:24
% DurationCPUTime: 0.63s
% Computational Cost: add. (558->173), mult. (935->217), div. (0->0), fcn. (762->6), ass. (0->63)
t108 = sin(qJ(4));
t112 = pkin(4) + pkin(5);
t147 = t108 * t112;
t146 = pkin(4) * t108;
t145 = pkin(8) - qJ(6);
t110 = cos(qJ(4));
t111 = cos(qJ(3));
t141 = t108 * t111;
t109 = sin(qJ(3));
t107 = cos(pkin(9));
t99 = -pkin(1) * t107 - pkin(2);
t86 = -pkin(3) * t111 - pkin(8) * t109 + t99;
t106 = sin(pkin(9));
t98 = pkin(1) * t106 + pkin(7);
t80 = t110 * t86 - t98 * t141;
t81 = t110 * t111 * t98 + t108 * t86;
t140 = t109 * t110;
t96 = qJ(5) * t140;
t79 = t96 + (-t98 - t147) * t109;
t144 = MDP(26) * t79;
t143 = qJ(5) * t110;
t142 = t108 * t109;
t139 = t111 * qJ(5);
t127 = qJ(5) * t108 + pkin(3);
t87 = t112 * t110 + t127;
t138 = t87 * MDP(26);
t91 = -pkin(4) * t110 - t127;
t137 = t91 * MDP(22);
t102 = t108 ^ 2;
t104 = t110 ^ 2;
t136 = t102 + t104;
t135 = MDP(13) * t110;
t134 = MDP(18) - MDP(21);
t133 = MDP(19) + MDP(23);
t132 = MDP(20) - MDP(25);
t131 = MDP(21) + MDP(24);
t130 = MDP(22) + MDP(26);
t129 = -0.2e1 * t139 + t81;
t101 = t111 * pkin(4);
t78 = t101 - t80;
t128 = t136 * pkin(8);
t126 = -pkin(4) * MDP(22) - MDP(19);
t125 = MDP(22) * pkin(8) + MDP(20);
t124 = t109 * t132;
t77 = -t139 + t81;
t123 = t80 * MDP(17) - t81 * MDP(18);
t122 = -MDP(26) * t112 - MDP(23) + t126;
t95 = qJ(6) * t142;
t76 = t77 + t95;
t82 = -t96 + (t98 + t146) * t109;
t121 = -t82 * MDP(19) + t79 * MDP(23) - t76 * MDP(25);
t75 = pkin(5) * t111 - qJ(6) * t140 + t78;
t120 = -t82 * MDP(21) + t79 * MDP(24) - t75 * MDP(25);
t92 = t145 * t108;
t93 = t145 * t110;
t119 = -t108 * MDP(14) + t92 * MDP(23) - t93 * MDP(24);
t118 = pkin(3) * MDP(17) - t91 * MDP(19) + t87 * MDP(23) - t93 * MDP(25);
t117 = -pkin(3) * MDP(18) - t91 * MDP(21) + t87 * MDP(24) - t92 * MDP(25);
t114 = qJ(5) ^ 2;
t105 = t111 ^ 2;
t103 = t109 ^ 2;
t94 = MDP(24) * t140;
t1 = [MDP(1) - 0.2e1 * t99 * t111 * MDP(10) + t105 * MDP(16) + (t77 ^ 2 + t78 ^ 2 + t82 ^ 2) * MDP(22) + (t75 ^ 2 + t76 ^ 2 + t79 ^ 2) * MDP(26) + (t106 ^ 2 + t107 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * (t78 * MDP(19) - t77 * MDP(21) + t75 * MDP(23) - t76 * MDP(24) - t123) * t111 + (t104 * MDP(12) - 0.2e1 * t108 * t135 + MDP(5) + 0.2e1 * (t108 * MDP(17) + t110 * MDP(18)) * t98) * t103 + 0.2e1 * (t99 * MDP(11) + (-t110 * MDP(14) + t108 * MDP(15) + MDP(6)) * t111 + (t78 * MDP(20) + t120) * t110 + (-t77 * MDP(20) - t121) * t108) * t109; (-MDP(22) * t82 + t144) * t111 + ((t108 * t78 + t110 * t77) * MDP(22) + (t108 * t75 + t110 * t76) * MDP(26)) * t109; MDP(4) + t130 * (t136 * t103 + t105); t82 * t137 + (t75 * t92 + t76 * t93 + t79 * t87) * MDP(26) + (MDP(17) + MDP(19)) * pkin(8) * t141 + (t125 * t77 + t121) * t110 + (t125 * t78 + t120) * t108 + (-MDP(11) * t98 + MDP(8) + (t134 * pkin(8) - MDP(15)) * t110 + t119) * t111 + (MDP(7) - t98 * MDP(10) + (-t102 + t104) * MDP(13) + (-t98 * MDP(17) + t117) * t110 + (t110 * MDP(12) + t98 * MDP(18) - t118) * t108) * t109; t136 * t124 + (-MDP(11) + (t108 * t92 + t110 * t93) * MDP(26) + MDP(22) * t128) * t109 + (-t137 + t138 + MDP(10) + (MDP(17) + t133) * t110 + (-MDP(18) + t131) * t108) * t111; MDP(9) + t102 * MDP(12) + (t136 * pkin(8) ^ 2 + t91 ^ 2) * MDP(22) + (t87 ^ 2 + t92 ^ 2 + t93 ^ 2) * MDP(26) + 0.2e1 * MDP(20) * t128 + 0.2e1 * t118 * t110 + 0.2e1 * (t117 + t135) * t108; (-0.2e1 * t101 + t80) * MDP(19) + t129 * MDP(21) + (-pkin(4) * t78 + qJ(5) * t77) * MDP(22) - t78 * MDP(23) + (t95 + t129) * MDP(24) + (qJ(5) * t76 - t112 * t75) * MDP(26) + (-MDP(16) + (-pkin(5) - t112) * MDP(23)) * t111 + ((-t132 * qJ(5) - MDP(15)) * t108 + (-MDP(20) * pkin(4) + MDP(23) * qJ(6) + MDP(25) * t112 + MDP(14)) * t110) * t109 + t123; t94 + t130 * t96 + (-t134 * t110 + (-MDP(17) + t122) * t108) * t109; t110 * MDP(15) + (t143 - t146) * MDP(20) + (-t143 + t147) * MDP(25) + (qJ(5) * t93 - t112 * t92) * MDP(26) + ((MDP(22) * qJ(5) - t134) * t110 + (-MDP(17) + t126) * t108) * pkin(8) - t119; MDP(16) + 0.2e1 * pkin(4) * MDP(19) + (pkin(4) ^ 2 + t114) * MDP(22) + 0.2e1 * t112 * MDP(23) + (t112 ^ 2 + t114) * MDP(26) + 0.2e1 * t131 * qJ(5); MDP(22) * t78 + MDP(26) * t75 + t110 * t124 + t133 * t111; t130 * t142; t92 * MDP(26) + (-MDP(25) + t125) * t108; t122; t130; -MDP(23) * t142 + t144 + t94; t111 * MDP(26); t110 * MDP(23) + t108 * MDP(24) + t138; 0; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
