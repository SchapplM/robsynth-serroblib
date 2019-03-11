% Calculate Gravitation load on the joints for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:45:57
% EndTime: 2019-03-09 16:45:58
% DurationCPUTime: 0.47s
% Computational Cost: add. (370->97), mult. (479->127), div. (0->0), fcn. (438->8), ass. (0->51)
t149 = MDP(17) - MDP(20);
t147 = MDP(27) + MDP(29);
t146 = MDP(28) - MDP(31);
t148 = MDP(16) - MDP(19) + MDP(30);
t107 = qJ(2) + qJ(3);
t104 = sin(t107);
t105 = cos(t107);
t136 = t105 * pkin(3) + t104 * qJ(4);
t110 = sin(qJ(1));
t113 = cos(qJ(1));
t92 = g(1) * t113 + g(2) * t110;
t145 = pkin(3) + pkin(9);
t109 = sin(qJ(2));
t144 = pkin(2) * t109;
t143 = pkin(3) * t104;
t139 = g(3) * t105;
t100 = t105 * pkin(9);
t108 = sin(qJ(5));
t133 = t105 * t108;
t124 = pkin(5) * t133;
t93 = t110 * t105 * qJ(4);
t138 = t110 * t124 + t93;
t131 = t105 * t113;
t95 = qJ(4) * t131;
t137 = t113 * t124 + t95;
t111 = cos(qJ(5));
t135 = qJ(6) * t111;
t134 = t104 * t113;
t132 = t105 * t111;
t130 = t108 * t110;
t129 = t108 * t113;
t128 = t110 * t111;
t127 = t111 * t113;
t114 = -pkin(8) - pkin(7);
t126 = t113 * t114;
t112 = cos(qJ(2));
t106 = t112 * pkin(2);
t103 = t106 + pkin(1);
t125 = pkin(3) * t131 + qJ(4) * t134 + t113 * t103;
t123 = qJ(6) * t132;
t122 = t104 * t108 * pkin(5) + t100 + t136;
t120 = -t143 - t144;
t118 = -t103 - t136;
t85 = -t104 * t127 + t130;
t87 = t104 * t128 + t129;
t68 = g(1) * t85 - g(2) * t87 + g(3) * t132;
t81 = t92 * t104 - t139;
t116 = t148 * t81 + (-t108 * t147 - t111 * t146 + t149) * (g(3) * t104 + t92 * t105);
t88 = -t104 * t130 + t127;
t86 = t104 * t129 + t128;
t1 = [(g(1) * t126 - g(2) * t125 + (-g(1) * t118 + g(2) * t114) * t110) * MDP(21) + (-g(1) * (t113 * pkin(4) + t88 * pkin(5) + t87 * qJ(6) - t126) - g(2) * (t86 * pkin(5) + pkin(9) * t131 + t85 * qJ(6) + t125) + (-g(1) * (t118 - t100) - g(2) * (pkin(4) - t114)) * t110) * MDP(32) + (MDP(3) - MDP(18)) * t92 + t147 * (-g(1) * t88 - g(2) * t86) + t146 * (g(1) * t87 + g(2) * t85) + (-t109 * MDP(10) + t112 * MDP(9) - t149 * t104 + t148 * t105 + MDP(2)) * (g(1) * t110 - g(2) * t113); (-g(3) * t112 + t92 * t109) * MDP(9) + (g(3) * t109 + t92 * t112) * MDP(10) + (-g(1) * (t120 * t113 + t95) - g(2) * (t120 * t110 + t93) - g(3) * (t106 + t136)) * MDP(21) + (-g(1) * t137 - g(2) * t138 - g(3) * (-t104 * t135 + t106 + t122) + t92 * (t145 * t104 + t123 + t144)) * MDP(32) + t116; (-g(1) * (-pkin(3) * t134 + t95) - g(2) * (-t110 * t143 + t93) - g(3) * t136) * MDP(21) + (-g(1) * (-t113 * t123 + t137) - g(2) * (-t110 * t123 + t138) - g(3) * t122 + (g(3) * t135 + t92 * t145) * t104) * MDP(32) + t116; (-MDP(21) - MDP(32)) * t81; (-g(1) * (-pkin(5) * t85 + qJ(6) * t86) - g(2) * (pkin(5) * t87 - qJ(6) * t88) - (-pkin(5) * t111 - qJ(6) * t108) * t139) * MDP(32) + t147 * t68 - t146 * (-g(1) * t86 + g(2) * t88 + g(3) * t133); -t68 * MDP(32);];
taug  = t1;
