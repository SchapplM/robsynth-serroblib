% Calculate Gravitation load on the joints for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:32:08
% EndTime: 2019-03-09 09:32:10
% DurationCPUTime: 0.69s
% Computational Cost: add. (246->94), mult. (597->154), div. (0->0), fcn. (682->10), ass. (0->45)
t138 = MDP(10) - MDP(13) - MDP(16);
t137 = -MDP(12) + MDP(9) + MDP(17);
t103 = cos(qJ(6));
t100 = sin(qJ(5));
t104 = cos(qJ(5));
t106 = cos(qJ(1));
t98 = sin(pkin(6));
t122 = t106 * t98;
t101 = sin(qJ(2));
t102 = sin(qJ(1));
t105 = cos(qJ(2));
t121 = cos(pkin(6));
t115 = t106 * t121;
t87 = t101 * t115 + t102 * t105;
t110 = -t87 * t100 + t104 * t122;
t86 = t101 * t102 - t105 * t115;
t99 = sin(qJ(6));
t134 = -t103 * t86 + t110 * t99;
t133 = t103 * t110 + t86 * t99;
t132 = g(3) * t98;
t124 = t105 * t98;
t128 = t101 * t98;
t130 = pkin(2) * t124 + qJ(3) * t128;
t129 = t100 * t99;
t127 = t102 * t98;
t125 = t104 * t98;
t123 = t105 * t99;
t120 = t100 * t103;
t119 = t103 * t105;
t118 = -t86 * pkin(2) + qJ(3) * t87;
t116 = t102 * t121;
t88 = t106 * t101 + t105 * t116;
t89 = -t101 * t116 + t105 * t106;
t117 = -t88 * pkin(2) + qJ(3) * t89;
t111 = t106 * pkin(1) + t89 * pkin(2) + pkin(8) * t127 + qJ(3) * t88;
t76 = t100 * t122 + t87 * t104;
t108 = -t102 * pkin(1) - t87 * pkin(2) + pkin(8) * t122 - t86 * qJ(3);
t107 = g(1) * t89 + g(2) * t87 + g(3) * t128;
t69 = -g(1) * t88 - g(2) * t86 + g(3) * t124;
t85 = t100 * t128 + t121 * t104;
t75 = t100 * t89 + t102 * t125;
t74 = -t100 * t127 + t104 * t89;
t68 = t103 * t75 - t88 * t99;
t67 = -t103 * t88 - t75 * t99;
t1 = [(g(1) * t102 - g(2) * t106) * MDP(2) + (-g(1) * t108 - g(2) * t111) * MDP(14) + (-g(1) * (pkin(3) * t122 - t87 * qJ(4) + t108) - g(2) * (pkin(3) * t127 + qJ(4) * t89 + t111)) * MDP(18) + (-g(1) * t110 - g(2) * t75) * MDP(24) + (g(1) * t76 - g(2) * t74) * MDP(25) + (-g(1) * t133 - g(2) * t68) * MDP(31) + (g(1) * t134 - g(2) * t67) * MDP(32) - t138 * (g(1) * t86 - g(2) * t88) + t137 * (g(1) * t87 - g(2) * t89) + (MDP(3) - (MDP(11) + MDP(15)) * t98) * (g(1) * t106 + g(2) * t102); (-g(1) * t117 - g(2) * t118 - g(3) * t130) * MDP(14) + (-g(1) * (-qJ(4) * t88 + t117) - g(2) * (-qJ(4) * t86 + t118) - g(3) * (qJ(4) * t124 + t130)) * MDP(18) + (-g(1) * (-t88 * t120 - t89 * t99) - g(2) * (-t86 * t120 - t87 * t99) - (t100 * t119 - t101 * t99) * t132) * MDP(31) + (-g(1) * (-t103 * t89 + t88 * t129) - g(2) * (-t87 * t103 + t86 * t129) - (-t100 * t123 - t101 * t103) * t132) * MDP(32) + (-t100 * MDP(24) - MDP(25) * t104 - t137) * t69 + t138 * t107; (MDP(14) + MDP(18)) * t69; -t107 * MDP(18); (g(1) * t75 - g(2) * t110 + g(3) * t85) * MDP(25) + (-MDP(31) * t103 + MDP(32) * t99 - MDP(24)) * (g(1) * t74 + g(2) * t76 + g(3) * (-t121 * t100 + t101 * t125)); (-g(1) * t67 - g(2) * t134 - g(3) * (t98 * t119 - t85 * t99)) * MDP(31) + (g(1) * t68 - g(2) * t133 - g(3) * (-t103 * t85 - t98 * t123)) * MDP(32);];
taug  = t1;
