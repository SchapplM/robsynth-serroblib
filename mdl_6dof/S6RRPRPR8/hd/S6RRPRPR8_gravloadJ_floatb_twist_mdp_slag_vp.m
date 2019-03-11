% Calculate Gravitation load on the joints for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:53:21
% EndTime: 2019-03-09 10:53:23
% DurationCPUTime: 0.75s
% Computational Cost: add. (368->95), mult. (522->137), div. (0->0), fcn. (543->10), ass. (0->41)
t102 = sin(qJ(6));
t105 = cos(qJ(6));
t98 = pkin(10) + qJ(4);
t92 = sin(t98);
t93 = cos(t98);
t110 = t102 * t92 + t105 * t93;
t111 = t102 * t93 - t105 * t92;
t104 = sin(qJ(1));
t106 = cos(qJ(2));
t107 = cos(qJ(1));
t120 = t106 * t107;
t84 = -t104 * t93 + t92 * t120;
t85 = t104 * t92 + t93 * t120;
t112 = t102 * t85 - t105 * t84;
t121 = t104 * t106;
t82 = t107 * t93 + t92 * t121;
t83 = -t107 * t92 + t93 * t121;
t113 = t102 * t82 + t105 * t83;
t103 = sin(qJ(2));
t127 = g(3) * t103;
t138 = -t102 * t83 + t82 * t105;
t75 = t102 * t84 + t105 * t85;
t144 = -(-g(1) * t112 + g(2) * t138 - t111 * t127) * MDP(31) + (g(1) * t75 + g(2) * t113 + t110 * t127) * MDP(32);
t117 = g(1) * t107 + g(2) * t104;
t142 = g(3) * t106 - t117 * t103;
t136 = MDP(20) + MDP(22);
t135 = MDP(21) - MDP(24);
t139 = MDP(10) - MDP(13) - MDP(23);
t130 = g(1) * t104;
t125 = t107 * pkin(1) + t104 * pkin(7);
t99 = sin(pkin(10));
t123 = t107 * t99;
t101 = -pkin(8) - qJ(3);
t122 = t101 * t103;
t115 = pkin(2) * t106 + qJ(3) * t103;
t100 = cos(pkin(10));
t91 = pkin(3) * t100 + pkin(2);
t109 = pkin(4) * t93 + qJ(5) * t92 + t91;
t73 = g(1) * t84 + g(2) * t82 + t92 * t127;
t95 = t107 * pkin(7);
t1 = [t117 * MDP(3) + (-g(1) * (-t100 * t121 + t123) - g(2) * (t100 * t120 + t104 * t99)) * MDP(11) + (-g(1) * (t100 * t107 + t99 * t121) - g(2) * (t104 * t100 - t99 * t120)) * MDP(12) + (-g(1) * t95 - g(2) * (t115 * t107 + t125) - (-pkin(1) - t115) * t130) * MDP(14) + (-g(1) * (pkin(3) * t123 - t83 * pkin(4) - t82 * qJ(5) + t95) - g(2) * (t85 * pkin(4) + t84 * qJ(5) - t107 * t122 + t91 * t120 + t125) + (-g(1) * (-t106 * t91 - pkin(1) + t122) - g(2) * pkin(3) * t99) * t104) * MDP(25) + (g(1) * t113 - g(2) * t75) * MDP(31) + (g(1) * t138 + g(2) * t112) * MDP(32) + t136 * (g(1) * t83 - g(2) * t85) - t135 * (g(1) * t82 - g(2) * t84) + (t106 * MDP(9) - t139 * t103 + MDP(2)) * (-g(2) * t107 + t130); (-g(3) * t115 + t117 * (pkin(2) * t103 - qJ(3) * t106)) * MDP(14) + ((-g(3) * t109 + t117 * t101) * t106 + (g(3) * t101 + t109 * t117) * t103) * MDP(25) + t139 * (t117 * t106 + t127) + (-t100 * MDP(11) + t99 * MDP(12) - t110 * MDP(31) + t111 * MDP(32) + t135 * t92 - t136 * t93 - MDP(9)) * t142; -(-MDP(14) - MDP(25)) * t142; (-g(1) * (-pkin(4) * t84 + qJ(5) * t85) - g(2) * (-pkin(4) * t82 + qJ(5) * t83) - (-pkin(4) * t92 + qJ(5) * t93) * t127) * MDP(25) + t136 * t73 + t135 * (g(1) * t85 + g(2) * t83 + t93 * t127) - t144; -t73 * MDP(25); t144;];
taug  = t1;
