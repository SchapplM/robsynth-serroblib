% Calculate Gravitation load on the joints for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR10V2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR10V2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:51:50
% EndTime: 2019-04-11 14:51:52
% DurationCPUTime: 0.52s
% Computational Cost: add. (412->95), mult. (626->162), div. (0->0), fcn. (710->12), ass. (0->49)
t104 = cos(qJ(6));
t103 = sin(qJ(1));
t105 = cos(qJ(5));
t100 = sin(qJ(5));
t98 = qJ(2) + qJ(3);
t96 = sin(t98);
t131 = t100 * t96;
t101 = sin(qJ(4));
t108 = cos(qJ(1));
t118 = t108 * t101;
t106 = cos(qJ(4));
t122 = t103 * t106;
t97 = cos(t98);
t91 = t122 * t97 - t118;
t79 = t103 * t131 + t105 * t91;
t119 = t106 * t108;
t123 = t103 * t101;
t90 = t123 * t97 + t119;
t99 = sin(qJ(6));
t139 = -t90 * t104 + t79 * t99;
t138 = t79 * t104 + t90 * t99;
t137 = -g(1) * t108 - g(2) * t103;
t136 = g(3) * t96;
t130 = t101 * t99;
t129 = t105 * t96;
t128 = t105 * t99;
t127 = t108 * t96;
t125 = t100 * t106;
t124 = t101 * t104;
t121 = t104 * t105;
t120 = t105 * t106;
t117 = t96 * t130;
t116 = t96 * t124;
t78 = -t100 * t91 + t103 * t129;
t88 = -t100 * t97 + t120 * t96;
t113 = t105 * t97 + t125 * t96;
t85 = t88 * t103;
t86 = t88 * t108;
t89 = t120 * t97 + t131;
t111 = (-t137 * t97 + t136) * MDP(17) + (g(1) * t86 + g(2) * t85 - g(3) * t89) * MDP(30) + (-g(3) * (-t125 * t97 + t129) + t137 * t113) * MDP(31) + (-g(1) * (-t86 * t104 - t108 * t117) - g(2) * (-t103 * t117 - t104 * t85) - g(3) * (t104 * t89 + t130 * t97)) * MDP(37) + (-g(1) * (-t108 * t116 + t86 * t99) - g(2) * (-t103 * t116 + t85 * t99) - g(3) * (t124 * t97 - t89 * t99)) * MDP(38) + (t106 * MDP(23) - t101 * MDP(24) + MDP(16)) * (-g(3) * t97 - t137 * t96);
t107 = cos(qJ(2));
t102 = sin(qJ(2));
t93 = t119 * t97 + t123;
t92 = t118 * t97 - t122;
t82 = t100 * t127 + t93 * t105;
t81 = -t93 * t100 + t105 * t127;
t75 = t104 * t82 + t92 * t99;
t74 = t104 * t92 - t82 * t99;
t1 = [-t137 * MDP(3) + (g(1) * t91 - g(2) * t93) * MDP(23) + (-g(1) * t90 + g(2) * t92) * MDP(24) + (g(1) * t79 - g(2) * t82) * MDP(30) + (g(1) * t78 - g(2) * t81) * MDP(31) + (g(1) * t138 - g(2) * t75) * MDP(37) + (-g(1) * t139 - g(2) * t74) * MDP(38) + (-MDP(10) * t102 + MDP(16) * t97 - MDP(17) * t96 + MDP(9) * t107 + MDP(2)) * (g(1) * t103 - g(2) * t108); (-g(3) * t107 - t102 * t137) * MDP(9) + (g(3) * t102 - t107 * t137) * MDP(10) + t111; t111; (g(1) * t93 + g(2) * t91 + t106 * t136) * MDP(24) + (-g(1) * (-t121 * t92 + t93 * t99) - g(2) * (-t121 * t90 + t91 * t99) - (-t101 * t121 + t106 * t99) * t136) * MDP(37) + (-g(1) * (t104 * t93 + t128 * t92) - g(2) * (t104 * t91 + t128 * t90) - (t101 * t128 + t104 * t106) * t136) * MDP(38) + (MDP(30) * t105 - MDP(31) * t100 + MDP(23)) * (g(1) * t92 + g(2) * t90 + t101 * t136); (g(1) * t82 + g(2) * t79 + g(3) * t88) * MDP(31) + (-MDP(37) * t104 + MDP(38) * t99 - MDP(30)) * (g(1) * t81 + g(2) * t78 - g(3) * t113); (-g(1) * t74 + g(2) * t139 - g(3) * (-t88 * t99 + t116)) * MDP(37) + (g(1) * t75 + g(2) * t138 - g(3) * (-t88 * t104 - t117)) * MDP(38);];
taug  = t1;
