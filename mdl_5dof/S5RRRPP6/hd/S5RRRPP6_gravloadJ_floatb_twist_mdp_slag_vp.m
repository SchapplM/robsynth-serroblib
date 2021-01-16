% Calculate Gravitation load on the joints for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:38
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:37:09
% EndTime: 2021-01-15 22:37:12
% DurationCPUTime: 0.36s
% Computational Cost: add. (265->85), mult. (409->122), div. (0->0), fcn. (387->10), ass. (0->50)
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t112 = sin(pkin(8));
t113 = cos(pkin(8));
t98 = pkin(4) * t113 + qJ(5) * t112 + pkin(3);
t99 = -t112 * pkin(4) + qJ(5) * t113;
t80 = -t98 * t115 + t99 * t118;
t144 = -pkin(6) + t80;
t125 = t99 * t115 + t98 * t118;
t143 = MDP(18) + MDP(22);
t142 = MDP(19) - MDP(24);
t140 = MDP(10) - MDP(20) - MDP(23);
t116 = sin(qJ(2));
t119 = cos(qJ(2));
t117 = sin(qJ(1));
t120 = cos(qJ(1));
t123 = g(1) * t120 + g(2) * t117;
t85 = -g(3) * t119 + t116 * t123;
t137 = g(3) * t116;
t114 = qJ(4) + pkin(7);
t103 = t114 * t116;
t135 = pkin(1) + t103;
t132 = t119 * t117;
t131 = t119 * t120;
t111 = qJ(3) + pkin(8);
t107 = sin(t111);
t130 = t120 * t107;
t108 = cos(t111);
t129 = t120 * t108;
t128 = t120 * t115;
t127 = t120 * t118;
t91 = t117 * t118 - t119 * t128;
t89 = t115 * t132 + t127;
t106 = t118 * pkin(3) + pkin(2);
t105 = t115 * pkin(3) + pkin(6);
t104 = t114 * t119;
t100 = t106 * t119;
t97 = t112 * t118 + t113 * t115;
t96 = t112 * t115 - t113 * t118;
t92 = t117 * t115 + t119 * t127;
t90 = -t118 * t132 + t128;
t87 = t100 + t135;
t84 = t117 * t107 + t119 * t129;
t83 = -t117 * t108 + t119 * t130;
t82 = t108 * t132 - t130;
t81 = t107 * t132 + t129;
t79 = pkin(2) + t125;
t78 = t79 * t119;
t74 = t78 + t135;
t1 = [t123 * MDP(3) + (-g(1) * t90 - g(2) * t92) * MDP(16) + (-g(1) * t89 - g(2) * t91) * MDP(17) + (-g(1) * (t105 * t120 - t87 * t117) - g(2) * (t105 * t117 + t87 * t120)) * MDP(21) + (-g(1) * (-t74 * t117 - t144 * t120) - g(2) * (-t144 * t117 + t74 * t120)) * MDP(25) + t143 * (g(1) * t82 - g(2) * t84) - t142 * (g(1) * t81 - g(2) * t83) + (MDP(9) * t119 - t140 * t116 + MDP(2)) * (g(1) * t117 - g(2) * t120); (-g(3) * (t100 + t103) - t123 * (-t116 * t106 + t104)) * MDP(21) + (-g(3) * (t78 + t103) - t123 * (-t79 * t116 + t104)) * MDP(25) + t140 * (t119 * t123 + t137) + (MDP(16) * t118 - MDP(17) * t115 - t142 * t107 + t143 * t108 + MDP(9)) * t85; (g(1) * t92 - g(2) * t90 + t118 * t137) * MDP(17) + (-g(1) * (t125 * t117 + t80 * t131) - g(2) * (-t125 * t120 + t80 * t132) - t80 * t137) * MDP(25) + t143 * (g(1) * t83 + g(2) * t81 + t107 * t137) + t142 * (g(1) * t84 + g(2) * t82 + t108 * t137) + (pkin(3) * MDP(21) + MDP(16)) * (-g(1) * t91 + g(2) * t89 + t115 * t137); (-MDP(21) - MDP(25)) * t85; (-g(1) * (t96 * t117 + t97 * t131) - g(2) * (-t96 * t120 + t97 * t132) - t97 * t137) * MDP(25);];
taug = t1;
