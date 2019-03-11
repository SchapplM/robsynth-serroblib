% Calculate Gravitation load on the joints for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:09:36
% EndTime: 2019-03-10 01:09:37
% DurationCPUTime: 0.28s
% Computational Cost: add. (342->73), mult. (368->106), div. (0->0), fcn. (335->10), ass. (0->43)
t134 = MDP(17) - MDP(32);
t105 = cos(qJ(2));
t104 = cos(qJ(4));
t99 = qJ(4) + qJ(5);
t94 = cos(t99);
t89 = t104 * pkin(4) + pkin(5) * t94;
t87 = pkin(3) + t89;
t100 = qJ(2) + qJ(3);
t93 = sin(t100);
t95 = cos(t100);
t98 = -qJ(6) - pkin(10) - pkin(9);
t115 = t95 * t87 - t93 * t98;
t133 = t105 * pkin(2) + t115;
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t114 = g(1) * t106 + g(2) * t103;
t129 = g(3) * t93;
t123 = t103 * t95;
t92 = sin(t99);
t77 = t106 * t94 + t123 * t92;
t122 = t106 * t95;
t79 = t103 * t94 - t122 * t92;
t132 = -g(1) * t79 + g(2) * t77 + t92 * t129;
t75 = -g(3) * t95 + t114 * t93;
t78 = t106 * t92 - t123 * t94;
t80 = t103 * t92 + t122 * t94;
t124 = t132 * MDP(30) + (g(1) * t80 - g(2) * t78 + t129 * t94) * MDP(31);
t101 = sin(qJ(4));
t88 = t101 * pkin(4) + pkin(5) * t92;
t121 = -pkin(8) - pkin(7) - t88;
t120 = t103 * t101;
t119 = t103 * t104;
t118 = t106 * t101;
t117 = t106 * t104;
t112 = t87 * t93 + t95 * t98;
t111 = pkin(1) + t133;
t110 = t134 * (t114 * t95 + t129) + (MDP(23) * t104 - MDP(24) * t101 + MDP(30) * t94 - MDP(31) * t92 + MDP(16)) * t75;
t102 = sin(qJ(2));
t85 = t117 * t95 + t120;
t84 = -t118 * t95 + t119;
t83 = -t119 * t95 + t118;
t82 = t120 * t95 + t117;
t1 = [t114 * MDP(3) + (-g(1) * t83 - g(2) * t85) * MDP(23) + (-g(1) * t82 - g(2) * t84) * MDP(24) + (-g(1) * t78 - g(2) * t80) * MDP(30) + (-g(1) * t77 - g(2) * t79) * MDP(31) + ((g(1) * t121 - g(2) * t111) * t106 + (g(1) * t111 + g(2) * t121) * t103) * MDP(33) + (-t102 * MDP(10) + t95 * MDP(16) + t105 * MDP(9) - t134 * t93 + MDP(2)) * (g(1) * t103 - g(2) * t106); (-g(3) * t105 + t102 * t114) * MDP(9) + (g(3) * t102 + t105 * t114) * MDP(10) + (-g(3) * t133 + t114 * (pkin(2) * t102 + t112)) * MDP(33) + t110; (-g(3) * t115 + t114 * t112) * MDP(33) + t110; (-g(1) * t84 + g(2) * t82 + t101 * t129) * MDP(23) + (g(1) * t85 - g(2) * t83 + t104 * t129) * MDP(24) + (-g(1) * (t103 * t89 - t122 * t88) - g(2) * (-t106 * t89 - t123 * t88) + t88 * t129) * MDP(33) + t124; t132 * MDP(33) * pkin(5) + t124; -t75 * MDP(33);];
taug  = t1;
