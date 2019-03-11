% Calculate Gravitation load on the joints for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:43
% EndTime: 2019-03-09 06:32:44
% DurationCPUTime: 0.48s
% Computational Cost: add. (296->81), mult. (419->112), div. (0->0), fcn. (412->8), ass. (0->41)
t102 = sin(qJ(1));
t105 = cos(qJ(1));
t137 = -g(1) * t102 + g(2) * t105;
t135 = MDP(26) + MDP(28);
t134 = MDP(27) - MDP(30);
t136 = MDP(13) - MDP(29);
t100 = sin(qJ(4));
t104 = cos(qJ(3));
t126 = g(3) * t104;
t103 = cos(qJ(4));
t118 = t103 * t105;
t101 = sin(qJ(3));
t122 = t101 * t102;
t83 = -t100 * t122 + t118;
t119 = t102 * t103;
t121 = t101 * t105;
t85 = t100 * t121 + t119;
t133 = -g(1) * t83 - g(2) * t85 + t100 * t126;
t125 = t105 * pkin(1) + t102 * qJ(2);
t99 = qJ(4) + qJ(5);
t93 = sin(t99);
t124 = t104 * t93;
t94 = cos(t99);
t123 = t104 * t94;
t106 = -pkin(9) - pkin(8);
t117 = t104 * t106;
t115 = pkin(4) * t100 + pkin(7);
t78 = -t105 * t94 + t93 * t122;
t80 = t102 * t94 + t93 * t121;
t72 = g(1) * t78 - g(2) * t80 + g(3) * t124;
t79 = t105 * t93 + t94 * t122;
t81 = -t102 * t93 + t94 * t121;
t114 = t135 * t72 + t134 * (g(1) * t79 - g(2) * t81 + g(3) * t123);
t92 = pkin(4) * t103 + pkin(3);
t112 = t101 * t92 + t117;
t111 = pkin(5) * t94 + qJ(6) * t93 + t92;
t107 = -g(1) * (-t78 * pkin(5) + qJ(6) * t79) - g(2) * (t80 * pkin(5) - qJ(6) * t81) - g(3) * (-pkin(5) * t124 + qJ(6) * t123);
t96 = t105 * qJ(2);
t86 = -t100 * t102 + t101 * t118;
t84 = t100 * t105 + t101 * t119;
t1 = [(-g(1) * (-t102 * pkin(1) + t96) - g(2) * t125) * MDP(6) + (-g(1) * t86 - g(2) * t84) * MDP(19) + (g(1) * t85 - g(2) * t83) * MDP(20) + (-g(1) * (t81 * pkin(5) + t80 * qJ(6) + t96) - g(2) * (t79 * pkin(5) + t78 * qJ(6) + t125) + (-g(1) * t112 - g(2) * t115) * t105 + (-g(1) * (-pkin(1) - t115) - g(2) * t112) * t102) * MDP(31) - (MDP(2) - MDP(4)) * t137 + t135 * (-g(1) * t81 - g(2) * t79) + t134 * (g(1) * t80 + g(2) * t78) + (-t101 * MDP(12) - t136 * t104 + MDP(3) - MDP(5)) * (g(1) * t105 + g(2) * t102); -(-MDP(31) - MDP(6)) * t137; (-g(3) * (-t111 * t101 - t117) + t137 * (-t101 * t106 + t111 * t104)) * MDP(31) + t136 * (-t101 * t137 + t126) + (-t103 * MDP(19) + t100 * MDP(20) + t134 * t93 - t135 * t94 - MDP(12)) * (-g(3) * t101 - t137 * t104); t133 * MDP(19) + (g(1) * t84 - g(2) * t86 + t103 * t126) * MDP(20) + (t133 * pkin(4) + t107) * MDP(31) + t114; t107 * MDP(31) + t114; -t72 * MDP(31);];
taug  = t1;
