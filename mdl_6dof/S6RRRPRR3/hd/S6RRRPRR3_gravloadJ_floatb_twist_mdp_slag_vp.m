% Calculate Gravitation load on the joints for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:13:23
% EndTime: 2019-03-09 18:13:24
% DurationCPUTime: 0.29s
% Computational Cost: add. (331->65), mult. (431->96), div. (0->0), fcn. (440->10), ass. (0->38)
t113 = cos(qJ(5));
t89 = qJ(2) + qJ(3);
t86 = sin(t89);
t87 = cos(t89);
t91 = sin(qJ(5));
t101 = t87 * t113 + t86 * t91;
t116 = t86 * t113 - t87 * t91;
t114 = g(3) * t116;
t93 = sin(qJ(1));
t68 = t116 * t93;
t69 = t101 * t93;
t96 = cos(qJ(1));
t70 = t116 * t96;
t71 = t101 * t96;
t90 = sin(qJ(6));
t94 = cos(qJ(6));
t122 = (MDP(34) * t94 - MDP(35) * t90 + MDP(27)) * (-g(1) * t70 - g(2) * t68 + g(3) * t101) + (g(1) * t71 + g(2) * t69 + t114) * MDP(28);
t111 = t87 * pkin(3) + t86 * qJ(4);
t95 = cos(qJ(2));
t121 = t95 * pkin(2) + t111;
t120 = MDP(16) + MDP(18);
t119 = MDP(17) - MDP(20);
t115 = pkin(3) * t86;
t110 = qJ(4) * t87;
t92 = sin(qJ(2));
t106 = -pkin(2) * t92 - t115;
t78 = g(1) * t96 + g(2) * t93;
t104 = t69 * t94 + t96 * t90;
t103 = t69 * t90 - t96 * t94;
t102 = pkin(1) + t121;
t66 = -g(3) * t87 + t78 * t86;
t98 = t119 * (g(3) * t86 + t78 * t87) + t120 * t66 - t122;
t97 = -pkin(8) - pkin(7);
t80 = t96 * t110;
t79 = t93 * t110;
t61 = t71 * t94 - t93 * t90;
t60 = -t71 * t90 - t93 * t94;
t1 = [((g(1) * t97 - g(2) * t102) * t96 + (g(1) * t102 + g(2) * t97) * t93) * MDP(21) + (g(1) * t69 - g(2) * t71) * MDP(27) + (g(1) * t68 - g(2) * t70) * MDP(28) + (g(1) * t104 - g(2) * t61) * MDP(34) + (-g(1) * t103 - g(2) * t60) * MDP(35) + (MDP(3) - MDP(19)) * t78 + (-t92 * MDP(10) + t95 * MDP(9) - t119 * t86 + t120 * t87 + MDP(2)) * (g(1) * t93 - g(2) * t96); (-g(3) * t95 + t78 * t92) * MDP(9) + (g(3) * t92 + t78 * t95) * MDP(10) + (-g(1) * (t106 * t96 + t80) - g(2) * (t106 * t93 + t79) - g(3) * t121) * MDP(21) + t98; (-g(1) * (-t96 * t115 + t80) - g(2) * (-t93 * t115 + t79) - g(3) * t111) * MDP(21) + t98; -t66 * MDP(21); t122; (-g(1) * t60 + g(2) * t103 + t90 * t114) * MDP(34) + (g(1) * t61 + g(2) * t104 + t94 * t114) * MDP(35);];
taug  = t1;
