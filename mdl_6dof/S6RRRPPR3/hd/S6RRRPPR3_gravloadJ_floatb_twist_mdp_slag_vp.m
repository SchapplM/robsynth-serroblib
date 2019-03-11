% Calculate Gravitation load on the joints for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:30:13
% EndTime: 2019-03-09 15:30:14
% DurationCPUTime: 0.38s
% Computational Cost: add. (278->74), mult. (343->93), div. (0->0), fcn. (292->8), ass. (0->42)
t87 = qJ(2) + qJ(3);
t84 = sin(t87);
t85 = cos(t87);
t103 = t85 * pkin(3) + t84 * qJ(4);
t92 = cos(qJ(2));
t100 = t92 * pkin(2) + t103;
t119 = pkin(1) + t100;
t93 = cos(qJ(1));
t118 = g(2) * t93;
t117 = MDP(16) + MDP(18) - MDP(23);
t116 = -MDP(17) + MDP(20) + MDP(22);
t111 = g(1) * t93;
t90 = sin(qJ(1));
t73 = g(2) * t90 + t111;
t115 = t73 * t84;
t89 = sin(qJ(2));
t113 = pkin(2) * t89;
t112 = pkin(3) * t84;
t108 = g(3) * t85;
t80 = t85 * pkin(4);
t88 = sin(qJ(6));
t107 = t90 * t88;
t91 = cos(qJ(6));
t106 = t90 * t91;
t105 = t93 * t88;
t104 = t93 * t91;
t102 = qJ(4) * t85;
t94 = -pkin(8) - pkin(7);
t101 = qJ(5) + t94;
t99 = t119 * t118;
t98 = -t112 - t113;
t72 = g(1) * t90 - t118;
t64 = -t108 + t115;
t96 = t117 * t64 + (-MDP(31) * t91 + MDP(32) * t88 - t116) * (g(3) * t84 + t73 * t85);
t95 = (pkin(3) + pkin(4)) * t115;
t76 = t93 * t102;
t74 = t90 * t102;
t71 = t84 * t104 - t107;
t70 = -t84 * t105 - t106;
t69 = -t84 * t106 - t105;
t68 = t84 * t107 - t104;
t1 = [(t94 * t111 - t99 + (g(1) * t119 + g(2) * t94) * t90) * MDP(21) + (-t99 + (g(1) * t101 - g(2) * t80) * t93 + (-g(1) * (-t119 - t80) + g(2) * t101) * t90) * MDP(25) + (-g(1) * t69 - g(2) * t71) * MDP(31) + (-g(1) * t68 - g(2) * t70) * MDP(32) + (MDP(3) - MDP(19) + MDP(24)) * t73 + (-t89 * MDP(10) + t92 * MDP(9) + t116 * t84 + t117 * t85 + MDP(2)) * t72; (-g(3) * t92 + t73 * t89) * MDP(9) + (g(3) * t89 + t73 * t92) * MDP(10) + (-g(1) * (t98 * t93 + t76) - g(2) * (t98 * t90 + t74) - g(3) * t100) * MDP(21) + (-g(1) * (-t93 * t113 + t76) - g(2) * (-t90 * t113 + t74) - g(3) * (t80 + t100) + t95) * MDP(25) + t96; (-g(1) * (-t93 * t112 + t76) - g(2) * (-t90 * t112 + t74) - g(3) * t103) * MDP(21) + (-g(1) * t76 - g(2) * t74 - g(3) * (t80 + t103) + t95) * MDP(25) + t96; (-MDP(21) - MDP(25)) * t64; t72 * MDP(25); (-g(1) * t70 + g(2) * t68 - t88 * t108) * MDP(31) + (g(1) * t71 - g(2) * t69 - t91 * t108) * MDP(32);];
taug  = t1;
