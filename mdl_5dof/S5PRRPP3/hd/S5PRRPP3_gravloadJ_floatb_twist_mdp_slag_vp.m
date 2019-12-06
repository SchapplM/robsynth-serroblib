% Calculate Gravitation load on the joints for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:33
% EndTime: 2019-12-05 16:13:35
% DurationCPUTime: 0.47s
% Computational Cost: add. (164->69), mult. (426->104), div. (0->0), fcn. (441->8), ass. (0->45)
t94 = sin(qJ(3));
t128 = qJ(4) * t94 + pkin(2);
t127 = MDP(11) - MDP(14) - MDP(17);
t91 = sin(pkin(7));
t93 = cos(pkin(7));
t107 = g(1) * t93 + g(2) * t91;
t95 = sin(qJ(2));
t126 = t107 * t95;
t125 = MDP(12) + MDP(16);
t124 = MDP(13) - MDP(18);
t97 = cos(qJ(2));
t122 = pkin(6) * t97;
t119 = t94 * t95;
t118 = t94 * t97;
t92 = cos(pkin(8));
t117 = t95 * t92;
t96 = cos(qJ(3));
t116 = t95 * t96;
t115 = t96 * t97;
t114 = t97 * t92;
t112 = MDP(15) + MDP(19);
t111 = g(3) * t119;
t74 = t91 * t118 + t93 * t96;
t75 = t91 * t115 - t93 * t94;
t110 = -t74 * pkin(3) + qJ(4) * t75;
t76 = t93 * t118 - t91 * t96;
t77 = t93 * t115 + t91 * t94;
t109 = -t76 * pkin(3) + qJ(4) * t77;
t108 = pkin(3) * t115 + t95 * pkin(6) + t128 * t97;
t90 = sin(pkin(8));
t106 = -pkin(4) * t92 - qJ(5) * t90;
t104 = -t92 * t116 + t90 * t97;
t103 = t90 * t116 + t114;
t101 = g(1) * t76 + g(2) * t74 + t111;
t98 = (pkin(3) * t96 + t128) * t126;
t82 = t93 * t122;
t81 = qJ(4) * t116;
t80 = t91 * t122;
t79 = t96 * t114 + t95 * t90;
t78 = t90 * t115 - t117;
t71 = t104 * t93;
t70 = t103 * t93;
t69 = t104 * t91;
t68 = t103 * t91;
t1 = [(-MDP(1) - t112) * g(3); (g(3) * t95 + t107 * t97) * MDP(4) + (-g(1) * t82 - g(2) * t80 - g(3) * t108 + t98) * MDP(15) + (-g(1) * (pkin(4) * t71 - qJ(5) * t70 + t82) - g(2) * (pkin(4) * t69 - qJ(5) * t68 + t80) - g(3) * (pkin(4) * t79 + qJ(5) * t78 + t108) + t98) * MDP(19) + t125 * (-g(1) * t71 - g(2) * t69 - g(3) * t79) - t124 * (g(1) * t70 + g(2) * t68 - g(3) * t78) + (t96 * MDP(10) - t127 * t94 + MDP(3)) * (-g(3) * t97 + t126); (-g(1) * t109 - g(2) * t110 - g(3) * (-pkin(3) * t119 + t81)) * MDP(15) + (-g(1) * (t106 * t76 + t109) - g(2) * (t106 * t74 + t110) - g(3) * t81 - (-pkin(3) + t106) * t111) * MDP(19) + t127 * (g(1) * t77 + g(2) * t75 + g(3) * t116) + (-t124 * t90 + t125 * t92 + MDP(10)) * t101; -t112 * t101; (-g(1) * (-t93 * t117 + t77 * t90) - g(2) * (-t91 * t117 + t75 * t90) - g(3) * t103) * MDP(19);];
taug = t1;
