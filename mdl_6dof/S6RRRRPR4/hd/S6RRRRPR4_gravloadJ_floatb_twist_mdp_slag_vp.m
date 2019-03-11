% Calculate Gravitation load on the joints for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:09:20
% EndTime: 2019-03-09 22:09:21
% DurationCPUTime: 0.29s
% Computational Cost: add. (331->67), mult. (337->93), div. (0->0), fcn. (308->10), ass. (0->42)
t128 = MDP(17) - MDP(25);
t96 = cos(qJ(4));
t85 = t96 * pkin(4) + pkin(3);
t91 = qJ(2) + qJ(3);
t88 = sin(t91);
t89 = cos(t91);
t92 = -qJ(5) - pkin(9);
t107 = t89 * t85 - t88 * t92;
t97 = cos(qJ(2));
t127 = t97 * pkin(2) + t107;
t95 = sin(qJ(1));
t98 = cos(qJ(1));
t106 = g(1) * t98 + g(2) * t95;
t71 = -g(3) * t89 + t106 * t88;
t121 = g(3) * t88;
t87 = qJ(4) + pkin(11) + qJ(6);
t83 = sin(t87);
t118 = t95 * t83;
t84 = cos(t87);
t117 = t95 * t84;
t93 = sin(qJ(4));
t116 = t95 * t93;
t115 = t95 * t96;
t114 = t98 * t83;
t113 = t98 * t84;
t112 = t98 * t93;
t111 = t98 * t96;
t73 = t89 * t118 + t113;
t74 = -t89 * t117 + t114;
t75 = -t89 * t114 + t117;
t76 = t89 * t113 + t118;
t110 = (-g(1) * t75 + g(2) * t73 + t83 * t121) * MDP(32) + (g(1) * t76 - g(2) * t74 + t84 * t121) * MDP(33);
t108 = pkin(4) * t93 + pkin(7) + pkin(8);
t104 = t85 * t88 + t89 * t92;
t80 = -t89 * t112 + t115;
t78 = t89 * t116 + t111;
t103 = pkin(1) + t127;
t102 = t128 * (t106 * t89 + t121) + (MDP(23) * t96 - MDP(24) * t93 + MDP(32) * t84 - MDP(33) * t83 + MDP(16)) * t71;
t94 = sin(qJ(2));
t81 = t89 * t111 + t116;
t79 = -t89 * t115 + t112;
t1 = [t106 * MDP(3) + (-g(1) * t79 - g(2) * t81) * MDP(23) + (-g(1) * t78 - g(2) * t80) * MDP(24) + ((-g(1) * t108 - g(2) * t103) * t98 + (g(1) * t103 - g(2) * t108) * t95) * MDP(26) + (-g(1) * t74 - g(2) * t76) * MDP(32) + (-g(1) * t73 - g(2) * t75) * MDP(33) + (-t94 * MDP(10) + t89 * MDP(16) + t97 * MDP(9) - t128 * t88 + MDP(2)) * (g(1) * t95 - g(2) * t98); (-g(3) * t97 + t106 * t94) * MDP(9) + (g(3) * t94 + t106 * t97) * MDP(10) + (-g(3) * t127 + t106 * (pkin(2) * t94 + t104)) * MDP(26) + t102; (-g(3) * t107 + t106 * t104) * MDP(26) + t102; (g(1) * t81 - g(2) * t79 + t96 * t121) * MDP(24) + t110 + (pkin(4) * MDP(26) + MDP(23)) * (-g(1) * t80 + g(2) * t78 + t93 * t121); -t71 * MDP(26); t110;];
taug  = t1;
