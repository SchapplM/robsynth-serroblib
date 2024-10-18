% Calculate Gravitation load on the joints for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:21
% EndTime: 2024-09-28 18:08:21
% DurationCPUTime: 0.25s
% Computational Cost: add. (457->80), mult. (411->137), div. (0->0), fcn. (446->22), ass. (0->57)
t97 = sin(pkin(5));
t120 = g(3) * t97;
t101 = sin(qJ(5));
t119 = t101 * t97;
t99 = cos(pkin(6));
t118 = t101 * t99;
t103 = cos(qJ(5));
t117 = t103 * t97;
t116 = t103 * t99;
t100 = cos(pkin(5));
t115 = t100 * t101;
t102 = sin(qJ(2));
t114 = t100 * t102;
t113 = t100 * t103;
t104 = cos(qJ(2));
t112 = t100 * t104;
t94 = qJ(2) + qJ(3);
t96 = sin(pkin(6));
t111 = t96 * t119;
t110 = t96 * t117;
t109 = t99 * t115;
t108 = t99 * t113;
t90 = pkin(5) - t94;
t89 = pkin(5) + t94;
t95 = sin(pkin(11));
t98 = cos(pkin(11));
t65 = -t103 * t98 + t95 * t109;
t66 = -t101 * t98 - t95 * t108;
t67 = t103 * t95 + t98 * t109;
t68 = -t101 * t95 + t98 * t108;
t69 = -t98 * t113 + t95 * t118;
t70 = -t98 * t115 - t95 * t116;
t71 = t95 * t113 + t98 * t118;
t72 = t95 * t115 - t98 * t116;
t86 = -qJ(4) + t90;
t77 = sin(t86) / 0.2e1;
t85 = qJ(4) + t89;
t81 = sin(t85);
t73 = t77 - t81 / 0.2e1;
t78 = cos(t85) / 0.2e1;
t82 = cos(t86);
t74 = t82 / 0.2e1 + t78;
t93 = qJ(4) + t94;
t87 = sin(t93);
t88 = cos(t93);
t107 = (-g(1) * (-t74 * t95 - t87 * t98) - g(2) * (t74 * t98 - t87 * t95) - g(3) * (t81 / 0.2e1 + t77)) * MDP(9) + (-g(1) * (-t73 * t95 - t88 * t98) - g(2) * (t73 * t98 - t88 * t95) - g(3) * (t78 - t82 / 0.2e1)) * MDP(10) + (-g(1) * (t65 * t87 - t71 * t88) - g(2) * (-t67 * t87 - t69 * t88) - (t103 * t88 - t87 * t118) * t120) * MDP(16) + (-g(1) * (-t66 * t87 + t72 * t88) - g(2) * (-t68 * t87 + t70 * t88) - (-t101 * t88 - t87 * t116) * t120) * MDP(17);
t79 = sin(t90) / 0.2e1;
t83 = sin(t89);
t75 = t79 - t83 / 0.2e1;
t80 = cos(t89) / 0.2e1;
t84 = cos(t90);
t76 = t84 / 0.2e1 + t80;
t91 = sin(t94);
t92 = cos(t94);
t106 = (-g(1) * (-t76 * t95 - t91 * t98) - g(2) * (t76 * t98 - t91 * t95) - g(3) * (t83 / 0.2e1 + t79)) * MDP(6) + (-g(1) * (-t75 * t95 - t92 * t98) - g(2) * (t75 * t98 - t92 * t95) - g(3) * (t80 - t84 / 0.2e1)) * MDP(7) + t107;
t105 = t88 * t97 * t99 + t100 * t96;
t1 = [-g(3) * MDP(1); (-g(1) * (-t98 * t102 - t95 * t112) - g(2) * (-t95 * t102 + t98 * t112) - t104 * t120) * MDP(3) + (-g(1) * (-t104 * t98 + t95 * t114) - g(2) * (-t104 * t95 - t98 * t114) + t102 * t120) * MDP(4) + t106; t106; t107; (-g(1) * (t95 * t110 + t66 * t88 + t72 * t87) - g(2) * (-t98 * t110 + t68 * t88 + t70 * t87) - g(3) * (t105 * t103 - t87 * t119)) * MDP(16) + (-g(1) * (-t95 * t111 + t65 * t88 + t71 * t87) - g(2) * (t98 * t111 - t67 * t88 + t69 * t87) - g(3) * (-t105 * t101 - t87 * t117)) * MDP(17);];
taug = t1;
