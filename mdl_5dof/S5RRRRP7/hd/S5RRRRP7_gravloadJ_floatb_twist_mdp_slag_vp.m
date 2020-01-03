% Calculate Gravitation load on the joints for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:44
% EndTime: 2019-12-31 21:57:45
% DurationCPUTime: 0.31s
% Computational Cost: add. (264->65), mult. (364->91), div. (0->0), fcn. (343->8), ass. (0->37)
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t119 = pkin(4) * t89 + qJ(5) * t86;
t118 = MDP(17) - MDP(26);
t115 = MDP(23) + MDP(25);
t114 = MDP(24) - MDP(27);
t85 = qJ(2) + qJ(3);
t82 = sin(t85);
t83 = cos(t85);
t117 = t83 * pkin(3) + t82 * pkin(8);
t88 = sin(qJ(1));
t91 = cos(qJ(1));
t100 = g(1) * t91 + g(2) * t88;
t116 = t100 * t82;
t87 = sin(qJ(2));
t113 = pkin(2) * t87;
t111 = pkin(8) * t83;
t108 = g(3) * t82;
t107 = t88 * t86;
t106 = t88 * t89;
t105 = t91 * t86;
t104 = t91 * t89;
t102 = t119 * t83 + t117;
t90 = cos(qJ(2));
t84 = t90 * pkin(2);
t98 = t84 + pkin(1) + t117;
t96 = t118 * (t100 * t83 + t108) + (-t114 * t86 + t115 * t89 + MDP(16)) * (-g(3) * t83 + t116);
t69 = t83 * t107 + t104;
t71 = t83 * t105 - t106;
t56 = g(1) * t71 + g(2) * t69 + t86 * t108;
t93 = (pkin(3) + t119) * t116;
t92 = -pkin(7) - pkin(6);
t76 = t91 * t111;
t74 = t88 * t111;
t72 = t83 * t104 + t107;
t70 = t83 * t106 - t105;
t1 = [t100 * MDP(3) + (-g(1) * (-t70 * pkin(4) - t69 * qJ(5)) - g(2) * (t72 * pkin(4) + t71 * qJ(5)) + (g(1) * t92 - g(2) * t98) * t91 + (g(1) * t98 + g(2) * t92) * t88) * MDP(28) + t115 * (g(1) * t70 - g(2) * t72) - t114 * (g(1) * t69 - g(2) * t71) + (-t87 * MDP(10) + t83 * MDP(16) + t90 * MDP(9) - t118 * t82 + MDP(2)) * (g(1) * t88 - g(2) * t91); (-g(3) * t90 + t100 * t87) * MDP(9) + (g(3) * t87 + t100 * t90) * MDP(10) + (-g(1) * (-t91 * t113 + t76) - g(2) * (-t88 * t113 + t74) - g(3) * (t84 + t102) + t93) * MDP(28) + t96; (-g(1) * t76 - g(2) * t74 - g(3) * t102 + t93) * MDP(28) + t96; (-g(1) * (-t71 * pkin(4) + t72 * qJ(5)) - g(2) * (-t69 * pkin(4) + t70 * qJ(5)) - (-pkin(4) * t86 + qJ(5) * t89) * t108) * MDP(28) + t114 * (g(1) * t72 + g(2) * t70 + t89 * t108) + t115 * t56; -t56 * MDP(28);];
taug = t1;
