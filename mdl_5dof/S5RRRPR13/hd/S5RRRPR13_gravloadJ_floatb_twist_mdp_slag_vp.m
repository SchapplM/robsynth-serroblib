% Calculate Gravitation load on the joints for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:34
% EndTime: 2019-12-31 21:46:36
% DurationCPUTime: 0.70s
% Computational Cost: add. (252->89), mult. (632->147), div. (0->0), fcn. (743->10), ass. (0->41)
t120 = MDP(16) - MDP(19);
t119 = MDP(10) - MDP(18);
t115 = MDP(17) - MDP(20);
t111 = cos(qJ(1));
t81 = sin(pkin(5));
t100 = t81 * t111;
t84 = sin(qJ(2));
t85 = sin(qJ(1));
t88 = cos(qJ(2));
t101 = cos(pkin(5));
t94 = t101 * t111;
t73 = t84 * t94 + t85 * t88;
t83 = sin(qJ(3));
t87 = cos(qJ(3));
t64 = t87 * t100 + t73 * t83;
t72 = t85 * t84 - t88 * t94;
t82 = sin(qJ(5));
t86 = cos(qJ(5));
t117 = t64 * t82 + t72 * t86;
t116 = t64 * t86 - t72 * t82;
t112 = g(3) * t81;
t108 = t81 * t84;
t107 = t81 * t85;
t106 = t81 * t87;
t105 = t82 * t83;
t104 = t82 * t88;
t103 = t83 * t86;
t102 = t86 * t88;
t65 = -t83 * t100 + t73 * t87;
t99 = t85 * t101;
t75 = t111 * t88 - t84 * t99;
t95 = -g(1) * t75 - g(2) * t73;
t68 = -t85 * t106 + t75 * t83;
t70 = -t101 * t87 + t83 * t108;
t92 = g(1) * t68 + g(2) * t64 + g(3) * t70;
t74 = t111 * t84 + t88 * t99;
t71 = t101 * t83 + t84 * t106;
t69 = t83 * t107 + t75 * t87;
t61 = t68 * t82 + t74 * t86;
t60 = t68 * t86 - t74 * t82;
t1 = [(g(1) * t85 - g(2) * t111) * MDP(2) + (g(1) * t111 + g(2) * t85) * MDP(3) + (g(1) * t73 - g(2) * t75) * MDP(9) + (-g(1) * (-t85 * pkin(1) - t73 * pkin(2) - pkin(3) * t65 + pkin(7) * t100 - t72 * pkin(8) - qJ(4) * t64) - g(2) * (t111 * pkin(1) + t75 * pkin(2) + t69 * pkin(3) + pkin(7) * t107 + t74 * pkin(8) + t68 * qJ(4))) * MDP(21) + (g(1) * t117 - g(2) * t61) * MDP(27) + (g(1) * t116 - g(2) * t60) * MDP(28) + t115 * (-g(1) * t64 + g(2) * t68) - t120 * (-g(1) * t65 + g(2) * t69) - t119 * (g(1) * t72 - g(2) * t74); (-t84 * t112 + t95) * MDP(21) * pkin(8) + (-g(1) * (-t74 * t105 + t75 * t86) - g(2) * (-t72 * t105 + t73 * t86) - (t83 * t104 + t84 * t86) * t112) * MDP(27) + (-g(1) * (-t74 * t103 - t75 * t82) - g(2) * (-t72 * t103 - t73 * t82) - (t83 * t102 - t82 * t84) * t112) * MDP(28) + t119 * (g(3) * t108 - t95) + (-MDP(9) + (-pkin(3) * t87 - qJ(4) * t83 - pkin(2)) * MDP(21) - t120 * t87 + t115 * t83) * (-g(1) * t74 - g(2) * t72 + t88 * t112); (-g(1) * (-t68 * pkin(3) + t69 * qJ(4)) - g(2) * (-t64 * pkin(3) + t65 * qJ(4)) - g(3) * (-t70 * pkin(3) + t71 * qJ(4))) * MDP(21) + (-MDP(27) * t82 - MDP(28) * t86 + t115) * (g(1) * t69 + g(2) * t65 + g(3) * t71) + t120 * t92; -t92 * MDP(21); (-g(1) * t60 - g(2) * t116 - g(3) * (t81 * t104 + t70 * t86)) * MDP(27) + (g(1) * t61 + g(2) * t117 - g(3) * (t81 * t102 - t70 * t82)) * MDP(28);];
taug = t1;
