% Calculate Gravitation load on the joints for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:15
% EndTime: 2019-03-09 05:07:17
% DurationCPUTime: 0.46s
% Computational Cost: add. (322->74), mult. (434->111), div. (0->0), fcn. (461->10), ass. (0->34)
t81 = sin(qJ(4));
t86 = cos(qJ(3));
t102 = t81 * t86;
t79 = qJ(1) + pkin(10);
t77 = sin(t79);
t78 = cos(t79);
t85 = cos(qJ(4));
t68 = t102 * t77 + t78 * t85;
t101 = t85 * t86;
t69 = t101 * t77 - t78 * t81;
t80 = sin(qJ(6));
t84 = cos(qJ(6));
t100 = -t68 * t84 + t69 * t80;
t82 = sin(qJ(3));
t106 = g(3) * t82;
t70 = t102 * t78 - t77 * t85;
t71 = t101 * t78 + t77 * t81;
t61 = t70 * t84 - t71 * t80;
t62 = t70 * t80 + t71 * t84;
t93 = t80 * t81 + t84 * t85;
t94 = t80 * t85 - t81 * t84;
t95 = t68 * t80 + t69 * t84;
t123 = (-g(1) * t61 + g(2) * t100 + t94 * t106) * MDP(28) + (g(1) * t62 + g(2) * t95 + t93 * t106) * MDP(29);
t119 = MDP(11) - MDP(20);
t115 = MDP(17) + MDP(19);
t114 = MDP(18) - MDP(21);
t98 = g(1) * t78 + g(2) * t77;
t104 = t82 * pkin(8);
t92 = pkin(3) * t86 + pkin(2) + t104;
t91 = pkin(4) * t85 + qJ(5) * t81 + pkin(3);
t60 = g(1) * t70 + g(2) * t68 + t106 * t81;
t87 = cos(qJ(1));
t83 = sin(qJ(1));
t1 = [(g(1) * t87 + g(2) * t83) * MDP(3) + (-g(1) * (-pkin(1) * t83 - pkin(4) * t69 - qJ(5) * t68) - g(2) * (pkin(1) * t87 + t71 * pkin(4) + t70 * qJ(5)) + (-g(1) * pkin(7) - g(2) * t92) * t78 + (-g(2) * pkin(7) + g(1) * t92) * t77) * MDP(22) + (g(1) * t95 - g(2) * t62) * MDP(28) + (-g(1) * t100 - g(2) * t61) * MDP(29) - t114 * (g(1) * t68 - g(2) * t70) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t83 - g(2) * t87) + t115 * (g(1) * t69 - g(2) * t71) + (t86 * MDP(10) - t119 * t82) * (g(1) * t77 - g(2) * t78); (-MDP(22) - MDP(4)) * g(3); (-g(3) * (t86 * t91 + t104) - t98 * (pkin(8) * t86 - t82 * t91)) * MDP(22) + t119 * (t86 * t98 + t106) + (-t93 * MDP(28) + t94 * MDP(29) + t114 * t81 - t115 * t85 - MDP(10)) * (g(3) * t86 - t82 * t98); (-g(1) * (-pkin(4) * t70 + qJ(5) * t71) - g(2) * (-pkin(4) * t68 + qJ(5) * t69) - (-pkin(4) * t81 + qJ(5) * t85) * t106) * MDP(22) + t114 * (g(1) * t71 + g(2) * t69 + t106 * t85) + t115 * t60 - t123; -t60 * MDP(22); t123;];
taug  = t1;
