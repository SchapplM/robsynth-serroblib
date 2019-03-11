% Calculate Gravitation load on the joints for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:40
% EndTime: 2019-03-08 19:37:42
% DurationCPUTime: 0.66s
% Computational Cost: add. (316->85), mult. (817->148), div. (0->0), fcn. (1013->12), ass. (0->48)
t135 = MDP(11) - MDP(14);
t134 = MDP(12) - MDP(15);
t103 = sin(qJ(2));
t100 = cos(pkin(6));
t106 = cos(qJ(2));
t120 = t100 * t106;
t97 = sin(pkin(10));
t99 = cos(pkin(10));
t133 = -t99 * t103 - t97 * t120;
t122 = cos(pkin(11));
t114 = t106 * t122;
t96 = sin(pkin(11));
t112 = -t103 * t96 + t114;
t108 = t100 * t112;
t89 = -t103 * t122 - t106 * t96;
t75 = t108 * t99 + t97 * t89;
t78 = -t108 * t97 + t89 * t99;
t98 = sin(pkin(6));
t125 = t98 * t103;
t85 = -t114 * t98 + t125 * t96;
t132 = -g(1) * t78 - g(2) * t75 + g(3) * t85;
t102 = sin(qJ(4));
t128 = t102 * t98;
t105 = cos(qJ(4));
t127 = t105 * t98;
t126 = t97 * t103;
t124 = t98 * t106;
t121 = t100 * t103;
t101 = sin(qJ(6));
t119 = t101 * t102;
t104 = cos(qJ(6));
t118 = t102 * t104;
t117 = MDP(16) + MDP(5);
t115 = t99 * t120;
t87 = t89 * t100;
t76 = t112 * t97 - t87 * t99;
t77 = -t112 * t99 - t97 * t87;
t70 = t102 * t76 + t127 * t99;
t72 = -t102 * t77 - t127 * t97;
t86 = t89 * t98;
t80 = -t100 * t105 - t102 * t86;
t111 = g(1) * t72 + g(2) * t70 + g(3) * t80;
t107 = -g(1) * t133 - g(3) * t124;
t90 = pkin(2) * t115;
t81 = t100 * t102 - t105 * t86;
t73 = -t105 * t77 + t128 * t97;
t71 = t105 * t76 - t128 * t99;
t1 = [(-MDP(1) - t117) * g(3); (-g(2) * (t115 - t126) + t107) * MDP(3) + (-g(1) * (-t106 * t99 + t121 * t97) - g(2) * (-t106 * t97 - t121 * t99) + g(3) * t125) * MDP(4) + (-g(2) * t90 + (g(2) * t126 + t107) * pkin(2)) * MDP(5) + (g(1) * t77 - g(2) * t76 + g(3) * t86) * MDP(13) + (-g(1) * (t133 * pkin(2) - t77 * pkin(8)) - g(2) * (-pkin(2) * t126 + pkin(8) * t76 + t90) - g(3) * (pkin(2) * t124 - t86 * pkin(8)) + t132 * (pkin(4) * t105 + qJ(5) * t102 + pkin(3))) * MDP(16) + (-g(1) * (-t104 * t77 + t119 * t78) - g(2) * (t104 * t76 + t119 * t75) - g(3) * (-t104 * t86 - t119 * t85)) * MDP(22) + (-g(1) * (t101 * t77 + t118 * t78) - g(2) * (-t101 * t76 + t118 * t75) - g(3) * (t101 * t86 - t118 * t85)) * MDP(23) + (-t134 * t102 + t135 * t105) * t132; t117 * (-g(3) * t100 + (-g(1) * t97 + g(2) * t99) * t98); (-g(1) * (-pkin(4) * t72 + qJ(5) * t73) - g(2) * (-pkin(4) * t70 + qJ(5) * t71) - g(3) * (-pkin(4) * t80 + qJ(5) * t81)) * MDP(16) + (-MDP(22) * t101 - MDP(23) * t104 + t134) * (g(1) * t73 + g(2) * t71 + g(3) * t81) + t135 * t111; -t111 * MDP(16); (-g(1) * (t101 * t78 + t104 * t72) - g(2) * (t101 * t75 + t104 * t70) - g(3) * (-t101 * t85 + t104 * t80)) * MDP(22) + (-g(1) * (-t101 * t72 + t104 * t78) - g(2) * (-t101 * t70 + t104 * t75) - g(3) * (-t101 * t80 - t104 * t85)) * MDP(23);];
taug  = t1;
