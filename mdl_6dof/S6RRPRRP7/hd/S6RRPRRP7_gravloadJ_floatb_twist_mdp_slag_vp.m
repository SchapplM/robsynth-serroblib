% Calculate Gravitation load on the joints for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:20:07
% EndTime: 2019-03-09 12:20:09
% DurationCPUTime: 0.70s
% Computational Cost: add. (305->94), mult. (750->127), div. (0->0), fcn. (824->8), ass. (0->43)
t114 = sin(qJ(5));
t118 = cos(qJ(5));
t117 = sin(qJ(1));
t115 = sin(qJ(4));
t119 = cos(qJ(2));
t116 = sin(qJ(2));
t138 = cos(qJ(4));
t129 = t116 * t138;
t98 = -t119 * t115 + t129;
t91 = t98 * t117;
t120 = cos(qJ(1));
t131 = t119 * t120;
t93 = t115 * t131 - t120 * t129;
t97 = t116 * t115 + t119 * t138;
t122 = g(1) * t93 - g(2) * t91 + g(3) * t97;
t144 = MDP(28) - MDP(31);
t145 = MDP(27) + MDP(29);
t139 = g(3) * t98;
t92 = t97 * t117;
t94 = t97 * t120;
t82 = g(1) * t94 + g(2) * t92 + t139;
t154 = -t82 * MDP(30) + (-t114 * t144 + t118 * t145 + MDP(20)) * t122;
t151 = t122 * (pkin(5) * t118 + qJ(6) * t114 + pkin(4));
t150 = MDP(9) + MDP(11);
t149 = MDP(10) - MDP(13);
t99 = g(1) * t120 + g(2) * t117;
t147 = t99 * t116;
t130 = t119 * pkin(2) + t116 * qJ(3);
t137 = pkin(3) * t119;
t136 = g(1) * t117;
t133 = t116 * t120;
t128 = t120 * pkin(1) + pkin(2) * t131 + t117 * pkin(7) + qJ(3) * t133;
t83 = t92 * t114 - t118 * t120;
t84 = t114 * t120 + t92 * t118;
t124 = -pkin(1) - t130;
t87 = t114 * t94 + t117 * t118;
t75 = g(1) * t87 + g(2) * t83 + t114 * t139;
t111 = t120 * pkin(7);
t104 = qJ(3) * t131;
t102 = t117 * t119 * qJ(3);
t89 = -g(3) * t119 + t147;
t88 = -t117 * t114 + t118 * t94;
t1 = [(-g(1) * t111 - g(2) * t128 - t124 * t136) * MDP(14) + (g(1) * t92 - g(2) * t94) * MDP(20) + (-g(1) * (-t92 * pkin(4) - pkin(5) * t84 - pkin(8) * t120 + t91 * pkin(9) - qJ(6) * t83 + t111) - g(2) * (pkin(3) * t131 + t94 * pkin(4) + t88 * pkin(5) + t93 * pkin(9) + t87 * qJ(6) + t128) + (-g(1) * (t124 - t137) + g(2) * pkin(8)) * t117) * MDP(32) + (MDP(3) - MDP(12)) * t99 + t145 * (g(1) * t84 - g(2) * t88) + t144 * (-g(1) * t83 + g(2) * t87) + (-MDP(21) + MDP(30)) * (-g(1) * t91 - g(2) * t93) + (-t149 * t116 + t150 * t119 + MDP(2)) * (-g(2) * t120 + t136); (-g(1) * (-pkin(2) * t133 + t104) - g(2) * (-pkin(2) * t116 * t117 + t102) - g(3) * t130) * MDP(14) - t82 * MDP(21) + (-g(1) * (-t94 * pkin(9) + t104) - g(2) * (-pkin(9) * t92 + t102) - g(3) * (-pkin(9) * t98 + t130 + t137) + (pkin(2) + pkin(3)) * t147 - t151) * MDP(32) + t149 * (g(3) * t116 + t99 * t119) + t150 * t89 - t154; (-MDP(14) - MDP(32)) * t89; (-pkin(9) * MDP(32) + MDP(21)) * t82 + MDP(32) * t151 + t154; (-g(1) * (-pkin(5) * t87 + qJ(6) * t88) - g(2) * (-pkin(5) * t83 + qJ(6) * t84) - (-pkin(5) * t114 + qJ(6) * t118) * t139) * MDP(32) + t145 * t75 + t144 * (g(1) * t88 + g(2) * t84 + t118 * t139); -t75 * MDP(32);];
taug  = t1;
