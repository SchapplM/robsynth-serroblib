% Calculate Gravitation load on the joints for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:11:00
% EndTime: 2019-03-09 09:11:02
% DurationCPUTime: 0.58s
% Computational Cost: add. (244->94), mult. (592->155), div. (0->0), fcn. (675->10), ass. (0->44)
t134 = MDP(10) - MDP(13) - MDP(16);
t101 = cos(qJ(6));
t102 = cos(qJ(5));
t104 = cos(qJ(1));
t96 = sin(pkin(6));
t120 = t104 * t96;
t100 = sin(qJ(1));
t103 = cos(qJ(2));
t119 = cos(pkin(6));
t113 = t104 * t119;
t99 = sin(qJ(2));
t85 = t100 * t103 + t99 * t113;
t98 = sin(qJ(5));
t74 = t85 * t102 + t98 * t120;
t84 = t100 * t99 - t103 * t113;
t97 = sin(qJ(6));
t133 = t101 * t84 + t74 * t97;
t132 = t101 * t74 - t84 * t97;
t130 = MDP(9) + MDP(11) + MDP(15);
t129 = g(3) * t96;
t127 = t96 * t99;
t121 = t103 * t96;
t126 = pkin(2) * t121 + qJ(3) * t127;
t125 = t100 * t96;
t123 = t102 * t96;
t122 = t102 * t97;
t118 = t101 * t102;
t117 = t102 * t103;
t116 = -t84 * pkin(2) + qJ(3) * t85;
t114 = t100 * t119;
t86 = t103 * t114 + t104 * t99;
t87 = t103 * t104 - t99 * t114;
t115 = -t86 * pkin(2) + qJ(3) * t87;
t110 = g(1) * t100 - g(2) * t104;
t109 = t104 * pkin(1) + t87 * pkin(2) + pkin(8) * t125 + qJ(3) * t86;
t73 = t102 * t120 - t85 * t98;
t107 = -t100 * pkin(1) - t85 * pkin(2) + pkin(8) * t120 - t84 * qJ(3);
t105 = -g(1) * t86 - g(2) * t84 + g(3) * t121;
t83 = -t119 * t98 + t99 * t123;
t77 = t102 * t87 - t98 * t125;
t76 = -t100 * t123 - t87 * t98;
t67 = t101 * t77 - t86 * t97;
t66 = -t101 * t86 - t77 * t97;
t1 = [t110 * MDP(2) + (-g(1) * t107 - g(2) * t109) * MDP(14) + (-g(1) * (-t85 * pkin(3) - qJ(4) * t120 + t107) - g(2) * (pkin(3) * t87 - qJ(4) * t125 + t109)) * MDP(18) + (g(1) * t74 - g(2) * t77) * MDP(24) + (g(1) * t73 - g(2) * t76) * MDP(25) + (g(1) * t132 - g(2) * t67) * MDP(31) + (-g(1) * t133 - g(2) * t66) * MDP(32) + t130 * (g(1) * t85 - g(2) * t87) - t134 * (g(1) * t84 - g(2) * t86) + (MDP(3) + (-MDP(12) + MDP(17)) * t96) * (g(1) * t104 + g(2) * t100); (-g(1) * t115 - g(2) * t116 - g(3) * t126) * MDP(14) + (-g(1) * (-pkin(3) * t86 + t115) - g(2) * (-pkin(3) * t84 + t116) - g(3) * (pkin(3) * t121 + t126)) * MDP(18) + (-g(1) * (-t86 * t118 - t87 * t97) - g(2) * (-t84 * t118 - t85 * t97) - (t101 * t117 - t97 * t99) * t129) * MDP(31) + (-g(1) * (-t101 * t87 + t86 * t122) - g(2) * (-t101 * t85 + t84 * t122) - (-t101 * t99 - t97 * t117) * t129) * MDP(32) + t134 * (g(1) * t87 + g(2) * t85 + g(3) * t127) + (-t102 * MDP(24) + MDP(25) * t98 - t130) * t105; (MDP(14) + MDP(18)) * t105; (g(3) * t119 + t110 * t96) * MDP(18); (g(1) * t77 + g(2) * t74 + g(3) * t83) * MDP(25) + (-MDP(31) * t101 + MDP(32) * t97 - MDP(24)) * (g(1) * t76 + g(2) * t73 + g(3) * (-t119 * t102 - t98 * t127)); (-g(1) * t66 + g(2) * t133 - g(3) * (t101 * t121 - t83 * t97)) * MDP(31) + (g(1) * t67 + g(2) * t132 - g(3) * (-t101 * t83 - t97 * t121)) * MDP(32);];
taug  = t1;
