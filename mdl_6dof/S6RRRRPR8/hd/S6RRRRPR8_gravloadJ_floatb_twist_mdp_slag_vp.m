% Calculate Gravitation load on the joints for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:44:49
% EndTime: 2019-03-09 22:44:51
% DurationCPUTime: 0.55s
% Computational Cost: add. (458->89), mult. (620->132), div. (0->0), fcn. (670->10), ass. (0->48)
t113 = qJ(3) + qJ(4);
t111 = sin(t113);
t112 = cos(t113);
t114 = sin(qJ(6));
t118 = cos(qJ(6));
t128 = t111 * t114 + t112 * t118;
t129 = t111 * t118 - t112 * t114;
t117 = sin(qJ(1));
t120 = cos(qJ(2));
t121 = cos(qJ(1));
t141 = t121 * t112;
t100 = t117 * t111 + t120 * t141;
t142 = t121 * t111;
t99 = -t117 * t112 + t120 * t142;
t131 = t100 * t114 - t99 * t118;
t143 = t117 * t120;
t97 = t111 * t143 + t141;
t98 = t112 * t143 - t142;
t132 = t97 * t114 + t98 * t118;
t136 = t98 * t114 - t97 * t118;
t116 = sin(qJ(2));
t148 = g(3) * t116;
t88 = t100 * t118 + t99 * t114;
t166 = (g(1) * t131 + g(2) * t136 - t129 * t148) * MDP(34) + (g(1) * t88 + g(2) * t132 + t128 * t148) * MDP(35);
t164 = MDP(10) - MDP(26);
t163 = MDP(23) + MDP(25);
t162 = MDP(24) - MDP(27);
t134 = g(1) * t121 + g(2) * t117;
t115 = sin(qJ(3));
t119 = cos(qJ(3));
t139 = t121 * t119;
t102 = t115 * t143 + t139;
t140 = t121 * t115;
t104 = t117 * t119 - t120 * t140;
t158 = -g(1) * t104 + g(2) * t102 + t115 * t148;
t145 = t111 * t116;
t144 = t112 * t116;
t137 = pkin(3) * t115 + pkin(7);
t86 = g(1) * t99 + g(2) * t97 + g(3) * t145;
t130 = t163 * t86 + t162 * (g(1) * t100 + g(2) * t98 + g(3) * t144) - t166;
t110 = t119 * pkin(3) + pkin(2);
t122 = -pkin(9) - pkin(8);
t127 = t120 * t110 - t116 * t122 + pkin(1);
t126 = pkin(4) * t112 + qJ(5) * t111 + t110;
t123 = -g(1) * (-t99 * pkin(4) + t100 * qJ(5)) - g(2) * (-t97 * pkin(4) + t98 * qJ(5)) - g(3) * (-pkin(4) * t145 + qJ(5) * t144);
t105 = t117 * t115 + t120 * t139;
t103 = -t119 * t143 + t140;
t1 = [t134 * MDP(3) + (-g(1) * t103 - g(2) * t105) * MDP(16) + (-g(1) * t102 - g(2) * t104) * MDP(17) + (-g(1) * (-t98 * pkin(4) - t97 * qJ(5)) - g(2) * (t100 * pkin(4) + t99 * qJ(5)) + (-g(1) * t137 - g(2) * t127) * t121 + (g(1) * t127 - g(2) * t137) * t117) * MDP(28) + (g(1) * t132 - g(2) * t88) * MDP(34) + (-g(1) * t136 + g(2) * t131) * MDP(35) + t163 * (g(1) * t98 - g(2) * t100) - t162 * (g(1) * t97 - g(2) * t99) + (t120 * MDP(9) - t116 * t164 + MDP(2)) * (g(1) * t117 - g(2) * t121); ((-g(3) * t126 + t134 * t122) * t120 + (g(3) * t122 + t134 * t126) * t116) * MDP(28) + t164 * (t134 * t120 + t148) + (t119 * MDP(16) - t115 * MDP(17) + t128 * MDP(34) + t129 * MDP(35) - t111 * t162 + t112 * t163 + MDP(9)) * (-g(3) * t120 + t134 * t116); t158 * MDP(16) + (g(1) * t105 - g(2) * t103 + t119 * t148) * MDP(17) + (t158 * pkin(3) + t123) * MDP(28) + t130; t123 * MDP(28) + t130; -t86 * MDP(28); t166;];
taug  = t1;
