% Calculate Gravitation load on the joints for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:34:16
% EndTime: 2019-03-09 19:34:18
% DurationCPUTime: 0.96s
% Computational Cost: add. (445->113), mult. (1148->185), div. (0->0), fcn. (1415->12), ass. (0->52)
t121 = sin(qJ(2));
t122 = sin(qJ(1));
t126 = cos(qJ(2));
t152 = cos(pkin(6));
t155 = cos(qJ(1));
t139 = t152 * t155;
t108 = t121 * t122 - t126 * t139;
t118 = sin(qJ(6));
t123 = cos(qJ(6));
t109 = t121 * t139 + t122 * t126;
t120 = sin(qJ(3));
t125 = cos(qJ(3));
t117 = sin(pkin(6));
t141 = t117 * t155;
t100 = t109 * t125 - t120 * t141;
t119 = sin(qJ(5));
t124 = cos(qJ(5));
t99 = t109 * t120 + t125 * t141;
t86 = t100 * t124 + t119 * t99;
t167 = t108 * t123 + t118 * t86;
t166 = -t108 * t118 + t123 * t86;
t165 = t100 * t119 - t124 * t99;
t149 = t117 * t121;
t106 = t120 * t149 - t125 * t152;
t147 = t117 * t125;
t107 = t120 * t152 + t121 * t147;
t140 = t122 * t152;
t111 = -t121 * t140 + t126 * t155;
t103 = t111 * t120 - t122 * t147;
t148 = t117 * t122;
t104 = t111 * t125 + t120 * t148;
t89 = t103 * t124 - t104 * t119;
t90 = t103 * t119 + t104 * t124;
t98 = t106 * t119 + t107 * t124;
t164 = (g(1) * t90 + g(2) * t86 + g(3) * t98) * MDP(28) + (-MDP(34) * t123 + MDP(35) * t118 - MDP(27)) * (g(1) * t89 - g(2) * t165 + g(3) * (t106 * t124 - t107 * t119));
t163 = MDP(10) - MDP(19);
t143 = MDP(16) + MDP(18);
t162 = MDP(17) - MDP(20);
t110 = t121 * t155 + t126 * t140;
t154 = g(1) * t110;
t153 = g(2) * t108;
t146 = t117 * t126;
t136 = -g(1) * t111 - g(2) * t109;
t134 = t119 * t120 + t124 * t125;
t133 = pkin(3) * t125 + qJ(4) * t120 + pkin(2);
t130 = g(1) * t103 + g(2) * t99 + g(3) * t106;
t105 = t134 * t146;
t95 = t134 * t110;
t94 = t134 * t108;
t82 = -t110 * t118 + t123 * t90;
t81 = -t110 * t123 - t118 * t90;
t1 = [(g(1) * t122 - g(2) * t155) * MDP(2) + (g(1) * t155 + g(2) * t122) * MDP(3) + (g(1) * t109 - g(2) * t111) * MDP(9) + (-g(1) * (-pkin(1) * t122 - pkin(2) * t109 - pkin(3) * t100 + pkin(8) * t141 - pkin(9) * t108 - qJ(4) * t99) - g(2) * (pkin(1) * t155 + pkin(2) * t111 + pkin(3) * t104 + pkin(8) * t148 + pkin(9) * t110 + qJ(4) * t103)) * MDP(21) + (g(1) * t86 - g(2) * t90) * MDP(27) + (-g(1) * t165 - g(2) * t89) * MDP(28) + (g(1) * t166 - g(2) * t82) * MDP(34) + (-g(1) * t167 - g(2) * t81) * MDP(35) + t143 * (g(1) * t100 - g(2) * t104) + t162 * (-g(1) * t99 + g(2) * t103) - t163 * (g(1) * t108 - g(2) * t110); (t136 * pkin(9) + t133 * t154 + t133 * t153 - g(3) * (pkin(9) * t121 + t126 * t133) * t117) * MDP(21) + (g(1) * t95 + g(2) * t94 - g(3) * t105) * MDP(27) + (-g(1) * (-t111 * t118 - t123 * t95) - g(2) * (-t109 * t118 - t123 * t94) - g(3) * (t105 * t123 - t118 * t149)) * MDP(34) + (-g(1) * (-t111 * t123 + t118 * t95) - g(2) * (-t109 * t123 + t118 * t94) - g(3) * (-t105 * t118 - t123 * t149)) * MDP(35) + t163 * (g(3) * t149 - t136) + (-MDP(9) - t143 * t125 + t162 * t120 + (t119 * t125 - t120 * t124) * MDP(28)) * (g(3) * t146 - t153 - t154); (-g(1) * (-pkin(3) * t103 + qJ(4) * t104) - g(2) * (-pkin(3) * t99 + qJ(4) * t100) - g(3) * (-pkin(3) * t106 + qJ(4) * t107)) * MDP(21) + t143 * t130 + t162 * (g(1) * t104 + g(2) * t100 + g(3) * t107) - t164; -t130 * MDP(21); t164; (-g(1) * t81 + g(2) * t167 - g(3) * (-t118 * t98 + t123 * t146)) * MDP(34) + (g(1) * t82 + g(2) * t166 - g(3) * (-t118 * t146 - t123 * t98)) * MDP(35);];
taug  = t1;
