% Calculate Gravitation load on the joints for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:36:01
% EndTime: 2019-03-09 15:36:03
% DurationCPUTime: 0.67s
% Computational Cost: add. (352->100), mult. (523->146), div. (0->0), fcn. (540->10), ass. (0->50)
t112 = qJ(3) + pkin(10);
t106 = sin(t112);
t107 = cos(t112);
t114 = sin(qJ(6));
t118 = cos(qJ(6));
t126 = t106 * t114 + t107 * t118;
t127 = t106 * t118 - t107 * t114;
t117 = sin(qJ(1));
t120 = cos(qJ(2));
t121 = cos(qJ(1));
t139 = t121 * t106;
t88 = -t117 * t107 + t120 * t139;
t138 = t121 * t107;
t89 = t117 * t106 + t120 * t138;
t129 = t89 * t114 - t88 * t118;
t140 = t117 * t120;
t86 = t106 * t140 + t138;
t87 = t107 * t140 - t139;
t130 = t86 * t114 + t87 * t118;
t134 = t87 * t114 - t86 * t118;
t116 = sin(qJ(2));
t146 = g(3) * t116;
t82 = t88 * t114 + t89 * t118;
t161 = (g(1) * t129 + g(2) * t134 - t127 * t146) * MDP(29) + (g(1) * t82 + g(2) * t130 + t126 * t146) * MDP(30);
t113 = -qJ(4) - pkin(8);
t119 = cos(qJ(3));
t105 = t119 * pkin(3) + pkin(2);
t98 = t120 * t105;
t160 = -t116 * t113 + t98;
t133 = g(1) * t121 + g(2) * t117;
t90 = -g(3) * t120 + t116 * t133;
t158 = -pkin(1) - t160;
t157 = MDP(10) - MDP(18) - MDP(21);
t115 = sin(qJ(3));
t142 = t117 * t115;
t141 = t117 * t119;
t137 = t121 * t115;
t136 = t121 * t119;
t135 = t120 * t137;
t131 = pkin(4) * t107 + qJ(5) * t106;
t92 = t115 * t140 + t136;
t124 = pkin(3) * t142 + t117 * pkin(7) - t158 * t121;
t123 = g(1) * t88 + g(2) * t86 + t106 * t146;
t122 = pkin(3) * t137 + t121 * pkin(7) + t158 * t117;
t91 = t120 * t133 + t146;
t103 = pkin(3) * t141;
t95 = t120 * t136 + t142;
t94 = -t135 + t141;
t93 = -t119 * t140 + t137;
t1 = [t133 * MDP(3) + (-g(1) * t93 - g(2) * t95) * MDP(16) + (-g(1) * t92 - g(2) * t94) * MDP(17) + (-g(1) * t122 - g(2) * t124) * MDP(19) + (g(1) * t87 - g(2) * t89) * MDP(20) + (g(1) * t86 - g(2) * t88) * MDP(22) + (-g(1) * (-t87 * pkin(4) - t86 * qJ(5) + t122) - g(2) * (t89 * pkin(4) + t88 * qJ(5) + t124)) * MDP(23) + (g(1) * t130 - g(2) * t82) * MDP(29) + (-g(1) * t134 + g(2) * t129) * MDP(30) + (t120 * MDP(9) - t157 * t116 + MDP(2)) * (g(1) * t117 - g(2) * t121); (-g(3) * t160 + t133 * (t105 * t116 + t113 * t120)) * MDP(19) + (-g(3) * t98 + (-g(3) * t131 + t113 * t133) * t120 + (g(3) * t113 + t133 * (t105 + t131)) * t116) * MDP(23) + t157 * t91 + (t119 * MDP(16) - t115 * MDP(17) + t107 * MDP(20) + t106 * MDP(22) + t126 * MDP(29) + t127 * MDP(30) + MDP(9)) * t90; (-g(1) * t94 + g(2) * t92 + t115 * t146) * MDP(16) + (g(1) * t95 - g(2) * t93 + t119 * t146) * MDP(17) + (-g(1) * t103 + (g(2) * t136 + t115 * t91) * pkin(3)) * MDP(19) + t123 * MDP(20) + (-g(1) * t89 - g(2) * t87 - t107 * t146) * MDP(22) + (-g(1) * (-pkin(3) * t135 - t88 * pkin(4) + t89 * qJ(5) + t103) - g(2) * (-pkin(3) * t92 - t86 * pkin(4) + t87 * qJ(5)) - (-pkin(3) * t115 - pkin(4) * t106 + qJ(5) * t107) * t146) * MDP(23) - t161; (-MDP(19) - MDP(23)) * t90; -t123 * MDP(23); t161;];
taug  = t1;
