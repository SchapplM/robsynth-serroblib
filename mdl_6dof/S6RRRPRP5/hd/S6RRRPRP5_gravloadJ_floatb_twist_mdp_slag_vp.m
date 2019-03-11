% Calculate Gravitation load on the joints for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:51:44
% EndTime: 2019-03-09 16:51:45
% DurationCPUTime: 0.60s
% Computational Cost: add. (460->100), mult. (487->137), div. (0->0), fcn. (472->10), ass. (0->49)
t148 = MDP(25) + MDP(27);
t149 = MDP(26) - MDP(29);
t115 = sin(qJ(1));
t118 = cos(qJ(1));
t124 = g(1) * t118 + g(2) * t115;
t145 = MDP(10) - MDP(18) - MDP(28);
t114 = sin(qJ(2));
t117 = cos(qJ(2));
t87 = -g(3) * t117 + t124 * t114;
t140 = g(3) * t114;
t113 = sin(qJ(3));
t138 = t113 * pkin(3);
t112 = -qJ(4) - pkin(8);
t111 = qJ(3) + pkin(10);
t104 = qJ(5) + t111;
t101 = sin(t104);
t137 = t101 * t114;
t102 = cos(t104);
t136 = t102 * t114;
t135 = t113 * t118;
t134 = t114 * t118;
t133 = t115 * t117;
t132 = t117 * t118;
t131 = t118 * t101;
t116 = cos(qJ(3));
t106 = t116 * pkin(3);
t96 = pkin(4) * cos(t111) + t106;
t130 = t118 * pkin(1) + t115 * pkin(7);
t83 = t101 * t133 + t102 * t118;
t84 = t102 * t133 - t131;
t128 = -t83 * pkin(5) + t84 * qJ(6);
t85 = -t115 * t102 + t117 * t131;
t86 = t115 * t101 + t102 * t132;
t127 = -t85 * pkin(5) + t86 * qJ(6);
t77 = g(1) * t85 + g(2) * t83 + g(3) * t137;
t126 = t148 * t77 + t149 * (g(1) * t86 + g(2) * t84 + g(3) * t136);
t103 = t106 + pkin(2);
t122 = t103 * t117 - t112 * t114;
t94 = pkin(2) + t96;
t120 = pkin(5) * t102 + qJ(6) * t101 + t94;
t91 = -t113 * t132 + t115 * t116;
t89 = t113 * t133 + t116 * t118;
t110 = -pkin(9) + t112;
t107 = t118 * pkin(7);
t97 = qJ(6) * t136;
t95 = pkin(4) * sin(t111) + t138;
t92 = t115 * t113 + t116 * t132;
t90 = -t116 * t133 + t135;
t1 = [t124 * MDP(3) + (-g(1) * t90 - g(2) * t92) * MDP(16) + (-g(1) * t89 - g(2) * t91) * MDP(17) + (-g(1) * (pkin(3) * t135 + t107) - g(2) * (t103 * t132 - t112 * t134 + t130) + (-g(1) * (-pkin(1) - t122) - g(2) * t138) * t115) * MDP(19) + (-g(1) * (-t84 * pkin(5) - t83 * qJ(6) + t118 * t95 + t107) - g(2) * (t86 * pkin(5) + t85 * qJ(6) - t110 * t134 + t94 * t132 + t130) + (-g(1) * (t110 * t114 - t117 * t94 - pkin(1)) - g(2) * t95) * t115) * MDP(30) + t148 * (g(1) * t84 - g(2) * t86) - t149 * (g(1) * t83 - g(2) * t85) + (t117 * MDP(9) - t145 * t114 + MDP(2)) * (g(1) * t115 - g(2) * t118); (-g(3) * t122 + t124 * (t103 * t114 + t112 * t117)) * MDP(19) + ((-g(3) * t120 + t124 * t110) * t117 + (g(3) * t110 + t124 * t120) * t114) * MDP(30) + t145 * (t124 * t117 + t140) + (t116 * MDP(16) - t113 * MDP(17) - t101 * t149 + t148 * t102 + MDP(9)) * t87; (g(1) * t92 - g(2) * t90 + t116 * t140) * MDP(17) + (-g(1) * (t115 * t96 - t95 * t132 + t127) - g(2) * (-t118 * t96 - t95 * t133 + t128) - g(3) * (t97 + (-pkin(5) * t101 - t95) * t114)) * MDP(30) + t126 + (pkin(3) * MDP(19) + MDP(16)) * (-g(1) * t91 + g(2) * t89 + t113 * t140); (-MDP(19) - MDP(30)) * t87; (-g(1) * t127 - g(2) * t128 - g(3) * (-pkin(5) * t137 + t97)) * MDP(30) + t126; -t77 * MDP(30);];
taug  = t1;
