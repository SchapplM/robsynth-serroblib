% Calculate Gravitation load on the joints for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:31:28
% EndTime: 2019-03-10 01:31:30
% DurationCPUTime: 0.52s
% Computational Cost: add. (586->99), mult. (570->141), div. (0->0), fcn. (578->10), ass. (0->50)
t148 = MDP(30) + MDP(32);
t150 = MDP(31) - MDP(34);
t116 = sin(qJ(1));
t119 = cos(qJ(1));
t128 = g(1) * t119 + g(2) * t116;
t149 = MDP(10) - MDP(33);
t113 = qJ(3) + qJ(4);
t108 = sin(t113);
t115 = sin(qJ(2));
t141 = g(3) * t115;
t109 = cos(t113);
t118 = cos(qJ(2));
t136 = t116 * t118;
t89 = t108 * t136 + t109 * t119;
t135 = t118 * t119;
t91 = -t108 * t135 + t109 * t116;
t146 = -g(1) * t91 + g(2) * t89 + t108 * t141;
t114 = sin(qJ(3));
t100 = pkin(3) * t114 + pkin(4) * t108;
t139 = pkin(7) + t100;
t110 = qJ(5) + t113;
t106 = sin(t110);
t138 = t106 * t115;
t107 = cos(t110);
t137 = t107 * t115;
t134 = t119 * t106;
t117 = cos(qJ(3));
t101 = pkin(3) * t117 + pkin(4) * t109;
t85 = t106 * t136 + t107 * t119;
t86 = t107 * t136 - t134;
t132 = -pkin(5) * t85 + t86 * qJ(6);
t87 = -t107 * t116 + t118 * t134;
t88 = t106 * t116 + t107 * t135;
t131 = -pkin(5) * t87 + t88 * qJ(6);
t77 = g(1) * t87 + g(2) * t85 + g(3) * t138;
t130 = t148 * t77 + t150 * (g(1) * t88 + g(2) * t86 + g(3) * t137);
t90 = t108 * t119 - t109 * t136;
t92 = t108 * t116 + t109 * t135;
t126 = t146 * MDP(23) + (g(1) * t92 - g(2) * t90 + t109 * t141) * MDP(24) + t130;
t112 = -pkin(10) - pkin(9) - pkin(8);
t99 = pkin(2) + t101;
t125 = t112 * t115 - t118 * t99 - pkin(1);
t124 = pkin(5) * t107 + qJ(6) * t106 + t99;
t102 = qJ(6) * t137;
t120 = -g(1) * t131 - g(2) * t132 - g(3) * (-pkin(5) * t138 + t102);
t97 = t114 * t116 + t117 * t135;
t96 = -t114 * t135 + t116 * t117;
t95 = t114 * t119 - t117 * t136;
t94 = t114 * t136 + t117 * t119;
t1 = [t128 * MDP(3) + (-g(1) * t95 - g(2) * t97) * MDP(16) + (-g(1) * t94 - g(2) * t96) * MDP(17) + (-g(1) * t90 - g(2) * t92) * MDP(23) + (-g(1) * t89 - g(2) * t91) * MDP(24) + (-g(1) * (-t86 * pkin(5) - t85 * qJ(6)) - g(2) * (t88 * pkin(5) + t87 * qJ(6)) + (-g(1) * t139 + g(2) * t125) * t119 + (-g(1) * t125 - g(2) * t139) * t116) * MDP(35) + t148 * (g(1) * t86 - g(2) * t88) - t150 * (g(1) * t85 - g(2) * t87) + (t118 * MDP(9) - t115 * t149 + MDP(2)) * (g(1) * t116 - g(2) * t119); ((-g(3) * t124 + t112 * t128) * t118 + (g(3) * t112 + t124 * t128) * t115) * MDP(35) + t149 * (t118 * t128 + t141) + (t117 * MDP(16) - MDP(17) * t114 + MDP(23) * t109 - MDP(24) * t108 - t106 * t150 + t107 * t148 + MDP(9)) * (-g(3) * t118 + t115 * t128); (-g(1) * t96 + g(2) * t94 + t114 * t141) * MDP(16) + (g(1) * t97 - g(2) * t95 + t117 * t141) * MDP(17) + (-g(1) * (-t100 * t135 + t101 * t116 + t131) - g(2) * (-t100 * t136 - t101 * t119 + t132) - g(3) * (t102 + (-pkin(5) * t106 - t100) * t115)) * MDP(35) + t126; (pkin(4) * t146 + t120) * MDP(35) + t126; MDP(35) * t120 + t130; -t77 * MDP(35);];
taug  = t1;
