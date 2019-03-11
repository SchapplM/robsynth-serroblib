% Calculate Gravitation load on the joints for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:27
% EndTime: 2019-03-09 06:05:28
% DurationCPUTime: 0.40s
% Computational Cost: add. (426->78), mult. (407->115), div. (0->0), fcn. (402->10), ass. (0->42)
t127 = MDP(11) - MDP(27);
t126 = MDP(24) + MDP(26);
t125 = MDP(25) - MDP(28);
t89 = qJ(1) + pkin(10);
t85 = sin(t89);
t86 = cos(t89);
t106 = g(1) * t86 + g(2) * t85;
t92 = sin(qJ(3));
t119 = g(3) * t92;
t91 = sin(qJ(4));
t95 = cos(qJ(3));
t113 = t91 * t95;
t94 = cos(qJ(4));
t76 = t85 * t113 + t86 * t94;
t78 = -t86 * t113 + t85 * t94;
t124 = -g(1) * t78 + g(2) * t76 + t91 * t119;
t90 = qJ(4) + qJ(5);
t87 = sin(t90);
t117 = t87 * t92;
t116 = t87 * t95;
t88 = cos(t90);
t115 = t88 * t92;
t114 = t88 * t95;
t97 = -pkin(9) - pkin(8);
t112 = t92 * t97;
t111 = t94 * t95;
t109 = pkin(4) * t91 + pkin(7);
t71 = t85 * t116 + t86 * t88;
t73 = t86 * t116 - t85 * t88;
t65 = g(1) * t73 + g(2) * t71 + g(3) * t117;
t72 = t85 * t114 - t86 * t87;
t74 = t86 * t114 + t85 * t87;
t108 = t126 * t65 + t125 * (g(1) * t74 + g(2) * t72 + g(3) * t115);
t84 = t94 * pkin(4) + pkin(3);
t103 = t95 * t84 + pkin(2) - t112;
t102 = pkin(5) * t88 + qJ(6) * t87 + t84;
t98 = -g(1) * (-t73 * pkin(5) + t74 * qJ(6)) - g(2) * (-t71 * pkin(5) + t72 * qJ(6)) - g(3) * (-pkin(5) * t117 + qJ(6) * t115);
t96 = cos(qJ(1));
t93 = sin(qJ(1));
t79 = t86 * t111 + t85 * t91;
t77 = -t85 * t111 + t86 * t91;
t1 = [(g(1) * t96 + g(2) * t93) * MDP(3) + (-g(1) * t77 - g(2) * t79) * MDP(17) + (-g(1) * t76 - g(2) * t78) * MDP(18) + (-g(1) * (-t93 * pkin(1) - t72 * pkin(5) - t71 * qJ(6)) - g(2) * (t96 * pkin(1) + t74 * pkin(5) + t73 * qJ(6)) + (-g(1) * t109 - g(2) * t103) * t86 + (g(1) * t103 - g(2) * t109) * t85) * MDP(29) + t126 * (g(1) * t72 - g(2) * t74) - t125 * (g(1) * t71 - g(2) * t73) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t93 - g(2) * t96) + (t95 * MDP(10) - t127 * t92) * (g(1) * t85 - g(2) * t86); (-MDP(29) - MDP(4)) * g(3); (-g(3) * (t102 * t95 - t112) + t106 * (t102 * t92 + t95 * t97)) * MDP(29) + t127 * (t106 * t95 + t119) + (t94 * MDP(17) - t91 * MDP(18) - t125 * t87 + t126 * t88 + MDP(10)) * (-g(3) * t95 + t106 * t92); t124 * MDP(17) + (g(1) * t79 - g(2) * t77 + t94 * t119) * MDP(18) + (t124 * pkin(4) + t98) * MDP(29) + t108; t98 * MDP(29) + t108; -t65 * MDP(29);];
taug  = t1;
