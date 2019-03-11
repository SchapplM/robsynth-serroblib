% Calculate Gravitation load on the joints for
% S6RPRRRP1
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
%   see S6RPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:57:07
% EndTime: 2019-03-09 05:57:08
% DurationCPUTime: 0.35s
% Computational Cost: add. (389->73), mult. (373->103), div. (0->0), fcn. (347->10), ass. (0->38)
t118 = MDP(18) - MDP(27);
t115 = MDP(24) + MDP(26);
t114 = MDP(25) - MDP(28);
t87 = qJ(3) + qJ(4);
t83 = sin(t87);
t84 = cos(t87);
t117 = t84 * pkin(4) + t83 * pkin(9);
t86 = qJ(1) + pkin(10);
t81 = sin(t86);
t82 = cos(t86);
t103 = g(1) * t82 + g(2) * t81;
t116 = t103 * t83;
t89 = sin(qJ(3));
t113 = pkin(3) * t89;
t112 = pkin(9) * t84;
t109 = g(3) * t83;
t88 = sin(qJ(5));
t108 = t84 * t88;
t91 = cos(qJ(5));
t107 = t84 * t91;
t106 = qJ(6) * t88;
t105 = pkin(5) * t107 + t84 * t106 + t117;
t92 = cos(qJ(3));
t85 = t92 * pkin(3);
t100 = t85 + pkin(2) + t117;
t98 = t118 * (t103 * t84 + t109) + (-t114 * t88 + t115 * t91 + MDP(17)) * (-g(3) * t84 + t116);
t68 = t81 * t108 + t82 * t91;
t70 = t82 * t108 - t81 * t91;
t55 = g(1) * t70 + g(2) * t68 + t88 * t109;
t95 = (pkin(5) * t91 + pkin(4) + t106) * t116;
t94 = -pkin(8) - pkin(7);
t93 = cos(qJ(1));
t90 = sin(qJ(1));
t73 = t82 * t112;
t72 = t81 * t112;
t71 = t82 * t107 + t81 * t88;
t69 = t81 * t107 - t82 * t88;
t1 = [(g(1) * t93 + g(2) * t90) * MDP(3) + (-g(1) * (-t90 * pkin(1) - t69 * pkin(5) - t68 * qJ(6)) - g(2) * (t93 * pkin(1) + t71 * pkin(5) + t70 * qJ(6)) + (g(1) * t94 - g(2) * t100) * t82 + (g(1) * t100 + g(2) * t94) * t81) * MDP(29) + t115 * (g(1) * t69 - g(2) * t71) - t114 * (g(1) * t68 - g(2) * t70) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t90 - g(2) * t93) + (t92 * MDP(10) - t89 * MDP(11) + t84 * MDP(17) - t118 * t83) * (g(1) * t81 - g(2) * t82); (-MDP(29) - MDP(4)) * g(3); (-g(3) * t92 + t103 * t89) * MDP(10) + (g(3) * t89 + t103 * t92) * MDP(11) + (-g(1) * (-t82 * t113 + t73) - g(2) * (-t81 * t113 + t72) - g(3) * (t85 + t105) + t95) * MDP(29) + t98; (-g(1) * t73 - g(2) * t72 - g(3) * t105 + t95) * MDP(29) + t98; (-g(1) * (-t70 * pkin(5) + t71 * qJ(6)) - g(2) * (-t68 * pkin(5) + t69 * qJ(6)) - (-pkin(5) * t88 + qJ(6) * t91) * t109) * MDP(29) + t114 * (g(1) * t71 + g(2) * t69 + t91 * t109) + t115 * t55; -t55 * MDP(29);];
taug  = t1;
