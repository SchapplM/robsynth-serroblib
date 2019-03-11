% Calculate Gravitation load on the joints for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:31
% EndTime: 2019-03-09 06:28:32
% DurationCPUTime: 0.32s
% Computational Cost: add. (192->71), mult. (286->101), div. (0->0), fcn. (261->8), ass. (0->41)
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t115 = -g(1) * t87 + g(2) * t90;
t114 = MDP(13) - MDP(28);
t89 = cos(qJ(3));
t107 = g(3) * t89;
t84 = qJ(4) + qJ(5);
t76 = sin(t84);
t104 = t87 * t76;
t86 = sin(qJ(3));
t77 = cos(t84);
t99 = t90 * t77;
t59 = -t86 * t104 + t99;
t100 = t90 * t76;
t103 = t87 * t77;
t61 = t86 * t100 + t103;
t113 = -g(1) * t59 - g(2) * t61 + t76 * t107;
t64 = -g(3) * t86 - t115 * t89;
t85 = sin(qJ(4));
t71 = t85 * pkin(4) + pkin(5) * t76;
t106 = pkin(7) + t71;
t105 = t71 * t86;
t102 = t87 * t85;
t88 = cos(qJ(4));
t101 = t87 * t88;
t98 = t90 * t85;
t97 = t90 * t88;
t60 = t86 * t103 + t100;
t62 = t86 * t99 - t104;
t96 = t113 * MDP(26) + (g(1) * t60 - g(2) * t62 + t77 * t107) * MDP(27);
t72 = t88 * pkin(4) + pkin(5) * t77;
t94 = g(2) * (t90 * pkin(1) + t87 * qJ(2));
t70 = pkin(3) + t72;
t83 = -qJ(6) - pkin(9) - pkin(8);
t92 = t86 * t70 + t89 * t83;
t79 = t90 * qJ(2);
t68 = t86 * t97 - t102;
t67 = t86 * t98 + t101;
t66 = t86 * t101 + t98;
t65 = -t86 * t102 + t97;
t1 = [(-g(1) * (-t87 * pkin(1) + t79) - t94) * MDP(6) + (-g(1) * t68 - g(2) * t66) * MDP(19) + (g(1) * t67 - g(2) * t65) * MDP(20) + (-g(1) * t62 - g(2) * t60) * MDP(26) + (g(1) * t61 - g(2) * t59) * MDP(27) + (-g(1) * t79 - t94 + (-g(1) * t92 - g(2) * t106) * t90 + (-g(1) * (-pkin(1) - t106) - g(2) * t92) * t87) * MDP(29) - (MDP(2) - MDP(4)) * t115 + (-t86 * MDP(12) - t114 * t89 + MDP(3) - MDP(5)) * (g(1) * t90 + g(2) * t87); -(-MDP(29) - MDP(6)) * t115; (g(3) * t92 + t115 * (t70 * t89 - t83 * t86)) * MDP(29) + t114 * (-t115 * t86 + t107) + (-MDP(19) * t88 + MDP(20) * t85 - MDP(26) * t77 + MDP(27) * t76 - MDP(12)) * t64; (-g(1) * t65 - g(2) * t67 + t85 * t107) * MDP(19) + (g(1) * t66 - g(2) * t68 + t88 * t107) * MDP(20) + (-g(1) * (-t87 * t105 + t90 * t72) - g(2) * (t90 * t105 + t87 * t72) + t71 * t107) * MDP(29) + t96; t113 * MDP(29) * pkin(5) + t96; t64 * MDP(29);];
taug  = t1;
