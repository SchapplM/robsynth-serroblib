% Calculate Gravitation load on the joints for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:41:44
% EndTime: 2019-03-09 14:41:45
% DurationCPUTime: 0.46s
% Computational Cost: add. (341->98), mult. (637->165), div. (0->0), fcn. (744->12), ass. (0->45)
t121 = MDP(9) - MDP(12);
t120 = -MDP(10) + MDP(13);
t107 = cos(pkin(6));
t95 = cos(qJ(1));
t104 = t95 * t107;
t90 = sin(qJ(2));
t91 = sin(qJ(1));
t94 = cos(qJ(2));
t78 = t90 * t104 + t91 * t94;
t88 = sin(qJ(6));
t92 = cos(qJ(6));
t87 = sin(pkin(6));
t110 = t87 * t95;
t77 = -t94 * t104 + t91 * t90;
t86 = qJ(4) + qJ(5);
t84 = sin(t86);
t85 = cos(t86);
t99 = t85 * t110 - t77 * t84;
t119 = t78 * t92 + t88 * t99;
t118 = -t78 * t88 + t92 * t99;
t117 = g(3) * t87;
t114 = t84 * t88;
t113 = t84 * t92;
t112 = t87 * t91;
t111 = t87 * t94;
t109 = t88 * t90;
t108 = t90 * t92;
t89 = sin(qJ(4));
t93 = cos(qJ(4));
t106 = -t93 * t110 + t77 * t89;
t105 = t91 * t107;
t79 = t94 * t105 + t95 * t90;
t67 = -t84 * t112 + t79 * t85;
t68 = t85 * t112 + t79 * t84;
t69 = t84 * t110 + t77 * t85;
t76 = t107 * t85 - t84 * t111;
t103 = (g(1) * t68 - g(2) * t99 + g(3) * t76) * MDP(28) + (-t92 * MDP(34) + t88 * MDP(35) - MDP(27)) * (g(1) * t67 + g(2) * t69 + g(3) * (-t107 * t84 - t85 * t111));
t98 = t89 * t110 + t77 * t93;
t66 = -g(1) * t79 - g(2) * t77 + g(3) * t111;
t80 = -t90 * t105 + t95 * t94;
t73 = t93 * t112 + t79 * t89;
t72 = -t89 * t112 + t79 * t93;
t65 = t68 * t92 + t80 * t88;
t64 = -t68 * t88 + t80 * t92;
t1 = [(g(1) * t91 - g(2) * t95) * MDP(2) + (-g(1) * (-t91 * pkin(1) - t78 * pkin(2) + pkin(8) * t110 - t77 * qJ(3)) - g(2) * (t95 * pkin(1) + t80 * pkin(2) + pkin(8) * t112 + t79 * qJ(3))) * MDP(14) + (g(1) * t106 - g(2) * t73) * MDP(20) + (g(1) * t98 - g(2) * t72) * MDP(21) + (-g(1) * t99 - g(2) * t68) * MDP(27) + (g(1) * t69 - g(2) * t67) * MDP(28) + (-g(1) * t118 - g(2) * t65) * MDP(34) + (g(1) * t119 - g(2) * t64) * MDP(35) + t120 * (g(1) * t77 - g(2) * t79) + t121 * (g(1) * t78 - g(2) * t80) + (-t87 * MDP(11) + MDP(3)) * (g(1) * t95 + g(2) * t91); (-g(1) * (-t79 * pkin(2) + t80 * qJ(3)) - g(2) * (-t77 * pkin(2) + t78 * qJ(3)) - (pkin(2) * t94 + qJ(3) * t90) * t117) * MDP(14) + (-g(1) * (t80 * t113 - t79 * t88) - g(2) * (t78 * t113 - t77 * t88) - (t84 * t108 + t88 * t94) * t117) * MDP(34) + (-g(1) * (-t80 * t114 - t79 * t92) - g(2) * (-t78 * t114 - t77 * t92) - (-t84 * t109 + t92 * t94) * t117) * MDP(35) - t121 * t66 + (-t89 * MDP(20) - t93 * MDP(21) - t84 * MDP(27) - MDP(28) * t85 - t120) * (g(1) * t80 + g(2) * t78 + t90 * t117); t66 * MDP(14); (-g(1) * t72 - g(2) * t98 - g(3) * (-t107 * t89 - t93 * t111)) * MDP(20) + (g(1) * t73 + g(2) * t106 - g(3) * (-t107 * t93 + t89 * t111)) * MDP(21) + t103; t103; (-g(1) * t64 - g(2) * t119 - g(3) * (t87 * t108 - t76 * t88)) * MDP(34) + (g(1) * t65 - g(2) * t118 - g(3) * (-t87 * t109 - t76 * t92)) * MDP(35);];
taug  = t1;
