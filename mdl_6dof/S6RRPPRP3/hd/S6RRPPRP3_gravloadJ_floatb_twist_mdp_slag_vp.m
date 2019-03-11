% Calculate Gravitation load on the joints for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:34
% EndTime: 2019-03-09 08:35:35
% DurationCPUTime: 0.50s
% Computational Cost: add. (160->79), mult. (329->102), div. (0->0), fcn. (282->6), ass. (0->49)
t87 = sin(qJ(2));
t78 = t87 * qJ(3);
t124 = pkin(1) + t78;
t123 = -MDP(10) + MDP(13) + MDP(15);
t122 = MDP(9) + MDP(11) - MDP(16) + MDP(26);
t88 = sin(qJ(1));
t113 = g(2) * t88;
t91 = cos(qJ(1));
t71 = g(1) * t91 + t113;
t121 = t71 * t87;
t90 = cos(qJ(2));
t63 = g(3) * t87 + t71 * t90;
t119 = pkin(2) + pkin(3);
t82 = t91 * pkin(7);
t117 = g(1) * t82;
t116 = g(1) * t88;
t111 = g(3) * t90;
t81 = t90 * pkin(2);
t80 = t90 * pkin(3);
t89 = cos(qJ(5));
t77 = pkin(5) * t89 + pkin(4);
t110 = t77 * t87;
t85 = -qJ(6) - pkin(8);
t109 = t85 * t90;
t108 = t87 * t91;
t86 = sin(qJ(5));
t107 = t88 * t86;
t106 = t88 * t89;
t105 = t89 * t91;
t104 = t90 * t91;
t103 = t81 + t78;
t102 = qJ(3) * t90;
t101 = MDP(18) + MDP(27);
t100 = -t85 + t119;
t98 = t80 + t103;
t96 = pkin(5) * t86 + qJ(4);
t95 = pkin(2) * t104 + t88 * pkin(7) + t124 * t91;
t94 = g(1) * t100;
t70 = -g(2) * t91 + t116;
t93 = -t124 - t81;
t66 = -t86 * t108 - t106;
t64 = t87 * t107 - t105;
t92 = g(2) * (pkin(3) * t104 + t95);
t74 = t91 * t102;
t72 = t88 * t102;
t67 = t87 * t105 - t107;
t65 = -t87 * t106 - t86 * t91;
t62 = -t111 + t121;
t1 = [(-g(2) * t95 - t93 * t116 - t117) * MDP(14) + (-g(1) * (-qJ(4) * t91 + t82) - t92 + (-g(1) * (t93 - t80) + g(2) * qJ(4)) * t88) * MDP(18) + (-g(1) * t65 - g(2) * t67) * MDP(24) + (-g(1) * t64 - g(2) * t66) * MDP(25) + (-t117 - t92 + (g(1) * t96 - g(2) * (-t109 + t110)) * t91 + (-g(1) * (-t124 - t110) + g(2) * t96 + t90 * t94) * t88) * MDP(27) + (MDP(3) - MDP(12) + MDP(17)) * t71 + (t122 * t90 + t123 * t87 + MDP(2)) * t70; (-g(1) * (-pkin(2) * t108 + t74) - g(2) * (-pkin(2) * t87 * t88 + t72) - g(3) * t103) * MDP(14) + (-g(1) * t74 - g(2) * t72 - g(3) * t98 + t119 * t121) * MDP(18) + (-g(1) * (t77 * t104 + t74) - g(2) * (t77 * t88 * t90 + t72) - g(3) * (t98 - t109) + (-g(3) * t77 + t100 * t113 + t91 * t94) * t87) * MDP(27) + t122 * t62 + (-t89 * MDP(24) + MDP(25) * t86 - t123) * t63; (-MDP(14) - t101) * t62; t101 * t70; (g(1) * t67 - g(2) * t65 - t89 * t111) * MDP(25) + (pkin(5) * MDP(27) + MDP(24)) * (-g(1) * t66 + g(2) * t64 - t86 * t111); -t63 * MDP(27);];
taug  = t1;
