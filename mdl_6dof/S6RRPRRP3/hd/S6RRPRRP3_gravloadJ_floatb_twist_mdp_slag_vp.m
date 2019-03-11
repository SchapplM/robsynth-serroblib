% Calculate Gravitation load on the joints for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:50:35
% EndTime: 2019-03-09 11:50:36
% DurationCPUTime: 0.34s
% Computational Cost: add. (264->78), mult. (299->111), div. (0->0), fcn. (271->10), ass. (0->46)
t89 = sin(qJ(1));
t92 = cos(qJ(1));
t72 = g(1) * t92 + g(2) * t89;
t84 = qJ(2) + pkin(10);
t76 = sin(t84);
t111 = g(3) * t76;
t85 = qJ(4) + qJ(5);
t79 = cos(t85);
t103 = t92 * t79;
t78 = sin(t85);
t108 = t89 * t78;
t77 = cos(t84);
t60 = t77 * t108 + t103;
t104 = t92 * t78;
t107 = t89 * t79;
t62 = -t77 * t104 + t107;
t116 = -g(1) * t62 + g(2) * t60 + t78 * t111;
t95 = -g(3) * t77 + t72 * t76;
t87 = sin(qJ(4));
t69 = t87 * pkin(4) + pkin(5) * t78;
t109 = t69 * t77;
t106 = t89 * t87;
t90 = cos(qJ(4));
t105 = t89 * t90;
t102 = t92 * t87;
t101 = t92 * t90;
t61 = -t77 * t107 + t104;
t63 = t77 * t103 + t108;
t100 = t116 * MDP(25) + (g(1) * t63 - g(2) * t61 + t79 * t111) * MDP(26);
t86 = -qJ(3) - pkin(7);
t99 = t69 - t86;
t70 = t90 * pkin(4) + pkin(5) * t79;
t71 = g(1) * t89 - g(2) * t92;
t68 = pkin(3) + t70;
t83 = -qJ(6) - pkin(9) - pkin(8);
t97 = t77 * t68 - t76 * t83;
t91 = cos(qJ(2));
t88 = sin(qJ(2));
t81 = t91 * pkin(2);
t75 = t81 + pkin(1);
t73 = t92 * t75;
t67 = t77 * t101 + t106;
t66 = -t77 * t102 + t105;
t65 = -t77 * t105 + t102;
t64 = t77 * t106 + t101;
t1 = [(-g(1) * (-t89 * t75 - t92 * t86) - g(2) * (-t89 * t86 + t73)) * MDP(12) + (-g(1) * t65 - g(2) * t67) * MDP(18) + (-g(1) * t64 - g(2) * t66) * MDP(19) + (-g(1) * t61 - g(2) * t63) * MDP(25) + (-g(1) * t60 - g(2) * t62) * MDP(26) + (-g(2) * t73 + (-g(1) * t99 - g(2) * t97) * t92 + (-g(1) * (-t75 - t97) - g(2) * t99) * t89) * MDP(28) + (MDP(3) - MDP(11)) * t72 + (-t88 * MDP(10) + t76 * MDP(27) + t91 * MDP(9) + MDP(2)) * t71; (g(3) * t88 + t72 * t91) * MDP(10) + (-t72 * t77 - t111) * MDP(27) + (-g(3) * (t81 + t97) + t72 * (pkin(2) * t88 + t68 * t76 + t77 * t83)) * MDP(28) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t91 + t72 * t88) + (t90 * MDP(18) - t87 * MDP(19) + t79 * MDP(25) - t78 * MDP(26)) * t95; (-MDP(12) - MDP(28)) * t71; (-g(1) * t66 + g(2) * t64 + t87 * t111) * MDP(18) + (g(1) * t67 - g(2) * t65 + t90 * t111) * MDP(19) + (-g(1) * (-t92 * t109 + t89 * t70) - g(2) * (-t89 * t109 - t92 * t70) + t69 * t111) * MDP(28) + t100; t116 * MDP(28) * pkin(5) + t100; -t95 * MDP(28);];
taug  = t1;
