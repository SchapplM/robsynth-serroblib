% Calculate Gravitation load on the joints for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:03:23
% EndTime: 2019-03-08 23:03:24
% DurationCPUTime: 0.44s
% Computational Cost: add. (335->90), mult. (508->158), div. (0->0), fcn. (580->14), ass. (0->43)
t83 = sin(pkin(6));
t108 = g(3) * t83;
t81 = qJ(3) + qJ(4);
t76 = pkin(12) + t81;
t74 = cos(t76);
t86 = sin(qJ(6));
t107 = t74 * t86;
t89 = cos(qJ(6));
t106 = t74 * t89;
t82 = sin(pkin(11));
t105 = t82 * t83;
t84 = cos(pkin(11));
t104 = t83 * t84;
t87 = sin(qJ(3));
t103 = t83 * t87;
t88 = sin(qJ(2));
t102 = t83 * t88;
t90 = cos(qJ(3));
t101 = t83 * t90;
t85 = cos(pkin(6));
t100 = t85 * t88;
t91 = cos(qJ(2));
t99 = t85 * t91;
t98 = t86 * t91;
t97 = t89 * t91;
t78 = cos(t81);
t70 = t90 * pkin(3) + pkin(4) * t78;
t65 = t84 * t100 + t82 * t91;
t67 = -t82 * t100 + t84 * t91;
t73 = sin(t76);
t77 = sin(t81);
t92 = -g(1) * (t78 * t105 - t67 * t77) - g(2) * (-t78 * t104 - t65 * t77) - g(3) * (-t77 * t102 + t85 * t78);
t96 = t92 * MDP(17) + (-g(1) * (-t77 * t105 - t67 * t78) - g(2) * (t77 * t104 - t65 * t78) - g(3) * (-t78 * t102 - t85 * t77)) * MDP(18) + (-t89 * MDP(26) + t86 * MDP(27)) * (g(1) * (t74 * t105 - t67 * t73) + g(2) * (-t74 * t104 - t65 * t73) + g(3) * (-t73 * t102 + t85 * t74));
t64 = t82 * t88 - t84 * t99;
t66 = t82 * t99 + t84 * t88;
t94 = -g(1) * t66 - g(2) * t64 + t91 * t108;
t80 = -qJ(5) - pkin(9) - pkin(8);
t69 = -t87 * pkin(3) - pkin(4) * t77;
t68 = pkin(2) + t70;
t63 = t74 * t102 + t85 * t73;
t61 = t73 * t105 + t67 * t74;
t59 = -t73 * t104 + t65 * t74;
t1 = [(-MDP(1) - MDP(20)) * g(3); (-g(1) * (-t66 * t68 - t67 * t80) - g(2) * (-t64 * t68 - t65 * t80) - (t68 * t91 - t80 * t88) * t108) * MDP(20) + (-g(1) * (-t66 * t106 + t67 * t86) - g(2) * (-t64 * t106 + t65 * t86) - (t74 * t97 + t86 * t88) * t108) * MDP(26) + (-g(1) * (t66 * t107 + t67 * t89) - g(2) * (t64 * t107 + t65 * t89) - (-t74 * t98 + t88 * t89) * t108) * MDP(27) + (MDP(4) - MDP(19)) * (g(1) * t67 + g(2) * t65 + g(3) * t102) + (-t90 * MDP(10) + t87 * MDP(11) - MDP(17) * t78 + MDP(18) * t77 - MDP(3)) * t94; (-g(1) * (t82 * t101 - t67 * t87) - g(2) * (-t84 * t101 - t65 * t87) - g(3) * (-t87 * t102 + t85 * t90)) * MDP(10) + (-g(1) * (-t82 * t103 - t67 * t90) - g(2) * (t84 * t103 - t65 * t90) - g(3) * (-t88 * t101 - t85 * t87)) * MDP(11) + (-g(1) * (t70 * t105 + t67 * t69) - g(2) * (-t70 * t104 + t65 * t69) - g(3) * (t69 * t102 + t85 * t70)) * MDP(20) + t96; t92 * pkin(4) * MDP(20) + t96; t94 * MDP(20); (-g(1) * (-t61 * t86 + t66 * t89) - g(2) * (-t59 * t86 + t64 * t89) - g(3) * (-t63 * t86 - t83 * t97)) * MDP(26) + (-g(1) * (-t61 * t89 - t66 * t86) - g(2) * (-t59 * t89 - t64 * t86) - g(3) * (-t63 * t89 + t83 * t98)) * MDP(27);];
taug  = t1;
