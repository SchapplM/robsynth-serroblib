% Calculate Gravitation load on the joints for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:46:39
% EndTime: 2021-01-16 03:46:43
% DurationCPUTime: 0.66s
% Computational Cost: add. (374->101), mult. (618->173), div. (0->0), fcn. (732->14), ass. (0->43)
t77 = sin(pkin(6));
t108 = g(3) * t77;
t107 = cos(qJ(2));
t74 = qJ(3) + pkin(12);
t71 = cos(t74);
t75 = qJ(5) + qJ(6);
t72 = sin(t75);
t106 = t71 * t72;
t73 = cos(t75);
t105 = t71 * t73;
t79 = sin(qJ(5));
t104 = t71 * t79;
t82 = cos(qJ(5));
t103 = t71 * t82;
t76 = sin(pkin(11));
t102 = t76 * t77;
t81 = sin(qJ(2));
t101 = t77 * t81;
t83 = cos(qJ(3));
t100 = t77 * t83;
t97 = cos(pkin(11));
t98 = cos(pkin(6));
t90 = t98 * t97;
t65 = t76 * t107 + t81 * t90;
t70 = sin(t74);
t91 = t77 * t97;
t58 = t65 * t71 - t70 * t91;
t92 = t76 * t98;
t63 = -t97 * t107 + t81 * t92;
t60 = t70 * t102 - t63 * t71;
t62 = t71 * t101 + t98 * t70;
t64 = -t107 * t90 + t76 * t81;
t66 = t107 * t92 + t97 * t81;
t95 = t77 * t107;
t99 = (-g(1) * (-t60 * t72 + t66 * t73) - g(2) * (-t58 * t72 + t64 * t73) - g(3) * (-t62 * t72 - t73 * t95)) * MDP(28) + (-g(1) * (-t60 * t73 - t66 * t72) - g(2) * (-t58 * t73 - t64 * t72) - g(3) * (-t62 * t73 + t72 * t95)) * MDP(29);
t96 = t71 * t107;
t94 = t79 * t107;
t93 = t82 * t107;
t84 = g(1) * t66 + g(2) * t64 - g(3) * t95;
t80 = sin(qJ(3));
t78 = qJ(4) + pkin(8);
t69 = pkin(3) * t83 + pkin(2);
t1 = [(-MDP(1) - MDP(15)) * g(3); (-g(1) * (-t63 * t78 - t66 * t69) - g(2) * (-t64 * t69 + t65 * t78) - (t107 * t69 + t78 * t81) * t108) * MDP(15) + (-g(1) * (-t66 * t103 - t63 * t79) - g(2) * (-t64 * t103 + t65 * t79) - (t71 * t93 + t79 * t81) * t108) * MDP(21) + (-g(1) * (t66 * t104 - t63 * t82) - g(2) * (t64 * t104 + t65 * t82) - (-t71 * t94 + t81 * t82) * t108) * MDP(22) + (-g(1) * (-t66 * t105 - t63 * t72) - g(2) * (-t64 * t105 + t65 * t72) - (t72 * t81 + t73 * t96) * t108) * MDP(28) + (-g(1) * (t66 * t106 - t63 * t73) - g(2) * (t64 * t106 + t65 * t73) - (-t72 * t96 + t73 * t81) * t108) * MDP(29) + (MDP(4) - MDP(14)) * (-g(1) * t63 + g(2) * t65 + g(3) * t101) + (t83 * MDP(10) - t80 * MDP(11) + t71 * MDP(12) - t70 * MDP(13) + MDP(3)) * t84; (-g(1) * (-t80 * t102 + t63 * t83) - g(2) * (-t65 * t83 + t80 * t91) - g(3) * (-t81 * t100 - t98 * t80)) * MDP(11) + (g(1) * t60 + g(2) * t58 + g(3) * t62) * MDP(13) + (-MDP(21) * t82 + MDP(22) * t79 - MDP(28) * t73 + MDP(29) * t72 - MDP(12)) * (g(1) * (t71 * t102 + t63 * t70) + g(2) * (-t65 * t70 - t71 * t91) + g(3) * (-t70 * t101 + t98 * t71)) + (pkin(3) * MDP(15) + MDP(10)) * (-g(1) * (t76 * t100 + t63 * t80) - g(2) * (-t65 * t80 - t83 * t91) - g(3) * (-t80 * t101 + t98 * t83)); -t84 * MDP(15); (-g(1) * (-t60 * t79 + t66 * t82) - g(2) * (-t58 * t79 + t64 * t82) - g(3) * (-t62 * t79 - t77 * t93)) * MDP(21) + (-g(1) * (-t60 * t82 - t66 * t79) - g(2) * (-t58 * t82 - t64 * t79) - g(3) * (-t62 * t82 + t77 * t94)) * MDP(22) + t99; t99;];
taug = t1;
