% Calculate Gravitation load on the joints for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:24:51
% EndTime: 2019-03-09 05:24:52
% DurationCPUTime: 0.30s
% Computational Cost: add. (187->64), mult. (258->89), div. (0->0), fcn. (237->8), ass. (0->37)
t79 = sin(qJ(1));
t82 = cos(qJ(1));
t106 = -g(1) * t79 + g(2) * t82;
t105 = MDP(13) - MDP(21);
t78 = sin(qJ(3));
t81 = cos(qJ(3));
t60 = -g(3) * t78 - t106 * t81;
t98 = g(3) * t81;
t71 = qJ(4) + pkin(10) + qJ(6);
t68 = sin(t71);
t97 = t79 * t68;
t69 = cos(t71);
t96 = t79 * t69;
t77 = sin(qJ(4));
t95 = t79 * t77;
t80 = cos(qJ(4));
t94 = t79 * t80;
t93 = t82 * t68;
t92 = t82 * t69;
t91 = t82 * t77;
t90 = t82 * t80;
t55 = -t78 * t97 + t92;
t56 = t78 * t96 + t93;
t57 = t78 * t93 + t96;
t58 = t78 * t92 - t97;
t89 = (-g(1) * t55 - g(2) * t57 + t68 * t98) * MDP(28) + (g(1) * t56 - g(2) * t58 + t69 * t98) * MDP(29);
t87 = pkin(4) * t77 + pkin(7);
t86 = g(2) * (t82 * pkin(1) + t79 * qJ(2));
t70 = t80 * pkin(4) + pkin(3);
t76 = -qJ(5) - pkin(8);
t84 = t78 * t70 + t81 * t76;
t63 = t78 * t91 + t94;
t61 = -t78 * t95 + t90;
t73 = t82 * qJ(2);
t64 = t78 * t90 - t95;
t62 = t78 * t94 + t91;
t1 = [(-g(1) * (-t79 * pkin(1) + t73) - t86) * MDP(6) + (-g(1) * t64 - g(2) * t62) * MDP(19) + (g(1) * t63 - g(2) * t61) * MDP(20) + (-g(1) * t73 - t86 + (-g(1) * t84 - g(2) * t87) * t82 + (-g(1) * (-pkin(1) - t87) - g(2) * t84) * t79) * MDP(22) + (-g(1) * t58 - g(2) * t56) * MDP(28) + (g(1) * t57 - g(2) * t55) * MDP(29) - (MDP(2) - MDP(4)) * t106 + (-t78 * MDP(12) - t105 * t81 + MDP(3) - MDP(5)) * (g(1) * t82 + g(2) * t79); -(-MDP(22) - MDP(6)) * t106; (g(3) * t84 + t106 * (t70 * t81 - t76 * t78)) * MDP(22) + t105 * (-t106 * t78 + t98) + (-MDP(19) * t80 + MDP(20) * t77 - MDP(28) * t69 + MDP(29) * t68 - MDP(12)) * t60; (g(1) * t62 - g(2) * t64 + t80 * t98) * MDP(20) + t89 + (MDP(22) * pkin(4) + MDP(19)) * (-g(1) * t61 - g(2) * t63 + t77 * t98); t60 * MDP(22); t89;];
taug  = t1;
