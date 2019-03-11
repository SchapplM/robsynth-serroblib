% Calculate Gravitation load on the joints for
% S6RPRRRP4
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
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:08:45
% EndTime: 2019-03-09 06:08:46
% DurationCPUTime: 0.31s
% Computational Cost: add. (284->61), mult. (267->80), div. (0->0), fcn. (228->10), ass. (0->33)
t75 = pkin(10) + qJ(3);
t71 = cos(t75);
t72 = qJ(4) + t75;
t67 = sin(t72);
t68 = cos(t72);
t81 = cos(qJ(5));
t69 = t81 * pkin(5) + pkin(4);
t78 = -qJ(6) - pkin(9);
t87 = -t67 * t78 + t68 * t69;
t103 = pkin(3) * t71 + t87;
t102 = MDP(21) - MDP(29);
t80 = sin(qJ(1));
t82 = cos(qJ(1));
t65 = g(1) * t82 + g(2) * t80;
t55 = -g(3) * t68 + t65 * t67;
t96 = g(3) * t67;
t79 = sin(qJ(5));
t93 = t80 * t79;
t92 = t80 * t81;
t91 = t82 * t79;
t90 = t82 * t81;
t88 = pkin(5) * t79 + pkin(7) + pkin(8) + qJ(2);
t86 = t102 * (t65 * t68 + t96) + (MDP(27) * t81 - MDP(28) * t79 + MDP(20)) * t55;
t64 = g(1) * t80 - g(2) * t82;
t85 = t67 * t69 + t68 * t78;
t60 = -t68 * t91 + t92;
t58 = t68 * t93 + t90;
t77 = cos(pkin(10));
t84 = -t77 * pkin(2) - pkin(1) - t103;
t70 = sin(t75);
t61 = t68 * t90 + t93;
t59 = -t68 * t92 + t91;
t1 = [(-g(1) * (-t80 * pkin(1) + t82 * qJ(2)) - g(2) * (t82 * pkin(1) + t80 * qJ(2))) * MDP(7) + (-g(1) * t59 - g(2) * t61) * MDP(27) + (-g(1) * t58 - g(2) * t60) * MDP(28) + ((-g(1) * t88 + g(2) * t84) * t82 + (-g(1) * t84 - g(2) * t88) * t80) * MDP(30) + (MDP(3) - MDP(6)) * t65 + (-t102 * t67 + t71 * MDP(13) - t70 * MDP(14) + t68 * MDP(20) + MDP(4) * t77 - MDP(5) * sin(pkin(10)) + MDP(2)) * t64; (-MDP(30) - MDP(7)) * t64; (-g(3) * t71 + t65 * t70) * MDP(13) + (g(3) * t70 + t65 * t71) * MDP(14) + (-g(3) * t103 + t65 * (pkin(3) * t70 + t85)) * MDP(30) + t86; (-g(3) * t87 + t65 * t85) * MDP(30) + t86; (g(1) * t61 - g(2) * t59 + t81 * t96) * MDP(28) + (pkin(5) * MDP(30) + MDP(27)) * (-g(1) * t60 + g(2) * t58 + t79 * t96); -t55 * MDP(30);];
taug  = t1;
