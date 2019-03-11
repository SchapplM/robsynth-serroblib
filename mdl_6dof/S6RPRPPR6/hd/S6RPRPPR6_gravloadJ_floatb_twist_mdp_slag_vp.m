% Calculate Gravitation load on the joints for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:54:38
% EndTime: 2019-03-09 02:54:40
% DurationCPUTime: 0.36s
% Computational Cost: add. (178->67), mult. (232->93), div. (0->0), fcn. (200->10), ass. (0->37)
t73 = sin(qJ(1));
t75 = cos(qJ(1));
t99 = -g(1) * t73 + g(2) * t75;
t68 = qJ(3) + pkin(9);
t60 = sin(t68);
t62 = cos(t68);
t78 = -g(3) * t60 - t62 * t99;
t72 = sin(qJ(3));
t98 = pkin(3) * t72;
t94 = g(3) * t62;
t67 = pkin(10) + qJ(6);
t59 = sin(t67);
t93 = t59 * t75;
t61 = cos(t67);
t92 = t61 * t75;
t69 = sin(pkin(10));
t91 = t69 * t75;
t70 = cos(pkin(10));
t90 = t70 * t75;
t89 = t73 * t59;
t88 = t73 * t61;
t87 = t73 * t69;
t86 = t73 * t70;
t85 = t75 * pkin(1) + t73 * qJ(2);
t84 = -MDP(15) - MDP(19);
t83 = -t73 * pkin(1) + t75 * qJ(2);
t55 = g(1) * t75 + g(2) * t73;
t82 = pkin(4) * t60 - qJ(5) * t62;
t71 = -qJ(4) - pkin(7);
t81 = t73 * t71 + t75 * t98 + t83;
t80 = -t71 * t75 + t73 * t98 + t85;
t74 = cos(qJ(3));
t53 = t60 * t92 - t89;
t52 = t60 * t93 + t88;
t51 = t60 * t88 + t93;
t50 = -t60 * t89 + t92;
t1 = [(-g(1) * t83 - g(2) * t85) * MDP(6) + (-g(1) * t81 - g(2) * t80) * MDP(15) + (-g(1) * (t60 * t90 - t87) - g(2) * (t60 * t86 + t91)) * MDP(16) + (-g(1) * (-t60 * t91 - t86) - g(2) * (-t60 * t87 + t90)) * MDP(17) + (-g(1) * (t82 * t75 + t81) - g(2) * (t82 * t73 + t80)) * MDP(19) + (-g(1) * t53 - g(2) * t51) * MDP(25) + (g(1) * t52 - g(2) * t50) * MDP(26) - (MDP(2) - MDP(4) + MDP(14)) * t99 + (-t72 * MDP(12) - t74 * MDP(13) + t62 * MDP(18) + MDP(3) - MDP(5)) * t55; -(-MDP(6) + t84) * t99; (g(3) * t74 - t72 * t99) * MDP(13) + (t60 * t99 - t94) * MDP(18) + (-g(3) * (-t82 - t98) + t99 * (pkin(3) * t74 + pkin(4) * t62 + qJ(5) * t60)) * MDP(19) + (pkin(3) * MDP(15) + MDP(12)) * (g(3) * t72 + t99 * t74) + (-MDP(16) * t70 + MDP(17) * t69 - MDP(25) * t61 + MDP(26) * t59) * t78; t84 * t55; t78 * MDP(19); (-g(1) * t50 - g(2) * t52 + t59 * t94) * MDP(25) + (g(1) * t51 - g(2) * t53 + t61 * t94) * MDP(26);];
taug  = t1;
