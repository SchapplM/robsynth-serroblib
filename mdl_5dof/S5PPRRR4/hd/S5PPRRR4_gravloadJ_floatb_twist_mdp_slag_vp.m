% Calculate Gravitation load on the joints for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:44
% EndTime: 2019-12-05 15:19:45
% DurationCPUTime: 0.24s
% Computational Cost: add. (256->63), mult. (720->119), div. (0->0), fcn. (925->14), ass. (0->43)
t103 = sin(pkin(6));
t80 = sin(pkin(5));
t101 = t80 * t103;
t105 = cos(pkin(6));
t81 = cos(pkin(10));
t106 = cos(pkin(5));
t102 = sin(pkin(10));
t79 = sin(pkin(11));
t98 = t102 * t79;
t104 = cos(pkin(11));
t99 = t81 * t104;
t89 = -t106 * t99 + t98;
t112 = t81 * t101 + t89 * t105;
t100 = t80 * t102;
t109 = t81 * t79;
t94 = t102 * t104;
t90 = t106 * t94 + t109;
t111 = -t103 * t100 + t90 * t105;
t110 = cos(qJ(3));
t82 = sin(qJ(5));
t86 = cos(qJ(4));
t108 = t82 * t86;
t85 = cos(qJ(5));
t107 = t85 * t86;
t96 = t106 * t103;
t95 = t105 * t104;
t84 = sin(qJ(3));
t83 = sin(qJ(4));
t75 = -t106 * t98 + t99;
t74 = t106 * t109 + t94;
t73 = -t104 * t101 + t106 * t105;
t70 = t105 * t100 + t90 * t103;
t69 = -t81 * t80 * t105 + t89 * t103;
t68 = t84 * t96 + (t110 * t79 + t84 * t95) * t80;
t67 = t79 * t80 * t84 + (-t95 * t80 - t96) * t110;
t66 = t68 * t86 + t73 * t83;
t64 = t75 * t110 - t111 * t84;
t63 = t111 * t110 + t75 * t84;
t62 = t74 * t110 - t112 * t84;
t61 = t112 * t110 + t74 * t84;
t60 = t64 * t86 + t70 * t83;
t58 = t62 * t86 + t69 * t83;
t1 = [(-MDP(1) - MDP(2)) * g(3); (-g(3) * t106 + (-t102 * g(1) + g(2) * t81) * t80) * MDP(2); (g(1) * t64 + g(2) * t62 + g(3) * t68) * MDP(5) + (-g(1) * (-t63 * t107 + t64 * t82) - g(2) * (-t61 * t107 + t62 * t82) - g(3) * (-t67 * t107 + t68 * t82)) * MDP(18) + (-g(1) * (t63 * t108 + t64 * t85) - g(2) * (t61 * t108 + t62 * t85) - g(3) * (t67 * t108 + t68 * t85)) * MDP(19) + (t86 * MDP(11) - MDP(12) * t83 + MDP(4)) * (g(1) * t63 + g(2) * t61 + g(3) * t67); (g(1) * t60 + g(2) * t58 + g(3) * t66) * MDP(12) + (-MDP(18) * t85 + MDP(19) * t82 - MDP(11)) * (g(1) * (-t64 * t83 + t70 * t86) + g(2) * (-t62 * t83 + t69 * t86) + g(3) * (-t68 * t83 + t73 * t86)); (-g(1) * (-t60 * t82 + t63 * t85) - g(2) * (-t58 * t82 + t61 * t85) - g(3) * (-t66 * t82 + t67 * t85)) * MDP(18) + (-g(1) * (-t60 * t85 - t63 * t82) - g(2) * (-t58 * t85 - t61 * t82) - g(3) * (-t66 * t85 - t67 * t82)) * MDP(19);];
taug = t1;
