% Calculate Gravitation load on the joints for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:47
% EndTime: 2019-03-09 02:26:48
% DurationCPUTime: 0.26s
% Computational Cost: add. (183->57), mult. (321->90), div. (0->0), fcn. (368->10), ass. (0->34)
t66 = sin(qJ(4));
t92 = t66 * MDP(16);
t64 = qJ(5) + qJ(6);
t58 = sin(t64);
t59 = cos(t64);
t65 = sin(qJ(5));
t67 = cos(qJ(5));
t91 = t67 * MDP(22) - t65 * MDP(23) + t59 * MDP(29) - t58 * MDP(30) + MDP(15);
t90 = g(3) * t66;
t89 = cos(qJ(1));
t88 = sin(qJ(1));
t68 = cos(qJ(4));
t87 = t58 * t68;
t86 = t59 * t68;
t85 = t65 * t68;
t84 = t67 * t68;
t80 = sin(pkin(10));
t81 = cos(pkin(10));
t51 = -t88 * t80 - t89 * t81;
t52 = t89 * t80 - t88 * t81;
t47 = t51 * t87 + t52 * t59;
t48 = -t51 * t86 + t52 * t58;
t71 = -t51 * t59 + t52 * t87;
t72 = t51 * t58 + t52 * t86;
t83 = (-g(1) * t47 - g(2) * t71 - t58 * t90) * MDP(29) + (g(1) * t48 - g(2) * t72 - t59 * t90) * MDP(30);
t82 = t89 * pkin(1) + t88 * qJ(2);
t75 = -t88 * pkin(1) + t89 * qJ(2);
t73 = g(1) * t51 + g(2) * t52;
t70 = t51 * t65 + t52 * t84;
t69 = -t51 * t67 + t52 * t85;
t53 = g(1) * t88 - g(2) * t89;
t50 = -t51 * t84 + t52 * t65;
t49 = t51 * t85 + t52 * t67;
t1 = [(-g(1) * t75 - g(2) * t82) * MDP(6) + t73 * MDP(8) + (-g(1) * (-t88 * pkin(2) + t75) - g(2) * (t89 * pkin(2) + t82)) * MDP(9) + (-g(1) * t70 - g(2) * t50) * MDP(22) + (g(1) * t69 - g(2) * t49) * MDP(23) + (-g(1) * t72 - g(2) * t48) * MDP(29) + (g(1) * t71 - g(2) * t47) * MDP(30) + (MDP(3) - MDP(5)) * (g(1) * t89 + g(2) * t88) + (MDP(2) + MDP(4)) * t53 + (-t68 * MDP(15) - MDP(7) + t92) * (g(1) * t52 - g(2) * t51); (-MDP(6) - MDP(9)) * t53; g(3) * MDP(9); (t91 * t68 - t92) * g(3) + (-MDP(16) * t68 - t91 * t66) * t73; (-g(1) * t49 - g(2) * t69 - t65 * t90) * MDP(22) + (g(1) * t50 - g(2) * t70 - t67 * t90) * MDP(23) + t83; t83;];
taug  = t1;
