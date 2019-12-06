% Calculate Gravitation load on the joints for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:18
% EndTime: 2019-12-05 16:52:19
% DurationCPUTime: 0.20s
% Computational Cost: add. (205->49), mult. (288->74), div. (0->0), fcn. (284->8), ass. (0->30)
t90 = MDP(17) + MDP(19);
t89 = MDP(18) - MDP(21);
t64 = sin(pkin(8));
t65 = cos(pkin(8));
t77 = g(1) * t65 + g(2) * t64;
t67 = sin(qJ(2));
t86 = g(3) * t67;
t63 = qJ(3) + qJ(4);
t61 = sin(t63);
t84 = t61 * t67;
t62 = cos(t63);
t83 = t62 * t67;
t69 = cos(qJ(2));
t82 = t64 * t69;
t81 = t65 * t69;
t66 = sin(qJ(3));
t80 = t66 * t69;
t68 = cos(qJ(3));
t79 = t68 * t69;
t52 = t61 * t82 + t65 * t62;
t54 = t61 * t81 - t64 * t62;
t47 = g(1) * t54 + g(2) * t52 + g(3) * t84;
t53 = -t65 * t61 + t62 * t82;
t55 = t64 * t61 + t62 * t81;
t78 = t89 * (g(1) * t55 + g(2) * t53 + g(3) * t83) + t90 * t47;
t76 = t68 * pkin(3) + pkin(4) * t62 + qJ(5) * t61 + pkin(2);
t72 = -g(1) * (-t54 * pkin(4) + t55 * qJ(5)) - g(2) * (-t52 * pkin(4) + t53 * qJ(5)) - g(3) * (-pkin(4) * t84 + qJ(5) * t83);
t71 = -g(1) * (t64 * t68 - t65 * t80) - g(2) * (-t64 * t80 - t65 * t68) + t66 * t86;
t70 = -pkin(7) - pkin(6);
t1 = [(-MDP(1) - MDP(22)) * g(3); (-g(3) * (-t67 * t70 + t76 * t69) + t77 * (t76 * t67 + t69 * t70)) * MDP(22) + (MDP(4) - MDP(20)) * (t77 * t69 + t86) + (t68 * MDP(10) - t66 * MDP(11) - t89 * t61 + t90 * t62 + MDP(3)) * (-g(3) * t69 + t77 * t67); t71 * MDP(10) + (-g(1) * (-t64 * t66 - t65 * t79) - g(2) * (-t64 * t79 + t65 * t66) + t68 * t86) * MDP(11) + (t71 * pkin(3) + t72) * MDP(22) + t78; t72 * MDP(22) + t78; -t47 * MDP(22);];
taug = t1;
