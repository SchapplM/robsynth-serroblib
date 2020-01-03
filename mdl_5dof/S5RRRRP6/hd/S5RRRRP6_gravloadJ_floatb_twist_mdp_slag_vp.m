% Calculate Gravitation load on the joints for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:34
% EndTime: 2019-12-31 21:54:35
% DurationCPUTime: 0.19s
% Computational Cost: add. (184->50), mult. (239->69), div. (0->0), fcn. (208->8), ass. (0->30)
t70 = cos(qJ(2));
t69 = cos(qJ(4));
t59 = t69 * pkin(4) + pkin(3);
t64 = qJ(2) + qJ(3);
t61 = sin(t64);
t62 = cos(t64);
t65 = -qJ(5) - pkin(8);
t79 = t62 * t59 - t61 * t65;
t95 = t70 * pkin(2) + t79;
t94 = MDP(17) - MDP(25);
t68 = sin(qJ(1));
t71 = cos(qJ(1));
t77 = g(1) * t71 + g(2) * t68;
t51 = -g(3) * t62 + t77 * t61;
t88 = g(3) * t61;
t66 = sin(qJ(4));
t85 = t68 * t66;
t84 = t68 * t69;
t83 = t71 * t66;
t82 = t71 * t69;
t80 = pkin(4) * t66 + pkin(6) + pkin(7);
t78 = t94 * (t77 * t62 + t88) + (MDP(23) * t69 - MDP(24) * t66 + MDP(16)) * t51;
t75 = t59 * t61 + t62 * t65;
t56 = -t62 * t83 + t84;
t54 = t62 * t85 + t82;
t74 = pkin(1) + t95;
t67 = sin(qJ(2));
t57 = t62 * t82 + t85;
t55 = -t62 * t84 + t83;
t1 = [t77 * MDP(3) + (-g(1) * t55 - g(2) * t57) * MDP(23) + (-g(1) * t54 - g(2) * t56) * MDP(24) + ((-g(1) * t80 - g(2) * t74) * t71 + (g(1) * t74 - g(2) * t80) * t68) * MDP(26) + (-t67 * MDP(10) + t62 * MDP(16) + t70 * MDP(9) - t94 * t61 + MDP(2)) * (g(1) * t68 - g(2) * t71); (-g(3) * t70 + t77 * t67) * MDP(9) + (g(3) * t67 + t77 * t70) * MDP(10) + (-g(3) * t95 + t77 * (pkin(2) * t67 + t75)) * MDP(26) + t78; (-g(3) * t79 + t77 * t75) * MDP(26) + t78; (g(1) * t57 - g(2) * t55 + t69 * t88) * MDP(24) + (pkin(4) * MDP(26) + MDP(23)) * (-g(1) * t56 + g(2) * t54 + t66 * t88); -t51 * MDP(26);];
taug = t1;
