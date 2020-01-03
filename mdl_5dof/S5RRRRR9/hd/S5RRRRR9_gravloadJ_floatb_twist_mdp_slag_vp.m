% Calculate Gravitation load on the joints for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:53
% EndTime: 2019-12-31 22:29:53
% DurationCPUTime: 0.21s
% Computational Cost: add. (233->57), mult. (268->87), div. (0->0), fcn. (274->10), ass. (0->38)
t65 = sin(qJ(2));
t89 = MDP(10) * t65;
t63 = qJ(3) + qJ(4);
t62 = qJ(5) + t63;
t58 = sin(t62);
t59 = cos(t62);
t60 = sin(t63);
t61 = cos(t63);
t64 = sin(qJ(3));
t67 = cos(qJ(3));
t88 = t67 * MDP(16) - t64 * MDP(17) + t61 * MDP(23) - t60 * MDP(24) + t59 * MDP(30) - t58 * MDP(31) + MDP(9);
t87 = g(3) * t65;
t66 = sin(qJ(1));
t68 = cos(qJ(2));
t86 = t66 * t68;
t69 = cos(qJ(1));
t85 = t69 * t58;
t84 = t69 * t59;
t83 = t69 * t60;
t82 = t69 * t61;
t81 = t69 * t64;
t80 = t69 * t67;
t46 = t58 * t86 + t84;
t47 = -t59 * t86 + t85;
t48 = t66 * t59 - t68 * t85;
t49 = t66 * t58 + t68 * t84;
t79 = (-g(1) * t48 + g(2) * t46 + t58 * t87) * MDP(30) + (g(1) * t49 - g(2) * t47 + t59 * t87) * MDP(31);
t50 = t60 * t86 + t82;
t51 = -t61 * t86 + t83;
t52 = t66 * t61 - t68 * t83;
t53 = t66 * t60 + t68 * t82;
t72 = (-g(1) * t52 + g(2) * t50 + t60 * t87) * MDP(23) + (g(1) * t53 - g(2) * t51 + t61 * t87) * MDP(24) + t79;
t71 = g(1) * t69 + g(2) * t66;
t57 = t66 * t64 + t68 * t80;
t56 = t66 * t67 - t68 * t81;
t55 = -t67 * t86 + t81;
t54 = t64 * t86 + t80;
t1 = [t71 * MDP(3) + (-g(1) * t55 - g(2) * t57) * MDP(16) + (-g(1) * t54 - g(2) * t56) * MDP(17) + (-g(1) * t51 - g(2) * t53) * MDP(23) + (-g(1) * t50 - g(2) * t52) * MDP(24) + (-g(1) * t47 - g(2) * t49) * MDP(30) + (-g(1) * t46 - g(2) * t48) * MDP(31) + (MDP(9) * t68 + MDP(2) - t89) * (g(1) * t66 - g(2) * t69); (-t88 * t68 + t89) * g(3) + (MDP(10) * t68 + t88 * t65) * t71; (-g(1) * t56 + g(2) * t54 + t64 * t87) * MDP(16) + (g(1) * t57 - g(2) * t55 + t67 * t87) * MDP(17) + t72; t72; t79;];
taug = t1;
