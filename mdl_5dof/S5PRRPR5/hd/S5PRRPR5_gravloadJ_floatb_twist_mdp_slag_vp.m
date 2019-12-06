% Calculate Gravitation load on the joints for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:55
% EndTime: 2019-12-05 16:27:57
% DurationCPUTime: 0.36s
% Computational Cost: add. (159->67), mult. (326->122), div. (0->0), fcn. (374->12), ass. (0->36)
t61 = sin(pkin(5));
t86 = g(3) * t61;
t59 = qJ(3) + pkin(10);
t58 = cos(t59);
t65 = sin(qJ(5));
t85 = t58 * t65;
t68 = cos(qJ(5));
t84 = t58 * t68;
t60 = sin(pkin(9));
t83 = t60 * t61;
t62 = cos(pkin(9));
t82 = t61 * t62;
t66 = sin(qJ(3));
t81 = t61 * t66;
t67 = sin(qJ(2));
t80 = t61 * t67;
t69 = cos(qJ(3));
t79 = t61 * t69;
t63 = cos(pkin(5));
t78 = t63 * t67;
t70 = cos(qJ(2));
t77 = t63 * t70;
t76 = t65 * t70;
t75 = t68 * t70;
t50 = t60 * t67 - t62 * t77;
t52 = t60 * t77 + t62 * t67;
t73 = -g(1) * t52 - g(2) * t50 + t70 * t86;
t64 = -qJ(4) - pkin(7);
t57 = sin(t59);
t56 = pkin(3) * t69 + pkin(2);
t53 = -t60 * t78 + t62 * t70;
t51 = t60 * t70 + t62 * t78;
t49 = t57 * t63 + t58 * t80;
t47 = t53 * t58 + t57 * t83;
t45 = t51 * t58 - t57 * t82;
t1 = [(-MDP(1) - MDP(13)) * g(3); (-g(1) * (-t52 * t56 - t53 * t64) - g(2) * (-t50 * t56 - t51 * t64) - (t56 * t70 - t64 * t67) * t86) * MDP(13) + (-g(1) * (-t52 * t84 + t53 * t65) - g(2) * (-t50 * t84 + t51 * t65) - (t58 * t75 + t65 * t67) * t86) * MDP(19) + (-g(1) * (t52 * t85 + t53 * t68) - g(2) * (t50 * t85 + t51 * t68) - (-t58 * t76 + t67 * t68) * t86) * MDP(20) + (MDP(4) - MDP(12)) * (g(1) * t53 + g(2) * t51 + g(3) * t80) + (-t69 * MDP(10) + t66 * MDP(11) - MDP(3)) * t73; (-g(1) * (-t53 * t69 - t60 * t81) - g(2) * (-t51 * t69 + t62 * t81) - g(3) * (-t63 * t66 - t67 * t79)) * MDP(11) + (-MDP(19) * t68 + MDP(20) * t65) * (g(1) * (-t53 * t57 + t58 * t83) + g(2) * (-t51 * t57 - t58 * t82) + g(3) * (-t57 * t80 + t58 * t63)) + (MDP(13) * pkin(3) + MDP(10)) * (-g(1) * (-t53 * t66 + t60 * t79) - g(2) * (-t51 * t66 - t62 * t79) - g(3) * (t63 * t69 - t66 * t80)); t73 * MDP(13); (-g(1) * (-t47 * t65 + t52 * t68) - g(2) * (-t45 * t65 + t50 * t68) - g(3) * (-t49 * t65 - t61 * t75)) * MDP(19) + (-g(1) * (-t47 * t68 - t52 * t65) - g(2) * (-t45 * t68 - t50 * t65) - g(3) * (-t49 * t68 + t61 * t76)) * MDP(20);];
taug = t1;
