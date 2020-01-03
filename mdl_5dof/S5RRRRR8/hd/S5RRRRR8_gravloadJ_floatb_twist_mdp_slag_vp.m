% Calculate Gravitation load on the joints for
% S5RRRRR8
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
%   see S5RRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:47
% EndTime: 2019-12-31 22:25:48
% DurationCPUTime: 0.15s
% Computational Cost: add. (201->47), mult. (236->70), div. (0->0), fcn. (226->10), ass. (0->33)
t63 = qJ(2) + qJ(3);
t59 = sin(t63);
t85 = g(3) * t59;
t62 = qJ(4) + qJ(5);
t58 = sin(t62);
t66 = sin(qJ(1));
t83 = t66 * t58;
t60 = cos(t62);
t82 = t66 * t60;
t64 = sin(qJ(4));
t81 = t66 * t64;
t67 = cos(qJ(4));
t80 = t66 * t67;
t69 = cos(qJ(1));
t79 = t69 * t58;
t78 = t69 * t60;
t77 = t69 * t64;
t76 = t69 * t67;
t61 = cos(t63);
t50 = t61 * t83 + t78;
t51 = -t61 * t82 + t79;
t52 = -t61 * t79 + t82;
t53 = t61 * t78 + t83;
t75 = (-g(1) * t52 + g(2) * t50 + t58 * t85) * MDP(30) + (g(1) * t53 - g(2) * t51 + t60 * t85) * MDP(31);
t74 = g(1) * t69 + g(2) * t66;
t72 = (t61 * t74 + t85) * MDP(17) + (t67 * MDP(23) - t64 * MDP(24) + t60 * MDP(30) - t58 * MDP(31) + MDP(16)) * (-g(3) * t61 + t59 * t74);
t68 = cos(qJ(2));
t65 = sin(qJ(2));
t57 = t61 * t76 + t81;
t56 = -t61 * t77 + t80;
t55 = -t61 * t80 + t77;
t54 = t61 * t81 + t76;
t1 = [t74 * MDP(3) + (-g(1) * t55 - g(2) * t57) * MDP(23) + (-g(1) * t54 - g(2) * t56) * MDP(24) + (-g(1) * t51 - g(2) * t53) * MDP(30) + (-g(1) * t50 - g(2) * t52) * MDP(31) + (-MDP(10) * t65 + MDP(16) * t61 - MDP(17) * t59 + MDP(9) * t68 + MDP(2)) * (g(1) * t66 - g(2) * t69); (-g(3) * t68 + t65 * t74) * MDP(9) + (g(3) * t65 + t68 * t74) * MDP(10) + t72; t72; (-g(1) * t56 + g(2) * t54 + t64 * t85) * MDP(23) + (g(1) * t57 - g(2) * t55 + t67 * t85) * MDP(24) + t75; t75;];
taug = t1;
