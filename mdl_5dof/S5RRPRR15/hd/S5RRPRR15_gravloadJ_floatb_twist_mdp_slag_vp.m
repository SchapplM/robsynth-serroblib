% Calculate Gravitation load on the joints for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR15_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR15_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:09
% EndTime: 2019-12-31 20:43:10
% DurationCPUTime: 0.20s
% Computational Cost: add. (130->53), mult. (226->77), div. (0->0), fcn. (213->8), ass. (0->34)
t64 = sin(qJ(1));
t67 = cos(qJ(1));
t57 = g(1) * t67 + g(2) * t64;
t88 = MDP(9) - MDP(12);
t87 = MDP(10) - MDP(13);
t66 = cos(qJ(2));
t82 = g(3) * t66;
t61 = qJ(4) + qJ(5);
t58 = sin(t61);
t81 = t64 * t58;
t59 = cos(t61);
t80 = t64 * t59;
t62 = sin(qJ(4));
t79 = t64 * t62;
t65 = cos(qJ(4));
t78 = t64 * t65;
t77 = t67 * t58;
t76 = t67 * t59;
t75 = t67 * t62;
t74 = t67 * t65;
t63 = sin(qJ(2));
t45 = t63 * t76 - t81;
t46 = t63 * t77 + t80;
t47 = t63 * t80 + t77;
t48 = -t63 * t81 + t76;
t73 = (-g(1) * t45 - g(2) * t47 + t59 * t82) * MDP(27) + (g(1) * t46 - g(2) * t48 - t58 * t82) * MDP(28);
t71 = t66 * pkin(2) + t63 * qJ(3);
t69 = pkin(1) + t71;
t54 = -t63 * t79 + t74;
t53 = t63 * t78 + t75;
t52 = t63 * t75 + t78;
t51 = t63 * t74 - t79;
t49 = t57 * t63 - t82;
t1 = [((-g(1) * pkin(6) - g(2) * t69) * t67 + (-g(2) * pkin(6) + g(1) * t69) * t64) * MDP(14) + (-g(1) * t54 - g(2) * t52) * MDP(20) + (g(1) * t53 - g(2) * t51) * MDP(21) + (-g(1) * t48 - g(2) * t46) * MDP(27) + (g(1) * t47 - g(2) * t45) * MDP(28) + (MDP(3) - MDP(11)) * t57 + (-t87 * t63 + t88 * t66 + MDP(2)) * (g(1) * t64 - g(2) * t67); (-g(3) * t71 + t57 * (pkin(2) * t63 - qJ(3) * t66)) * MDP(14) + t88 * t49 + (-MDP(20) * t62 - MDP(21) * t65 - MDP(27) * t58 - MDP(28) * t59 + t87) * (g(3) * t63 + t57 * t66); -t49 * MDP(14); (-g(1) * t51 - g(2) * t53 + t65 * t82) * MDP(20) + (g(1) * t52 - g(2) * t54 - t62 * t82) * MDP(21) + t73; t73;];
taug = t1;
