% Calculate Gravitation load on the joints for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:33:53
% EndTime: 2019-12-31 20:33:54
% DurationCPUTime: 0.30s
% Computational Cost: add. (195->62), mult. (248->94), div. (0->0), fcn. (238->10), ass. (0->35)
t65 = sin(qJ(1));
t67 = cos(qJ(1));
t73 = g(1) * t67 + g(2) * t65;
t86 = MDP(10) - MDP(13);
t64 = sin(qJ(2));
t66 = cos(qJ(2));
t53 = -g(3) * t66 + t73 * t64;
t83 = g(3) * t64;
t81 = t65 * t66;
t61 = pkin(9) + qJ(4);
t60 = qJ(5) + t61;
t56 = sin(t60);
t80 = t67 * t56;
t57 = cos(t60);
t79 = t67 * t57;
t58 = sin(t61);
t78 = t67 * t58;
t59 = cos(t61);
t77 = t67 * t59;
t62 = sin(pkin(9));
t76 = t67 * t62;
t63 = cos(pkin(9));
t75 = t67 * t63;
t45 = t56 * t81 + t79;
t46 = -t57 * t81 + t80;
t47 = t65 * t57 - t66 * t80;
t48 = t65 * t56 + t66 * t79;
t74 = (-g(1) * t47 + g(2) * t45 + t56 * t83) * MDP(27) + (g(1) * t48 - g(2) * t46 + t57 * t83) * MDP(28);
t71 = t66 * pkin(2) + t64 * qJ(3);
t69 = pkin(1) + t71;
t52 = t65 * t58 + t66 * t77;
t51 = t65 * t59 - t66 * t78;
t50 = -t59 * t81 + t78;
t49 = t58 * t81 + t77;
t1 = [t73 * MDP(3) + (-g(1) * (-t63 * t81 + t76) - g(2) * (t65 * t62 + t66 * t75)) * MDP(11) + (-g(1) * (t62 * t81 + t75) - g(2) * (t65 * t63 - t66 * t76)) * MDP(12) + ((-g(1) * pkin(6) - g(2) * t69) * t67 + (-g(2) * pkin(6) + g(1) * t69) * t65) * MDP(14) + (-g(1) * t50 - g(2) * t52) * MDP(20) + (-g(1) * t49 - g(2) * t51) * MDP(21) + (-g(1) * t46 - g(2) * t48) * MDP(27) + (-g(1) * t45 - g(2) * t47) * MDP(28) + (t66 * MDP(9) - t86 * t64 + MDP(2)) * (g(1) * t65 - g(2) * t67); (-g(3) * t71 + t73 * (pkin(2) * t64 - qJ(3) * t66)) * MDP(14) + t86 * (t73 * t66 + t83) + (t63 * MDP(11) - t62 * MDP(12) + t59 * MDP(20) - t58 * MDP(21) + t57 * MDP(27) - t56 * MDP(28) + MDP(9)) * t53; -t53 * MDP(14); (-g(1) * t51 + g(2) * t49 + t58 * t83) * MDP(20) + (g(1) * t52 - g(2) * t50 + t59 * t83) * MDP(21) + t74; t74;];
taug = t1;
