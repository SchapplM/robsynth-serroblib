% Calculate Gravitation load on the joints for
% S5RRPRR6
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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:55
% EndTime: 2020-01-03 12:05:55
% DurationCPUTime: 0.14s
% Computational Cost: add. (215->47), mult. (201->76), div. (0->0), fcn. (196->10), ass. (0->31)
t72 = sin(pkin(9));
t86 = g(1) * t72;
t71 = qJ(1) + qJ(2);
t67 = sin(t71);
t73 = cos(pkin(9));
t85 = t67 * t73;
t69 = cos(t71);
t84 = t69 * t73;
t74 = sin(qJ(4));
t83 = t73 * t74;
t76 = cos(qJ(4));
t82 = t73 * t76;
t70 = qJ(4) + qJ(5);
t66 = sin(t70);
t68 = cos(t70);
t48 = -t66 * t85 - t68 * t69;
t49 = -t66 * t69 + t68 * t85;
t50 = t66 * t84 - t67 * t68;
t51 = t66 * t67 + t68 * t84;
t81 = (-g(2) * t48 - g(3) * t50 + t66 * t86) * MDP(23) + (g(2) * t49 - g(3) * t51 + t68 * t86) * MDP(24);
t80 = t69 * pkin(2) + t67 * qJ(3);
t79 = t67 * pkin(2) - qJ(3) * t69;
t62 = g(2) * t69 + g(3) * t67;
t54 = -t67 * t83 - t69 * t76;
t55 = t67 * t82 - t69 * t74;
t56 = -t67 * t76 + t69 * t83;
t57 = t67 * t74 + t69 * t82;
t78 = (g(2) * t50 - g(3) * t48) * MDP(24) + (-g(2) * t51 - g(3) * t49) * MDP(23) + (g(2) * t56 - g(3) * t54) * MDP(17) + (-g(2) * t57 - g(3) * t55) * MDP(16) + (MDP(6) - MDP(9)) * (g(2) * t67 - g(3) * t69) + (-t73 * MDP(7) + t72 * MDP(8) - MDP(5)) * t62;
t77 = cos(qJ(1));
t75 = sin(qJ(1));
t1 = [(-g(2) * t77 - g(3) * t75) * MDP(2) + (g(2) * t75 - g(3) * t77) * MDP(3) + (-g(2) * (pkin(1) * t77 + t80) - g(3) * (pkin(1) * t75 + t79)) * MDP(10) + t78; (-g(2) * t80 - g(3) * t79) * MDP(10) + t78; t62 * MDP(10); (-g(2) * t54 - g(3) * t56 + t74 * t86) * MDP(16) + (g(2) * t55 - g(3) * t57 + t76 * t86) * MDP(17) + t81; t81;];
taug = t1;
