% Calculate Gravitation load on the joints for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:58:45
% EndTime: 2019-03-09 04:58:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (209->52), mult. (185->77), div. (0->0), fcn. (155->12), ass. (0->32)
t61 = qJ(3) + qJ(4);
t55 = pkin(11) + t61;
t50 = sin(t55);
t79 = g(3) * t50;
t60 = qJ(1) + pkin(10);
t53 = sin(t60);
t62 = sin(qJ(6));
t77 = t53 * t62;
t65 = cos(qJ(6));
t76 = t53 * t65;
t54 = cos(t60);
t75 = t54 * t62;
t74 = t54 * t65;
t57 = cos(t61);
t66 = cos(qJ(3));
t73 = t66 * pkin(3) + pkin(4) * t57;
t51 = cos(t55);
t56 = sin(t61);
t71 = g(1) * t54 + g(2) * t53;
t68 = -g(3) * t57 + t71 * t56;
t72 = t68 * MDP(17) + (g(3) * t56 + t71 * t57) * MDP(18) + (t65 * MDP(26) - t62 * MDP(27)) * (-g(3) * t51 + t71 * t50);
t70 = g(1) * t53 - g(2) * t54;
t67 = cos(qJ(1));
t64 = sin(qJ(1));
t63 = sin(qJ(3));
t59 = -qJ(5) - pkin(8) - pkin(7);
t48 = pkin(2) + t73;
t47 = t51 * t74 + t77;
t46 = -t51 * t75 + t76;
t45 = -t51 * t76 + t75;
t44 = t51 * t77 + t74;
t1 = [(g(1) * t67 + g(2) * t64) * MDP(3) - t71 * MDP(19) + (-g(1) * (-t64 * pkin(1) - t53 * t48 - t54 * t59) - g(2) * (t67 * pkin(1) + t54 * t48 - t53 * t59)) * MDP(20) + (-g(1) * t45 - g(2) * t47) * MDP(26) + (-g(1) * t44 - g(2) * t46) * MDP(27) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t64 - g(2) * t67) + (t66 * MDP(10) - t63 * MDP(11) + MDP(17) * t57 - MDP(18) * t56) * t70; (-MDP(20) - MDP(4)) * g(3); (-g(3) * t66 + t71 * t63) * MDP(10) + (g(3) * t63 + t71 * t66) * MDP(11) + (-g(3) * t73 - t71 * (-t63 * pkin(3) - pkin(4) * t56)) * MDP(20) + t72; t68 * MDP(20) * pkin(4) + t72; -t70 * MDP(20); (-g(1) * t46 + g(2) * t44 + t62 * t79) * MDP(26) + (g(1) * t47 - g(2) * t45 + t65 * t79) * MDP(27);];
taug  = t1;
