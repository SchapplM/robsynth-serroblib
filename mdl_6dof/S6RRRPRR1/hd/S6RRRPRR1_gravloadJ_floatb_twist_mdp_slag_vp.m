% Calculate Gravitation load on the joints for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:04:06
% EndTime: 2019-03-09 18:04:07
% DurationCPUTime: 0.16s
% Computational Cost: add. (281->50), mult. (240->71), div. (0->0), fcn. (205->10), ass. (0->30)
t59 = qJ(2) + qJ(3);
t54 = pkin(11) + qJ(5) + t59;
t51 = sin(t54);
t77 = g(3) * t51;
t60 = sin(qJ(6));
t65 = cos(qJ(1));
t75 = t60 * t65;
t62 = sin(qJ(1));
t74 = t62 * t60;
t63 = cos(qJ(6));
t73 = t62 * t63;
t72 = t63 * t65;
t56 = cos(t59);
t64 = cos(qJ(2));
t71 = t64 * pkin(2) + pkin(3) * t56;
t52 = cos(t54);
t69 = g(1) * t65 + g(2) * t62;
t70 = (t69 * t52 + t77) * MDP(26) + (t63 * MDP(32) - t60 * MDP(33) + MDP(25)) * (-g(3) * t52 + t69 * t51);
t49 = g(1) * t62 - g(2) * t65;
t55 = sin(t59);
t66 = -g(3) * t56 + t69 * t55;
t68 = t66 * MDP(16) + (g(3) * t55 + t69 * t56) * MDP(17) + t70;
t61 = sin(qJ(2));
t58 = -qJ(4) - pkin(8) - pkin(7);
t47 = pkin(1) + t71;
t46 = t52 * t72 + t74;
t45 = -t52 * t75 + t73;
t44 = -t52 * t73 + t75;
t43 = t52 * t74 + t72;
t1 = [(-g(1) * (-t62 * t47 - t58 * t65) - g(2) * (t47 * t65 - t62 * t58)) * MDP(19) + (-g(1) * t44 - g(2) * t46) * MDP(32) + (-g(1) * t43 - g(2) * t45) * MDP(33) + (MDP(3) - MDP(18)) * t69 + (-t61 * MDP(10) + MDP(16) * t56 - MDP(17) * t55 + MDP(25) * t52 - MDP(26) * t51 + t64 * MDP(9) + MDP(2)) * t49; (-g(3) * t64 + t69 * t61) * MDP(9) + (g(3) * t61 + t69 * t64) * MDP(10) + (-g(3) * t71 - t69 * (-pkin(2) * t61 - pkin(3) * t55)) * MDP(19) + t68; t66 * MDP(19) * pkin(3) + t68; -t49 * MDP(19); t70; (-g(1) * t45 + g(2) * t43 + t60 * t77) * MDP(32) + (g(1) * t46 - g(2) * t44 + t63 * t77) * MDP(33);];
taug  = t1;
