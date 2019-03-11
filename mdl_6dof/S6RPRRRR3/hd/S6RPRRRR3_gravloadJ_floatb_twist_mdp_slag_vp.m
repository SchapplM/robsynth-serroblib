% Calculate Gravitation load on the joints for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:35
% EndTime: 2019-03-09 07:02:36
% DurationCPUTime: 0.23s
% Computational Cost: add. (328->62), mult. (274->99), div. (0->0), fcn. (276->12), ass. (0->39)
t71 = sin(qJ(3));
t95 = t71 * MDP(11);
t69 = qJ(4) + qJ(5);
t67 = qJ(6) + t69;
t61 = sin(t67);
t62 = cos(t67);
t65 = sin(t69);
t66 = cos(t69);
t70 = sin(qJ(4));
t73 = cos(qJ(4));
t94 = t73 * MDP(17) - t70 * MDP(18) + t66 * MDP(24) - t65 * MDP(25) + t62 * MDP(31) - t61 * MDP(32) + MDP(10);
t93 = g(3) * t71;
t68 = qJ(1) + pkin(11);
t63 = sin(t68);
t74 = cos(qJ(3));
t92 = t63 * t74;
t64 = cos(t68);
t91 = t64 * t74;
t90 = t65 * t74;
t89 = t66 * t74;
t88 = t70 * t74;
t87 = t73 * t74;
t49 = t61 * t92 + t64 * t62;
t50 = t64 * t61 - t62 * t92;
t51 = -t61 * t91 + t63 * t62;
t52 = t63 * t61 + t62 * t91;
t86 = (-g(1) * t51 + g(2) * t49 + t61 * t93) * MDP(31) + (g(1) * t52 - g(2) * t50 + t62 * t93) * MDP(32);
t53 = t63 * t90 + t64 * t66;
t54 = -t63 * t89 + t64 * t65;
t55 = t63 * t66 - t64 * t90;
t56 = t63 * t65 + t64 * t89;
t79 = (-g(1) * t55 + g(2) * t53 + t65 * t93) * MDP(24) + (g(1) * t56 - g(2) * t54 + t66 * t93) * MDP(25) + t86;
t75 = cos(qJ(1));
t72 = sin(qJ(1));
t60 = t63 * t70 + t64 * t87;
t59 = t63 * t73 - t64 * t88;
t58 = -t63 * t87 + t64 * t70;
t57 = t63 * t88 + t64 * t73;
t1 = [(g(1) * t75 + g(2) * t72) * MDP(3) + (-g(1) * t58 - g(2) * t60) * MDP(17) + (-g(1) * t57 - g(2) * t59) * MDP(18) + (-g(1) * t54 - g(2) * t56) * MDP(24) + (-g(1) * t53 - g(2) * t55) * MDP(25) + (-g(1) * t50 - g(2) * t52) * MDP(31) + (-g(1) * t49 - g(2) * t51) * MDP(32) + (t74 * MDP(10) - t95) * (g(1) * t63 - g(2) * t64) + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t72 - g(2) * t75); -g(3) * MDP(4); (-t94 * t74 + t95) * g(3) + (MDP(11) * t74 + t94 * t71) * (g(1) * t64 + g(2) * t63); (-g(1) * t59 + g(2) * t57 + t70 * t93) * MDP(17) + (g(1) * t60 - g(2) * t58 + t73 * t93) * MDP(18) + t79; t79; t86;];
taug  = t1;
