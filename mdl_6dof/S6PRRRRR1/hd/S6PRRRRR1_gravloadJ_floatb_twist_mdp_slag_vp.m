% Calculate Gravitation load on the joints for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:39:52
% EndTime: 2019-03-09 00:39:53
% DurationCPUTime: 0.34s
% Computational Cost: add. (448->78), mult. (585->137), div. (0->0), fcn. (694->14), ass. (0->38)
t72 = sin(pkin(6));
t96 = g(3) * t72;
t70 = qJ(3) + qJ(4);
t69 = qJ(5) + t70;
t66 = cos(t69);
t75 = sin(qJ(6));
t95 = t66 * t75;
t78 = cos(qJ(6));
t94 = t66 * t78;
t71 = sin(pkin(12));
t93 = t71 * t72;
t73 = cos(pkin(12));
t92 = t72 * t73;
t76 = sin(qJ(3));
t91 = t72 * t76;
t77 = sin(qJ(2));
t90 = t72 * t77;
t79 = cos(qJ(3));
t89 = t72 * t79;
t74 = cos(pkin(6));
t88 = t74 * t77;
t80 = cos(qJ(2));
t87 = t74 * t80;
t86 = t75 * t80;
t85 = t78 * t80;
t61 = t71 * t80 + t73 * t88;
t65 = sin(t69);
t55 = t61 * t66 - t65 * t92;
t63 = -t71 * t88 + t73 * t80;
t57 = t63 * t66 + t65 * t93;
t59 = t74 * t65 + t66 * t90;
t84 = (g(1) * t57 + g(2) * t55 + g(3) * t59) * MDP(25) + (-t78 * MDP(31) + t75 * MDP(32) - MDP(24)) * (g(1) * (-t63 * t65 + t66 * t93) + g(2) * (-t61 * t65 - t66 * t92) + g(3) * (-t65 * t90 + t74 * t66));
t67 = sin(t70);
t68 = cos(t70);
t83 = (-g(1) * (-t63 * t67 + t68 * t93) - g(2) * (-t61 * t67 - t68 * t92) - g(3) * (-t67 * t90 + t74 * t68)) * MDP(17) + (-g(1) * (-t63 * t68 - t67 * t93) - g(2) * (-t61 * t68 + t67 * t92) - g(3) * (-t74 * t67 - t68 * t90)) * MDP(18) + t84;
t62 = t71 * t87 + t73 * t77;
t60 = t71 * t77 - t73 * t87;
t1 = [-g(3) * MDP(1); (g(1) * t63 + g(2) * t61 + g(3) * t90) * MDP(4) + (-g(1) * (-t62 * t94 + t63 * t75) - g(2) * (-t60 * t94 + t61 * t75) - (t66 * t85 + t75 * t77) * t96) * MDP(31) + (-g(1) * (t62 * t95 + t63 * t78) - g(2) * (t60 * t95 + t61 * t78) - (-t66 * t86 + t77 * t78) * t96) * MDP(32) + (-MDP(10) * t79 + MDP(11) * t76 - MDP(17) * t68 + MDP(18) * t67 - t66 * MDP(24) + MDP(25) * t65 - MDP(3)) * (-g(1) * t62 - g(2) * t60 + t80 * t96); (-g(1) * (-t63 * t76 + t71 * t89) - g(2) * (-t61 * t76 - t73 * t89) - g(3) * (t74 * t79 - t76 * t90)) * MDP(10) + (-g(1) * (-t63 * t79 - t71 * t91) - g(2) * (-t61 * t79 + t73 * t91) - g(3) * (-t74 * t76 - t77 * t89)) * MDP(11) + t83; t83; t84; (-g(1) * (-t57 * t75 + t62 * t78) - g(2) * (-t55 * t75 + t60 * t78) - g(3) * (-t59 * t75 - t72 * t85)) * MDP(31) + (-g(1) * (-t57 * t78 - t62 * t75) - g(2) * (-t55 * t78 - t60 * t75) - g(3) * (-t59 * t78 + t72 * t86)) * MDP(32);];
taug  = t1;
