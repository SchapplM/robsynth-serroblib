% Calculate Gravitation load on the joints for
% S6RPRRPR2
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
%   see S6RPRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:02:34
% EndTime: 2019-03-09 05:02:35
% DurationCPUTime: 0.26s
% Computational Cost: add. (259->61), mult. (246->92), div. (0->0), fcn. (227->10), ass. (0->37)
t68 = qJ(1) + pkin(10);
t65 = sin(t68);
t66 = cos(t68);
t83 = g(1) * t66 + g(2) * t65;
t98 = MDP(11) - MDP(19);
t71 = sin(qJ(3));
t74 = cos(qJ(3));
t55 = -g(3) * t74 + t83 * t71;
t92 = g(3) * t71;
t90 = t65 * t74;
t89 = t66 * t74;
t70 = sin(qJ(4));
t88 = t70 * t74;
t73 = cos(qJ(4));
t87 = t73 * t74;
t67 = qJ(4) + pkin(11) + qJ(6);
t62 = sin(t67);
t63 = cos(t67);
t51 = t62 * t90 + t66 * t63;
t52 = t66 * t62 - t63 * t90;
t53 = -t62 * t89 + t65 * t63;
t54 = t65 * t62 + t63 * t89;
t86 = (-g(1) * t53 + g(2) * t51 + t62 * t92) * MDP(26) + (g(1) * t54 - g(2) * t52 + t63 * t92) * MDP(27);
t84 = pkin(4) * t70 + pkin(7);
t72 = sin(qJ(1));
t75 = cos(qJ(1));
t81 = g(1) * t72 - g(2) * t75;
t64 = t73 * pkin(4) + pkin(3);
t69 = -qJ(5) - pkin(8);
t80 = t74 * t64 - t71 * t69;
t78 = pkin(2) + t80;
t59 = t65 * t73 - t66 * t88;
t57 = t65 * t88 + t66 * t73;
t77 = t81 * pkin(1);
t60 = t65 * t70 + t66 * t87;
t58 = -t65 * t87 + t66 * t70;
t1 = [t81 * MDP(2) + (g(1) * t75 + g(2) * t72) * MDP(3) + MDP(4) * t77 + (-g(1) * t58 - g(2) * t60) * MDP(17) + (-g(1) * t57 - g(2) * t59) * MDP(18) + (t77 + (-g(1) * t84 - g(2) * t78) * t66 + (g(1) * t78 - g(2) * t84) * t65) * MDP(20) + (-g(1) * t52 - g(2) * t54) * MDP(26) + (-g(1) * t51 - g(2) * t53) * MDP(27) + (t74 * MDP(10) - t98 * t71) * (g(1) * t65 - g(2) * t66); (-MDP(20) - MDP(4)) * g(3); (-g(3) * t80 + t83 * (t64 * t71 + t69 * t74)) * MDP(20) + t98 * (t83 * t74 + t92) + (MDP(17) * t73 - MDP(18) * t70 + MDP(26) * t63 - MDP(27) * t62 + MDP(10)) * t55; (g(1) * t60 - g(2) * t58 + t73 * t92) * MDP(18) + t86 + (pkin(4) * MDP(20) + MDP(17)) * (-g(1) * t59 + g(2) * t57 + t70 * t92); -t55 * MDP(20); t86;];
taug  = t1;
