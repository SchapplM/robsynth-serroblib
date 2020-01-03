% Calculate Gravitation load on the joints for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:09
% EndTime: 2019-12-31 22:02:10
% DurationCPUTime: 0.28s
% Computational Cost: add. (179->62), mult. (265->93), div. (0->0), fcn. (247->8), ass. (0->36)
t72 = sin(qJ(1));
t75 = cos(qJ(1));
t81 = g(1) * t75 + g(2) * t72;
t97 = MDP(10) - MDP(25);
t69 = qJ(3) + qJ(4);
t65 = sin(t69);
t66 = cos(t69);
t86 = t75 * t66;
t74 = cos(qJ(2));
t88 = t72 * t74;
t50 = t65 * t88 + t86;
t87 = t75 * t65;
t52 = t66 * t72 - t74 * t87;
t71 = sin(qJ(2));
t91 = g(3) * t71;
t96 = -g(1) * t52 + g(2) * t50 + t65 * t91;
t54 = -g(3) * t74 + t71 * t81;
t70 = sin(qJ(3));
t62 = pkin(3) * t70 + pkin(4) * t65;
t89 = pkin(6) + t62;
t85 = t75 * t70;
t73 = cos(qJ(3));
t84 = t75 * t73;
t51 = -t66 * t88 + t87;
t53 = t65 * t72 + t74 * t86;
t83 = t96 * MDP(23) + (g(1) * t53 - g(2) * t51 + t66 * t91) * MDP(24);
t63 = t73 * pkin(3) + pkin(4) * t66;
t61 = pkin(2) + t63;
t68 = -qJ(5) - pkin(8) - pkin(7);
t79 = t61 * t74 - t68 * t71;
t77 = pkin(1) + t79;
t59 = t70 * t72 + t74 * t84;
t58 = t72 * t73 - t74 * t85;
t57 = -t73 * t88 + t85;
t56 = t70 * t88 + t84;
t1 = [t81 * MDP(3) + (-g(1) * t57 - g(2) * t59) * MDP(16) + (-g(1) * t56 - g(2) * t58) * MDP(17) + (-g(1) * t51 - g(2) * t53) * MDP(23) + (-g(1) * t50 - g(2) * t52) * MDP(24) + ((-g(1) * t89 - g(2) * t77) * t75 + (g(1) * t77 - g(2) * t89) * t72) * MDP(26) + (t74 * MDP(9) - t71 * t97 + MDP(2)) * (g(1) * t72 - g(2) * t75); (-g(3) * t79 + t81 * (t61 * t71 + t68 * t74)) * MDP(26) + t97 * (t74 * t81 + t91) + (MDP(16) * t73 - MDP(17) * t70 + MDP(23) * t66 - MDP(24) * t65 + MDP(9)) * t54; (-g(1) * t58 + g(2) * t56 + t70 * t91) * MDP(16) + (g(1) * t59 - g(2) * t57 + t73 * t91) * MDP(17) + (-g(1) * (-t62 * t74 * t75 + t63 * t72) - g(2) * (-t62 * t88 - t63 * t75) + t62 * t91) * MDP(26) + t83; MDP(26) * pkin(4) * t96 + t83; -t54 * MDP(26);];
taug = t1;
