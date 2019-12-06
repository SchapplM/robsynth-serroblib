% Calculate Gravitation load on the joints for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:48
% EndTime: 2019-12-05 16:32:51
% DurationCPUTime: 0.51s
% Computational Cost: add. (218->77), mult. (492->133), div. (0->0), fcn. (579->12), ass. (0->40)
t99 = MDP(11) - MDP(14);
t73 = sin(qJ(2));
t75 = cos(qJ(2));
t86 = cos(pkin(9));
t87 = cos(pkin(5));
t81 = t87 * t86;
t85 = sin(pkin(9));
t55 = t73 * t85 - t75 * t81;
t80 = t87 * t85;
t57 = t73 * t86 + t75 * t80;
t98 = -g(1) * t57 - g(2) * t55;
t70 = sin(pkin(5));
t95 = g(3) * t70;
t68 = pkin(10) + qJ(5);
t66 = sin(t68);
t74 = cos(qJ(3));
t94 = t66 * t74;
t67 = cos(t68);
t93 = t67 * t74;
t69 = sin(pkin(10));
t92 = t69 * t74;
t91 = t70 * t73;
t90 = t70 * t75;
t71 = cos(pkin(10));
t89 = t71 * t74;
t88 = t74 * t75;
t84 = t70 * t86;
t83 = t70 * t85;
t56 = t73 * t81 + t75 * t85;
t58 = -t73 * t80 + t75 * t86;
t82 = -g(1) * t58 - g(2) * t56;
t72 = sin(qJ(3));
t51 = t56 * t72 + t74 * t84;
t53 = t58 * t72 - t74 * t83;
t59 = t72 * t91 - t74 * t87;
t78 = g(1) * t53 + g(2) * t51 + g(3) * t59;
t60 = t72 * t87 + t74 * t91;
t54 = t58 * t74 + t72 * t83;
t52 = t56 * t74 - t72 * t84;
t1 = [(-MDP(1) - MDP(15)) * g(3); (g(3) * t91 - t82) * MDP(4) + (-g(1) * (-t57 * t89 + t58 * t69) - g(2) * (-t55 * t89 + t56 * t69) - (t69 * t73 + t71 * t88) * t95) * MDP(12) + (-g(1) * (t57 * t92 + t58 * t71) - g(2) * (t55 * t92 + t56 * t71) - (-t69 * t88 + t71 * t73) * t95) * MDP(13) + ((-t73 * t95 + t82) * pkin(7) + (-t75 * t95 - t98) * (pkin(3) * t74 + qJ(4) * t72 + pkin(2))) * MDP(15) + (-g(1) * (-t57 * t93 + t58 * t66) - g(2) * (-t55 * t93 + t56 * t66) - (t66 * t73 + t67 * t88) * t95) * MDP(21) + (-g(1) * (t57 * t94 + t58 * t67) - g(2) * (t55 * t94 + t56 * t67) - (-t66 * t88 + t67 * t73) * t95) * MDP(22) + (-t74 * MDP(10) + t99 * t72 - MDP(3)) * (g(3) * t90 + t98); (-g(1) * (-pkin(3) * t53 + qJ(4) * t54) - g(2) * (-pkin(3) * t51 + qJ(4) * t52) - g(3) * (-pkin(3) * t59 + qJ(4) * t60)) * MDP(15) + t99 * (g(1) * t54 + g(2) * t52 + g(3) * t60) + (MDP(12) * t71 - MDP(13) * t69 + MDP(21) * t67 - MDP(22) * t66 + MDP(10)) * t78; -t78 * MDP(15); (-g(1) * (-t54 * t66 + t57 * t67) - g(2) * (-t52 * t66 + t55 * t67) - g(3) * (-t60 * t66 - t67 * t90)) * MDP(21) + (-g(1) * (-t54 * t67 - t57 * t66) - g(2) * (-t52 * t67 - t55 * t66) - g(3) * (-t60 * t67 + t66 * t90)) * MDP(22);];
taug = t1;
