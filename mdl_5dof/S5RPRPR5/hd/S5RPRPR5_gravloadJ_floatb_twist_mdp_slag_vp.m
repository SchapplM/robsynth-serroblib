% Calculate Gravitation load on the joints for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:26:01
% EndTime: 2022-01-23 09:26:01
% DurationCPUTime: 0.18s
% Computational Cost: add. (173->65), mult. (213->98), div. (0->0), fcn. (209->10), ass. (0->43)
t74 = sin(pkin(8));
t95 = g(3) * t74;
t73 = qJ(3) + pkin(9);
t71 = qJ(5) + t73;
t66 = sin(t71);
t77 = sin(qJ(1));
t94 = t77 * t66;
t67 = cos(t71);
t93 = t77 * t67;
t69 = sin(t73);
t92 = t77 * t69;
t70 = cos(t73);
t91 = t77 * t70;
t76 = sin(qJ(3));
t90 = t77 * t76;
t78 = cos(qJ(3));
t89 = t77 * t78;
t79 = cos(qJ(1));
t88 = t79 * t66;
t87 = t79 * t67;
t86 = t79 * t69;
t85 = t79 * t70;
t84 = t79 * t76;
t83 = t79 * t78;
t75 = cos(pkin(8));
t51 = t75 * t94 + t87;
t52 = -t75 * t93 + t88;
t53 = -t75 * t88 + t93;
t54 = t75 * t87 + t94;
t82 = (-g(1) * t53 + g(2) * t51 + t66 * t95) * MDP(23) + (g(1) * t54 - g(2) * t52 + t67 * t95) * MDP(24);
t80 = g(1) * t79 + g(2) * t77;
t64 = g(1) * t77 - g(2) * t79;
t62 = -t75 * t84 + t89;
t60 = t75 * t90 + t83;
t68 = t76 * pkin(3) + qJ(2);
t63 = t75 * t83 + t90;
t61 = -t75 * t89 + t84;
t59 = (qJ(4) + pkin(6)) * t74 + pkin(1) + (pkin(3) * t78 + pkin(2)) * t75;
t58 = t75 * t85 + t92;
t57 = -t75 * t86 + t91;
t56 = -t75 * t91 + t86;
t55 = t75 * t92 + t85;
t1 = [(-g(1) * (-t77 * pkin(1) + t79 * qJ(2)) - g(2) * (t79 * pkin(1) + t77 * qJ(2))) * MDP(6) + (-g(1) * t61 - g(2) * t63) * MDP(12) + (-g(1) * t60 - g(2) * t62) * MDP(13) + (-g(1) * t56 - g(2) * t58) * MDP(14) + (-g(1) * t55 - g(2) * t57) * MDP(15) + (-g(1) * (-t59 * t77 + t68 * t79) - g(2) * (t59 * t79 + t68 * t77)) * MDP(17) + (-g(1) * t52 - g(2) * t54) * MDP(23) + (-g(1) * t51 - g(2) * t53) * MDP(24) + (MDP(3) - MDP(5)) * t80 + (MDP(16) * t74 + MDP(4) * t75 + MDP(2)) * t64; (-MDP(17) - MDP(6)) * t64; (g(1) * t63 - g(2) * t61 + t78 * t95) * MDP(13) + (-g(1) * t57 + g(2) * t55 + t69 * t95) * MDP(14) + (g(1) * t58 - g(2) * t56 + t70 * t95) * MDP(15) + t82 + (pkin(3) * MDP(17) + MDP(12)) * (-g(1) * t62 + g(2) * t60 + t76 * t95); (g(3) * t75 - t80 * t74) * MDP(17); t82;];
taug = t1;
