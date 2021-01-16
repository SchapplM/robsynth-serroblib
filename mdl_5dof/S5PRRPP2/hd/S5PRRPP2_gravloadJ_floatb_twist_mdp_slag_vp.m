% Calculate Gravitation load on the joints for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:32:34
% EndTime: 2021-01-15 15:32:36
% DurationCPUTime: 0.31s
% Computational Cost: add. (172->64), mult. (272->93), div. (0->0), fcn. (252->8), ass. (0->39)
t67 = qJ(3) + pkin(8);
t64 = sin(t67);
t65 = cos(t67);
t97 = pkin(4) * t65 + qJ(5) * t64;
t68 = sin(pkin(7));
t69 = cos(pkin(7));
t78 = g(1) * t69 + g(2) * t68;
t96 = MDP(12) + MDP(16);
t95 = MDP(13) - MDP(18);
t72 = sin(qJ(2));
t74 = cos(qJ(2));
t55 = -g(3) * t74 + t78 * t72;
t91 = g(3) * t72;
t73 = cos(qJ(3));
t63 = t73 * pkin(3) + pkin(2);
t89 = t63 * t72;
t88 = t68 * t73;
t87 = t68 * t74;
t86 = t69 * t73;
t85 = t69 * t74;
t70 = qJ(4) + pkin(6);
t84 = t70 * t74;
t71 = sin(qJ(3));
t83 = t71 * t74;
t82 = t73 * t74;
t80 = -MDP(15) - MDP(19);
t79 = t69 * t83;
t77 = -t68 * t83 - t86;
t51 = t64 * t87 + t69 * t65;
t53 = t64 * t85 - t68 * t65;
t48 = g(1) * t53 + g(2) * t51 + t64 * t91;
t56 = t78 * t74 + t91;
t62 = pkin(3) * t88;
t59 = t74 * t63;
t58 = t69 * t84;
t57 = t68 * t84;
t54 = t68 * t64 + t65 * t85;
t52 = -t69 * t64 + t65 * t87;
t1 = [(-MDP(1) + t80) * g(3); (-g(1) * (-t69 * t89 + t58) - g(2) * (-t68 * t89 + t57) - g(3) * (t72 * t70 + t59)) * MDP(15) + (-g(1) * t58 - g(2) * t57 - g(3) * (t97 * t74 + t59) + (-g(3) * t70 + t78 * (t63 + t97)) * t72) * MDP(19) + (MDP(4) - MDP(14) - MDP(17)) * t56 + (t73 * MDP(10) - t71 * MDP(11) - t95 * t64 + t96 * t65 + MDP(3)) * t55; (-g(1) * (-t79 + t88) - g(2) * t77 + t71 * t91) * MDP(10) + (-g(1) * (-t68 * t71 - t69 * t82) - g(2) * (-t68 * t82 + t69 * t71) + t73 * t91) * MDP(11) + (-g(1) * t62 + (g(2) * t86 + t56 * t71) * pkin(3)) * MDP(15) + (-g(1) * (-pkin(3) * t79 - t53 * pkin(4) + t54 * qJ(5) + t62) - g(2) * (pkin(3) * t77 - t51 * pkin(4) + t52 * qJ(5)) - (-pkin(3) * t71 - pkin(4) * t64 + qJ(5) * t65) * t91) * MDP(19) + t95 * (g(1) * t54 + g(2) * t52 + t65 * t91) + t96 * t48; t80 * t55; -t48 * MDP(19);];
taug = t1;
