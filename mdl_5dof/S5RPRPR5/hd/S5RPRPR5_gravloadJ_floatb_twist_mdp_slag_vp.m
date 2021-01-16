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
% Datum: 2021-01-15 11:56
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
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
% StartTime: 2021-01-15 11:55:14
% EndTime: 2021-01-15 11:55:15
% DurationCPUTime: 0.20s
% Computational Cost: add. (173->65), mult. (213->98), div. (0->0), fcn. (209->10), ass. (0->43)
t74 = sin(pkin(8));
t75 = cos(pkin(8));
t79 = cos(qJ(3));
t101 = (pkin(3) * t79 + pkin(2)) * t75 + (qJ(4) + pkin(6)) * t74 + pkin(1);
t100 = g(1) * t74;
t73 = qJ(3) + pkin(9);
t72 = qJ(5) + t73;
t67 = sin(t72);
t78 = sin(qJ(1));
t97 = t78 * t67;
t68 = cos(t72);
t96 = t78 * t68;
t70 = sin(t73);
t95 = t78 * t70;
t71 = cos(t73);
t94 = t78 * t71;
t77 = sin(qJ(3));
t93 = t78 * t77;
t92 = t78 * t79;
t80 = cos(qJ(1));
t91 = t80 * t67;
t90 = t80 * t68;
t89 = t80 * t70;
t88 = t80 * t71;
t87 = t80 * t77;
t86 = t80 * t79;
t53 = -t75 * t97 - t90;
t54 = t75 * t96 - t91;
t55 = t75 * t91 - t96;
t56 = t75 * t90 + t97;
t85 = (-g(2) * t53 - g(3) * t55 + t67 * t100) * MDP(23) + (g(2) * t54 - g(3) * t56 + t68 * t100) * MDP(24);
t66 = g(2) * t80 + g(3) * t78;
t81 = g(2) * t78 - g(3) * t80;
t63 = t75 * t87 - t92;
t61 = -t75 * t93 - t86;
t69 = -t77 * pkin(3) - qJ(2);
t64 = t75 * t86 + t93;
t62 = t75 * t92 - t87;
t60 = t75 * t88 + t95;
t59 = t75 * t89 - t94;
t58 = t75 * t94 - t89;
t57 = -t75 * t95 - t88;
t1 = [(-g(2) * (t80 * pkin(1) + t78 * qJ(2)) - g(3) * (t78 * pkin(1) - t80 * qJ(2))) * MDP(6) + (-g(2) * t64 - g(3) * t62) * MDP(12) + (g(2) * t63 - g(3) * t61) * MDP(13) + (-g(2) * t60 - g(3) * t58) * MDP(14) + (g(2) * t59 - g(3) * t57) * MDP(15) + (-g(2) * (t101 * t80 - t69 * t78) - g(3) * (t101 * t78 + t69 * t80)) * MDP(17) + (-g(2) * t56 - g(3) * t54) * MDP(23) + (g(2) * t55 - g(3) * t53) * MDP(24) + (MDP(3) - MDP(5)) * t81 + (-t74 * MDP(16) - t75 * MDP(4) - MDP(2)) * t66; (MDP(17) + MDP(6)) * t66; (g(2) * t62 - g(3) * t64 + t79 * t100) * MDP(13) + (-g(2) * t57 - g(3) * t59 + t70 * t100) * MDP(14) + (g(2) * t58 - g(3) * t60 + t71 * t100) * MDP(15) + t85 + (pkin(3) * MDP(17) + MDP(12)) * (-g(2) * t61 - g(3) * t63 + t77 * t100); (g(1) * t75 - t81 * t74) * MDP(17); t85;];
taug = t1;
