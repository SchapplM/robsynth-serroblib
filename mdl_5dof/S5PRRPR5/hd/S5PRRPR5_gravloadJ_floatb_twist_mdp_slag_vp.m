% Calculate Gravitation load on the joints for
% S5PRRPR5
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
%   see S5PRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:44
% EndTime: 2021-01-15 16:04:46
% DurationCPUTime: 0.41s
% Computational Cost: add. (203->74), mult. (392->125), div. (0->0), fcn. (450->12), ass. (0->37)
t65 = sin(pkin(5));
t93 = g(3) * t65;
t92 = cos(qJ(2));
t63 = qJ(3) + pkin(10);
t62 = cos(t63);
t68 = sin(qJ(5));
t91 = t62 * t68;
t71 = cos(qJ(5));
t90 = t62 * t71;
t64 = sin(pkin(9));
t89 = t64 * t65;
t70 = sin(qJ(2));
t88 = t64 * t70;
t87 = t65 * t70;
t72 = cos(qJ(3));
t86 = t65 * t72;
t85 = cos(pkin(9));
t84 = t64 * t92;
t83 = t68 * t92;
t82 = t71 * t92;
t81 = t65 * t85;
t80 = t85 * t70;
t79 = t85 * t92;
t66 = cos(pkin(5));
t54 = t66 * t88 - t79;
t61 = sin(t63);
t51 = -t54 * t62 + t61 * t89;
t55 = -t66 * t79 + t88;
t57 = t66 * t84 + t80;
t73 = g(1) * t57 + g(2) * t55 - t92 * t93;
t69 = sin(qJ(3));
t67 = qJ(4) + pkin(7);
t60 = t72 * pkin(3) + pkin(2);
t56 = t66 * t80 + t84;
t53 = t66 * t61 + t62 * t87;
t49 = t56 * t62 - t61 * t81;
t1 = [(-MDP(1) - MDP(15)) * g(3); (-g(1) * (-t54 * t67 - t57 * t60) - g(2) * (-t55 * t60 + t56 * t67) - (t92 * t60 + t67 * t70) * t93) * MDP(15) + (-g(1) * (-t54 * t68 - t57 * t90) - g(2) * (-t55 * t90 + t56 * t68) - (t62 * t82 + t68 * t70) * t93) * MDP(21) + (-g(1) * (-t54 * t71 + t57 * t91) - g(2) * (t55 * t91 + t56 * t71) - (-t62 * t83 + t70 * t71) * t93) * MDP(22) + (MDP(4) - MDP(14)) * (-g(1) * t54 + g(2) * t56 + g(3) * t87) + (t72 * MDP(10) - t69 * MDP(11) + t62 * MDP(12) - t61 * MDP(13) + MDP(3)) * t73; (-g(1) * (t54 * t72 - t69 * t89) - g(2) * (-t56 * t72 + t69 * t81) - g(3) * (-t66 * t69 - t70 * t86)) * MDP(11) + (g(1) * t51 + g(2) * t49 + g(3) * t53) * MDP(13) + (-MDP(21) * t71 + MDP(22) * t68 - MDP(12)) * (g(1) * (t54 * t61 + t62 * t89) + g(2) * (-t56 * t61 - t62 * t81) + g(3) * (-t61 * t87 + t66 * t62)) + (pkin(3) * MDP(15) + MDP(10)) * (-g(1) * (t54 * t69 + t64 * t86) - g(2) * (-t56 * t69 - t72 * t81) - g(3) * (t66 * t72 - t69 * t87)); -t73 * MDP(15); (-g(1) * (-t51 * t68 + t57 * t71) - g(2) * (-t49 * t68 + t55 * t71) - g(3) * (-t53 * t68 - t65 * t82)) * MDP(21) + (-g(1) * (-t51 * t71 - t57 * t68) - g(2) * (-t49 * t71 - t55 * t68) - g(3) * (-t53 * t71 + t65 * t83)) * MDP(22);];
taug = t1;
