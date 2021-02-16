% Calculate joint inertia matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:45
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR4_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:44:16
% EndTime: 2021-01-15 11:44:17
% DurationCPUTime: 0.19s
% Computational Cost: add. (262->67), mult. (462->108), div. (0->0), fcn. (477->8), ass. (0->37)
t67 = sin(pkin(9));
t69 = cos(pkin(9));
t72 = sin(qJ(3));
t74 = cos(qJ(3));
t61 = t67 * t74 + t69 * t72;
t81 = t61 * MDP(13);
t59 = t67 * t72 - t69 * t74;
t82 = t59 * MDP(12);
t71 = sin(qJ(5));
t73 = cos(qJ(5));
t50 = t73 * t59 + t71 * t61;
t46 = t50 * MDP(21);
t51 = -t71 * t59 + t73 * t61;
t86 = -t51 * MDP(22) - t46;
t90 = t81 + t82 - t86;
t89 = MDP(15) * pkin(3);
t70 = cos(pkin(8));
t66 = -t70 * pkin(1) - pkin(2);
t62 = -t74 * pkin(3) + t66;
t88 = 0.2e1 * t59 * pkin(4) + 0.2e1 * t62;
t87 = pkin(3) * t67;
t68 = sin(pkin(8));
t64 = t68 * pkin(1) + pkin(6);
t85 = qJ(4) + t64;
t65 = t69 * pkin(3) + pkin(4);
t84 = (t73 * t65 - t71 * t87) * MDP(21);
t83 = (-t71 * t65 - t73 * t87) * MDP(22);
t80 = t62 * MDP(15);
t79 = t74 * MDP(10);
t57 = t85 * t72;
t58 = t85 * t74;
t44 = -t69 * t57 - t67 * t58;
t42 = -t61 * pkin(7) + t44;
t45 = -t67 * t57 + t69 * t58;
t43 = -t59 * pkin(7) + t45;
t78 = t51 * MDP(18) - t50 * MDP(19) + (t73 * t42 - t71 * t43) * MDP(21) + (-t71 * t42 - t73 * t43) * MDP(22);
t1 = [MDP(1) - 0.2e1 * t66 * t79 + 0.2e1 * (-t44 * t61 - t45 * t59) * MDP(14) + (t44 ^ 2 + t45 ^ 2) * MDP(15) + t46 * t88 + (t68 ^ 2 + t70 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t80 + 0.2e1 * t81 + 0.2e1 * t82) * t62 + (0.2e1 * t66 * MDP(11) + MDP(5) * t72 + 0.2e1 * t74 * MDP(6)) * t72 + (MDP(16) * t51 - 0.2e1 * t50 * MDP(17) + MDP(22) * t88) * t51; (-t44 * t59 + t45 * t61) * MDP(15); MDP(4) + (t59 ^ 2 + t61 ^ 2) * MDP(15); t44 * MDP(12) - t45 * MDP(13) + t72 * MDP(7) + t74 * MDP(8) + (-t72 * MDP(10) - t74 * MDP(11)) * t64 + ((-t59 * t67 - t61 * t69) * MDP(14) + (t44 * t69 + t45 * t67) * MDP(15)) * pkin(3) + t78; t79 - t72 * MDP(11) + (-t59 * t69 + t61 * t67) * t89 - t90; MDP(9) + MDP(20) + 0.2e1 * t84 + 0.2e1 * t83 + (0.2e1 * t69 * MDP(12) - 0.2e1 * t67 * MDP(13) + (t67 ^ 2 + t69 ^ 2) * t89) * pkin(3); t80 + t90; 0; 0; MDP(15); t78; t86; MDP(20) + t83 + t84; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
