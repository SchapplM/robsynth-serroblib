% Calculate joint inertia matrix for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRP9_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:12
% EndTime: 2019-12-31 18:49:13
% DurationCPUTime: 0.23s
% Computational Cost: add. (407->83), mult. (730->116), div. (0->0), fcn. (779->6), ass. (0->38)
t79 = MDP(20) + MDP(22);
t90 = MDP(21) - MDP(24);
t68 = cos(pkin(8));
t61 = -t68 * pkin(2) - pkin(1);
t89 = 0.2e1 * t61;
t88 = 2 * MDP(22);
t87 = 2 * MDP(24);
t86 = pkin(1) * MDP(7);
t85 = pkin(6) + qJ(2);
t67 = sin(pkin(8));
t83 = t67 * MDP(5);
t82 = t68 * MDP(4);
t70 = sin(qJ(3));
t72 = cos(qJ(3));
t55 = t70 * t67 - t72 * t68;
t81 = t55 * MDP(13);
t69 = sin(qJ(4));
t80 = t69 * MDP(21);
t57 = t85 * t67;
t58 = t85 * t68;
t78 = -t72 * t57 - t70 * t58;
t77 = pkin(4) * t88 + MDP(19);
t76 = t70 * t57 - t72 * t58;
t52 = t55 * pkin(3) + t61;
t47 = -t55 * pkin(7) - t76;
t71 = cos(qJ(4));
t56 = t72 * t67 + t70 * t68;
t74 = -t56 * pkin(7) + t78;
t44 = t69 * t47 - t71 * t74;
t45 = t71 * t47 + t69 * t74;
t50 = t71 * t55 + t69 * t56;
t51 = -t69 * t55 + t71 * t56;
t75 = t51 * MDP(17) - t50 * MDP(18) - t79 * t44 - t90 * t45;
t64 = t69 * pkin(3);
t62 = t71 * pkin(3) + pkin(4);
t60 = t64 + qJ(5);
t43 = t50 * pkin(4) - t51 * qJ(5) + t52;
t1 = [MDP(1) + t81 * t89 + (t43 ^ 2 + t44 ^ 2 + t45 ^ 2) * MDP(25) + (0.2e1 * t82 - 0.2e1 * t83 + t86) * pkin(1) + 0.2e1 * (t52 * MDP(20) + t43 * MDP(22) - MDP(23) * t45) * t50 + (MDP(15) * t51 - 0.2e1 * t50 * MDP(16) + 0.2e1 * t52 * MDP(21) + 0.2e1 * MDP(23) * t44 - 0.2e1 * t43 * MDP(24)) * t51 + (MDP(14) * t89 + MDP(8) * t56 - 0.2e1 * t55 * MDP(9)) * t56 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t67 ^ 2 + t68 ^ 2) * qJ(2); t56 * MDP(14) + t43 * MDP(25) + t79 * t50 + t90 * t51 + t81 - t82 + t83 - t86; MDP(7) + MDP(25); t56 * MDP(10) - t55 * MDP(11) + t78 * MDP(13) + t76 * MDP(14) + (-t60 * t50 - t62 * t51) * MDP(23) + (-t44 * t62 + t45 * t60) * MDP(25) + t75; 0; MDP(12) + MDP(19) + (t60 ^ 2 + t62 ^ 2) * MDP(25) + 0.2e1 * (t71 * MDP(20) - t80) * pkin(3) + t62 * t88 + t60 * t87; (-pkin(4) * t51 - t50 * qJ(5)) * MDP(23) + (-t44 * pkin(4) + t45 * qJ(5)) * MDP(25) + t75; 0; (0.2e1 * qJ(5) + t64) * MDP(24) + (t62 * pkin(4) + t60 * qJ(5)) * MDP(25) + (t79 * t71 - t80) * pkin(3) + t77; qJ(5) * t87 + ((pkin(4) ^ 2) + qJ(5) ^ 2) * MDP(25) + t77; t51 * MDP(23) + t44 * MDP(25); 0; -t62 * MDP(25) - MDP(22); -MDP(25) * pkin(4) - MDP(22); MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
