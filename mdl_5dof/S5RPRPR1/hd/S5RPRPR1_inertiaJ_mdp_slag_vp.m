% Calculate joint inertia matrix for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR1_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:33:54
% EndTime: 2021-01-15 11:33:55
% DurationCPUTime: 0.16s
% Computational Cost: add. (258->73), mult. (417->107), div. (0->0), fcn. (429->6), ass. (0->36)
t93 = MDP(17) * pkin(3);
t70 = sin(pkin(8));
t71 = cos(pkin(8));
t73 = sin(qJ(3));
t75 = cos(qJ(3));
t63 = t70 * t75 + t71 * t73;
t67 = t73 * pkin(3) + qJ(2);
t92 = 0.2e1 * t63 * pkin(4) + 0.2e1 * t67;
t91 = pkin(3) * t70;
t90 = (pkin(1) * MDP(6));
t62 = -t70 * t73 + t71 * t75;
t72 = sin(qJ(5));
t74 = cos(qJ(5));
t52 = t72 * t62 + t74 * t63;
t53 = t74 * t62 - t72 * t63;
t89 = t53 * MDP(23) - t52 * MDP(24);
t76 = -pkin(1) - pkin(6);
t88 = -qJ(4) + t76;
t87 = t52 * MDP(23);
t68 = t71 * pkin(3) + pkin(4);
t86 = (t74 * t68 - t72 * t91) * MDP(23);
t85 = (-t72 * t68 - t74 * t91) * MDP(24);
t84 = t62 * MDP(15);
t83 = t63 * MDP(14);
t82 = t67 * MDP(17);
t81 = t62 ^ 2 + t63 ^ 2;
t64 = t88 * t73;
t65 = t88 * t75;
t54 = -t70 * t64 + t71 * t65;
t46 = -t62 * pkin(7) + t54;
t55 = t71 * t64 + t70 * t65;
t47 = -t63 * pkin(7) + t55;
t80 = t53 * MDP(20) - t52 * MDP(21) + (t74 * t46 - t72 * t47) * MDP(23) + (-t72 * t46 - t74 * t47) * MDP(24);
t79 = t54 * t62 + t55 * t63;
t78 = t62 * t71 + t63 * t70;
t1 = [MDP(1) - 0.2e1 * t79 * MDP(16) + (t54 ^ 2 + t55 ^ 2) * MDP(17) + t87 * t92 + (MDP(7) * t75 - 0.2e1 * t73 * MDP(8)) * t75 + (t82 + 0.2e1 * t83 + 0.2e1 * t84) * t67 + ((-2 * MDP(4) + t90) * pkin(1)) + (MDP(18) * t53 - 0.2e1 * t52 * MDP(19) + MDP(24) * t92) * t53 + (0.2e1 * t73 * MDP(12) + 0.2e1 * t75 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2); -t81 * MDP(16) + t79 * MDP(17) + MDP(4) - t90; t81 * MDP(17) + MDP(6); t54 * MDP(14) - t55 * MDP(15) + (t76 * MDP(12) + MDP(9)) * t75 + (-t76 * MDP(13) - MDP(10)) * t73 + (-t78 * MDP(16) + (t54 * t71 + t55 * t70) * MDP(17)) * pkin(3) + t80; t75 * MDP(12) - t73 * MDP(13) + t62 * MDP(14) - t63 * MDP(15) + t78 * t93 + t89; MDP(11) + MDP(22) + 0.2e1 * t86 + 0.2e1 * t85 + (0.2e1 * t71 * MDP(14) - 0.2e1 * t70 * MDP(15) + (t70 ^ 2 + t71 ^ 2) * t93) * pkin(3); t53 * MDP(24) + t82 + t83 + t84 + t87; 0; 0; MDP(17); t80; t89; MDP(22) + t85 + t86; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
