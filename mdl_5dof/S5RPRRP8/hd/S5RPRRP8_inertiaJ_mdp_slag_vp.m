% Calculate joint inertia matrix for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP8_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:30
% EndTime: 2019-12-31 18:47:31
% DurationCPUTime: 0.21s
% Computational Cost: add. (223->74), mult. (325->90), div. (0->0), fcn. (212->4), ass. (0->37)
t52 = sin(qJ(4));
t49 = t52 ^ 2;
t54 = cos(qJ(4));
t77 = t54 ^ 2 + t49;
t64 = t77 * MDP(20);
t84 = pkin(7) * t64;
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t56 = -pkin(1) - pkin(2);
t43 = t55 * qJ(2) + t53 * t56;
t41 = -pkin(7) + t43;
t62 = t41 * t64;
t83 = t77 * MDP(18);
t68 = MDP(16) - MDP(19);
t82 = -t68 * t52 + (MDP(15) + MDP(17)) * t54 + MDP(8);
t81 = 0.2e1 * MDP(17);
t80 = 2 * MDP(19);
t42 = t53 * qJ(2) - t55 * t56;
t40 = pkin(3) + t42;
t79 = pkin(3) + t40;
t44 = -t54 * pkin(4) - t52 * qJ(5) - pkin(3);
t37 = t42 - t44;
t78 = -t37 + t44;
t76 = t42 * MDP(8);
t75 = t43 * MDP(9);
t74 = t49 * MDP(10) + MDP(7);
t73 = t37 * MDP(20);
t72 = t53 * t83;
t71 = t44 * MDP(20);
t70 = t54 * MDP(11);
t67 = 0.2e1 * t52 * t70 + t74;
t63 = -MDP(20) * pkin(4) - MDP(17);
t61 = t52 * t80 + t54 * t81;
t60 = 0.2e1 * t54 * MDP(15) - 0.2e1 * t52 * MDP(16);
t59 = -t52 * MDP(12) - t54 * MDP(13) - (-t52 * pkin(4) + t54 * qJ(5)) * MDP(18);
t58 = (MDP(20) * qJ(5) - t68) * t54 + (-MDP(15) + t63) * t52;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + t40 * t60 + (t61 + t73) * t37 + 0.2e1 * t76 + 0.2e1 * t75 + t67 + (-0.2e1 * t83 + t62) * t41; -t72 - pkin(1) * MDP(6) - MDP(4) + (MDP(9) + t62) * t53 + (-t73 - t82) * t55; MDP(6) + (t77 * t53 ^ 2 + t55 ^ 2) * MDP(20); t37 * t71 - t76 - t75 + t41 * t83 + (-t83 + t62) * pkin(7) + (-t79 * MDP(15) + t78 * MDP(17)) * t54 + (t79 * MDP(16) + t78 * MDP(19) - 0.2e1 * t70) * t52 - t74; t72 + (-MDP(9) + t84) * t53 + (-t71 + t82) * t55; (-t61 + t71) * t44 + pkin(3) * t60 + t67 + (0.2e1 * t83 + t84) * pkin(7); t58 * t41 + t59; t58 * t53; t58 * pkin(7) - t59; MDP(14) + pkin(4) * t81 + qJ(5) * t80 + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(20); (t41 * MDP(20) - MDP(18)) * t52; t52 * t53 * MDP(20); (MDP(20) * pkin(7) + MDP(18)) * t52; t63; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
