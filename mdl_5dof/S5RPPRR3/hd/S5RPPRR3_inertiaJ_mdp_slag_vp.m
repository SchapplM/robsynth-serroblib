% Calculate joint inertia matrix for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRR3_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:28:17
% EndTime: 2020-01-03 11:28:19
% DurationCPUTime: 0.21s
% Computational Cost: add. (190->49), mult. (334->73), div. (0->0), fcn. (354->8), ass. (0->36)
t62 = sin(pkin(9));
t64 = cos(pkin(9));
t67 = sin(qJ(4));
t82 = cos(qJ(4));
t54 = t62 * t82 + t67 * t64;
t53 = t67 * t62 - t64 * t82;
t76 = t53 * MDP(14);
t66 = sin(qJ(5));
t68 = cos(qJ(5));
t46 = t68 * t53 + t66 * t54;
t42 = t46 * MDP(21);
t47 = -t66 * t53 + t68 * t54;
t81 = -t47 * MDP(22) - t42;
t86 = -t54 * MDP(15) - t76 + t81;
t65 = cos(pkin(8));
t59 = -t65 * pkin(1) - pkin(2);
t55 = -t64 * pkin(3) + t59;
t85 = 0.2e1 * t53 * pkin(4) + 0.2e1 * t55;
t84 = 0.2e1 * t55;
t63 = sin(pkin(8));
t57 = t63 * pkin(1) + qJ(3);
t83 = pkin(6) + t57;
t80 = t62 ^ 2 + t64 ^ 2;
t79 = t59 * MDP(8);
t78 = t62 * MDP(6);
t77 = t64 * MDP(5);
t51 = t83 * t62;
t52 = t83 * t64;
t75 = -t82 * t51 - t67 * t52;
t74 = t80 * MDP(8);
t40 = -t54 * pkin(7) + t75;
t71 = t67 * t51 - t52 * t82;
t41 = -t53 * pkin(7) - t71;
t73 = t47 * MDP(18) - t46 * MDP(19) + (t68 * t40 - t66 * t41) * MDP(21) + (-t66 * t40 - t68 * t41) * MDP(22);
t70 = (MDP(21) * t68 - MDP(22) * t66) * pkin(4);
t1 = [t76 * t84 + t42 * t85 + MDP(1) + (t63 ^ 2 + t65 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t77 + 0.2e1 * t78 + t79) * t59 + (-0.2e1 * t53 * MDP(10) + MDP(15) * t84 + MDP(9) * t54) * t54 + (MDP(16) * t47 - 0.2e1 * t46 * MDP(17) + MDP(22) * t85) * t47 + (0.2e1 * t80 * MDP(7) + t74 * t57) * t57; 0; MDP(4) + t74; -t77 + t78 + t79 - t86; 0; MDP(8); t54 * MDP(11) - t53 * MDP(12) + MDP(14) * t75 + MDP(15) * t71 + t73; t86; 0; MDP(13) + MDP(20) + 0.2e1 * t70; t73; t81; 0; MDP(20) + t70; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
