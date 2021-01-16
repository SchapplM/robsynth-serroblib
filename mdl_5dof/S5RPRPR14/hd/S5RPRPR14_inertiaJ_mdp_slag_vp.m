% Calculate joint inertia matrix for
% S5RPRPR14
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
%   see S5RPRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR14_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR14_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:06
% EndTime: 2021-01-15 12:17:07
% DurationCPUTime: 0.20s
% Computational Cost: add. (267->84), mult. (457->129), div. (0->0), fcn. (451->6), ass. (0->37)
t67 = sin(qJ(3));
t59 = t67 * pkin(3) + qJ(2);
t86 = 0.2e1 * t59;
t85 = -pkin(1) - pkin(6);
t84 = (pkin(1) * MDP(6));
t78 = -qJ(4) + t85;
t56 = t78 * t67;
t64 = sin(pkin(8));
t65 = cos(pkin(8));
t69 = cos(qJ(3));
t76 = t78 * t69;
t48 = t64 * t56 - t65 * t76;
t54 = -t64 * t67 + t65 * t69;
t83 = t48 * t54;
t66 = sin(qJ(5));
t68 = cos(qJ(5));
t82 = t66 * t68;
t81 = MDP(17) * pkin(3);
t80 = t66 * MDP(23);
t79 = t68 * MDP(24);
t77 = MDP(19) * t82;
t75 = MDP(20) * t68 - MDP(21) * t66;
t74 = t68 * MDP(23) - t66 * MDP(24);
t73 = t79 + t80;
t72 = MDP(14) + t74;
t71 = t66 * MDP(20) + t68 * MDP(21) - t73 * (t64 * pkin(3) + pkin(7));
t63 = t68 ^ 2;
t62 = t66 ^ 2;
t60 = -t65 * pkin(3) - pkin(4);
t55 = t64 * t69 + t65 * t67;
t53 = t55 ^ 2;
t52 = t54 ^ 2;
t50 = t65 * t56 + t64 * t76;
t47 = t55 * pkin(4) - t54 * pkin(7) + t59;
t46 = t66 * t47 + t68 * t50;
t45 = t68 * t47 - t66 * t50;
t1 = [MDP(1) + t55 * MDP(14) * t86 + (t48 ^ 2 + t50 ^ 2 + t59 ^ 2) * MDP(17) + t53 * MDP(22) + (MDP(7) * t69 - 0.2e1 * t67 * MDP(8)) * t69 + (t63 * MDP(18) - 0.2e1 * t77) * t52 + ((-2 * MDP(4) + t84) * pkin(1)) + (MDP(15) * t86 + 0.2e1 * t55 * t75) * t54 + (0.2e1 * t67 * MDP(12) + 0.2e1 * t69 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (-t50 * t55 + t83) * MDP(16) + 0.2e1 * (t45 * t55 + t66 * t83) * MDP(23) + 0.2e1 * (-t46 * t55 + t68 * t83) * MDP(24); -t53 * t80 - t84 + MDP(4) + (t50 * MDP(17) + (-MDP(16) - t79) * t55) * t55 + (-t48 * MDP(17) + (-MDP(16) - t73) * t54) * t54; MDP(6) + (t53 + t52) * MDP(17); -t50 * MDP(15) + (t85 * MDP(12) + MDP(9)) * t69 + (-t85 * MDP(13) - MDP(10)) * t67 - t72 * t48 + t71 * t55 + (MDP(18) * t82 + (-t62 + t63) * MDP(19) + t73 * t60) * t54 + ((-t54 * t65 - t55 * t64) * MDP(16) + (-t48 * t65 + t50 * t64) * MDP(17)) * pkin(3); t69 * MDP(12) - t67 * MDP(13) + (t64 * t81 - MDP(15)) * t55 + (t65 * t81 + t72) * t54; 0.2e1 * t77 + t62 * MDP(18) + MDP(11) + (t64 ^ 2 + t65 ^ 2) * MDP(17) * pkin(3) ^ 2 - 0.2e1 * t74 * t60 + 0.2e1 * (t65 * MDP(14) - t64 * MDP(15)) * pkin(3); t54 * MDP(15) + t59 * MDP(17) + t72 * t55; 0; 0; MDP(17); t55 * MDP(22) + t45 * MDP(23) - t46 * MDP(24) + t75 * t54; -t73 * t55; t71; t74; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
