% Calculate joint inertia matrix for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRPP4_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:03
% EndTime: 2021-01-15 11:24:04
% DurationCPUTime: 0.17s
% Computational Cost: add. (259->69), mult. (392->95), div. (0->0), fcn. (366->4), ass. (0->32)
t92 = 2 * MDP(14);
t63 = sin(pkin(7));
t64 = cos(pkin(7));
t65 = sin(qJ(3));
t66 = cos(qJ(3));
t53 = -t63 * t65 + t64 * t66;
t54 = t63 * t66 + t64 * t65;
t91 = t53 ^ 2 + t54 ^ 2;
t84 = -pkin(1) - pkin(6);
t75 = -qJ(4) + t84;
t56 = t75 * t65;
t71 = t75 * t66;
t48 = t63 * t56 - t64 * t71;
t50 = t64 * t56 + t63 * t71;
t90 = -t48 * t53 + t50 * t54;
t89 = MDP(17) * pkin(3);
t88 = MDP(16) + MDP(19);
t87 = MDP(17) + MDP(21);
t85 = 2 * MDP(18);
t83 = (pkin(1) * MDP(6));
t59 = t65 * pkin(3) + qJ(2);
t46 = t54 * pkin(4) - t53 * qJ(5) + t59;
t80 = t46 * MDP(21);
t79 = t59 * MDP(17);
t78 = MDP(15) - MDP(20);
t60 = t64 * pkin(3) + pkin(4);
t72 = -t60 * MDP(21) - MDP(18);
t70 = t53 * t64 + t54 * t63;
t69 = MDP(14) - t72;
t57 = t63 * pkin(3) + qJ(5);
t68 = MDP(21) * t57 - t78;
t1 = [MDP(1) + (MDP(7) * t66 - 0.2e1 * t65 * MDP(8)) * t66 + (0.2e1 * t53 * MDP(15) + t54 * t92 + t79) * t59 + (-0.2e1 * t53 * MDP(20) + t54 * t85 + t80) * t46 + ((-2 * MDP(4) + t83) * pkin(1)) + (0.2e1 * t65 * MDP(12) + 0.2e1 * t66 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + t87 * (t48 ^ 2 + t50 ^ 2) - 0.2e1 * t88 * t90; t87 * t90 - t88 * t91 + MDP(4) - t83; t87 * t91 + MDP(6); (-t60 * t53 - t54 * t57) * MDP(19) + (t84 * MDP(12) + MDP(9)) * t66 + (-t84 * MDP(13) - MDP(10)) * t65 + t68 * t50 - t69 * t48 + (-t70 * MDP(16) + (-t48 * t64 + t50 * t63) * MDP(17)) * pkin(3); t66 * MDP(12) - t65 * MDP(13) + t69 * t53 + t68 * t54 + t70 * t89; MDP(11) + (t57 ^ 2 + t60 ^ 2) * MDP(21) + t60 * t85 + 0.2e1 * t57 * MDP(20) + (t64 * t92 - 0.2e1 * t63 * MDP(15) + (t63 ^ 2 + t64 ^ 2) * t89) * pkin(3); t79 + t80 + (MDP(14) + MDP(18)) * t54 + t78 * t53; 0; 0; t87; t53 * MDP(19) + t48 * MDP(21); -t53 * MDP(21); t72; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
