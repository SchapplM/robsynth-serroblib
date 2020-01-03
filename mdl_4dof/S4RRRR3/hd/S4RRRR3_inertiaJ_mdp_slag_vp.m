% Calculate joint inertia matrix for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR3_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:37
% EndTime: 2019-12-31 17:24:37
% DurationCPUTime: 0.12s
% Computational Cost: add. (193->58), mult. (353->85), div. (0->0), fcn. (357->6), ass. (0->33)
t63 = sin(qJ(3));
t64 = sin(qJ(2));
t66 = cos(qJ(3));
t67 = cos(qJ(2));
t53 = t63 * t64 - t66 * t67;
t61 = -t67 * pkin(2) - pkin(1);
t82 = 0.2e1 * t53 * pkin(3) + 0.2e1 * t61;
t81 = 0.2e1 * t61;
t80 = 0.2e1 * t67;
t79 = pkin(5) + pkin(6);
t78 = pkin(2) * t63;
t60 = t66 * pkin(2) + pkin(3);
t65 = cos(qJ(4));
t58 = t65 * t60;
t62 = sin(qJ(4));
t77 = (-t62 * t78 + t58) * MDP(23);
t76 = (-t62 * t60 - t65 * t78) * MDP(24);
t75 = t62 * MDP(24);
t74 = t66 * MDP(16);
t73 = MDP(15) + MDP(22);
t56 = t79 * t64;
t57 = t79 * t67;
t72 = -t66 * t56 - t63 * t57;
t54 = t63 * t67 + t66 * t64;
t39 = -t54 * pkin(7) + t72;
t70 = t63 * t56 - t66 * t57;
t40 = -t53 * pkin(7) - t70;
t43 = t65 * t53 + t62 * t54;
t44 = -t62 * t53 + t65 * t54;
t71 = t44 * MDP(20) - t43 * MDP(21) + (t65 * t39 - t62 * t40) * MDP(23) + (-t62 * t39 - t65 * t40) * MDP(24);
t69 = (t65 * MDP(23) - t75) * pkin(3);
t68 = t54 * MDP(13) - t53 * MDP(14) + t72 * MDP(16) + t70 * MDP(17) + t71;
t1 = [t53 * MDP(16) * t81 + t43 * MDP(23) * t82 + pkin(1) * MDP(9) * t80 + MDP(1) + (-0.2e1 * MDP(10) * pkin(1) + MDP(4) * t64 + MDP(5) * t80) * t64 + (MDP(11) * t54 - 0.2e1 * MDP(12) * t53 + MDP(17) * t81) * t54 + (MDP(18) * t44 - 0.2e1 * MDP(19) * t43 + MDP(24) * t82) * t44; t64 * MDP(6) + t67 * MDP(7) + (-MDP(10) * t67 - MDP(9) * t64) * pkin(5) + t68; MDP(8) + 0.2e1 * (-t63 * MDP(17) + t74) * pkin(2) + 0.2e1 * t77 + 0.2e1 * t76 + t73; t68; (t65 * pkin(3) + t58) * MDP(23) + (-pkin(3) - t60) * t75 + (t74 + (-MDP(23) * t62 - MDP(24) * t65 - MDP(17)) * t63) * pkin(2) + t73; 0.2e1 * t69 + t73; t71; MDP(22) + t76 + t77; MDP(22) + t69; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
