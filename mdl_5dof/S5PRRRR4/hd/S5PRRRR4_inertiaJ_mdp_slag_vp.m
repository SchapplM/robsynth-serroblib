% Calculate joint inertia matrix for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRRRR4_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:59
% EndTime: 2019-12-05 17:08:00
% DurationCPUTime: 0.17s
% Computational Cost: add. (141->50), mult. (252->61), div. (0->0), fcn. (225->6), ass. (0->30)
t66 = sin(qJ(3));
t57 = pkin(2) * t66 + pkin(7);
t65 = sin(qJ(4));
t48 = (-pkin(8) - t57) * t65;
t68 = cos(qJ(4));
t63 = t68 * pkin(8);
t49 = t57 * t68 + t63;
t64 = sin(qJ(5));
t67 = cos(qJ(5));
t86 = (t48 * t67 - t49 * t64) * MDP(20) + (-t48 * t64 - t49 * t67) * MDP(21);
t53 = (-pkin(7) - pkin(8)) * t65;
t54 = pkin(7) * t68 + t63;
t85 = (t53 * t67 - t54 * t64) * MDP(20) + (-t53 * t64 - t54 * t67) * MDP(21);
t50 = t64 * t65 - t67 * t68;
t51 = t64 * t68 + t65 * t67;
t84 = t50 * MDP(20) + t51 * MDP(21);
t74 = MDP(13) * t68 - MDP(14) * t65;
t69 = cos(qJ(3));
t83 = pkin(2) * t69;
t80 = t51 * MDP(17) - t50 * MDP(18);
t59 = -t68 * pkin(4) - pkin(3);
t76 = t65 * MDP(10) + t68 * MDP(11) + t80;
t75 = MDP(5) + (MDP(8) * t65 + 0.2e1 * MDP(9) * t68) * t65 + (MDP(15) * t51 - 0.2e1 * MDP(16) * t50) * t51;
t73 = -MDP(13) * t65 - MDP(14) * t68;
t72 = (MDP(6) * t69 - MDP(7) * t66) * pkin(2);
t71 = (MDP(20) * t67 - MDP(21) * t64) * pkin(4);
t70 = -0.2e1 * t84;
t58 = -pkin(3) - t83;
t52 = t59 - t83;
t1 = [MDP(1); 0; -t52 * t70 - 0.2e1 * t58 * t74 + MDP(2) + 0.2e1 * t72 + t75; 0; t72 + t75 + t74 * (pkin(3) - t58) + t84 * (t52 + t59); 0.2e1 * pkin(3) * t74 - t59 * t70 + t75; t74 - t84; t73 * t57 + t76 + t86; t73 * pkin(7) + t76 + t85; MDP(12) + MDP(19) + 0.2e1 * t71; -t84; t80 + t86; t80 + t85; MDP(19) + t71; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
