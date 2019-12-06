% Calculate joint inertia matrix for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP2_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:59
% EndTime: 2019-12-05 16:42:00
% DurationCPUTime: 0.15s
% Computational Cost: add. (126->53), mult. (218->73), div. (0->0), fcn. (125->4), ass. (0->29)
t52 = sin(qJ(3));
t41 = t52 * pkin(2) + pkin(7);
t51 = sin(qJ(4));
t49 = t51 ^ 2;
t53 = cos(qJ(4));
t67 = t53 ^ 2 + t49;
t69 = t67 * t41;
t73 = 2 * MDP(16);
t54 = cos(qJ(3));
t72 = t54 * pkin(2);
t42 = -pkin(3) - t72;
t71 = pkin(3) - t42;
t36 = -t53 * pkin(4) - t51 * qJ(5) - pkin(3);
t34 = t36 - t72;
t70 = -t34 - t36;
t68 = t67 * pkin(7);
t66 = MDP(18) * t51;
t65 = t51 * MDP(10) + t53 * MDP(11) + (-t51 * pkin(4) + t53 * qJ(5)) * MDP(16);
t64 = 0.2e1 * t51 * t53 * MDP(9) + t49 * MDP(8) + MDP(5);
t63 = t67 * MDP(18);
t62 = -pkin(4) * MDP(18) - MDP(15);
t61 = MDP(13) - t62;
t60 = MDP(18) * qJ(5) - MDP(14) + MDP(17);
t59 = t53 * MDP(13) - t51 * MDP(14);
t58 = -0.2e1 * t53 * MDP(15) - 0.2e1 * t51 * MDP(17);
t57 = (t54 * MDP(6) - t52 * MDP(7)) * pkin(2);
t56 = -t51 * t61 + t53 * t60;
t43 = t51 * MDP(16);
t1 = [MDP(1) + t63; 0; MDP(2) + t69 * t73 + t41 ^ 2 * t63 + (t34 * MDP(18) + t58) * t34 + t64 - 0.2e1 * t59 * t42 + 0.2e1 * t57; 0; (t68 + t69) * MDP(16) + (pkin(7) * t69 + t34 * t36) * MDP(18) + t57 + (t71 * MDP(13) + t70 * MDP(15)) * t53 + (-t71 * MDP(14) + t70 * MDP(17)) * t51 + t64; t68 * t73 + pkin(7) ^ 2 * t63 + (MDP(18) * t36 + t58) * t36 + 0.2e1 * t59 * pkin(3) + t64; t51 * t60 + t53 * t61; t41 * t56 + t65; pkin(7) * t56 + t65; MDP(12) + 0.2e1 * pkin(4) * MDP(15) + 0.2e1 * qJ(5) * MDP(17) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(18); -t53 * MDP(18); t41 * t66 + t43; pkin(7) * t66 + t43; t62; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
