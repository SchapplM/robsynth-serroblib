% Calculate joint inertia matrix for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRPP2_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:10
% EndTime: 2019-12-05 16:10:11
% DurationCPUTime: 0.21s
% Computational Cost: add. (204->67), mult. (408->102), div. (0->0), fcn. (398->6), ass. (0->29)
t81 = MDP(12) + MDP(15);
t75 = MDP(13) + MDP(17);
t80 = -qJ(4) - pkin(6);
t60 = sin(pkin(8));
t61 = cos(pkin(8));
t62 = sin(qJ(3));
t64 = cos(qJ(3));
t50 = t60 * t62 - t61 * t64;
t51 = t60 * t64 + t61 * t62;
t58 = -t64 * pkin(3) - pkin(2);
t40 = t50 * pkin(4) - t51 * qJ(5) + t58;
t79 = t40 * MDP(17);
t78 = t50 * MDP(14);
t77 = t51 * MDP(16);
t76 = t64 * MDP(10);
t53 = t80 * t64;
t73 = t80 * t62;
t43 = -t60 * t53 - t61 * t73;
t45 = -t61 * t53 + t60 * t73;
t74 = t43 ^ 2 + t45 ^ 2;
t68 = -t62 * MDP(10) - t64 * MDP(11);
t67 = t58 * MDP(13) - t77 + t78 + t79;
t65 = cos(qJ(2));
t63 = sin(qJ(2));
t56 = t61 * pkin(3) + pkin(4);
t54 = t60 * pkin(3) + qJ(5);
t49 = t50 * t63;
t47 = t51 * t63;
t1 = [MDP(1) + t75 * (t47 ^ 2 + t49 ^ 2 + t65 ^ 2); -t63 * MDP(4) + (-t62 * MDP(11) + MDP(3) - t67 + t76) * t65 + t75 * (t47 * t43 - t49 * t45) + t81 * (t47 * t51 + t49 * t50); MDP(2) + 0.2e1 * pkin(2) * t76 + (t58 ^ 2 + t74) * MDP(13) + t74 * MDP(17) + (-0.2e1 * t77 + 0.2e1 * t78 + t79) * t40 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t62 + 0.2e1 * t64 * MDP(6)) * t62 + 0.2e1 * t81 * (t43 * t51 - t45 * t50); -t47 * MDP(14) - t49 * MDP(16) + (-t47 * t56 - t49 * t54) * MDP(17) + t68 * t63 + (-t47 * t61 - t49 * t60) * MDP(13) * pkin(3); t62 * MDP(7) + t64 * MDP(8) - t43 * MDP(14) + (-t54 * t50 - t56 * t51) * MDP(15) + t45 * MDP(16) + (-t43 * t56 + t45 * t54) * MDP(17) + t68 * pkin(6) + ((-t50 * t60 - t51 * t61) * MDP(12) + (-t43 * t61 + t45 * t60) * MDP(13)) * pkin(3); MDP(9) + (t54 ^ 2 + t56 ^ 2) * MDP(17) + (t60 ^ 2 + t61 ^ 2) * MDP(13) * pkin(3) ^ 2 + 0.2e1 * t56 * MDP(14) + 0.2e1 * t54 * MDP(16); -t75 * t65; t67; 0; t75; t47 * MDP(17); t51 * MDP(15) + t43 * MDP(17); -t56 * MDP(17) - MDP(14); 0; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
