% Calculate joint inertia matrix for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR3_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:42
% EndTime: 2019-12-05 16:19:42
% DurationCPUTime: 0.16s
% Computational Cost: add. (190->56), mult. (365->92), div. (0->0), fcn. (389->6), ass. (0->28)
t58 = sin(pkin(9));
t59 = cos(pkin(9));
t61 = sin(qJ(3));
t63 = cos(qJ(3));
t49 = -t58 * t61 + t59 * t63;
t57 = -t63 * pkin(3) - pkin(2);
t72 = -0.2e1 * t49 * pkin(4) + 0.2e1 * t57;
t71 = pkin(3) * t58;
t70 = -qJ(4) - pkin(6);
t50 = t58 * t63 + t59 * t61;
t60 = sin(qJ(5));
t62 = cos(qJ(5));
t41 = -t62 * t49 + t60 * t50;
t37 = t41 * MDP(19);
t42 = t60 * t49 + t62 * t50;
t69 = -t42 * MDP(20) - t37;
t54 = t70 * t61;
t55 = t70 * t63;
t44 = t58 * t54 - t59 * t55;
t56 = t59 * pkin(3) + pkin(4);
t68 = (t62 * t56 - t60 * t71) * MDP(19);
t67 = (-t60 * t56 - t62 * t71) * MDP(20);
t66 = t63 * MDP(10);
t43 = t59 * t54 + t58 * t55;
t35 = -t50 * pkin(7) + t43;
t36 = t49 * pkin(7) + t44;
t65 = t42 * MDP(16) - t41 * MDP(17) + (t62 * t35 - t60 * t36) * MDP(19) + (-t60 * t35 - t62 * t36) * MDP(20);
t1 = [MDP(1) + (t49 ^ 2 + t50 ^ 2) * MDP(13); (t49 * t43 + t50 * t44) * MDP(13); MDP(2) + 0.2e1 * pkin(2) * t66 + 0.2e1 * (-t43 * t50 + t44 * t49) * MDP(12) + (t43 ^ 2 + t44 ^ 2 + t57 ^ 2) * MDP(13) + t37 * t72 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t61 + 0.2e1 * t63 * MDP(6)) * t61 + (MDP(14) * t42 - 0.2e1 * t41 * MDP(15) + MDP(20) * t72) * t42; t66 - t61 * MDP(11) + (t49 * t59 + t50 * t58) * MDP(13) * pkin(3) + t69; t61 * MDP(7) + t63 * MDP(8) + (-t61 * MDP(10) - t63 * MDP(11)) * pkin(6) + ((t49 * t58 - t50 * t59) * MDP(12) + (t43 * t59 + t44 * t58) * MDP(13)) * pkin(3) + t65; MDP(9) + MDP(18) + (t58 ^ 2 + t59 ^ 2) * MDP(13) * pkin(3) ^ 2 + 0.2e1 * t68 + 0.2e1 * t67; 0; t57 * MDP(13) - t69; 0; MDP(13); t69; t65; MDP(18) + t67 + t68; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
