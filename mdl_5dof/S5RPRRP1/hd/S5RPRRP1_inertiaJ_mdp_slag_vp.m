% Calculate joint inertia matrix for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP1_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:22
% EndTime: 2021-01-15 12:27:23
% DurationCPUTime: 0.19s
% Computational Cost: add. (244->73), mult. (374->98), div. (0->0), fcn. (352->4), ass. (0->35)
t82 = -MDP(20) - MDP(22);
t58 = sin(qJ(4));
t59 = sin(qJ(3));
t60 = cos(qJ(4));
t61 = cos(qJ(3));
t50 = -t58 * t59 + t60 * t61;
t48 = t50 ^ 2;
t49 = t58 * t61 + t60 * t59;
t81 = t49 ^ 2 + t48;
t79 = 2 * MDP(21);
t78 = pkin(3) * t58;
t62 = -pkin(1) - pkin(6);
t77 = -pkin(7) + t62;
t76 = (pkin(1) * MDP(6));
t75 = MDP(24) * pkin(4);
t74 = t59 * pkin(3) + qJ(2);
t40 = t49 * pkin(4) + t74;
t73 = t40 * MDP(24);
t72 = t49 * MDP(21);
t71 = t50 * MDP(22);
t57 = t60 * pkin(3);
t55 = t57 + pkin(4);
t70 = t55 * MDP(24);
t69 = t60 * MDP(19);
t51 = t77 * t59;
t52 = t77 * t61;
t68 = -t58 * t51 + t60 * t52;
t67 = (MDP(19) + MDP(21)) * t50 + t82 * t49;
t36 = -t50 * qJ(5) + t68;
t65 = -t60 * t51 - t58 * t52;
t37 = -t49 * qJ(5) - t65;
t66 = t36 * t50 + t37 * t49;
t64 = t50 * MDP(16) - t49 * MDP(17) + t68 * MDP(19) + t65 * MDP(20) + t36 * MDP(21) - t37 * MDP(22);
t63 = t49 * t78 + t50 * t55;
t1 = [MDP(1) + t48 * MDP(14) - 0.2e1 * t49 * t50 * MDP(15) - 0.2e1 * t66 * MDP(23) + (t36 ^ 2 + t37 ^ 2) * MDP(24) + (MDP(7) * t61 - 0.2e1 * t59 * MDP(8)) * t61 + 0.2e1 * (t49 * MDP(19) + t50 * MDP(20)) * t74 + (0.2e1 * t71 + 0.2e1 * t72 + t73) * t40 + ((-2 * MDP(4) + t76) * pkin(1)) + (0.2e1 * t59 * MDP(12) + 0.2e1 * t61 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2); -t81 * MDP(23) + t66 * MDP(24) + MDP(4) - t76; t81 * MDP(24) + MDP(6); -t63 * MDP(23) + (t36 * t55 + t37 * t78) * MDP(24) + (t62 * MDP(12) + MDP(9)) * t61 + (-t62 * MDP(13) - MDP(10)) * t59 + t64; t61 * MDP(12) - t59 * MDP(13) + MDP(24) * t63 + t67; MDP(11) + MDP(18) + (t79 + t70) * t55 + (0.2e1 * t69 + (MDP(24) * t78 - 0.2e1 * MDP(20) - 0.2e1 * MDP(22)) * t58) * pkin(3); (-MDP(23) * t50 + t36 * MDP(24)) * pkin(4) + t64; t50 * t75 + t67; MDP(18) + (0.2e1 * pkin(4) + t57) * MDP(21) + pkin(4) * t70 + (t82 * t58 + t69) * pkin(3); MDP(18) + (t79 + t75) * pkin(4); t71 + t72 + t73; 0; 0; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
