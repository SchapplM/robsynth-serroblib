% Calculate joint inertia matrix for
% S5PRRPR4
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
%   see S5PRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR4_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:23
% EndTime: 2019-12-05 16:23:24
% DurationCPUTime: 0.18s
% Computational Cost: add. (206->67), mult. (425->108), div. (0->0), fcn. (455->8), ass. (0->34)
t66 = sin(pkin(9));
t67 = cos(pkin(9));
t69 = sin(qJ(3));
t72 = cos(qJ(3));
t57 = -t66 * t69 + t67 * t72;
t65 = -t72 * pkin(3) - pkin(2);
t85 = -0.2e1 * t57 * pkin(4) + 0.2e1 * t65;
t84 = pkin(3) * t66;
t83 = -qJ(4) - pkin(6);
t58 = t66 * t72 + t67 * t69;
t70 = sin(qJ(2));
t52 = t58 * t70;
t53 = t57 * t70;
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t82 = (-t71 * t52 - t68 * t53) * MDP(19) + (t68 * t52 - t71 * t53) * MDP(20);
t62 = t83 * t69;
t63 = t83 * t72;
t50 = t66 * t62 - t67 * t63;
t47 = -t71 * t57 + t68 * t58;
t81 = t47 * MDP(19);
t64 = t67 * pkin(3) + pkin(4);
t80 = (t71 * t64 - t68 * t84) * MDP(19);
t79 = (-t68 * t64 - t71 * t84) * MDP(20);
t78 = t72 * MDP(10);
t49 = t67 * t62 + t66 * t63;
t43 = -t58 * pkin(7) + t49;
t44 = t57 * pkin(7) + t50;
t48 = t68 * t57 + t71 * t58;
t77 = t48 * MDP(16) - t47 * MDP(17) + (t71 * t43 - t68 * t44) * MDP(19) + (-t68 * t43 - t71 * t44) * MDP(20);
t76 = -t69 * MDP(10) - t72 * MDP(11);
t75 = t65 * MDP(13) + t48 * MDP(20) + t81;
t73 = cos(qJ(2));
t1 = [MDP(1) + (t52 ^ 2 + t53 ^ 2 + t73 ^ 2) * MDP(13); -t70 * MDP(4) + (t52 * t58 + t53 * t57) * MDP(12) + (-t52 * t49 + t53 * t50) * MDP(13) + (-t69 * MDP(11) + MDP(3) - t75 + t78) * t73; MDP(2) + 0.2e1 * pkin(2) * t78 + 0.2e1 * (-t49 * t58 + t50 * t57) * MDP(12) + (t49 ^ 2 + t50 ^ 2 + t65 ^ 2) * MDP(13) + t81 * t85 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t69 + 0.2e1 * t72 * MDP(6)) * t69 + (MDP(14) * t48 - 0.2e1 * t47 * MDP(15) + MDP(20) * t85) * t48; t76 * t70 + (-t52 * t67 + t53 * t66) * MDP(13) * pkin(3) + t82; t69 * MDP(7) + t72 * MDP(8) + t76 * pkin(6) + ((t57 * t66 - t58 * t67) * MDP(12) + (t49 * t67 + t50 * t66) * MDP(13)) * pkin(3) + t77; MDP(9) + MDP(18) + (t66 ^ 2 + t67 ^ 2) * MDP(13) * pkin(3) ^ 2 + 0.2e1 * t80 + 0.2e1 * t79; -t73 * MDP(13); t75; 0; MDP(13); t82; t77; MDP(18) + t79 + t80; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
