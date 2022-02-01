% Calculate joint inertia matrix for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR1_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:43
% EndTime: 2022-01-20 09:51:44
% DurationCPUTime: 0.23s
% Computational Cost: add. (200->61), mult. (342->85), div. (0->0), fcn. (291->8), ass. (0->39)
t84 = sin(pkin(9));
t86 = cos(pkin(9));
t100 = t84 ^ 2 + t86 ^ 2;
t85 = sin(pkin(8));
t77 = t85 * pkin(2) + qJ(4);
t101 = t100 * t77;
t88 = sin(qJ(5));
t90 = cos(qJ(5));
t68 = t88 * t84 - t90 * t86;
t69 = t90 * t84 + t88 * t86;
t95 = t68 * MDP(16) + t69 * MDP(17);
t107 = 2 * MDP(9);
t89 = sin(qJ(2));
t106 = pkin(1) * t89;
t105 = t86 * pkin(4);
t91 = cos(qJ(2));
t80 = t91 * pkin(1) + pkin(2);
t87 = cos(pkin(8));
t60 = t87 * t106 + t85 * t80;
t57 = qJ(4) + t60;
t103 = t100 * t57;
t102 = t69 * MDP(13) - t68 * MDP(14);
t99 = t86 * MDP(8);
t98 = -0.2e1 * t99;
t97 = MDP(4) + (MDP(11) * t69 - 0.2e1 * MDP(12) * t68) * t69;
t79 = -t87 * pkin(2) - pkin(3);
t59 = -t85 * t106 + t87 * t80;
t58 = -pkin(3) - t59;
t96 = -t99 + t95;
t94 = (t91 * MDP(5) - t89 * MDP(6)) * pkin(1);
t93 = 0.2e1 * t95;
t81 = t86 * pkin(7);
t70 = t79 - t105;
t62 = t86 * t77 + t81;
t61 = (-pkin(7) - t77) * t84;
t53 = t58 - t105;
t52 = t86 * t57 + t81;
t51 = (-pkin(7) - t57) * t84;
t1 = [MDP(1) + (t59 ^ 2 + t60 ^ 2) * MDP(7) + t58 * t98 + t103 * t107 + (t100 * t57 ^ 2 + t58 ^ 2) * MDP(10) + t53 * t93 + 0.2e1 * t94 + t97; (t101 + t103) * MDP(9) + (t57 * t101 + t58 * t79) * MDP(10) + (-t58 - t79) * t99 + (t59 * t87 + t60 * t85) * MDP(7) * pkin(2) + t94 + t97 + t95 * (t53 + t70); t79 * t98 + t101 * t107 + (t100 * t77 ^ 2 + t79 ^ 2) * MDP(10) + (t85 ^ 2 + t87 ^ 2) * MDP(7) * pkin(2) ^ 2 + t70 * t93 + t97; 0; 0; t100 * MDP(10) + MDP(7); t58 * MDP(10) + t96; t79 * MDP(10) + t96; 0; MDP(10); (t90 * t51 - t88 * t52) * MDP(16) + (-t88 * t51 - t90 * t52) * MDP(17) + t102; (t90 * t61 - t88 * t62) * MDP(16) + (-t88 * t61 - t90 * t62) * MDP(17) + t102; -t95; 0; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
