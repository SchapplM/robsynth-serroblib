% Calculate joint inertia matrix for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRRPR2_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:53
% EndTime: 2022-01-20 11:30:54
% DurationCPUTime: 0.13s
% Computational Cost: add. (220->69), mult. (372->83), div. (0->0), fcn. (294->8), ass. (0->41)
t73 = sin(qJ(2));
t92 = pkin(1) * t73;
t72 = sin(qJ(3));
t91 = t72 * pkin(2);
t76 = cos(qJ(2));
t64 = t76 * pkin(1) + pkin(2);
t75 = cos(qJ(3));
t58 = t75 * t64;
t54 = -t72 * t92 + t58;
t51 = pkin(3) + t54;
t55 = t72 * t64 + t75 * t92;
t69 = sin(pkin(9));
t70 = cos(pkin(9));
t45 = t69 * t51 + t70 * t55;
t68 = t75 * pkin(2);
t63 = t68 + pkin(3);
t53 = t69 * t63 + t70 * t91;
t71 = sin(qJ(5));
t74 = cos(qJ(5));
t90 = t71 * MDP(13) + t74 * MDP(14);
t89 = MDP(10) * pkin(3);
t88 = t54 * MDP(8);
t87 = t55 * MDP(9);
t86 = t72 * MDP(9);
t85 = t76 * MDP(5);
t84 = t74 * MDP(16);
t83 = MDP(7) + (MDP(11) * t71 + 0.2e1 * MDP(12) * t74) * t71;
t82 = MDP(4) + t83;
t44 = t70 * t51 - t69 * t55;
t81 = -t71 * MDP(17) + t84;
t80 = -MDP(16) * t71 - MDP(17) * t74;
t52 = t70 * t63 - t69 * t91;
t79 = (t75 * MDP(8) - t86) * pkin(2);
t78 = 0.2e1 * t81;
t62 = -t70 * pkin(3) - pkin(4);
t57 = t62 * t71;
t49 = -pkin(4) - t52;
t46 = t49 * t71;
t42 = -pkin(4) - t44;
t41 = t42 * t71;
t1 = [MDP(1) + (t44 ^ 2 + t45 ^ 2) * MDP(10) - t42 * t78 + 0.2e1 * (-t73 * MDP(6) + t85) * pkin(1) + 0.2e1 * t88 - 0.2e1 * t87 + t82; (t58 + t68) * MDP(8) + (t44 * t52 + t45 * t53) * MDP(10) + (t46 + t41) * MDP(17) + (-t42 - t49) * t84 + (-pkin(2) - t64) * t86 + (t85 + (-MDP(8) * t72 - MDP(9) * t75 - MDP(6)) * t73) * pkin(1) + t82; (t52 ^ 2 + t53 ^ 2) * MDP(10) - t49 * t78 + 0.2e1 * t79 + t82; t88 - t87 + (t57 + t41) * MDP(17) + (-t42 - t62) * t84 + (t44 * t70 + t45 * t69) * t89 + t83; (t57 + t46) * MDP(17) + (-t49 - t62) * t84 + (t52 * t70 + t53 * t69) * t89 + t79 + t83; (t69 ^ 2 + t70 ^ 2) * MDP(10) * pkin(3) ^ 2 - t62 * t78 + t83; 0; 0; 0; MDP(10); t80 * (pkin(8) + t45) + t90; t80 * (pkin(8) + t53) + t90; t80 * (t69 * pkin(3) + pkin(8)) + t90; t81; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
