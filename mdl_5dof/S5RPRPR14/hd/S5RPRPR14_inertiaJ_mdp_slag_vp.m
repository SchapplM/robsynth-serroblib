% Calculate joint inertia matrix for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR14_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR14_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:17
% EndTime: 2019-12-31 18:35:18
% DurationCPUTime: 0.25s
% Computational Cost: add. (239->75), mult. (413->116), div. (0->0), fcn. (411->6), ass. (0->36)
t67 = sin(pkin(8));
t68 = cos(pkin(8));
t71 = cos(qJ(3));
t87 = sin(qJ(3));
t54 = -t67 * t87 + t68 * t71;
t69 = sin(qJ(5));
t70 = cos(qJ(5));
t90 = (MDP(18) * t70 - MDP(19) * t69) * t54;
t72 = -pkin(1) - pkin(6);
t89 = -qJ(4) + t72;
t82 = t70 * MDP(22);
t83 = t69 * MDP(21);
t74 = t82 + t83;
t88 = t74 * (t67 * pkin(3) + pkin(7)) - t69 * MDP(18) - t70 * MDP(19);
t86 = (pkin(1) * MDP(6));
t57 = t89 * t87;
t80 = t89 * t71;
t48 = t67 * t57 - t68 * t80;
t85 = t48 * t54;
t84 = t69 * t70;
t63 = t87 * pkin(3) + qJ(2);
t81 = MDP(17) * t84;
t79 = t87 * MDP(13);
t55 = -t67 * t71 - t68 * t87;
t78 = t54 * t68 - t55 * t67;
t75 = t70 * MDP(21) - t69 * MDP(22);
t66 = t70 ^ 2;
t65 = t69 ^ 2;
t62 = -t68 * pkin(3) - pkin(4);
t53 = t55 ^ 2;
t52 = t54 ^ 2;
t50 = t68 * t57 + t67 * t80;
t47 = -t55 * pkin(4) - t54 * pkin(7) + t63;
t46 = t69 * t47 + t70 * t50;
t45 = t70 * t47 - t69 * t50;
t1 = [MDP(1) + (t48 ^ 2 + t50 ^ 2 + t63 ^ 2) * MDP(15) + t53 * MDP(20) + (MDP(7) * t71 - 0.2e1 * t87 * MDP(8)) * t71 - 0.2e1 * t55 * t90 + (t66 * MDP(16) - 0.2e1 * t81) * t52 + ((-2 * MDP(4) + t86) * pkin(1)) + (0.2e1 * t87 * MDP(12) + 0.2e1 * t71 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (t50 * t55 + t85) * MDP(14) + 0.2e1 * (-t45 * t55 + t69 * t85) * MDP(21) + 0.2e1 * (t46 * t55 + t70 * t85) * MDP(22); -t53 * t83 - t86 + MDP(4) + (-t50 * MDP(15) + (-MDP(14) - t82) * t55) * t55 + (-t48 * MDP(15) + (-MDP(14) - t74) * t54) * t54; MDP(6) + (t53 + t52) * MDP(15); -t72 * t79 - t87 * MDP(10) + (t72 * MDP(12) + MDP(9)) * t71 - t75 * t48 + t88 * t55 + (MDP(16) * t84 + (-t65 + t66) * MDP(17) + t74 * t62) * t54 + (-t78 * MDP(14) + (-t48 * t68 + t50 * t67) * MDP(15)) * pkin(3); t78 * MDP(15) * pkin(3) + t71 * MDP(12) + t75 * t54 - t79; 0.2e1 * t81 + t65 * MDP(16) + MDP(11) + (t67 ^ 2 + t68 ^ 2) * MDP(15) * pkin(3) ^ 2 - 0.2e1 * t75 * t62; t63 * MDP(15) - t75 * t55; 0; 0; MDP(15); -t55 * MDP(20) + t45 * MDP(21) - t46 * MDP(22) + t90; t74 * t55; -t88; t75; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
