% Calculate joint inertia matrix for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPPRR5_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:43
% EndTime: 2019-03-09 01:37:43
% DurationCPUTime: 0.20s
% Computational Cost: add. (202->84), mult. (297->118), div. (0->0), fcn. (238->6), ass. (0->39)
t66 = sin(qJ(6));
t68 = cos(qJ(6));
t75 = MDP(25) * t66 + MDP(26) * t68;
t90 = t66 * MDP(22) + t68 * MDP(23) - t75 * pkin(8);
t61 = sin(pkin(9));
t62 = cos(pkin(9));
t63 = pkin(3) + qJ(2);
t64 = pkin(1) + qJ(3);
t55 = t61 * t63 - t62 * t64;
t53 = pkin(7) + t55;
t89 = t53 * t66;
t88 = t53 * t68;
t87 = t66 * t68;
t69 = cos(qJ(5));
t86 = t66 * t69;
t85 = t68 * t69;
t84 = MDP(9) + (t61 ^ 2 + t62 ^ 2) * MDP(12);
t67 = sin(qJ(5));
t83 = t67 * MDP(19);
t82 = t69 * MDP(24);
t81 = MDP(21) * t87;
t54 = t61 * t64 + t62 * t63;
t52 = -pkin(4) - t54;
t80 = MDP(22) * t68 - MDP(23) * t66;
t78 = (-t61 * t86 - t68 * t62) * MDP(25) - (t61 * t85 - t66 * t62) * MDP(26);
t77 = (t68 * t61 - t62 * t86) * MDP(25) - (t66 * t61 + t62 * t85) * MDP(26);
t76 = t68 * MDP(25) - t66 * MDP(26);
t74 = MDP(18) + t76;
t73 = t54 * MDP(12) + t69 * MDP(18) + MDP(10) - t83;
t72 = -t69 * MDP(19) - t74 * t67;
t59 = t67 ^ 2;
t71 = t55 * MDP(12) + t75 * t59 - MDP(11);
t70 = (qJ(2) ^ 2);
t60 = t68 ^ 2;
t58 = t66 ^ 2;
t47 = -t69 * pkin(5) - t67 * pkin(8) + t52;
t46 = t66 * t47 + t53 * t85;
t45 = t68 * t47 - t53 * t86;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2 + t70) * MDP(6)) + (t64 ^ 2 + t70) * MDP(9) + (t54 ^ 2 + t55 ^ 2) * MDP(12) + 0.2e1 * t52 * t83 + (t60 * MDP(20) + MDP(13) - 0.2e1 * t81) * t59 + (-0.2e1 * t52 * MDP(18) + t82) * t69 + 0.2e1 * t64 * MDP(8) + 0.2e1 * t54 * MDP(10) - 0.2e1 * t55 * MDP(11) + 0.2e1 * (-t45 * t69 + t59 * t89) * MDP(25) + 0.2e1 * (t46 * t69 + t59 * t88) * MDP(26) + (2 * (MDP(5) + MDP(7)) * qJ(2)) + 0.2e1 * (MDP(14) - t80) * t67 * t69; -(pkin(1) * MDP(6)) - t64 * MDP(9) - t73 * t61 + t71 * t62 - t77 * t69 + MDP(4) - MDP(8); MDP(6) + t84; qJ(2) * MDP(9) + t71 * t61 + t73 * t62 - t78 * t69 + MDP(7); 0; t84; 0; 0; 0; MDP(12); (-t53 * MDP(19) + MDP(16) - t90) * t69 + (MDP(15) - t53 * MDP(18) + MDP(20) * t87 + (-t58 + t60) * MDP(21) + (-pkin(5) * t66 - t88) * MDP(25) + (-pkin(5) * t68 + t89) * MDP(26)) * t67; t72 * t62; t72 * t61; t74 * t69 - t83; t58 * MDP(20) + 0.2e1 * pkin(5) * t76 + MDP(17) + 0.2e1 * t81; t45 * MDP(25) - t46 * MDP(26) + t80 * t67 - t82; t77; t78; -t75 * t67; t90; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
