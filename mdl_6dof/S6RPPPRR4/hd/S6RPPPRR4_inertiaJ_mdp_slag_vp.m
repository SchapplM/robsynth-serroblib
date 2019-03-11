% Calculate joint inertia matrix for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPPRR4_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:53
% EndTime: 2019-03-09 01:35:53
% DurationCPUTime: 0.23s
% Computational Cost: add. (195->82), mult. (292->111), div. (0->0), fcn. (225->6), ass. (0->37)
t61 = sin(qJ(6));
t63 = cos(qJ(6));
t67 = t61 * MDP(25) + t63 * MDP(26);
t88 = t61 * MDP(22) + t63 * MDP(23) - pkin(8) * t67;
t62 = sin(qJ(5));
t64 = cos(qJ(5));
t68 = t63 * MDP(25) - t61 * MDP(26);
t66 = MDP(18) + t68;
t87 = -t62 * MDP(19) + t64 * t66;
t75 = t64 * MDP(19);
t86 = t62 * MDP(18) + MDP(11) + t75;
t59 = sin(pkin(9));
t60 = cos(pkin(9));
t65 = -pkin(1) - pkin(2);
t50 = t59 * qJ(2) - t60 * t65;
t49 = pkin(3) + t50;
t46 = pkin(7) + t49;
t84 = t46 * t61;
t83 = t46 * t63;
t82 = t61 * t62;
t81 = t61 * t63;
t80 = t62 * t63;
t74 = MDP(9) + MDP(12);
t73 = MDP(21) * t81;
t72 = t49 * MDP(12) - MDP(10);
t52 = t60 * qJ(2) + t59 * t65;
t71 = -MDP(22) * t63 + MDP(23) * t61;
t69 = (t63 * t59 + t60 * t82) * MDP(25) + (-t61 * t59 + t60 * t80) * MDP(26);
t47 = -qJ(4) + t52;
t58 = t64 ^ 2;
t57 = t63 ^ 2;
t56 = t62 ^ 2;
t55 = t61 ^ 2;
t43 = -t62 * pkin(5) + t64 * pkin(8) + t47;
t42 = t61 * t43 + t46 * t80;
t41 = t63 * t43 - t46 * t82;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t50 ^ 2 + t52 ^ 2) * MDP(9) + (t47 ^ 2 + t49 ^ 2) * MDP(12) + t56 * MDP(24) + (t57 * MDP(20) + MDP(13) - 0.2e1 * t73) * t58 + 0.2e1 * (-MDP(14) - t71) * t64 * t62 + 0.2e1 * t50 * MDP(7) + 0.2e1 * t52 * MDP(8) - 0.2e1 * t49 * MDP(10) + 0.2e1 * (-t41 * t62 + t58 * t84) * MDP(25) + 0.2e1 * (t42 * t62 + t58 * t83) * MDP(26) - 0.2e1 * t86 * t47; -pkin(1) * MDP(6) - MDP(4) - t69 * t62 + (-t50 * MDP(9) - t58 * t67 - MDP(7) - t72) * t60 + (t47 * MDP(12) + t52 * MDP(9) + MDP(8) - t86) * t59; MDP(6) + t74 * (t59 ^ 2 + t60 ^ 2); 0; 0; t74; t72 + t67 * (t56 + t58); -t60 * MDP(12); 0; MDP(12); (-t46 * MDP(19) + MDP(16) - t88) * t62 + (-MDP(15) + t46 * MDP(18) - MDP(20) * t81 + (t55 - t57) * MDP(21) + (pkin(5) * t61 + t83) * MDP(25) + (pkin(5) * t63 - t84) * MDP(26)) * t64; -t87 * t60; -t62 * t66 - t75; t87; t55 * MDP(20) + 0.2e1 * pkin(5) * t68 + MDP(17) + 0.2e1 * t73; -t62 * MDP(24) + t41 * MDP(25) - t42 * MDP(26) + t64 * t71; t69; -t67 * t64; -t67 * t62; t88; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
