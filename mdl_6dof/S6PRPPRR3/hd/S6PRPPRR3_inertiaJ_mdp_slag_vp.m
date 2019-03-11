% Calculate joint inertia matrix for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRPPRR3_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:03
% EndTime: 2019-03-08 19:23:04
% DurationCPUTime: 0.21s
% Computational Cost: add. (213->95), mult. (421->147), div. (0->0), fcn. (405->10), ass. (0->47)
t78 = sin(pkin(11));
t79 = sin(pkin(6));
t80 = cos(pkin(11));
t84 = sin(qJ(2));
t87 = cos(qJ(2));
t62 = (-t78 * t87 + t80 * t84) * t79;
t81 = cos(pkin(6));
t83 = sin(qJ(5));
t86 = cos(qJ(5));
t57 = t62 * t83 + t81 * t86;
t106 = t57 * t83;
t88 = -pkin(2) - pkin(3);
t69 = t80 * qJ(3) + t78 * t88;
t66 = -pkin(8) + t69;
t82 = sin(qJ(6));
t105 = t66 * t82;
t85 = cos(qJ(6));
t104 = t66 * t85;
t103 = t82 * t85;
t102 = t82 * t86;
t101 = t85 * t86;
t100 = MDP(22) * t86;
t99 = t83 * MDP(17);
t98 = MDP(19) * t103;
t97 = -pkin(2) * MDP(7) - MDP(5);
t67 = qJ(3) * t78 - t80 * t88;
t96 = t69 * MDP(10) + MDP(9);
t65 = pkin(4) + t67;
t95 = -MDP(20) * t85 + MDP(21) * t82;
t94 = MDP(23) * (-t78 * t102 - t80 * t85) - MDP(24) * (t78 * t101 - t80 * t82);
t93 = t85 * MDP(23) - t82 * MDP(24);
t92 = -MDP(23) * t82 - MDP(24) * t85;
t91 = MDP(16) + t93;
t90 = t67 * MDP(10) + t86 * MDP(16) + MDP(8) - t99;
t89 = t82 * MDP(20) + t85 * MDP(21) + t92 * pkin(9);
t77 = t85 ^ 2;
t76 = t83 ^ 2;
t75 = t82 ^ 2;
t74 = t81 ^ 2;
t60 = (-t78 * t84 - t80 * t87) * t79;
t59 = pkin(5) * t86 + pkin(9) * t83 + t65;
t58 = t62 * t86 - t81 * t83;
t56 = t66 * t101 + t59 * t82;
t55 = -t66 * t102 + t59 * t85;
t54 = t58 * t85 - t60 * t82;
t53 = -t58 * t82 - t60 * t85;
t1 = [MDP(1) + (t74 + (t84 ^ 2 + t87 ^ 2) * t79 ^ 2) * MDP(7) + (t60 ^ 2 + t62 ^ 2 + t74) * MDP(10); (-t82 * t106 + t53 * t86) * MDP(23) + (-t85 * t106 - t54 * t86) * MDP(24) + t96 * t62 - t90 * t60 + ((MDP(3) - t97) * t87 + (qJ(3) * MDP(7) - MDP(4) + MDP(6)) * t84) * t79; MDP(2) + (2 * pkin(2) * MDP(5)) + 0.2e1 * qJ(3) * MDP(6) + ((pkin(2) ^ 2) + qJ(3) ^ 2) * MDP(7) + (t67 ^ 2 + t69 ^ 2) * MDP(10) - 0.2e1 * t65 * t99 + (t77 * MDP(18) + MDP(11) - 0.2e1 * t98) * t76 + (0.2e1 * t65 * MDP(16) + t100 + 0.2e1 * (MDP(12) + t95) * t83) * t86 + 0.2e1 * t67 * MDP(8) + 0.2e1 * t69 * MDP(9) + 0.2e1 * (-t76 * t105 + t55 * t86) * MDP(23) + 0.2e1 * (-t76 * t104 - t56 * t86) * MDP(24); -t79 * t87 * MDP(7) + (t60 * t80 + t62 * t78) * MDP(10); t94 * t86 + (t92 * t76 + t96) * t78 - t90 * t80 + t97; MDP(7) + (t78 ^ 2 + t80 ^ 2) * MDP(10); -t81 * MDP(10); 0; 0; MDP(10); -MDP(17) * t58 - t91 * t57; (-t66 * MDP(17) - MDP(14) + t89) * t86 + (-MDP(13) - t66 * MDP(16) - MDP(18) * t103 + (t75 - t77) * MDP(19) + (pkin(5) * t82 - t104) * MDP(23) + (pkin(5) * t85 + t105) * MDP(24)) * t83; (-MDP(17) * t86 - t91 * t83) * t78; t91 * t86 - t99; MDP(18) * t75 + 0.2e1 * pkin(5) * t93 + MDP(15) + 0.2e1 * t98; t53 * MDP(23) - t54 * MDP(24); t55 * MDP(23) - t56 * MDP(24) + t95 * t83 + t100; t94; t92 * t83; t89; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
