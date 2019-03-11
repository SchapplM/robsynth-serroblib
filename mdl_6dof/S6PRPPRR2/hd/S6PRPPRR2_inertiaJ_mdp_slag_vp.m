% Calculate joint inertia matrix for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6PRPPRR2_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:11
% EndTime: 2019-03-08 19:20:11
% DurationCPUTime: 0.21s
% Computational Cost: add. (168->77), mult. (371->119), div. (0->0), fcn. (357->10), ass. (0->42)
t75 = sin(qJ(6));
t78 = cos(qJ(6));
t102 = t75 * MDP(21) + t78 * MDP(22);
t90 = MDP(5) + MDP(8);
t76 = sin(qJ(5));
t79 = cos(qJ(5));
t91 = t79 * MDP(15);
t101 = t76 * MDP(14) + MDP(7) + t91;
t99 = MDP(5) * pkin(2);
t71 = sin(pkin(11));
t72 = sin(pkin(6));
t73 = cos(pkin(11));
t77 = sin(qJ(2));
t80 = cos(qJ(2));
t56 = (t71 * t77 - t73 * t80) * t72;
t74 = cos(pkin(6));
t54 = -t56 * t79 + t74 * t76;
t98 = t54 * t79;
t65 = -t73 * pkin(2) - pkin(3);
t62 = -pkin(8) + t65;
t97 = t62 * t75;
t96 = t62 * t78;
t95 = t75 * t78;
t89 = MDP(17) * t95;
t63 = t71 * pkin(2) + qJ(4);
t88 = t65 * MDP(8) + MDP(6);
t86 = MDP(18) * t78 - MDP(19) * t75;
t85 = t78 * MDP(21) - t75 * MDP(22);
t83 = MDP(14) + t85;
t82 = t75 * MDP(18) + t78 * MDP(19) - pkin(9) * t102;
t70 = t79 ^ 2;
t69 = t78 ^ 2;
t68 = t76 ^ 2;
t67 = t75 ^ 2;
t60 = t76 * pkin(5) - t79 * pkin(9) + t63;
t58 = (-t71 * t80 - t73 * t77) * t72;
t55 = t56 * t76 + t74 * t79;
t53 = t75 * t60 + t76 * t96;
t52 = t78 * t60 - t76 * t97;
t51 = t55 * t78 - t58 * t75;
t50 = -t55 * t75 - t58 * t78;
t1 = [MDP(1) + t90 * (t56 ^ 2 + t58 ^ 2 + t74 ^ 2); (t50 * t76 + t75 * t98) * MDP(21) + (-t51 * t76 + t78 * t98) * MDP(22) + (t80 * MDP(3) - t77 * MDP(4)) * t72 + (-t73 * t99 + t88) * t56 - (t63 * MDP(8) + t71 * t99 + t101) * t58; MDP(2) + (t63 ^ 2 + t65 ^ 2) * MDP(8) + t68 * MDP(20) + (t71 ^ 2 + t73 ^ 2) * MDP(5) * pkin(2) ^ 2 + (t69 * MDP(16) + MDP(9) - 0.2e1 * t89) * t70 + 0.2e1 * (-MDP(10) + t86) * t79 * t76 + 0.2e1 * t65 * MDP(6) + 0.2e1 * (t52 * t76 - t70 * t97) * MDP(21) + 0.2e1 * (-t53 * t76 - t70 * t96) * MDP(22) + 0.2e1 * t101 * t63; t90 * t74; 0; t90; t56 * MDP(8); t88 + t102 * (-t68 - t70); 0; MDP(8); -t55 * MDP(15) - t83 * t54; (-t62 * MDP(15) - MDP(12) + t82) * t76 + (MDP(11) + t62 * MDP(14) + MDP(16) * t95 + (-t67 + t69) * MDP(17) + (-pkin(5) * t75 + t96) * MDP(21) + (-pkin(5) * t78 - t97) * MDP(22)) * t79; -t83 * t76 - t91; -t76 * MDP(15) + t83 * t79; t67 * MDP(16) + 0.2e1 * pkin(5) * t85 + MDP(13) + 0.2e1 * t89; t50 * MDP(21) - t51 * MDP(22); t76 * MDP(20) + t52 * MDP(21) - t53 * MDP(22) + t86 * t79; -t102 * t79; -t102 * t76; t82; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
