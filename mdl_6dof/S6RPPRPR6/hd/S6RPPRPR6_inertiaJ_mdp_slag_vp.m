% Calculate joint inertia matrix for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRPR6_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:34
% EndTime: 2019-03-09 01:51:35
% DurationCPUTime: 0.22s
% Computational Cost: add. (171->88), mult. (256->104), div. (0->0), fcn. (165->4), ass. (0->43)
t58 = sin(qJ(6));
t60 = cos(qJ(6));
t67 = t60 * MDP(26) - t58 * MDP(27);
t88 = -MDP(17) - t67;
t86 = -2 * MDP(18);
t85 = pkin(4) + pkin(8);
t55 = -pkin(7) + qJ(2);
t84 = pkin(5) - t55;
t59 = sin(qJ(4));
t46 = t84 * t59;
t83 = t46 * t59;
t82 = t58 * t60;
t56 = pkin(1) + qJ(3);
t52 = t59 ^ 2;
t61 = cos(qJ(4));
t54 = t61 ^ 2;
t49 = t52 + t54;
t81 = (MDP(20) * pkin(4));
t69 = -t61 * qJ(5) + t56;
t45 = t59 * pkin(4) + t69;
t80 = t45 * MDP(20);
t79 = t55 ^ 2 * MDP(20);
t78 = t58 * MDP(26);
t75 = t60 * MDP(27);
t74 = MDP(20) * qJ(5);
t73 = MDP(15) - MDP(18);
t72 = MDP(19) - MDP(16);
t71 = MDP(22) * t82;
t70 = MDP(18) - t81;
t68 = MDP(23) * t58 + MDP(24) * t60;
t66 = t75 + t78;
t65 = t66 + t72;
t64 = (-MDP(26) * t85 + MDP(23)) * t60 + (MDP(27) * t85 - MDP(24)) * t58;
t63 = (qJ(2) ^ 2);
t53 = t60 ^ 2;
t51 = t58 ^ 2;
t48 = t61 * pkin(4) + t59 * qJ(5);
t47 = t84 * t61;
t44 = t49 * t55;
t43 = t85 * t59 + t69;
t42 = t60 * t43 + t58 * t47;
t41 = -t58 * t43 + t60 * t47;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2 + t63) * MDP(6)) + (t56 ^ 2 + t63) * MDP(9) + (MDP(10) + MDP(25) + t79) * t54 + (t51 * MDP(21) + 0.2e1 * t71 + t79) * t52 - 0.2e1 * t44 * MDP(17) + 0.2e1 * (t41 * t61 + t60 * t83) * MDP(26) + 0.2e1 * (-t42 * t61 - t58 * t83) * MDP(27) + (-0.2e1 * t61 * MDP(19) + t86 * t59 + t80) * t45 + (2 * (MDP(5) + MDP(7)) * qJ(2)) + 0.2e1 * (-MDP(11) + t68) * t61 * t59 + 0.2e1 * (t59 * MDP(15) + t61 * MDP(16) + MDP(8)) * t56; -(pkin(1) * MDP(6)) - t56 * MDP(9) - t73 * t59 + t65 * t61 + MDP(4) - MDP(8) - t80; MDP(6) + MDP(9) + MDP(20); t44 * MDP(20) + qJ(2) * MDP(9) + t88 * t49 + MDP(7); 0; t49 * MDP(20) + MDP(9); -t48 * MDP(17) - t66 * t46 + (MDP(12) + t64) * t61 + (-MDP(13) + MDP(21) * t82 + (-t51 + t53) * MDP(22) - t67 * qJ(5)) * t59 + ((MDP(15) - t70) * t61 + (t72 + t74) * t59) * t55; 0; t48 * MDP(20) + t59 * t65 + t61 * t73; -0.2e1 * t71 + t53 * MDP(21) + MDP(14) + (t86 + t81) * pkin(4) + (0.2e1 * MDP(19) + t74 + 0.2e1 * t75 + 0.2e1 * t78) * qJ(5); (-MDP(20) * t55 - t88) * t61; 0; -t61 * MDP(20); t70; MDP(20); t61 * MDP(25) + t41 * MDP(26) - t42 * MDP(27) + t59 * t68; t66; -t67 * t61; t64; t67; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
