% Calculate joint inertia matrix for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR3_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:56
% EndTime: 2022-01-23 09:20:56
% DurationCPUTime: 0.17s
% Computational Cost: add. (208->66), mult. (401->100), div. (0->0), fcn. (300->8), ass. (0->46)
t84 = sin(pkin(9));
t88 = sin(qJ(5));
t90 = cos(qJ(5));
t116 = (-MDP(13) * t90 + MDP(14) * t88) * t84;
t115 = t90 * MDP(16) - t88 * MDP(17);
t114 = 2 * MDP(9);
t113 = 0.2e1 * MDP(16);
t112 = 0.2e1 * MDP(17);
t85 = sin(pkin(8));
t111 = pkin(1) * t85;
t87 = cos(pkin(8));
t78 = t87 * pkin(1) + pkin(2);
t89 = sin(qJ(3));
t91 = cos(qJ(3));
t67 = -t91 * t111 - t89 * t78;
t64 = qJ(4) - t67;
t86 = cos(pkin(9));
t110 = t64 * t86;
t82 = t84 ^ 2;
t61 = t82 * t64;
t109 = t82 * t90;
t66 = -t89 * t111 + t91 * t78;
t71 = -t86 * pkin(4) - t84 * pkin(7) - pkin(3);
t55 = -t66 + t71;
t54 = t90 * t110 + t88 * t55;
t108 = t64 * t109 + t54 * t86;
t103 = qJ(4) * t86;
t60 = t90 * t103 + t88 * t71;
t80 = t82 * qJ(4);
t107 = t60 * t86 + t90 * t80;
t83 = t86 ^ 2;
t106 = t83 * t64 + t61;
t105 = t83 * qJ(4) + t80;
t104 = t82 + t83;
t102 = t66 * MDP(6);
t101 = t67 * MDP(7);
t100 = t86 * MDP(8);
t96 = t90 ^ 2 * t82 * MDP(11) - 0.2e1 * t88 * MDP(12) * t109 + t83 * MDP(15) + 0.2e1 * t116 * t86 + MDP(5);
t95 = (-MDP(8) - t115) * t86;
t94 = -t86 * MDP(15) - t116;
t75 = t88 * t80;
t65 = -pkin(3) - t66;
t59 = -t88 * t103 + t90 * t71;
t57 = t88 * t61;
t53 = -t88 * t110 + t90 * t55;
t1 = [MDP(1) - 0.2e1 * t65 * t100 + (t104 * t64 ^ 2 + t65 ^ 2) * MDP(10) + (t85 ^ 2 + t87 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * t102 + 0.2e1 * t101 + t106 * t114 + (-t53 * t86 + t57) * t113 + t108 * t112 + t96; 0; t104 * MDP(10) + MDP(4); t102 + t101 + (t105 + t106) * MDP(9) + (t104 * t64 * qJ(4) - t65 * pkin(3)) * MDP(10) + (t57 + t75) * MDP(16) + (t107 + t108) * MDP(17) + ((pkin(3) - t65) * MDP(8) + (-t53 - t59) * MDP(16)) * t86 + t96; 0; 0.2e1 * pkin(3) * t100 + (t104 * qJ(4) ^ 2 + pkin(3) ^ 2) * MDP(10) + t105 * t114 + (-t59 * t86 + t75) * t113 + t107 * t112 + t96; t65 * MDP(10) + t95; 0; -pkin(3) * MDP(10) + t95; MDP(10); t53 * MDP(16) - t54 * MDP(17) + t94; (-MDP(16) * t88 - MDP(17) * t90) * t84; t59 * MDP(16) - t60 * MDP(17) + t94; t115; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
