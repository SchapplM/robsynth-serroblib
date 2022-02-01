% Calculate joint inertia matrix for
% S5RPRPR2
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
%   see S5RPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR2_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:07
% EndTime: 2022-01-23 09:19:08
% DurationCPUTime: 0.20s
% Computational Cost: add. (175->56), mult. (311->77), div. (0->0), fcn. (264->8), ass. (0->37)
t73 = sin(pkin(9));
t75 = cos(pkin(9));
t90 = t73 ^ 2 + t75 ^ 2;
t91 = t90 * qJ(4);
t77 = sin(qJ(5));
t79 = cos(qJ(5));
t60 = t77 * t73 - t79 * t75;
t61 = t79 * t73 + t77 * t75;
t84 = t60 * MDP(16) + t61 * MDP(17);
t97 = 2 * MDP(9);
t74 = sin(pkin(8));
t96 = pkin(1) * t74;
t95 = t75 * pkin(4);
t76 = cos(pkin(8));
t66 = t76 * pkin(1) + pkin(2);
t78 = sin(qJ(3));
t80 = cos(qJ(3));
t54 = -t78 * t66 - t80 * t96;
t51 = qJ(4) - t54;
t93 = t90 * t51;
t92 = t61 * MDP(13) - t60 * MDP(14);
t53 = t80 * t66 - t78 * t96;
t89 = t53 * MDP(6);
t88 = t54 * MDP(7);
t87 = t75 * MDP(8);
t86 = MDP(5) + (MDP(11) * t61 - 0.2e1 * MDP(12) * t60) * t61;
t52 = -pkin(3) - t53;
t85 = -t87 + t84;
t83 = 0.2e1 * t84;
t70 = t75 * pkin(7);
t67 = -pkin(3) - t95;
t63 = t75 * qJ(4) + t70;
t62 = (-pkin(7) - qJ(4)) * t73;
t47 = t52 - t95;
t46 = t75 * t51 + t70;
t45 = (-pkin(7) - t51) * t73;
t1 = [MDP(1) - 0.2e1 * t52 * t87 + (t90 * t51 ^ 2 + t52 ^ 2) * MDP(10) + (t74 ^ 2 + t76 ^ 2) * MDP(4) * pkin(1) ^ 2 + t47 * t83 + 0.2e1 * t89 + 0.2e1 * t88 + t93 * t97 + t86; 0; t90 * MDP(10) + MDP(4); t89 + t88 + (t91 + t93) * MDP(9) + (-t52 * pkin(3) + t51 * t91) * MDP(10) + (pkin(3) - t52) * t87 + t86 + t84 * (t47 + t67); 0; 0.2e1 * pkin(3) * t87 + t91 * t97 + (t90 * qJ(4) ^ 2 + pkin(3) ^ 2) * MDP(10) + t67 * t83 + t86; t52 * MDP(10) + t85; 0; -pkin(3) * MDP(10) + t85; MDP(10); (t79 * t45 - t77 * t46) * MDP(16) + (-t77 * t45 - t79 * t46) * MDP(17) + t92; -t84; (t79 * t62 - t77 * t63) * MDP(16) + (-t77 * t62 - t79 * t63) * MDP(17) + t92; 0; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
