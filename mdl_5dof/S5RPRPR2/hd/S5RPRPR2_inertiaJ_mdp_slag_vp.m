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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRPR2_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:39
% EndTime: 2019-12-05 17:49:39
% DurationCPUTime: 0.17s
% Computational Cost: add. (187->57), mult. (330->73), div. (0->0), fcn. (278->8), ass. (0->41)
t75 = sin(pkin(9));
t77 = cos(pkin(9));
t96 = t75 ^ 2 + t77 ^ 2;
t97 = t96 * qJ(4);
t79 = sin(qJ(5));
t81 = cos(qJ(5));
t61 = t75 * t79 - t81 * t77;
t62 = t75 * t81 + t77 * t79;
t86 = t61 * MDP(17) + t62 * MDP(18);
t105 = MDP(8) * t77 - t75 * MDP(9);
t104 = 2 * MDP(10);
t76 = sin(pkin(8));
t103 = pkin(1) * t76;
t102 = pkin(4) * t77;
t78 = cos(pkin(8));
t67 = pkin(1) * t78 + pkin(2);
t80 = sin(qJ(3));
t82 = cos(qJ(3));
t55 = -t82 * t103 - t67 * t80;
t52 = qJ(4) - t55;
t99 = t96 * t52;
t98 = t62 * MDP(14) - t61 * MDP(15);
t94 = pkin(3) * MDP(11);
t54 = -t80 * t103 + t67 * t82;
t93 = t54 * MDP(6);
t92 = t55 * MDP(7);
t53 = -pkin(3) - t54;
t91 = MDP(11) * t53;
t90 = MDP(5) + (MDP(12) * t62 - 0.2e1 * MDP(13) * t61) * t62;
t89 = t96 * MDP(11);
t88 = 0.2e1 * t105;
t87 = -t105 + t86;
t85 = 0.2e1 * t86;
t72 = t77 * pkin(7);
t68 = -pkin(3) - t102;
t64 = qJ(4) * t77 + t72;
t63 = (-pkin(7) - qJ(4)) * t75;
t48 = t53 - t102;
t47 = t52 * t77 + t72;
t46 = (-pkin(7) - t52) * t75;
t1 = [MDP(1) + (t76 ^ 2 + t78 ^ 2) * MDP(4) * pkin(1) ^ 2 + t52 ^ 2 * t89 + (-t88 + t91) * t53 + t48 * t85 + 0.2e1 * t93 + 0.2e1 * t92 + t99 * t104 + t90; 0; MDP(4) + t89; t93 + t92 + (t97 + t99) * MDP(10) + (-t53 * pkin(3) + t52 * t97) * MDP(11) + t90 + t105 * (pkin(3) - t53) + t86 * (t48 + t68); 0; t97 * t104 + qJ(4) ^ 2 * t89 + t68 * t85 + (t88 + t94) * pkin(3) + t90; t87 + t91; 0; t87 - t94; MDP(11); (t46 * t81 - t47 * t79) * MDP(17) + (-t46 * t79 - t47 * t81) * MDP(18) + t98; -t86; (t63 * t81 - t64 * t79) * MDP(17) + (-t63 * t79 - t64 * t81) * MDP(18) + t98; 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
