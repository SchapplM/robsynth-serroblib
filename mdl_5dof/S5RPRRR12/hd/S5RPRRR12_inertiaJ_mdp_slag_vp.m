% Calculate joint inertia matrix for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR12_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:10
% EndTime: 2019-12-31 19:13:11
% DurationCPUTime: 0.27s
% Computational Cost: add. (262->81), mult. (461->113), div. (0->0), fcn. (455->6), ass. (0->40)
t84 = sin(qJ(5));
t87 = cos(qJ(5));
t104 = t84 * MDP(23) + t87 * MDP(24);
t113 = t84 * MDP(26) + t87 * MDP(27);
t94 = t87 * MDP(26) - t84 * MDP(27);
t108 = cos(qJ(4));
t85 = sin(qJ(4));
t86 = sin(qJ(3));
t88 = cos(qJ(3));
t68 = -t108 * t88 + t85 * t86;
t111 = t68 ^ 2;
t69 = t108 * t86 + t85 * t88;
t66 = t69 ^ 2;
t89 = -pkin(1) - pkin(6);
t109 = -pkin(7) + t89;
t107 = (pkin(1) * MDP(6));
t71 = t109 * t86;
t72 = t109 * t88;
t56 = -t108 * t72 + t85 * t71;
t53 = t56 * t84;
t106 = t56 * t87;
t105 = t84 * t87;
t74 = t86 * pkin(3) + qJ(2);
t99 = MDP(22) * t105;
t82 = t84 ^ 2;
t98 = t82 * MDP(21) + MDP(18) + 0.2e1 * t99;
t97 = pkin(4) * t68 - pkin(8) * t69;
t76 = t85 * pkin(3) + pkin(8);
t77 = -t108 * pkin(3) - pkin(4);
t96 = -t68 * t77 - t69 * t76;
t95 = -MDP(23) * t87 + MDP(24) * t84;
t92 = -t69 * MDP(20) + (-MDP(19) - t94) * t68;
t57 = t108 * t71 + t85 * t72;
t83 = t87 ^ 2;
t91 = -t56 * MDP(19) - t57 * MDP(20) + (-MDP(17) + t104) * t69 + ((t82 - t83) * MDP(22) - MDP(21) * t105 - MDP(16)) * t68;
t90 = (t108 * MDP(19) - t85 * MDP(20)) * pkin(3);
t52 = t69 * pkin(4) + t68 * pkin(8) + t74;
t50 = t84 * t52 + t87 * t57;
t49 = t87 * t52 - t84 * t57;
t1 = [0.2e1 * t74 * t69 * MDP(19) + t66 * MDP(25) + MDP(1) + (MDP(7) * t88 - 0.2e1 * t86 * MDP(8)) * t88 + ((-2 * MDP(4) + t107) * pkin(1)) + 0.2e1 * (-t74 * MDP(20) + (MDP(15) + t95) * t69) * t68 + 0.2e1 * (t49 * t69 - t68 * t53) * MDP(26) + 0.2e1 * (-t68 * t106 - t50 * t69) * MDP(27) + (t83 * MDP(21) + MDP(14) - 0.2e1 * t99) * t111 + (0.2e1 * t86 * MDP(12) + 0.2e1 * t88 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2); MDP(4) - t107 + t113 * (-t66 - t111); MDP(6); (t96 * t84 - t106) * MDP(26) + (t96 * t87 + t53) * MDP(27) + (t89 * MDP(12) + MDP(9)) * t88 + (-t89 * MDP(13) - MDP(10)) * t86 + t91; t88 * MDP(12) - t86 * MDP(13) + t92; -0.2e1 * t77 * t94 + MDP(11) + 0.2e1 * t90 + t98; (t97 * t84 - t106) * MDP(26) + (t97 * t87 + t53) * MDP(27) + t91; t92; t90 + t98 + t94 * (pkin(4) - t77); 0.2e1 * pkin(4) * t94 + t98; t69 * MDP(25) + t49 * MDP(26) - t50 * MDP(27) + t95 * t68; -t113 * t69; -t113 * t76 + t104; -pkin(8) * t113 + t104; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
