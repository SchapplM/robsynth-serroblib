% Calculate joint inertia matrix for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRPR5_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:57:14
% EndTime: 2019-12-05 17:57:15
% DurationCPUTime: 0.32s
% Computational Cost: add. (415->84), mult. (854->133), div. (0->0), fcn. (864->8), ass. (0->42)
t102 = sin(qJ(5));
t104 = cos(qJ(5));
t101 = cos(pkin(8));
t100 = cos(pkin(9));
t105 = cos(qJ(3));
t103 = sin(qJ(3));
t121 = qJ(2) * t103;
t99 = sin(pkin(8));
t124 = qJ(4) * t99;
t91 = -t101 * pkin(2) - t99 * pkin(6) - pkin(1);
t87 = t105 * t91;
t74 = -t105 * t124 + t87 + (-pkin(3) - t121) * t101;
t114 = t105 * t101 * qJ(2);
t79 = t114 + (t91 - t124) * t103;
t98 = sin(pkin(9));
t64 = t100 * t74 - t98 * t79;
t88 = t100 * t105 - t98 * t103;
t84 = t88 * t99;
t62 = -t101 * pkin(4) - t84 * pkin(7) + t64;
t65 = t100 * t79 + t98 * t74;
t89 = t100 * t103 + t98 * t105;
t83 = t89 * t99;
t63 = -t83 * pkin(7) + t65;
t68 = t102 * t84 + t104 * t83;
t69 = -t102 * t83 + t104 * t84;
t132 = t69 * MDP(19) - t68 * MDP(20) - (t102 * t62 + t104 * t63) * MDP(23) + (-t102 * t63 + t104 * t62) * MDP(22);
t131 = (t105 * MDP(10) - t103 * MDP(11)) * t99 + (-t101 * t121 + t87) * MDP(13) - (t103 * t91 + t114) * MDP(14) + t132;
t125 = (-t102 * t89 + t104 * t88) * MDP(22) - (t102 * t88 + t104 * t89) * MDP(23);
t129 = -t105 * MDP(13) + t103 * MDP(14) - t125;
t127 = pkin(3) * t98;
t90 = (pkin(3) * t103 + qJ(2)) * t99;
t122 = t99 * MDP(5);
t94 = t100 * pkin(3) + pkin(4);
t117 = (-t102 * t127 + t104 * t94) * MDP(22);
t116 = (t102 * t94 + t104 * t127) * MDP(23);
t115 = MDP(12) + MDP(21);
t111 = t68 * MDP(22) + t69 * MDP(23);
t108 = MDP(21) - t116 + t117;
t106 = qJ(2) ^ 2;
t97 = t101 ^ 2;
t96 = t99 ^ 2;
t1 = [MDP(1) - 0.2e1 * pkin(1) * t122 + (pkin(1) ^ 2 + t96 * t106) * MDP(7) + (t64 ^ 2 + t65 ^ 2 + t90 ^ 2) * MDP(16) + (t106 * MDP(7) + t115) * t97 + (MDP(17) * t69 - 0.2e1 * t68 * MDP(18)) * t69 + (MDP(8) * t105 - 0.2e1 * MDP(9) * t103) * t96 * t105 + 0.2e1 * (-t64 * t84 - t65 * t83) * MDP(15) + 0.2e1 * t111 * (t83 * pkin(4) + t90) + 0.2e1 * (t97 * MDP(6) + (t103 * MDP(13) + t105 * MDP(14) + MDP(6)) * t96) * qJ(2) + 0.2e1 * (pkin(1) * MDP(4) - t131) * t101; t122 - pkin(1) * MDP(7) + (-t89 * t83 - t88 * t84) * MDP(15) + (t64 * t88 + t65 * t89) * MDP(16) + (-MDP(4) + t129) * t101; MDP(7) + (t88 ^ 2 + t89 ^ 2) * MDP(16); (-MDP(12) - t108) * t101 + ((-t100 * t84 - t83 * t98) * MDP(15) + (t100 * t64 + t65 * t98) * MDP(16)) * pkin(3) + t131; (t100 * t88 + t89 * t98) * MDP(16) * pkin(3) - t129; (t100 ^ 2 + t98 ^ 2) * MDP(16) * pkin(3) ^ 2 + 0.2e1 * t117 - 0.2e1 * t116 + t115; t90 * MDP(16) + t111; 0; 0; MDP(16); -t101 * MDP(21) + t132; t125; t108; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
