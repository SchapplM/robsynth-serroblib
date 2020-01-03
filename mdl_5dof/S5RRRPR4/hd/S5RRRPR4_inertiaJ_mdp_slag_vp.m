% Calculate joint inertia matrix for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRPR4_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:40
% EndTime: 2019-12-31 21:11:41
% DurationCPUTime: 0.27s
% Computational Cost: add. (247->87), mult. (410->118), div. (0->0), fcn. (312->6), ass. (0->47)
t103 = sin(qJ(3));
t100 = t103 ^ 2;
t106 = cos(qJ(3));
t121 = t106 ^ 2 + t100;
t104 = sin(qJ(2));
t89 = t104 * pkin(1) + pkin(7);
t129 = t121 * t89;
t139 = 0.2e1 * t106;
t138 = t106 * pkin(3) + t103 * qJ(4);
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t77 = t103 * t102 + t106 * t105;
t78 = -t106 * t102 + t103 * t105;
t137 = t77 * MDP(23) + t78 * MDP(24);
t131 = t78 * MDP(20) - t77 * MDP(21);
t136 = 2 * MDP(15);
t135 = pkin(7) - pkin(8);
t107 = cos(qJ(2));
t90 = -t107 * pkin(1) - pkin(2);
t134 = pkin(2) - t90;
t133 = -pkin(8) + t89;
t119 = pkin(2) + t138;
t73 = t90 - t138;
t130 = -t73 + t119;
t128 = t121 * pkin(7);
t108 = -pkin(3) - pkin(4);
t123 = (t102 * qJ(4) - t105 * t108) * MDP(23);
t122 = (t105 * qJ(4) + t102 * t108) * MDP(24);
t120 = MDP(17) * t103;
t118 = -pkin(3) * MDP(17) - MDP(14);
t117 = (-t103 * pkin(3) + t106 * qJ(4)) * MDP(15) + t106 * MDP(10) + t103 * MDP(9) - t131;
t116 = t103 * MDP(8) * t139 + t100 * MDP(7) + MDP(4) + (MDP(18) * t78 - 0.2e1 * MDP(19) * t77) * t78;
t75 = t133 * t103;
t76 = t133 * t106;
t115 = (t102 * t76 - t105 * t75) * MDP(23) + (t102 * t75 + t105 * t76) * MDP(24);
t83 = t135 * t103;
t84 = t135 * t106;
t114 = (t102 * t84 - t105 * t83) * MDP(23) + (t102 * t83 + t105 * t84) * MDP(24);
t113 = t105 * MDP(23) - t102 * MDP(24);
t112 = (t107 * MDP(5) - t104 * MDP(6)) * pkin(1);
t111 = 0.2e1 * t137;
t110 = (MDP(17) * qJ(4) - MDP(13) + MDP(16)) * t106 + (-MDP(12) + t118) * t103;
t98 = t106 * pkin(4);
t91 = t103 * MDP(15);
t74 = t98 + t119;
t67 = -t73 + t98;
t1 = [MDP(1) + t129 * t136 + (t121 * t89 ^ 2 + t73 ^ 2) * MDP(17) + t67 * t111 + t116 + (-t90 * MDP(12) - t73 * MDP(14)) * t139 + 0.2e1 * (t90 * MDP(13) - t73 * MDP(16)) * t103 + 0.2e1 * t112; (t128 + t129) * MDP(15) + (pkin(7) * t129 - t119 * t73) * MDP(17) + t112 + (t134 * MDP(12) + t130 * MDP(14)) * t106 + (-t134 * MDP(13) + t130 * MDP(16)) * t103 + t116 + t137 * (t67 + t74); t128 * t136 + t121 * MDP(17) * pkin(7) ^ 2 - (-0.2e1 * t106 * MDP(14) - 0.2e1 * t103 * MDP(16) - MDP(17) * t119) * t119 + t74 * t111 + 0.2e1 * (t106 * MDP(12) - t103 * MDP(13)) * pkin(2) + t116; t110 * t89 + t115 + t117; t110 * pkin(7) + t114 + t117; MDP(11) + 0.2e1 * pkin(3) * MDP(14) + 0.2e1 * qJ(4) * MDP(16) + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(17) + MDP(22) + 0.2e1 * t123 + 0.2e1 * t122; t89 * t120 + t91; pkin(7) * t120 + t91; -t113 + t118; MDP(17); -t115 + t131; -t114 + t131; -MDP(22) - t122 - t123; t113; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
