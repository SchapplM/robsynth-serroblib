% Calculate Gravitation load on the joints for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:05
% EndTime: 2019-12-31 22:06:06
% DurationCPUTime: 0.37s
% Computational Cost: add. (283->71), mult. (398->104), div. (0->0), fcn. (398->8), ass. (0->39)
t116 = MDP(23) + MDP(25);
t118 = MDP(24) - MDP(27);
t85 = sin(qJ(1));
t88 = cos(qJ(1));
t96 = g(1) * t88 + g(2) * t85;
t117 = MDP(10) - MDP(26);
t84 = sin(qJ(2));
t109 = g(3) * t84;
t86 = cos(qJ(3));
t101 = t88 * t86;
t87 = cos(qJ(2));
t105 = t85 * t87;
t83 = sin(qJ(3));
t71 = t83 * t105 + t101;
t102 = t88 * t83;
t73 = -t87 * t102 + t85 * t86;
t114 = -g(1) * t73 + g(2) * t71 + t83 * t109;
t82 = qJ(3) + qJ(4);
t80 = sin(t82);
t107 = t80 * t84;
t81 = cos(t82);
t106 = t81 * t84;
t104 = t88 * t80;
t103 = t88 * t81;
t99 = pkin(3) * t83 + pkin(6);
t66 = t80 * t105 + t103;
t68 = t87 * t104 - t85 * t81;
t60 = g(1) * t68 + g(2) * t66 + g(3) * t107;
t67 = t81 * t105 - t104;
t69 = t87 * t103 + t85 * t80;
t98 = t118 * (g(1) * t69 + g(2) * t67 + g(3) * t106) + t116 * t60;
t79 = t86 * pkin(3) + pkin(2);
t89 = -pkin(8) - pkin(7);
t94 = t87 * t79 - t84 * t89 + pkin(1);
t93 = pkin(4) * t81 + qJ(5) * t80 + t79;
t90 = -g(1) * (-t68 * pkin(4) + t69 * qJ(5)) - g(2) * (-t66 * pkin(4) + t67 * qJ(5)) - g(3) * (-pkin(4) * t107 + qJ(5) * t106);
t74 = t87 * t101 + t85 * t83;
t72 = -t86 * t105 + t102;
t1 = [t96 * MDP(3) + (-g(1) * t72 - g(2) * t74) * MDP(16) + (-g(1) * t71 - g(2) * t73) * MDP(17) + (-g(1) * (-t67 * pkin(4) - t66 * qJ(5)) - g(2) * (t69 * pkin(4) + t68 * qJ(5)) + (-g(1) * t99 - g(2) * t94) * t88 + (g(1) * t94 - g(2) * t99) * t85) * MDP(28) - t118 * (g(1) * t66 - g(2) * t68) + t116 * (g(1) * t67 - g(2) * t69) + (t87 * MDP(9) - t117 * t84 + MDP(2)) * (g(1) * t85 - g(2) * t88); ((-g(3) * t93 + t96 * t89) * t87 + (g(3) * t89 + t96 * t93) * t84) * MDP(28) + t117 * (t96 * t87 + t109) + (t86 * MDP(16) - t83 * MDP(17) + t116 * t81 - t118 * t80 + MDP(9)) * (-g(3) * t87 + t96 * t84); t114 * MDP(16) + (g(1) * t74 - g(2) * t72 + t86 * t109) * MDP(17) + (t114 * pkin(3) + t90) * MDP(28) + t98; MDP(28) * t90 + t98; -t60 * MDP(28);];
taug = t1;
