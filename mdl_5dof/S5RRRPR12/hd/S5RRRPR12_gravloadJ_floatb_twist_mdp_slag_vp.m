% Calculate Gravitation load on the joints for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:36
% EndTime: 2019-12-31 21:40:39
% DurationCPUTime: 0.75s
% Computational Cost: add. (318->111), mult. (712->186), div. (0->0), fcn. (847->12), ass. (0->45)
t116 = MDP(17) - MDP(20);
t87 = sin(qJ(2));
t88 = sin(qJ(1));
t90 = cos(qJ(2));
t100 = cos(pkin(5));
t112 = cos(qJ(1));
t95 = t100 * t112;
t72 = t87 * t95 + t88 * t90;
t86 = sin(qJ(3));
t89 = cos(qJ(3));
t84 = sin(pkin(5));
t99 = t84 * t112;
t64 = t72 * t89 - t86 * t99;
t71 = t87 * t88 - t90 * t95;
t82 = pkin(10) + qJ(5);
t80 = sin(t82);
t81 = cos(t82);
t119 = t64 * t80 - t71 * t81;
t118 = t64 * t81 + t71 * t80;
t98 = t88 * t100;
t73 = t112 * t87 + t90 * t98;
t117 = -g(1) * t73 - g(2) * t71;
t113 = g(3) * t84;
t109 = t80 * t89;
t108 = t81 * t89;
t83 = sin(pkin(10));
t107 = t83 * t89;
t106 = t84 * t87;
t105 = t84 * t88;
t104 = t84 * t89;
t103 = t84 * t90;
t85 = cos(pkin(10));
t102 = t85 * t89;
t101 = t89 * t90;
t74 = t112 * t90 - t87 * t98;
t96 = -g(1) * t74 - g(2) * t72;
t63 = t72 * t86 + t89 * t99;
t67 = -t88 * t104 + t74 * t86;
t69 = -t100 * t89 + t86 * t106;
t93 = g(1) * t67 + g(2) * t63 + g(3) * t69;
t70 = t100 * t86 + t87 * t104;
t68 = t86 * t105 + t74 * t89;
t61 = t68 * t81 + t73 * t80;
t60 = -t68 * t80 + t73 * t81;
t1 = [(g(1) * t88 - g(2) * t112) * MDP(2) + (g(1) * t112 + g(2) * t88) * MDP(3) + (g(1) * t72 - g(2) * t74) * MDP(9) + (-g(1) * t71 + g(2) * t73) * MDP(10) + (g(1) * t64 - g(2) * t68) * MDP(16) + (-g(1) * (-t64 * t85 - t71 * t83) - g(2) * (t68 * t85 + t73 * t83)) * MDP(18) + (-g(1) * (t64 * t83 - t71 * t85) - g(2) * (-t68 * t83 + t73 * t85)) * MDP(19) + (-g(1) * (-t88 * pkin(1) - t72 * pkin(2) - pkin(3) * t64 + pkin(7) * t99 - t71 * pkin(8) - qJ(4) * t63) - g(2) * (t112 * pkin(1) + t74 * pkin(2) + t68 * pkin(3) + pkin(7) * t105 + t73 * pkin(8) + t67 * qJ(4))) * MDP(21) + (g(1) * t118 - g(2) * t61) * MDP(27) + (-g(1) * t119 - g(2) * t60) * MDP(28) + t116 * (-g(1) * t63 + g(2) * t67); (g(3) * t106 - t96) * MDP(10) + (-g(1) * (-t73 * t102 + t74 * t83) - g(2) * (-t71 * t102 + t72 * t83) - (t85 * t101 + t83 * t87) * t113) * MDP(18) + (-g(1) * (t73 * t107 + t74 * t85) - g(2) * (t71 * t107 + t72 * t85) - (-t83 * t101 + t85 * t87) * t113) * MDP(19) + ((-t87 * t113 + t96) * pkin(8) + (-t90 * t113 - t117) * (pkin(3) * t89 + qJ(4) * t86 + pkin(2))) * MDP(21) + (-g(1) * (-t73 * t108 + t74 * t80) - g(2) * (-t71 * t108 + t72 * t80) - (t81 * t101 + t80 * t87) * t113) * MDP(27) + (-g(1) * (t73 * t109 + t74 * t81) - g(2) * (t71 * t109 + t72 * t81) - (-t80 * t101 + t81 * t87) * t113) * MDP(28) + (-t89 * MDP(16) + t116 * t86 - MDP(9)) * (g(3) * t103 + t117); (-g(1) * (-pkin(3) * t67 + qJ(4) * t68) - g(2) * (-pkin(3) * t63 + qJ(4) * t64) - g(3) * (-pkin(3) * t69 + qJ(4) * t70)) * MDP(21) + t116 * (g(1) * t68 + g(2) * t64 + g(3) * t70) + (t85 * MDP(18) - t83 * MDP(19) + t81 * MDP(27) - t80 * MDP(28) + MDP(16)) * t93; -t93 * MDP(21); (-g(1) * t60 + g(2) * t119 - g(3) * (-t81 * t103 - t70 * t80)) * MDP(27) + (g(1) * t61 + g(2) * t118 - g(3) * (t80 * t103 - t70 * t81)) * MDP(28);];
taug = t1;
