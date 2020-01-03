% Calculate Gravitation load on the joints for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:16
% EndTime: 2019-12-31 19:39:17
% DurationCPUTime: 0.35s
% Computational Cost: add. (150->65), mult. (268->97), div. (0->0), fcn. (255->8), ass. (0->41)
t111 = MDP(9) + MDP(11);
t110 = MDP(10) - MDP(13);
t88 = sin(qJ(1));
t90 = cos(qJ(1));
t71 = g(1) * t90 + g(2) * t88;
t87 = sin(qJ(2));
t109 = t71 * t87;
t78 = t87 * qJ(3);
t89 = cos(qJ(2));
t100 = t89 * pkin(2) + t78;
t84 = pkin(8) + qJ(5);
t76 = sin(t84);
t77 = cos(t84);
t97 = t89 * t76 - t87 * t77;
t55 = t97 * t88;
t96 = t87 * t76 + t89 * t77;
t56 = t96 * t88;
t101 = t89 * t90;
t102 = t87 * t90;
t57 = t76 * t101 - t77 * t102;
t58 = t96 * t90;
t108 = (g(1) * t57 + g(2) * t55 + g(3) * t96) * MDP(24) + (g(1) * t58 + g(2) * t56 - g(3) * t97) * MDP(25);
t106 = g(1) * t88;
t103 = t89 * pkin(3);
t99 = qJ(3) * t89;
t98 = pkin(2) * t101 + t88 * pkin(6) + (pkin(1) + t78) * t90;
t70 = -g(2) * t90 + t106;
t85 = sin(pkin(8));
t86 = cos(pkin(8));
t95 = t89 * t85 - t87 * t86;
t94 = t87 * t85 + t89 * t86;
t93 = -pkin(1) - t100;
t81 = t90 * pkin(6);
t74 = t90 * t99;
t72 = t88 * t99;
t64 = t94 * t90;
t63 = t95 * t90;
t62 = t94 * t88;
t61 = t95 * t88;
t59 = -g(3) * t89 + t109;
t1 = [(-g(1) * t81 - g(2) * t98 - t93 * t106) * MDP(14) + (g(1) * t62 - g(2) * t64) * MDP(15) + (-g(1) * t61 + g(2) * t63) * MDP(16) + (-g(1) * (-t90 * qJ(4) + t81) - g(2) * (pkin(3) * t101 + t98) + (-g(1) * (t93 - t103) + g(2) * qJ(4)) * t88) * MDP(18) + (g(1) * t56 - g(2) * t58) * MDP(24) + (-g(1) * t55 + g(2) * t57) * MDP(25) + (MDP(3) - MDP(12) + MDP(17)) * t71 + (-t110 * t87 + t111 * t89 + MDP(2)) * t70; (-g(1) * (-pkin(2) * t102 + t74) - g(2) * (-t88 * t87 * pkin(2) + t72) - g(3) * t100) * MDP(14) + (-g(1) * t63 - g(2) * t61 - g(3) * t94) * MDP(15) + (-g(1) * t64 - g(2) * t62 + g(3) * t95) * MDP(16) + (-g(1) * t74 - g(2) * t72 - g(3) * (t100 + t103) + (pkin(2) + pkin(3)) * t109) * MDP(18) - t108 + t110 * (g(3) * t87 + t71 * t89) + t111 * t59; (-MDP(14) - MDP(18)) * t59; t70 * MDP(18); t108;];
taug = t1;
