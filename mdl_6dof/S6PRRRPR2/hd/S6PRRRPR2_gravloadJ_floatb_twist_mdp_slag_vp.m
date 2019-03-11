% Calculate Gravitation load on the joints for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:08:59
% EndTime: 2019-03-08 23:09:00
% DurationCPUTime: 0.60s
% Computational Cost: add. (488->96), mult. (767->159), div. (0->0), fcn. (900->14), ass. (0->47)
t134 = MDP(18) - MDP(21);
t100 = sin(pkin(11));
t104 = sin(qJ(2));
t106 = cos(qJ(2));
t123 = cos(pkin(11));
t124 = cos(pkin(6));
t115 = t124 * t123;
t83 = t100 * t104 - t106 * t115;
t118 = t100 * t124;
t85 = t104 * t123 + t106 * t118;
t133 = -g(1) * t85 - g(2) * t83;
t101 = sin(pkin(6));
t130 = g(3) * t101;
t97 = pkin(12) + qJ(6);
t93 = sin(t97);
t98 = qJ(3) + qJ(4);
t96 = cos(t98);
t129 = t93 * t96;
t94 = cos(t97);
t128 = t94 * t96;
t99 = sin(pkin(12));
t127 = t96 * t99;
t102 = cos(pkin(12));
t126 = t102 * t96;
t125 = t106 * t96;
t122 = t100 * t101;
t121 = t101 * t104;
t105 = cos(qJ(3));
t120 = t101 * t105;
t119 = t101 * t106;
t117 = t101 * t123;
t84 = t100 * t106 + t104 * t115;
t86 = -t104 * t118 + t106 * t123;
t116 = g(1) * t86 + g(2) * t84;
t95 = sin(t98);
t76 = t117 * t96 + t84 * t95;
t78 = -t122 * t96 + t86 * t95;
t81 = t121 * t95 - t124 * t96;
t112 = g(1) * t78 + g(2) * t76 + g(3) * t81;
t77 = -t117 * t95 + t84 * t96;
t79 = t122 * t95 + t86 * t96;
t82 = t121 * t96 + t124 * t95;
t113 = t134 * (g(1) * t79 + g(2) * t77 + g(3) * t82) + (MDP(19) * t102 - MDP(20) * t99 + MDP(28) * t94 - MDP(29) * t93 + MDP(17)) * t112;
t109 = -g(1) * (-t78 * pkin(4) + qJ(5) * t79) - g(2) * (-t76 * pkin(4) + t77 * qJ(5)) - g(3) * (-t81 * pkin(4) + t82 * qJ(5));
t103 = sin(qJ(3));
t108 = -g(1) * (t100 * t120 - t103 * t86) - g(2) * (-t84 * t103 - t105 * t117) - g(3) * (-t103 * t121 + t105 * t124);
t1 = [(-MDP(1) - MDP(22)) * g(3); (g(3) * t121 + t116) * MDP(4) + (-g(1) * (-t126 * t85 + t86 * t99) - g(2) * (-t83 * t126 + t84 * t99) - (t102 * t125 + t104 * t99) * t130) * MDP(19) + (-g(1) * (t102 * t86 + t127 * t85) - g(2) * (t102 * t84 + t127 * t83) - (t102 * t104 - t125 * t99) * t130) * MDP(20) + ((t104 * t130 + t116) * (-pkin(9) - pkin(8)) + (-t106 * t130 - t133) * (pkin(3) * t105 + pkin(4) * t96 + qJ(5) * t95 + pkin(2))) * MDP(22) + (-g(1) * (-t128 * t85 + t86 * t93) - g(2) * (-t128 * t83 + t84 * t93) - (t104 * t93 + t125 * t94) * t130) * MDP(28) + (-g(1) * (t129 * t85 + t86 * t94) - g(2) * (t129 * t83 + t84 * t94) - (t104 * t94 - t125 * t93) * t130) * MDP(29) + (-t105 * MDP(10) + t103 * MDP(11) - t96 * MDP(17) + t134 * t95 - MDP(3)) * (g(3) * t119 + t133); t108 * MDP(10) + (-g(1) * (-t103 * t122 - t105 * t86) - g(2) * (t103 * t117 - t84 * t105) - g(3) * (-t103 * t124 - t104 * t120)) * MDP(11) + (pkin(3) * t108 + t109) * MDP(22) + t113; t109 * MDP(22) + t113; -t112 * MDP(22); (-g(1) * (-t79 * t93 + t85 * t94) - g(2) * (-t77 * t93 + t83 * t94) - g(3) * (-t119 * t94 - t82 * t93)) * MDP(28) + (-g(1) * (-t79 * t94 - t85 * t93) - g(2) * (-t77 * t94 - t83 * t93) - g(3) * (t119 * t93 - t82 * t94)) * MDP(29);];
taug  = t1;
