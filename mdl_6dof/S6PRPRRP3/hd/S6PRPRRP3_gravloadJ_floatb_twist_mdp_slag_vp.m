% Calculate Gravitation load on the joints for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:38:33
% EndTime: 2021-01-16 01:38:36
% DurationCPUTime: 0.69s
% Computational Cost: add. (401->93), mult. (690->150), div. (0->0), fcn. (795->11), ass. (0->54)
t122 = cos(pkin(10));
t98 = sin(pkin(6));
t118 = t98 * t122;
t100 = cos(pkin(6));
t104 = sin(qJ(2));
t114 = t122 * t104;
t106 = cos(qJ(2));
t97 = sin(pkin(10));
t123 = t97 * t106;
t85 = t100 * t114 + t123;
t96 = pkin(11) + qJ(4);
t94 = sin(t96);
t95 = cos(t96);
t77 = t95 * t118 + t85 * t94;
t133 = t97 * t98;
t113 = t122 * t106;
t124 = t97 * t104;
t83 = t100 * t124 - t113;
t79 = t95 * t133 + t83 * t94;
t130 = t104 * t98;
t81 = -t100 * t95 + t94 * t130;
t139 = g(1) * t79 - g(2) * t77 - g(3) * t81;
t138 = -MDP(14) + MDP(24);
t137 = MDP(20) + MDP(22);
t136 = MDP(21) + MDP(23);
t134 = g(3) * t98;
t132 = t100 * t97;
t103 = sin(qJ(5));
t131 = t103 * t95;
t105 = cos(qJ(5));
t129 = t105 * t95;
t128 = t83 * t103;
t127 = t85 * t103;
t126 = t97 * qJ(3);
t102 = qJ(3) + pkin(8);
t125 = t97 * t102;
t121 = t103 * t106;
t120 = t105 * t106;
t119 = MDP(25) + MDP(7);
t76 = -t94 * t133 + t83 * t95;
t117 = t100 * t122;
t116 = t122 * qJ(3);
t115 = t122 * t102;
t101 = -qJ(6) - pkin(9);
t93 = pkin(5) * t105 + pkin(4);
t111 = t94 * t101 - t95 * t93;
t84 = -t100 * t113 + t124;
t86 = t100 * t123 + t114;
t107 = -g(1) * t86 - g(2) * t84 + t106 * t134;
t99 = cos(pkin(11));
t92 = pkin(3) * t99 + pkin(2);
t82 = t100 * t94 + t95 * t130;
t78 = -t94 * t118 + t85 * t95;
t1 = [(-MDP(1) - t119) * g(3); (-g(1) * (-(t122 * pkin(2) + t100 * t126) * t104 + (-pkin(2) * t132 + t116) * t106) - g(2) * (-(t97 * pkin(2) - t100 * t116) * t104 + (pkin(2) * t117 + t126) * t106) - (pkin(2) * t106 + qJ(3) * t104) * t134) * MDP(7) + (-g(1) * (-pkin(5) * t128 - (t100 * t125 + t122 * t92) * t104 + (-t92 * t132 + t115) * t106 + t111 * t86) - g(2) * (pkin(5) * t127 - (-t100 * t115 + t97 * t92) * t104 + (t92 * t117 + t125) * t106 + t111 * t84) - ((t103 * pkin(5) + t102) * t104 + (-t111 + t92) * t106) * t134) * MDP(25) + t137 * (-g(1) * (-t86 * t129 - t128) - g(2) * (-t84 * t129 + t127) - (t103 * t104 + t95 * t120) * t134) + t136 * (-g(1) * (-t105 * t83 + t86 * t131) - g(2) * (t105 * t85 + t84 * t131) - (t104 * t105 - t95 * t121) * t134) + (MDP(4) - MDP(6)) * (-g(1) * t83 + g(2) * t85 + g(3) * t130) + (-t95 * MDP(13) - t99 * MDP(5) - t138 * t94 - MDP(3)) * t107; t119 * t107; (-g(1) * (t101 * t76 + t79 * t93) - g(2) * (-t101 * t78 - t77 * t93) - g(3) * (-t101 * t82 - t81 * t93)) * MDP(25) + t138 * (g(1) * t76 - g(2) * t78 - g(3) * t82) + (t103 * t136 - t105 * t137 - MDP(13)) * t139; t136 * (-g(1) * (-t103 * t86 + t105 * t76) - g(2) * (-t103 * t84 - t105 * t78) - g(3) * (-t82 * t105 + t98 * t121)) + (MDP(25) * pkin(5) + t137) * (-g(1) * (t103 * t76 + t105 * t86) - g(2) * (-t103 * t78 + t105 * t84) - g(3) * (-t82 * t103 - t98 * t120)); t139 * MDP(25);];
taug = t1;
