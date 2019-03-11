% Calculate Gravitation load on the joints for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:01:10
% EndTime: 2019-03-09 11:01:12
% DurationCPUTime: 1.08s
% Computational Cost: add. (527->140), mult. (872->222), div. (0->0), fcn. (1020->14), ass. (0->55)
t147 = MDP(10) - MDP(13);
t143 = MDP(21) - MDP(24);
t103 = pkin(11) + qJ(4);
t100 = cos(t103);
t106 = sin(pkin(6));
t140 = cos(qJ(1));
t125 = t106 * t140;
t110 = sin(qJ(2));
t111 = sin(qJ(1));
t112 = cos(qJ(2));
t133 = cos(pkin(6));
t118 = t133 * t140;
t86 = t110 * t118 + t111 * t112;
t98 = sin(t103);
t78 = t86 * t100 - t98 * t125;
t85 = t111 * t110 - t112 * t118;
t102 = pkin(12) + qJ(6);
t97 = sin(t102);
t99 = cos(t102);
t146 = t78 * t97 - t85 * t99;
t145 = t78 * t99 + t85 * t97;
t123 = t111 * t133;
t87 = t140 * t110 + t112 * t123;
t144 = -g(1) * t87 - g(2) * t85;
t139 = g(3) * t106;
t136 = t100 * t97;
t135 = t100 * t99;
t128 = t106 * t111;
t134 = t140 * pkin(1) + pkin(8) * t128;
t104 = sin(pkin(12));
t132 = t100 * t104;
t107 = cos(pkin(12));
t131 = t100 * t107;
t130 = t100 * t112;
t129 = t106 * t110;
t127 = t106 * t112;
t105 = sin(pkin(11));
t126 = t105 * t128;
t124 = -t111 * pkin(1) + pkin(8) * t125;
t122 = t105 * t125;
t88 = -t110 * t123 + t140 * t112;
t119 = g(1) * t88 + g(2) * t86;
t77 = t100 * t125 + t86 * t98;
t81 = -t100 * t128 + t88 * t98;
t83 = -t133 * t100 + t98 * t129;
t116 = g(1) * t81 + g(2) * t77 + g(3) * t83;
t113 = g(3) * t127 + t144;
t109 = -pkin(9) - qJ(3);
t108 = cos(pkin(11));
t96 = t108 * pkin(3) + pkin(2);
t84 = t100 * t129 + t133 * t98;
t82 = t88 * t100 + t98 * t128;
t74 = t82 * t99 + t87 * t97;
t73 = -t82 * t97 + t87 * t99;
t1 = [(g(1) * t111 - g(2) * t140) * MDP(2) + (g(1) * t140 + g(2) * t111) * MDP(3) + (g(1) * t86 - g(2) * t88) * MDP(9) + (-g(1) * (-t86 * t108 + t122) - g(2) * (t88 * t108 + t126)) * MDP(11) + (-g(1) * (t86 * t105 + t108 * t125) - g(2) * (-t88 * t105 + t108 * t128)) * MDP(12) + (-g(1) * (-t86 * pkin(2) - t85 * qJ(3) + t124) - g(2) * (t88 * pkin(2) + t87 * qJ(3) + t134)) * MDP(14) + (g(1) * t78 - g(2) * t82) * MDP(20) + (-g(1) * (-t85 * t104 - t107 * t78) - g(2) * (t87 * t104 + t82 * t107)) * MDP(22) + (-g(1) * (t104 * t78 - t85 * t107) - g(2) * (-t82 * t104 + t87 * t107)) * MDP(23) + (-g(1) * (pkin(3) * t122 - pkin(4) * t78 - qJ(5) * t77 + t85 * t109 - t86 * t96 + t124) - g(2) * (pkin(3) * t126 + t82 * pkin(4) + t81 * qJ(5) - t87 * t109 + t88 * t96 + t134)) * MDP(25) + (g(1) * t145 - g(2) * t74) * MDP(31) + (-g(1) * t146 - g(2) * t73) * MDP(32) + t143 * (-g(1) * t77 + g(2) * t81) - t147 * (g(1) * t85 - g(2) * t87); (-g(1) * (-t87 * pkin(2) + t88 * qJ(3)) - g(2) * (-t85 * pkin(2) + t86 * qJ(3)) - (pkin(2) * t112 + qJ(3) * t110) * t139) * MDP(14) + (-g(1) * (t88 * t104 - t87 * t131) - g(2) * (t86 * t104 - t85 * t131) - (t104 * t110 + t107 * t130) * t139) * MDP(22) + (-g(1) * (t88 * t107 + t87 * t132) - g(2) * (t86 * t107 + t85 * t132) - (-t104 * t130 + t107 * t110) * t139) * MDP(23) + ((t110 * t139 + t119) * t109 + (-t112 * t139 - t144) * (pkin(4) * t100 + qJ(5) * t98 + t96)) * MDP(25) + (-g(1) * (-t87 * t135 + t88 * t97) - g(2) * (-t85 * t135 + t86 * t97) - (t110 * t97 + t99 * t130) * t139) * MDP(31) + (-g(1) * (t87 * t136 + t88 * t99) - g(2) * (t85 * t136 + t86 * t99) - (t110 * t99 - t97 * t130) * t139) * MDP(32) + t147 * (g(3) * t129 + t119) + (-t108 * MDP(11) + t105 * MDP(12) - t100 * MDP(20) + t143 * t98 - MDP(9)) * t113; (MDP(14) + MDP(25)) * t113; (-g(1) * (-t81 * pkin(4) + t82 * qJ(5)) - g(2) * (-t77 * pkin(4) + t78 * qJ(5)) - g(3) * (-t83 * pkin(4) + t84 * qJ(5))) * MDP(25) + t143 * (g(1) * t82 + g(2) * t78 + g(3) * t84) + (MDP(22) * t107 - MDP(23) * t104 + MDP(31) * t99 - MDP(32) * t97 + MDP(20)) * t116; -t116 * MDP(25); (-g(1) * t73 + g(2) * t146 - g(3) * (-t99 * t127 - t84 * t97)) * MDP(31) + (g(1) * t74 + g(2) * t145 - g(3) * (t97 * t127 - t84 * t99)) * MDP(32);];
taug  = t1;
