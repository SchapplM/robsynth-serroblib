% Calculate Gravitation load on the joints for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:32:55
% EndTime: 2019-03-08 19:32:58
% DurationCPUTime: 0.82s
% Computational Cost: add. (379->97), mult. (910->170), div. (0->0), fcn. (1135->14), ass. (0->52)
t108 = sin(pkin(11));
t114 = sin(qJ(2));
t116 = cos(qJ(2));
t144 = cos(pkin(11));
t95 = -t108 * t116 - t114 * t144;
t152 = MDP(12) - MDP(15);
t109 = sin(pkin(10));
t111 = cos(pkin(10));
t112 = cos(pkin(6));
t135 = t112 * t116;
t151 = -t109 * t135 - t111 * t114;
t119 = -t108 * t114 + t116 * t144;
t118 = t119 * t112;
t82 = t109 * t95 + t111 * t118;
t85 = -t109 * t118 + t111 * t95;
t143 = sin(pkin(6));
t124 = t144 * t143;
t127 = t114 * t143;
t92 = t108 * t127 - t116 * t124;
t150 = -g(1) * t85 - g(2) * t82 + g(3) * t92;
t133 = t95 * t112;
t86 = t109 * t133 + t111 * t119;
t81 = -t109 * t119 + t111 * t133;
t106 = pkin(12) + qJ(6);
t104 = sin(t106);
t115 = cos(qJ(4));
t142 = t104 * t115;
t105 = cos(t106);
t141 = t105 * t115;
t107 = sin(pkin(12));
t140 = t107 * t115;
t139 = t109 * t114;
t110 = cos(pkin(12));
t138 = t110 * t115;
t136 = t112 * t114;
t132 = MDP(16) + MDP(5);
t130 = t111 * t135;
t113 = sin(qJ(4));
t129 = t113 * t143;
t126 = t115 * t143;
t125 = t116 * t143;
t77 = t111 * t126 - t113 * t81;
t79 = -t109 * t126 + t113 * t86;
t93 = t108 * t125 + t114 * t124;
t87 = -t112 * t115 + t113 * t93;
t122 = g(1) * t79 + g(2) * t77 + g(3) * t87;
t117 = -g(1) * t151 - g(3) * t125;
t96 = pkin(2) * t130;
t88 = t112 * t113 + t115 * t93;
t80 = t109 * t129 + t115 * t86;
t78 = -t111 * t129 - t115 * t81;
t1 = [(-MDP(1) - t132) * g(3); (-g(2) * (t130 - t139) + t117) * MDP(3) + (-g(1) * (t109 * t136 - t111 * t116) - g(2) * (-t109 * t116 - t111 * t136) + g(3) * t127) * MDP(4) + (-g(2) * t96 + (g(2) * t139 + t117) * pkin(2)) * MDP(5) + (-g(1) * (t107 * t86 + t138 * t85) - g(2) * (-t107 * t81 + t138 * t82) - g(3) * (t107 * t93 - t138 * t92)) * MDP(13) + (-g(1) * (t110 * t86 - t140 * t85) - g(2) * (-t110 * t81 - t140 * t82) - g(3) * (t110 * t93 + t140 * t92)) * MDP(14) + (-g(1) * (pkin(2) * t151 + pkin(8) * t86) - g(2) * (-pkin(2) * t139 - t81 * pkin(8) + t96) - g(3) * (pkin(2) * t125 + t93 * pkin(8)) + t150 * (pkin(4) * t115 + qJ(5) * t113 + pkin(3))) * MDP(16) + (-g(1) * (t104 * t86 + t141 * t85) - g(2) * (-t104 * t81 + t141 * t82) - g(3) * (t104 * t93 - t141 * t92)) * MDP(22) + (-g(1) * (t105 * t86 - t142 * t85) - g(2) * (-t105 * t81 - t142 * t82) - g(3) * (t105 * t93 + t142 * t92)) * MDP(23) + (t115 * MDP(11) - t152 * t113) * t150; t132 * (-g(3) * t112 + (-g(1) * t109 + g(2) * t111) * t143); (-g(1) * (-pkin(4) * t79 + qJ(5) * t80) - g(2) * (-pkin(4) * t77 + qJ(5) * t78) - g(3) * (-pkin(4) * t87 + qJ(5) * t88)) * MDP(16) + t152 * (g(1) * t80 + g(2) * t78 + g(3) * t88) + (MDP(13) * t110 - MDP(14) * t107 + MDP(22) * t105 - MDP(23) * t104 + MDP(11)) * t122; -t122 * MDP(16); (-g(1) * (-t104 * t80 - t105 * t85) - g(2) * (-t104 * t78 - t105 * t82) - g(3) * (-t104 * t88 + t105 * t92)) * MDP(22) + (-g(1) * (t104 * t85 - t105 * t80) - g(2) * (t104 * t82 - t105 * t78) - g(3) * (-t104 * t92 - t105 * t88)) * MDP(23);];
taug  = t1;
