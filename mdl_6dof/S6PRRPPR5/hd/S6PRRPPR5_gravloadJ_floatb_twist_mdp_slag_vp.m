% Calculate Gravitation load on the joints for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:21:16
% EndTime: 2019-03-08 21:21:18
% DurationCPUTime: 0.78s
% Computational Cost: add. (319->101), mult. (739->163), div. (0->0), fcn. (863->12), ass. (0->51)
t150 = -MDP(10) + MDP(13) - MDP(18);
t149 = MDP(11) - MDP(14);
t113 = sin(qJ(3));
t115 = cos(qJ(3));
t148 = -pkin(3) * t115 - qJ(4) * t113 - pkin(2);
t146 = pkin(4) + pkin(8);
t111 = sin(pkin(6));
t144 = g(3) * t111;
t143 = cos(pkin(6));
t142 = cos(pkin(10));
t140 = qJ(5) * t115;
t108 = pkin(11) + qJ(6);
t106 = sin(t108);
t139 = t106 * t113;
t107 = cos(t108);
t138 = t107 * t113;
t109 = sin(pkin(11));
t137 = t109 * t113;
t114 = sin(qJ(2));
t136 = t111 * t114;
t135 = t111 * t115;
t116 = cos(qJ(2));
t134 = t111 * t116;
t112 = cos(pkin(11));
t133 = t112 * t113;
t132 = t113 * t116;
t131 = t115 * t116;
t130 = MDP(15) + MDP(19);
t110 = sin(pkin(10));
t121 = t143 * t142;
t91 = t110 * t114 - t116 * t121;
t129 = t148 * t91;
t124 = t110 * t143;
t93 = t142 * t114 + t116 * t124;
t128 = t148 * t93;
t123 = t111 * t142;
t92 = t110 * t116 + t114 * t121;
t80 = t92 * t113 + t115 * t123;
t81 = -t113 * t123 + t115 * t92;
t127 = -t80 * pkin(3) + qJ(4) * t81;
t94 = -t114 * t124 + t142 * t116;
t82 = -t110 * t135 + t113 * t94;
t83 = t110 * t111 * t113 + t115 * t94;
t126 = -t82 * pkin(3) + qJ(4) * t83;
t95 = t113 * t136 - t143 * t115;
t96 = t143 * t113 + t114 * t135;
t125 = -t95 * pkin(3) + qJ(4) * t96;
t122 = pkin(2) * t134 + pkin(8) * t136 + (pkin(3) * t131 + qJ(4) * t132) * t111;
t120 = g(1) * t82 + g(2) * t80 + g(3) * t95;
t119 = g(1) * t83 + g(2) * t81 + g(3) * t96;
t1 = [(-MDP(1) - t130) * g(3); (-g(1) * (pkin(8) * t94 + t128) - g(2) * (pkin(8) * t92 + t129) - g(3) * t122) * MDP(15) + (-g(1) * (t112 * t94 - t93 * t137) - g(2) * (t112 * t92 - t91 * t137) - (t109 * t132 + t112 * t114) * t144) * MDP(16) + (-g(1) * (-t109 * t94 - t93 * t133) - g(2) * (-t109 * t92 - t91 * t133) - (-t109 * t114 + t112 * t132) * t144) * MDP(17) + (-g(1) * (-t93 * t140 + t146 * t94 + t128) - g(2) * (-t91 * t140 + t146 * t92 + t129) - g(3) * ((pkin(4) * t114 + qJ(5) * t131) * t111 + t122)) * MDP(19) + (-g(1) * (t107 * t94 - t93 * t139) - g(2) * (t107 * t92 - t91 * t139) - (t106 * t132 + t107 * t114) * t144) * MDP(25) + (-g(1) * (-t106 * t94 - t93 * t138) - g(2) * (-t106 * t92 - t91 * t138) - (-t106 * t114 + t107 * t132) * t144) * MDP(26) + (MDP(4) - MDP(12)) * (g(1) * t94 + g(2) * t92 + g(3) * t136) + (t149 * t113 + t150 * t115 - MDP(3)) * (-g(1) * t93 - g(2) * t91 + g(3) * t134); (-g(1) * t126 - g(2) * t127 - g(3) * t125) * MDP(15) + (-g(1) * (-qJ(5) * t82 + t126) - g(2) * (-t80 * qJ(5) + t127) - g(3) * (-qJ(5) * t95 + t125)) * MDP(19) - t150 * t120 + (-MDP(16) * t109 - MDP(17) * t112 - MDP(25) * t106 - MDP(26) * t107 + t149) * t119; -t130 * t120; -t119 * MDP(19); (-g(1) * (-t106 * t93 + t107 * t82) - g(2) * (-t106 * t91 + t107 * t80) - g(3) * (t106 * t134 + t95 * t107)) * MDP(25) + (-g(1) * (-t106 * t82 - t107 * t93) - g(2) * (-t106 * t80 - t107 * t91) - g(3) * (-t95 * t106 + t107 * t134)) * MDP(26);];
taug  = t1;
