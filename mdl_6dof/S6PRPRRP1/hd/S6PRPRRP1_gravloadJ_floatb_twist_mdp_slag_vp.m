% Calculate Gravitation load on the joints for
% S6PRPRRP1
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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:24:20
% EndTime: 2021-01-16 01:24:23
% DurationCPUTime: 0.60s
% Computational Cost: add. (560->101), mult. (672->163), div. (0->0), fcn. (761->18), ass. (0->69)
t125 = qJ(2) + pkin(11);
t123 = sin(t125);
t131 = cos(pkin(6));
t137 = cos(qJ(4));
t128 = sin(pkin(6));
t134 = sin(qJ(4));
t154 = t128 * t134;
t103 = t123 * t154 - t131 * t137;
t124 = cos(t125);
t127 = sin(pkin(10));
t113 = t127 * t124;
t130 = cos(pkin(10));
t151 = t130 * t131;
t101 = t123 * t151 + t113;
t153 = t128 * t137;
t90 = t101 * t134 + t130 * t153;
t114 = t130 * t124;
t155 = t127 * t131;
t98 = t123 * t155 - t114;
t92 = t127 * t153 + t134 * t98;
t167 = g(1) * t92 - g(2) * t90 - g(3) * t103;
t166 = -MDP(12) + MDP(22);
t165 = MDP(18) + MDP(20);
t164 = MDP(19) + MDP(21);
t162 = t127 / 0.2e1;
t161 = -t130 / 0.2e1;
t160 = g(3) * t128;
t133 = sin(qJ(5));
t121 = pkin(6) + t125;
t115 = sin(t121);
t122 = pkin(6) - t125;
t116 = sin(t122);
t107 = -t115 + t116;
t95 = t107 * t161 + t113;
t159 = t95 * t133;
t97 = t107 * t162 + t114;
t158 = t97 * t133;
t117 = cos(t121);
t118 = cos(t122);
t106 = t118 / 0.2e1 - t117 / 0.2e1;
t157 = t106 * t133;
t156 = t127 * t123;
t152 = t130 * t123;
t135 = sin(qJ(2));
t150 = t131 * t135;
t138 = cos(qJ(2));
t149 = t131 * t138;
t148 = t133 * t137;
t136 = cos(qJ(5));
t147 = t136 * t137;
t146 = MDP(23) + MDP(5);
t145 = t124 * t153;
t89 = -t127 * t154 + t137 * t98;
t120 = pkin(5) * t136 + pkin(4);
t132 = -qJ(6) - pkin(9);
t143 = t137 * t120 - t134 * t132;
t129 = cos(pkin(11));
t126 = sin(pkin(11));
t110 = -pkin(3) * t126 + pkin(8) * t129;
t109 = pkin(3) * t129 + pkin(8) * t126 + pkin(2);
t108 = t117 + t118;
t105 = -t116 / 0.2e1 - t115 / 0.2e1;
t104 = t123 * t153 + t131 * t134;
t100 = t124 * t151 - t156;
t99 = t124 * t155 + t152;
t96 = t108 * t162 + t152;
t94 = t108 * t161 + t156;
t91 = t101 * t137 - t130 * t154;
t1 = [(-MDP(1) - t146) * g(3); (-g(1) * (t127 * t150 - t130 * t138) - g(2) * (-t127 * t138 - t130 * t150) + t135 * t160) * MDP(4) + (-g(1) * (pkin(5) * t158 - (t130 * t109 + t110 * t155) * t135 + (-t109 * t155 + t130 * t110) * t138 - t143 * t99) - g(2) * (pkin(5) * t159 - (t109 * t127 - t110 * t151) * t135 + (t109 * t151 + t110 * t127) * t138 + t143 * t100) - g(3) * pkin(5) * t157 - (t109 * t138 + t110 * t135 + t143 * t124) * t160) * MDP(23) + t165 * (-g(1) * (-t99 * t147 + t158) - g(2) * (t100 * t147 + t159) - g(3) * (t136 * t145 + t157)) + t164 * (-g(1) * (t136 * t97 + t99 * t148) - g(2) * (-t100 * t148 + t136 * t95) - g(3) * (t106 * t136 - t133 * t145)) + (pkin(2) * MDP(5) + MDP(3)) * (-g(1) * (-t127 * t149 - t130 * t135) - g(2) * (-t127 * t135 + t130 * t149) - t138 * t160) + (t137 * MDP(11) + t166 * t134) * (g(1) * t99 - g(2) * t100 - t124 * t160); t146 * (-g(3) * t131 + (-g(1) * t127 + g(2) * t130) * t128); (-g(1) * (t120 * t92 + t132 * t89) - g(2) * (-t120 * t90 - t132 * t91) - g(3) * (-t103 * t120 - t104 * t132)) * MDP(23) + t166 * (g(1) * t89 - g(2) * t91 - g(3) * t104) + (t164 * t133 - t165 * t136 - MDP(11)) * t167; t164 * (-g(1) * (-t133 * t96 + t136 * t89) - g(2) * (-t133 * t94 - t136 * t91) - g(3) * (-t104 * t136 - t105 * t133)) + (MDP(23) * pkin(5) + t165) * (-g(1) * (t133 * t89 + t136 * t96) - g(2) * (-t133 * t91 + t136 * t94) - g(3) * (-t104 * t133 + t105 * t136)); t167 * MDP(23);];
taug = t1;
