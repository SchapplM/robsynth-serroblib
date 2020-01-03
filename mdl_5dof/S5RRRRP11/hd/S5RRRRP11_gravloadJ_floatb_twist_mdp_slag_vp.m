% Calculate Gravitation load on the joints for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:33
% EndTime: 2019-12-31 22:18:36
% DurationCPUTime: 0.74s
% Computational Cost: add. (393->105), mult. (997->167), div. (0->0), fcn. (1211->10), ass. (0->49)
t169 = MDP(17) - MDP(26);
t151 = MDP(23) + MDP(25);
t168 = MDP(24) - MDP(27);
t131 = sin(qJ(2));
t132 = sin(qJ(1));
t135 = cos(qJ(2));
t161 = cos(pkin(5));
t167 = cos(qJ(1));
t145 = t161 * t167;
t118 = t131 * t145 + t132 * t135;
t130 = sin(qJ(3));
t134 = cos(qJ(3));
t128 = sin(pkin(5));
t149 = t128 * t167;
t103 = t118 * t134 - t130 * t149;
t117 = t131 * t132 - t135 * t145;
t129 = sin(qJ(4));
t133 = cos(qJ(4));
t89 = t103 * t129 - t117 * t133;
t90 = t103 * t133 + t117 * t129;
t148 = t132 * t161;
t119 = t167 * t131 + t135 * t148;
t170 = g(1) * t119 + g(2) * t117;
t158 = t128 * t131;
t157 = t128 * t132;
t156 = t128 * t134;
t155 = t128 * t135;
t154 = t129 * t134;
t153 = t133 * t134;
t152 = t133 * t135;
t150 = t129 * t155;
t147 = -t118 * t130 - t134 * t149;
t143 = pkin(3) * t134 + pkin(9) * t130 + pkin(2);
t116 = t161 * t130 + t131 * t156;
t100 = t116 * t129 + t128 * t152;
t120 = -t131 * t148 + t167 * t135;
t107 = t120 * t134 + t130 * t157;
t93 = t107 * t129 - t119 * t133;
t141 = g(1) * t93 + g(2) * t89 + g(3) * t100;
t109 = (t129 * t131 + t134 * t152) * t128;
t108 = -t133 * t158 + t134 * t150;
t106 = t120 * t130 - t132 * t156;
t101 = t116 * t133 - t150;
t99 = -t119 * t153 + t120 * t129;
t98 = -t119 * t154 - t120 * t133;
t97 = -t117 * t153 + t118 * t129;
t96 = -t117 * t154 - t118 * t133;
t94 = t107 * t133 + t119 * t129;
t1 = [(g(1) * t132 - g(2) * t167) * MDP(2) + (g(1) * t167 + g(2) * t132) * MDP(3) + (g(1) * t118 - g(2) * t120) * MDP(9) + (-g(1) * t117 + g(2) * t119) * MDP(10) + (g(1) * t103 - g(2) * t107) * MDP(16) + (-g(1) * (-t132 * pkin(1) - t118 * pkin(2) - pkin(3) * t103 - pkin(4) * t90 + pkin(7) * t149 - t117 * pkin(8) + pkin(9) * t147 - qJ(5) * t89) - g(2) * (t167 * pkin(1) + t120 * pkin(2) + t107 * pkin(3) + t94 * pkin(4) + pkin(7) * t157 + t119 * pkin(8) + t106 * pkin(9) + t93 * qJ(5))) * MDP(28) + t151 * (g(1) * t90 - g(2) * t94) + t168 * (-g(1) * t89 + g(2) * t93) + t169 * (g(1) * t147 + g(2) * t106); (g(1) * t120 + g(2) * t118 + g(3) * t158) * MDP(10) + (-g(1) * (pkin(4) * t99 + pkin(8) * t120 + qJ(5) * t98) - g(2) * (pkin(4) * t97 + pkin(8) * t118 + qJ(5) * t96) + t170 * t143 + (-t109 * pkin(4) - t108 * qJ(5) - (pkin(8) * t131 + t143 * t135) * t128) * g(3)) * MDP(28) + t151 * (-g(1) * t99 - g(2) * t97 - g(3) * t109) + t168 * (g(1) * t98 + g(2) * t96 + g(3) * t108) + (-t134 * MDP(16) + t130 * t169 - MDP(9)) * (g(3) * t155 - t170); (-pkin(9) * MDP(28) + t169) * (g(1) * t107 + g(2) * t103 + g(3) * t116) + (MDP(16) + t151 * t133 - t168 * t129 + MDP(28) * (pkin(4) * t133 + qJ(5) * t129 + pkin(3))) * (-g(3) * (-t130 * t158 + t161 * t134) - g(2) * t147 + g(1) * t106); (-g(1) * (-pkin(4) * t93 + qJ(5) * t94) - g(2) * (-pkin(4) * t89 + qJ(5) * t90) - g(3) * (-pkin(4) * t100 + qJ(5) * t101)) * MDP(28) + t151 * t141 + t168 * (g(1) * t94 + g(2) * t90 + g(3) * t101); -t141 * MDP(28);];
taug = t1;
