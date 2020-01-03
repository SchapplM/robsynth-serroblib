% Calculate Gravitation load on the joints for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:31
% EndTime: 2019-12-31 19:24:33
% DurationCPUTime: 0.51s
% Computational Cost: add. (256->95), mult. (672->137), div. (0->0), fcn. (727->8), ass. (0->59)
t129 = sin(pkin(5));
t132 = sin(qJ(2));
t168 = t129 * t132;
t119 = qJ(3) * t168;
t134 = cos(qJ(2));
t155 = t134 * pkin(2) + t119;
t175 = -MDP(12) + MDP(17) + MDP(20);
t174 = MDP(13) + MDP(15) + MDP(19);
t173 = -MDP(16) + MDP(11) + MDP(21);
t133 = sin(qJ(1));
t171 = g(1) * t133;
t128 = sin(pkin(8));
t169 = t128 * t132;
t167 = t129 * t133;
t166 = t129 * t134;
t135 = cos(qJ(1));
t165 = t129 * t135;
t131 = cos(pkin(5));
t164 = t131 * t134;
t130 = cos(pkin(8));
t163 = t132 * t130;
t162 = t132 * t133;
t161 = t132 * t135;
t160 = t133 * t131;
t159 = t133 * t134;
t158 = t134 * t135;
t157 = t135 * t131;
t156 = t135 * pkin(7) + qJ(3) * t157;
t154 = MDP(18) + MDP(22);
t153 = pkin(2) * t162;
t152 = pkin(2) * t161;
t151 = qJ(3) * t166;
t150 = t131 * t169;
t149 = g(1) * t135 + g(2) * t133;
t141 = -t130 * t164 + t169;
t102 = t141 * t133;
t140 = t128 * t164 + t163;
t103 = t140 * t133;
t113 = t133 * t151;
t147 = -t103 * pkin(3) - qJ(4) * t102 + t113;
t104 = t141 * t135;
t105 = t140 * t135;
t115 = t135 * t151;
t146 = -t105 * pkin(3) - t104 * qJ(4) + t115;
t145 = pkin(2) * t158 + t133 * pkin(7) + qJ(3) * t160 + (pkin(1) + t119) * t135;
t95 = t128 * t159 + (t132 * t160 + t165) * t130;
t96 = -t128 * t165 + t130 * t159 - t133 * t150;
t144 = -t96 * pkin(3) - qJ(4) * t95 + t156;
t143 = -pkin(2) * t132 + pkin(4) * t166;
t107 = t128 * t134 + t131 * t163;
t108 = t130 * t134 - t150;
t142 = t108 * pkin(3) + qJ(4) * t107 + t155;
t97 = t107 * t135 - t130 * t167;
t98 = t128 * t167 + t130 * t158 - t135 * t150;
t137 = t98 * pkin(3) + qJ(4) * t97 + t145;
t136 = (-pkin(1) - t155) * t171;
t110 = t129 * t161 + t160;
t109 = t129 * t162 - t157;
t1 = [t149 * MDP(3) + (-g(1) * t156 - g(2) * t145 - t136) * MDP(14) + (-g(1) * t144 - g(2) * t137 - t136) * MDP(18) + (-g(1) * (-pkin(4) * t109 - qJ(5) * t96 + t144) - g(2) * (pkin(4) * t110 + qJ(5) * t98 + t137) - t136) * MDP(22) + t174 * (g(1) * t109 - g(2) * t110) + t173 * (g(1) * t96 - g(2) * t98) + t175 * (g(1) * t95 - g(2) * t97) + (-t132 * MDP(10) + t134 * MDP(9) + MDP(2)) * (-g(2) * t135 + t171); (-g(3) * t134 + t149 * t132) * MDP(9) + (-g(1) * (t115 - t152) - g(2) * (t113 - t153) - g(3) * t155) * MDP(14) + (-g(1) * (t146 - t152) - g(2) * (t147 - t153) - g(3) * t142) * MDP(18) + (-g(1) * (-t105 * qJ(5) + t143 * t135 + t146) - g(2) * (-qJ(5) * t103 + t143 * t133 + t147) - g(3) * (pkin(4) * t168 + qJ(5) * t108 + t142)) * MDP(22) + t175 * (g(1) * t104 + g(2) * t102 - g(3) * t107) + t173 * (g(1) * t105 + g(2) * t103 - g(3) * t108) + (-t174 * t129 + MDP(10)) * (g(3) * t132 + t149 * t134); (MDP(14) + t154) * (-g(1) * t110 - g(2) * t109 + g(3) * t166); t154 * (-g(1) * t97 - g(2) * t95 - g(3) * t141); (-g(1) * t98 - g(2) * t96 - g(3) * t140) * MDP(22);];
taug = t1;
