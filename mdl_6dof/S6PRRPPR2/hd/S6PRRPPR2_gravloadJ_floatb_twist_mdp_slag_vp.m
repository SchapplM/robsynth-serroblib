% Calculate Gravitation load on the joints for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:19
% EndTime: 2021-01-16 02:21:24
% DurationCPUTime: 0.87s
% Computational Cost: add. (339->89), mult. (622->141), div. (0->0), fcn. (707->12), ass. (0->53)
t119 = sin(pkin(10));
t124 = sin(qJ(2));
t160 = cos(pkin(10));
t161 = cos(pkin(6));
t140 = t161 * t160;
t166 = cos(qJ(2));
t104 = t119 * t166 + t124 * t140;
t118 = qJ(3) + pkin(11);
t116 = sin(t118);
t117 = cos(t118);
t120 = sin(pkin(6));
t145 = t120 * t160;
t131 = -t104 * t116 - t117 * t145;
t146 = t119 * t161;
t102 = t124 * t146 - t160 * t166;
t155 = t119 * t120;
t137 = t102 * t116 + t117 * t155;
t154 = t120 * t124;
t99 = t116 * t154 - t117 * t161;
t179 = g(1) * t137 + g(2) * t131 - g(3) * t99;
t123 = sin(qJ(3));
t126 = cos(qJ(3));
t178 = -t123 * t154 + t161 * t126;
t153 = t120 * t126;
t177 = t102 * t123 + t119 * t153;
t172 = MDP(12) - MDP(17);
t171 = -MDP(13) + MDP(18);
t165 = pkin(4) * t117;
t163 = g(3) * t120;
t105 = t124 * t160 + t146 * t166;
t115 = pkin(3) * t126 + pkin(2);
t121 = qJ(4) + pkin(8);
t162 = -t102 * t121 - t105 * t115;
t122 = sin(qJ(6));
t157 = t116 * t122;
t125 = cos(qJ(6));
t156 = t116 * t125;
t147 = t120 * t166;
t152 = t115 * t147 + t121 * t154;
t151 = MDP(15) + MDP(19);
t148 = t116 * t166;
t103 = t119 * t124 - t140 * t166;
t143 = -t103 * t115 + t104 * t121;
t142 = t177 * pkin(3);
t139 = -qJ(5) * t116 - t165;
t138 = t178 * pkin(3);
t136 = -t102 * t117 + t116 * t155;
t132 = -t104 * t123 - t126 * t145;
t95 = t104 * t117 - t116 * t145;
t128 = t132 * pkin(3);
t127 = g(1) * t105 + g(2) * t103 - g(3) * t147;
t100 = t116 * t161 + t117 * t154;
t1 = [(-MDP(1) - t151) * g(3); (-g(1) * t162 - g(2) * t143 - g(3) * t152) * MDP(15) + (-g(1) * (t105 * t139 + t162) - g(2) * (t103 * t139 + t143) - g(3) * ((qJ(5) * t148 + t165 * t166) * t120 + t152)) * MDP(19) + (-g(1) * (-t102 * t125 - t105 * t157) - g(2) * (-t103 * t157 + t104 * t125) - (t122 * t148 + t124 * t125) * t163) * MDP(25) + (-g(1) * (t102 * t122 - t105 * t156) - g(2) * (-t103 * t156 - t104 * t122) - (-t122 * t124 + t125 * t148) * t163) * MDP(26) + (MDP(4) - MDP(14) - MDP(16)) * (-g(1) * t102 + g(2) * t104 + g(3) * t154) + (t126 * MDP(10) - t123 * MDP(11) + t171 * t116 + t172 * t117 + MDP(3)) * t127; (-g(1) * t177 - g(2) * t132 - g(3) * t178) * MDP(10) + (-g(1) * (t102 * t126 - t123 * t155) - g(2) * (-t104 * t126 + t123 * t145) - g(3) * (-t123 * t161 - t124 * t153)) * MDP(11) + (-g(1) * t142 - g(2) * t128 - g(3) * t138) * MDP(15) + (-g(1) * (pkin(4) * t137 + qJ(5) * t136 + t142) - g(2) * (pkin(4) * t131 + t95 * qJ(5) + t128) - g(3) * (-pkin(4) * t99 + qJ(5) * t100 + t138)) * MDP(19) - t172 * t179 + (t122 * MDP(25) + t125 * MDP(26) + t171) * (-g(1) * t136 - g(2) * t95 - g(3) * t100); -t151 * t127; t179 * MDP(19); (t127 * t122 + t125 * t179) * MDP(25) + (-t122 * t179 + t127 * t125) * MDP(26);];
taug = t1;
