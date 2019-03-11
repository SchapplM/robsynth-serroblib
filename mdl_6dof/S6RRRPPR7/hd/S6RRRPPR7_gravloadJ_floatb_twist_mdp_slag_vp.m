% Calculate Gravitation load on the joints for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:01:07
% EndTime: 2019-03-09 16:01:10
% DurationCPUTime: 1.02s
% Computational Cost: add. (332->119), mult. (699->169), div. (0->0), fcn. (756->10), ass. (0->60)
t146 = sin(qJ(1));
t147 = cos(qJ(3));
t148 = cos(qJ(2));
t144 = sin(qJ(3));
t149 = cos(qJ(1));
t174 = t149 * t144;
t121 = -t146 * t147 + t148 * t174;
t175 = t148 * t149;
t122 = t146 * t144 + t147 * t175;
t141 = pkin(10) + qJ(6);
t134 = sin(t141);
t135 = cos(t141);
t101 = t121 * t135 - t122 * t134;
t102 = t121 * t134 + t122 * t135;
t157 = t134 * t144 + t135 * t147;
t158 = t134 * t147 - t135 * t144;
t177 = t146 * t148;
t119 = t144 * t177 + t147 * t149;
t176 = t147 * t148;
t120 = t146 * t176 - t174;
t161 = t119 * t134 + t120 * t135;
t145 = sin(qJ(2));
t184 = g(3) * t145;
t199 = t119 * t135 - t120 * t134;
t204 = -(g(1) * t101 + g(2) * t199 - t158 * t184) * MDP(31) + (g(1) * t102 + g(2) * t161 + t157 * t184) * MDP(32);
t164 = g(1) * t149 + g(2) * t146;
t195 = t145 * t164;
t202 = g(3) * t148 - t195;
t197 = MDP(16) + MDP(18);
t196 = MDP(17) - MDP(20);
t200 = MDP(10) - MDP(19) + MDP(24);
t192 = -pkin(3) - pkin(4);
t189 = g(1) * t146;
t136 = t145 * pkin(8);
t138 = t148 * pkin(2);
t182 = qJ(4) * t144;
t180 = t144 * t145;
t179 = t145 * t147;
t178 = t145 * t149;
t173 = -pkin(1) - t138;
t172 = -pkin(2) - t182;
t171 = -t119 * pkin(3) + qJ(4) * t120;
t170 = -t121 * pkin(3) + qJ(4) * t122;
t142 = sin(pkin(10));
t143 = cos(pkin(10));
t168 = -t119 * t143 + t120 * t142;
t167 = t121 * t142 + t122 * t143;
t166 = pkin(3) * t176 + t148 * t182 + t136 + t138;
t162 = -t120 * pkin(3) + t149 * pkin(7) - qJ(4) * t119;
t160 = t119 * t142 + t120 * t143;
t159 = t121 * t143 - t122 * t142;
t156 = t142 * t147 - t143 * t144;
t155 = t142 * t144 + t143 * t147;
t153 = g(3) * t155;
t151 = t149 * pkin(1) + pkin(2) * t175 + t122 * pkin(3) + t146 * pkin(7) + pkin(8) * t178 + t121 * qJ(4);
t104 = g(1) * t121 + g(2) * t119 + g(3) * t180;
t129 = pkin(8) * t175;
t126 = pkin(8) * t177;
t124 = qJ(4) * t179;
t1 = [t164 * MDP(3) + (-g(1) * t162 - g(2) * t151 - (t173 - t136) * t189) * MDP(21) + (g(1) * t160 - g(2) * t167) * MDP(22) + (-g(1) * t168 - g(2) * t159) * MDP(23) + (-g(1) * (-pkin(4) * t120 + t162) - g(2) * (t122 * pkin(4) - qJ(5) * t178 + t151) - ((-pkin(8) + qJ(5)) * t145 + t173) * t189) * MDP(25) + (g(1) * t161 - g(2) * t102) * MDP(31) + (g(1) * t199 - g(2) * t101) * MDP(32) - t196 * (g(1) * t119 - g(2) * t121) + (t148 * MDP(9) + MDP(2)) * (-g(2) * t149 + t189) + t197 * (g(1) * t120 - g(2) * t122) - t200 * (-g(2) * t178 + t145 * t189); (-g(1) * t129 - g(2) * t126 - g(3) * t166 + (pkin(3) * t147 - t172) * t195) * MDP(21) + (-t148 * t153 + t155 * t195) * MDP(22) + (-g(1) * (-qJ(5) * t175 + t129) - g(2) * (-qJ(5) * t177 + t126) - g(3) * (pkin(4) * t176 + t166) + (g(3) * qJ(5) + t164 * (-t147 * t192 - t172)) * t145) * MDP(25) + t200 * (t148 * t164 + t184) + (t156 * MDP(23) - t157 * MDP(31) + t158 * MDP(32) + t196 * t144 - t197 * t147 - MDP(9)) * t202; (-g(1) * t170 - g(2) * t171 - g(3) * (-pkin(3) * t180 + t124)) * MDP(21) + (g(1) * t159 - g(2) * t168 - t156 * t184) * MDP(22) + (-g(1) * t167 - g(2) * t160 - t145 * t153) * MDP(23) + (-g(1) * (-pkin(4) * t121 + t170) - g(2) * (-pkin(4) * t119 + t171) - g(3) * (t180 * t192 + t124)) * MDP(25) + t196 * (g(1) * t122 + g(2) * t120 + g(3) * t179) + t197 * t104 - t204; -(MDP(21) + MDP(25)) * t104; -t202 * MDP(25); t204;];
taug  = t1;
