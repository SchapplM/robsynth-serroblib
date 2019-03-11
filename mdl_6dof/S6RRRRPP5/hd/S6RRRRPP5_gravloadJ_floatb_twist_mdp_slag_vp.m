% Calculate Gravitation load on the joints for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:09:09
% EndTime: 2019-03-09 21:09:11
% DurationCPUTime: 0.72s
% Computational Cost: add. (506->111), mult. (681->147), div. (0->0), fcn. (686->8), ass. (0->59)
t191 = MDP(24) - MDP(27) - MDP(30);
t140 = sin(qJ(3));
t142 = sin(qJ(1));
t143 = cos(qJ(3));
t144 = cos(qJ(2));
t145 = cos(qJ(1));
t168 = t144 * t145;
t118 = -t140 * t168 + t142 * t143;
t187 = MDP(23) + MDP(25) + MDP(29);
t180 = g(2) * t142;
t157 = g(1) * t145 + t180;
t188 = MDP(10) - MDP(26) + MDP(31);
t141 = sin(qJ(2));
t114 = -g(3) * t144 + t157 * t141;
t185 = -pkin(4) - pkin(5);
t184 = pkin(3) * t140;
t139 = qJ(3) + qJ(4);
t135 = cos(t139);
t183 = pkin(4) * t135;
t182 = g(1) * t142;
t178 = g(3) * t141;
t134 = sin(t139);
t176 = qJ(5) * t134;
t175 = t134 * t141;
t174 = t135 * t141;
t173 = t140 * t142;
t172 = t140 * t145;
t146 = -pkin(9) - pkin(8);
t171 = t141 * t146;
t169 = t142 * t144;
t133 = pkin(3) * t143 + pkin(2);
t124 = t144 * t133;
t167 = t145 * t134;
t166 = qJ(6) + t146;
t165 = t185 * t134;
t163 = -pkin(1) - t124;
t162 = t166 * t145;
t110 = t134 * t169 + t135 * t145;
t111 = t135 * t169 - t167;
t161 = -t110 * pkin(4) + qJ(5) * t111;
t112 = -t142 * t135 + t144 * t167;
t113 = t134 * t142 + t135 * t168;
t160 = -t112 * pkin(4) + qJ(5) * t113;
t159 = -t133 - t176;
t158 = g(3) * (t124 + (t176 + t183) * t144);
t98 = g(1) * t112 + g(2) * t110 + g(3) * t175;
t155 = t187 * t98 + t191 * (g(1) * t113 + g(2) * t111 + g(3) * t174);
t153 = t157 * t144;
t116 = t140 * t169 + t143 * t145;
t152 = pkin(3) * t172 - t111 * pkin(4) + t145 * pkin(7) - qJ(5) * t110 + t142 * t171;
t150 = t145 * pkin(1) + pkin(3) * t173 + t113 * pkin(4) + t142 * pkin(7) + t112 * qJ(5) + t133 * t168;
t148 = t118 * pkin(3) + t160;
t147 = -t116 * pkin(3) + t161;
t122 = qJ(5) * t174;
t119 = t143 * t168 + t173;
t117 = -t143 * t169 + t172;
t107 = t112 * pkin(5);
t104 = t110 * pkin(5);
t1 = [t157 * MDP(3) + (-g(1) * t117 - g(2) * t119) * MDP(16) + (-g(1) * t116 - g(2) * t118) * MDP(17) + (-g(1) * (t163 * t142 + t152) - g(2) * (-t145 * t171 + t150)) * MDP(28) + (-g(1) * (-pkin(5) * t111 + t152) - g(2) * (t113 * pkin(5) - t141 * t162 + t150) - (t141 * qJ(6) + t163) * t182) * MDP(32) + t187 * (g(1) * t111 - g(2) * t113) - t191 * (g(1) * t110 - g(2) * t112) + (t144 * MDP(9) - t188 * t141 + MDP(2)) * (-g(2) * t145 + t182); (-t158 + t146 * t153 + (g(3) * t146 + t157 * (-t159 + t183)) * t141) * MDP(28) + (-t158 + (-g(3) * pkin(5) * t135 + g(1) * t162 + t166 * t180) * t144 + (g(3) * t166 + t157 * (-t185 * t135 - t159)) * t141) * MDP(32) + t188 * (t153 + t178) + (t143 * MDP(16) - t140 * MDP(17) - t134 * t191 + t187 * t135 + MDP(9)) * t114; (-g(1) * t118 + g(2) * t116 + t140 * t178) * MDP(16) + (g(1) * t119 - g(2) * t117 + t143 * t178) * MDP(17) + (-g(1) * t148 - g(2) * t147 - g(3) * (t122 + (-pkin(4) * t134 - t184) * t141)) * MDP(28) + (-g(1) * (-t107 + t148) - g(2) * (-t104 + t147) - g(3) * t122 - (t165 - t184) * t178) * MDP(32) + t155; (-g(1) * t160 - g(2) * t161 - g(3) * (-pkin(4) * t175 + t122)) * MDP(28) + (-g(1) * (-t107 + t160) - g(2) * (-t104 + t161) - g(3) * (t141 * t165 + t122)) * MDP(32) + t155; -(MDP(28) + MDP(32)) * t98; t114 * MDP(32);];
taug  = t1;
