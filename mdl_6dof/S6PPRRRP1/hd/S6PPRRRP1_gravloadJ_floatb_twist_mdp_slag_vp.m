% Calculate Gravitation load on the joints for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:51
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:50:15
% EndTime: 2021-01-16 00:50:20
% DurationCPUTime: 1.08s
% Computational Cost: add. (629->105), mult. (1631->181), div. (0->0), fcn. (2097->14), ass. (0->68)
t162 = cos(pkin(12));
t164 = cos(pkin(7));
t165 = cos(pkin(6));
t185 = t164 * t165;
t160 = sin(pkin(7));
t161 = sin(pkin(6));
t187 = t161 * t160;
t147 = t162 * t185 - t187;
t159 = sin(pkin(11));
t163 = cos(pkin(11));
t158 = sin(pkin(12));
t192 = t158 * t164;
t137 = -t147 * t163 + t159 * t192;
t191 = t158 * t165;
t149 = t159 * t162 + t163 * t191;
t169 = sin(qJ(3));
t172 = cos(qJ(3));
t124 = t137 * t169 - t149 * t172;
t188 = t160 * t165;
t145 = t161 * t164 + t162 * t188;
t193 = t158 * t160;
t133 = t145 * t163 - t159 * t193;
t168 = sin(qJ(4));
t171 = cos(qJ(4));
t113 = t124 * t168 - t133 * t171;
t134 = t147 * t159 + t163 * t192;
t150 = -t159 * t191 + t162 * t163;
t118 = t134 * t169 - t150 * t172;
t132 = t145 * t159 + t163 * t193;
t115 = t118 * t168 + t132 * t171;
t186 = t162 * t164;
t148 = t161 * t186 + t188;
t189 = t158 * t172;
t138 = t148 * t169 + t161 * t189;
t146 = t162 * t187 - t185;
t205 = t138 * t168 + t146 * t171;
t211 = g(1) * t115 + g(2) * t113 - g(3) * t205;
t112 = t118 * t171 - t132 * t168;
t117 = t124 * t171 + t133 * t168;
t119 = t134 * t172 + t150 * t169;
t125 = t137 * t172 + t149 * t169;
t178 = g(1) * t119 + g(2) * t125;
t127 = t138 * t171 - t146 * t168;
t204 = -MDP(12) + MDP(22);
t180 = MDP(18) + MDP(20);
t179 = MDP(19) + MDP(21);
t167 = sin(qJ(5));
t203 = t118 * t167;
t202 = t124 * t167;
t131 = t169 * t188 + (t169 * t186 + t189) * t161;
t201 = t131 * t167;
t190 = t158 * t169;
t166 = -qJ(6) - pkin(10);
t184 = t166 * t168;
t183 = t167 * t171;
t170 = cos(qJ(5));
t182 = t170 * t171;
t181 = MDP(2) + MDP(23);
t154 = pkin(3) * t186 + pkin(9) * t158;
t177 = pkin(3) * t187 - t154 * t165;
t153 = -t158 * pkin(3) + pkin(9) * t186;
t176 = pkin(9) * t187 - t153 * t165;
t141 = -t148 * t172 + t161 * t190;
t156 = pkin(5) * t170 + pkin(4);
t152 = pkin(3) * t192 - pkin(9) * t162;
t151 = pkin(3) * t162 + pkin(9) * t192;
t130 = -t172 * t188 + (-t172 * t186 + t190) * t161;
t1 = [(-MDP(1) - t181) * g(3); t181 * (-g(3) * t165 + (-g(1) * t159 + g(2) * t163) * t161); (g(3) * t130 + t178) * MDP(4) + (-g(1) * t118 - g(2) * t124 + g(3) * t131) * MDP(5) + (-g(1) * (t119 * t184 - pkin(5) * t203 - (t163 * t151 - t176 * t159) * t169 + (-t163 * t152 + t177 * t159) * t172) - g(2) * (t125 * t184 - pkin(5) * t202 - (t159 * t151 + t176 * t163) * t169 + (-t159 * t152 - t177 * t163) * t172) - g(3) * (t141 * t184 + pkin(5) * t201 - (-pkin(9) * t188 - t153 * t161) * t169 + (pkin(3) * t188 + t154 * t161) * t172)) * MDP(23) + t180 * (-g(1) * (-t119 * t182 - t203) - g(2) * (-t125 * t182 - t202) - g(3) * (-t141 * t182 + t201)) + t179 * (-g(1) * (-t118 * t170 + t119 * t183) - g(2) * (-t124 * t170 + t125 * t183) - g(3) * (t131 * t170 + t141 * t183)) + (t204 * t168 + (MDP(23) * t156 + MDP(11)) * t171) * (g(3) * t141 + t178); (-g(1) * (t112 * t166 + t115 * t156) - g(2) * (t113 * t156 + t117 * t166) - g(3) * (-t127 * t166 - t156 * t205)) * MDP(23) + t204 * (g(1) * t112 + g(2) * t117 - g(3) * t127) + (t179 * t167 - t180 * t170 - MDP(11)) * t211; t179 * (-g(1) * (t112 * t170 - t119 * t167) - g(2) * (t117 * t170 - t125 * t167) - g(3) * (-t127 * t170 - t130 * t167)) + (MDP(23) * pkin(5) + t180) * (-g(1) * (t112 * t167 + t119 * t170) - g(2) * (t117 * t167 + t125 * t170) - g(3) * (-t127 * t167 + t130 * t170)); t211 * MDP(23);];
taug = t1;
