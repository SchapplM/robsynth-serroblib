% Calculate Gravitation load on the joints for
% S6PPRRRP2
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
%   see S6PPRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:09
% EndTime: 2019-03-08 18:58:11
% DurationCPUTime: 0.53s
% Computational Cost: add. (745->89), mult. (2080->147), div. (0->0), fcn. (2685->14), ass. (0->61)
t195 = -MDP(12) + MDP(21);
t191 = MDP(18) + MDP(20);
t190 = MDP(19) - MDP(22);
t175 = sin(pkin(12));
t176 = sin(pkin(11));
t163 = t176 * t175;
t179 = cos(pkin(12));
t180 = cos(pkin(11));
t171 = t180 * t179;
t182 = cos(pkin(6));
t150 = -t182 * t171 + t163;
t177 = sin(pkin(7));
t178 = sin(pkin(6));
t168 = t178 * t177;
t181 = cos(pkin(7));
t194 = t150 * t181 + t180 * t168;
t164 = t176 * t179;
t169 = t180 * t175;
t151 = t182 * t164 + t169;
t167 = t178 * t176;
t193 = t151 * t181 - t177 * t167;
t192 = t179 * t181 * t178 + t177 * t182;
t189 = cos(qJ(3));
t143 = sin(qJ(5));
t147 = cos(qJ(4));
t174 = t143 * t147;
t146 = cos(qJ(5));
t173 = t146 * t147;
t172 = MDP(2) + MDP(23);
t170 = t180 * t178;
t166 = t178 * t175;
t137 = t182 * t169 + t164;
t145 = sin(qJ(3));
t122 = t137 * t189 - t145 * t194;
t131 = t150 * t177 - t181 * t170;
t144 = sin(qJ(4));
t110 = t122 * t147 + t131 * t144;
t121 = t137 * t145 + t189 * t194;
t100 = t110 * t143 - t121 * t146;
t138 = -t182 * t163 + t171;
t124 = t138 * t189 - t145 * t193;
t132 = t151 * t177 + t181 * t167;
t112 = t124 * t147 + t132 * t144;
t123 = t138 * t145 + t189 * t193;
t102 = t112 * t143 - t123 * t146;
t130 = t145 * t192 + t189 * t166;
t136 = -t179 * t168 + t182 * t181;
t126 = t130 * t147 + t136 * t144;
t129 = t145 * t166 - t192 * t189;
t113 = t126 * t143 - t129 * t146;
t160 = g(1) * t102 + g(2) * t100 + g(3) * t113;
t116 = -t129 * t173 + t130 * t143;
t115 = -t129 * t174 - t130 * t146;
t114 = t126 * t146 + t129 * t143;
t108 = -t123 * t173 + t124 * t143;
t107 = -t123 * t174 - t124 * t146;
t106 = -t121 * t173 + t122 * t143;
t105 = -t121 * t174 - t122 * t146;
t103 = t112 * t146 + t123 * t143;
t101 = t110 * t146 + t121 * t143;
t1 = [(-MDP(1) - t172) * g(3); t172 * (-g(1) * t167 + g(2) * t170 - g(3) * t182); (g(1) * t124 + g(2) * t122 + g(3) * t130) * MDP(5) + (-g(1) * (t108 * pkin(5) + t124 * pkin(9) + t107 * qJ(6)) - g(2) * (t106 * pkin(5) + t122 * pkin(9) + t105 * qJ(6)) - g(3) * (t116 * pkin(5) + t130 * pkin(9) + t115 * qJ(6))) * MDP(23) + t191 * (-g(1) * t108 - g(2) * t106 - g(3) * t116) + t190 * (g(1) * t107 + g(2) * t105 + g(3) * t115) + (MDP(4) + t147 * MDP(11) + (pkin(4) * t147 + pkin(10) * t144 + pkin(3)) * MDP(23) + t195 * t144) * (g(1) * t123 + g(2) * t121 + g(3) * t129); (-MDP(23) * pkin(10) - t195) * (g(1) * t112 + g(2) * t110 + g(3) * t126) + (-MDP(11) - t191 * t146 + t190 * t143 - MDP(23) * (pkin(5) * t146 + qJ(6) * t143 + pkin(4))) * (g(3) * (-t130 * t144 + t136 * t147) + g(2) * (-t122 * t144 + t131 * t147) + g(1) * (-t124 * t144 + t132 * t147)); (-g(1) * (-pkin(5) * t102 + qJ(6) * t103) - g(2) * (-pkin(5) * t100 + qJ(6) * t101) - g(3) * (-pkin(5) * t113 + qJ(6) * t114)) * MDP(23) + t191 * t160 + t190 * (g(1) * t103 + g(2) * t101 + g(3) * t114); -t160 * MDP(23);];
taug  = t1;
