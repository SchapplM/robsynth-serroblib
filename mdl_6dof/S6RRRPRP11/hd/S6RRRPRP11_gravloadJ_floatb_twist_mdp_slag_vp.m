% Calculate Gravitation load on the joints for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:49:04
% EndTime: 2019-03-09 17:49:07
% DurationCPUTime: 0.88s
% Computational Cost: add. (419->124), mult. (1020->188), div. (0->0), fcn. (1195->10), ass. (0->58)
t193 = -MDP(16) + MDP(19) - MDP(29);
t191 = MDP(17) - MDP(20);
t149 = cos(qJ(3));
t199 = pkin(3) * t149 + pkin(2);
t198 = MDP(30) * pkin(5) + MDP(27);
t145 = sin(qJ(3));
t197 = -qJ(4) * t145 - t199;
t146 = sin(qJ(2));
t147 = sin(qJ(1));
t150 = cos(qJ(2));
t184 = cos(pkin(6));
t189 = cos(qJ(1));
t162 = t184 * t189;
t125 = t146 * t162 + t147 * t150;
t142 = sin(pkin(6));
t169 = t142 * t189;
t107 = t125 * t145 + t149 * t169;
t124 = t146 * t147 - t150 * t162;
t144 = sin(qJ(5));
t148 = cos(qJ(5));
t196 = t107 * t144 + t124 * t148;
t195 = t107 * t148 - t124 * t144;
t192 = MDP(10) - MDP(18);
t187 = g(3) * t142;
t140 = pkin(5) * t148 + pkin(4);
t186 = pkin(9) + t140;
t180 = t142 * t146;
t179 = t142 * t147;
t178 = t142 * t149;
t177 = t142 * t150;
t176 = t144 * t145;
t175 = t145 * t148;
t174 = t145 * t150;
t173 = t148 * t150;
t171 = t197 * t124;
t167 = t147 * t184;
t126 = t146 * t189 + t150 * t167;
t170 = t197 * t126;
t168 = pkin(5) * t144 + qJ(4);
t108 = t125 * t149 - t145 * t169;
t166 = -t107 * pkin(3) + qJ(4) * t108;
t127 = -t146 * t167 + t150 * t189;
t111 = t127 * t145 - t147 * t178;
t112 = t127 * t149 + t145 * t179;
t165 = -t111 * pkin(3) + qJ(4) * t112;
t122 = t145 * t180 - t149 * t184;
t123 = t145 * t184 + t146 * t178;
t164 = -t122 * pkin(3) + qJ(4) * t123;
t163 = t189 * pkin(1) + t127 * pkin(2) + t112 * pkin(3) + pkin(8) * t179;
t99 = t111 * t148 - t126 * t144;
t158 = g(3) * (t142 * qJ(4) * t174 + pkin(9) * t180 + t199 * t177);
t157 = -pkin(1) * t147 - t125 * pkin(2) - pkin(3) * t108 + pkin(8) * t169;
t143 = -qJ(6) - pkin(10);
t156 = pkin(5) * t176 - t143 * t149;
t155 = g(1) * t111 + g(2) * t107 + g(3) * t122;
t154 = g(1) * t112 + g(2) * t108 + g(3) * t123;
t100 = t111 * t144 + t126 * t148;
t1 = [(g(1) * t147 - g(2) * t189) * MDP(2) + (g(1) * t189 + g(2) * t147) * MDP(3) + (g(1) * t125 - g(2) * t127) * MDP(9) + (-g(1) * (-pkin(9) * t124 - qJ(4) * t107 + t157) - g(2) * (pkin(9) * t126 + qJ(4) * t111 + t163)) * MDP(21) + (g(1) * t196 - g(2) * t100) * MDP(27) + (g(1) * t195 - g(2) * t99) * MDP(28) + (-g(1) * (-t107 * t168 + t108 * t143 - t124 * t186 + t157) - g(2) * (t111 * t168 - t112 * t143 + t126 * t186 + t163)) * MDP(30) + t191 * (-g(1) * t107 + g(2) * t111) + t193 * (-g(1) * t108 + g(2) * t112) - t192 * (g(1) * t124 - g(2) * t126); (-g(1) * (t127 * pkin(9) + t170) - g(2) * (t125 * pkin(9) + t171) - t158) * MDP(21) + (-g(1) * (-t126 * t176 + t127 * t148) - g(2) * (-t124 * t176 + t125 * t148) - (t144 * t174 + t146 * t148) * t187) * MDP(27) + (-g(1) * (-t126 * t175 - t127 * t144) - g(2) * (-t124 * t175 - t125 * t144) - (-t144 * t146 + t145 * t173) * t187) * MDP(28) + (-g(1) * (-t126 * t156 + t127 * t186 + t170) - g(2) * (-t124 * t156 + t125 * t186 + t171) - t158 - (t140 * t146 + t150 * t156) * t187) * MDP(30) + t192 * (g(1) * t127 + g(2) * t125 + g(3) * t180) + (t145 * t191 + t149 * t193 - MDP(9)) * (-g(1) * t126 - g(2) * t124 + g(3) * t177); (-g(1) * t165 - g(2) * t166 - g(3) * t164) * MDP(21) + (-g(1) * (t111 * t143 + t165) - g(2) * (t107 * t143 + t166) - g(3) * (t122 * t143 + t164)) * MDP(30) - t193 * t155 + (-t148 * MDP(28) - t198 * t144 + t191) * t154; -(MDP(21) + MDP(30)) * t155; (g(1) * t100 + g(2) * t196 - g(3) * (-t122 * t144 + t142 * t173)) * MDP(28) + t198 * (-g(2) * t195 - g(3) * (t122 * t148 + t144 * t177) - g(1) * t99); -t154 * MDP(30);];
taug  = t1;
