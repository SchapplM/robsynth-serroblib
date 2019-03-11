% Calculate Gravitation load on the joints for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:55:51
% EndTime: 2019-03-10 01:55:54
% DurationCPUTime: 0.98s
% Computational Cost: add. (789->128), mult. (1352->200), div. (0->0), fcn. (1628->12), ass. (0->64)
t184 = sin(qJ(2));
t188 = cos(qJ(2));
t189 = cos(qJ(1));
t185 = sin(qJ(1));
t221 = cos(pkin(6));
t207 = t185 * t221;
t168 = -t184 * t207 + t188 * t189;
t180 = qJ(3) + qJ(4);
t178 = sin(t180);
t179 = cos(t180);
t181 = sin(pkin(6));
t214 = t181 * t185;
t149 = t168 * t178 - t179 * t214;
t206 = t189 * t221;
t166 = t184 * t206 + t185 * t188;
t211 = t181 * t189;
t205 = -t166 * t178 - t179 * t211;
t215 = t181 * t184;
t239 = -g(3) * (-t178 * t215 + t179 * t221) - g(2) * t205 + g(1) * t149;
t235 = MDP(24) - MDP(33);
t234 = MDP(30) + MDP(32);
t233 = MDP(31) - MDP(34);
t146 = t166 * t179 - t178 * t211;
t165 = t184 * t185 - t188 * t206;
t182 = sin(qJ(5));
t186 = cos(qJ(5));
t127 = t146 * t182 - t165 * t186;
t128 = t146 * t186 + t165 * t182;
t150 = t168 * t179 + t178 * t214;
t159 = t178 * t221 + t179 * t215;
t231 = g(1) * t150 + g(2) * t146 + g(3) * t159;
t183 = sin(qJ(3));
t187 = cos(qJ(3));
t213 = t181 * t187;
t151 = -t168 * t183 + t185 * t213;
t200 = t166 * t183 + t187 * t211;
t230 = g(2) * t200 - g(3) * (-t183 * t215 + t187 * t221) - g(1) * t151;
t229 = g(1) * t189 + g(2) * t185;
t167 = t189 * t184 + t188 * t207;
t228 = g(1) * t167 + g(2) * t165;
t217 = t179 * t182;
t216 = t179 * t186;
t212 = t181 * t188;
t210 = t186 * t188;
t209 = t182 * t212;
t204 = t166 * t187 - t183 * t211;
t177 = pkin(3) * t187 + pkin(2);
t201 = pkin(4) * t179 + pkin(11) * t178 + t177;
t131 = t150 * t182 - t167 * t186;
t143 = t159 * t182 + t181 * t210;
t199 = g(1) * t131 + g(2) * t127 + g(3) * t143;
t194 = t235 * t231 + (-t182 * t233 + t186 * t234 + MDP(23)) * t239;
t191 = -t231 * pkin(11) + t239 * (pkin(5) * t186 + qJ(6) * t182 + pkin(4));
t190 = -pkin(10) - pkin(9);
t156 = (t179 * t210 + t182 * t184) * t181;
t155 = t179 * t209 - t186 * t215;
t152 = t168 * t187 + t183 * t214;
t144 = t159 * t186 - t209;
t136 = -t167 * t216 + t168 * t182;
t135 = -t167 * t217 - t168 * t186;
t134 = -t165 * t216 + t166 * t182;
t133 = -t165 * t217 - t166 * t186;
t132 = t150 * t186 + t167 * t182;
t1 = [(g(1) * t185 - g(2) * t189) * MDP(2) + t229 * MDP(3) + (g(1) * t166 - g(2) * t168) * MDP(9) + (-g(1) * t165 + g(2) * t167) * MDP(10) + (g(1) * t204 - g(2) * t152) * MDP(16) + (-g(1) * t200 - g(2) * t151) * MDP(17) + (g(1) * t146 - g(2) * t150) * MDP(23) + (-g(1) * (-t185 * pkin(1) - pkin(4) * t146 - pkin(5) * t128 + pkin(11) * t205 - qJ(6) * t127 + t165 * t190 - t166 * t177) - g(2) * (t189 * pkin(1) + t150 * pkin(4) + t132 * pkin(5) + t149 * pkin(11) + t131 * qJ(6) - t167 * t190 + t168 * t177) - t229 * t181 * (pkin(3) * t183 + pkin(8))) * MDP(35) + t233 * (-g(1) * t127 + g(2) * t131) + t235 * (g(1) * t205 + g(2) * t149) + t234 * (g(1) * t128 - g(2) * t132); (g(1) * t168 + g(2) * t166 + g(3) * t215) * MDP(10) + (-g(1) * (t136 * pkin(5) + t135 * qJ(6) - t168 * t190) - g(2) * (t134 * pkin(5) + t133 * qJ(6) - t166 * t190) + t228 * t201 + (-t156 * pkin(5) - t155 * qJ(6) - (-t184 * t190 + t188 * t201) * t181) * g(3)) * MDP(35) + t233 * (g(1) * t135 + g(2) * t133 + g(3) * t155) + t234 * (-g(1) * t136 - g(2) * t134 - g(3) * t156) + (-MDP(16) * t187 + MDP(17) * t183 - MDP(23) * t179 + t235 * t178 - MDP(9)) * (g(3) * t212 - t228); t230 * MDP(16) + (g(1) * t152 + g(2) * t204 - g(3) * (-t183 * t221 - t184 * t213)) * MDP(17) + (pkin(3) * t230 + t191) * MDP(35) + t194; MDP(35) * t191 + t194; (-g(1) * (-pkin(5) * t131 + qJ(6) * t132) - g(2) * (-pkin(5) * t127 + qJ(6) * t128) - g(3) * (-pkin(5) * t143 + qJ(6) * t144)) * MDP(35) + t234 * t199 + t233 * (g(1) * t132 + g(2) * t128 + g(3) * t144); -t199 * MDP(35);];
taug  = t1;
