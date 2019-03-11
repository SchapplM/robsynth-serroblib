% Calculate Gravitation load on the joints for
% S6RRRPRP12
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
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:59:19
% EndTime: 2019-03-09 17:59:22
% DurationCPUTime: 1.02s
% Computational Cost: add. (543->140), mult. (1357->209), div. (0->0), fcn. (1631->10), ass. (0->64)
t236 = -MDP(16) + MDP(19) - MDP(30);
t186 = sin(qJ(3));
t238 = -qJ(4) * t186 - pkin(2);
t234 = MDP(17) - MDP(20);
t212 = MDP(27) + MDP(29);
t233 = MDP(28) - MDP(31);
t187 = sin(qJ(2));
t188 = sin(qJ(1));
t191 = cos(qJ(2));
t227 = cos(pkin(6));
t231 = cos(qJ(1));
t205 = t227 * t231;
t167 = t187 * t205 + t188 * t191;
t190 = cos(qJ(3));
t184 = sin(pkin(6));
t208 = t184 * t231;
t145 = t167 * t186 + t190 * t208;
t166 = t187 * t188 - t191 * t205;
t185 = sin(qJ(5));
t189 = cos(qJ(5));
t125 = t145 * t185 + t166 * t189;
t126 = t145 * t189 - t166 * t185;
t235 = MDP(10) - MDP(18);
t232 = pkin(4) + pkin(9);
t223 = t166 * t190;
t207 = t188 * t227;
t168 = t187 * t231 + t191 * t207;
t222 = t168 * t190;
t221 = t184 * t187;
t220 = t184 * t188;
t219 = t184 * t190;
t218 = t184 * t191;
t217 = t185 * t186;
t216 = t186 * t189;
t215 = t186 * t191;
t214 = t190 * t191;
t211 = t189 * t218;
t210 = -pkin(3) * t223 + t238 * t166;
t209 = -pkin(3) * t222 + t238 * t168;
t146 = t167 * t190 - t186 * t208;
t206 = pkin(2) * t218 + pkin(9) * t221 + (pkin(3) * t214 + qJ(4) * t215) * t184;
t169 = -t187 * t207 + t191 * t231;
t149 = t169 * t186 - t188 * t219;
t150 = t169 * t190 + t186 * t220;
t199 = t231 * pkin(1) + t169 * pkin(2) + t150 * pkin(3) + pkin(8) * t220 + qJ(4) * t149;
t128 = -t149 * t189 + t168 * t185;
t164 = t186 * t221 - t190 * t227;
t143 = t164 * t189 + t185 * t218;
t198 = g(1) * t128 - g(2) * t126 - g(3) * t143;
t122 = g(1) * t149 + g(2) * t145 + g(3) * t164;
t192 = -pkin(1) * t188 - t167 * pkin(2) - pkin(3) * t146 + pkin(8) * t208 - qJ(4) * t145;
t165 = t186 * t227 + t187 * t219;
t159 = t164 * pkin(3);
t154 = (t185 * t215 + t187 * t189) * t184;
t153 = t185 * t221 - t186 * t211;
t144 = -t164 * t185 + t211;
t141 = t149 * pkin(3);
t139 = t145 * pkin(3);
t135 = -t168 * t217 + t169 * t189;
t134 = t168 * t216 + t169 * t185;
t133 = -t166 * t217 + t167 * t189;
t132 = t166 * t216 + t167 * t185;
t129 = t149 * t185 + t168 * t189;
t1 = [(g(1) * t188 - g(2) * t231) * MDP(2) + (g(1) * t231 + g(2) * t188) * MDP(3) + (g(1) * t167 - g(2) * t169) * MDP(9) + (-g(1) * (-pkin(9) * t166 + t192) - g(2) * (pkin(9) * t168 + t199)) * MDP(21) + (-g(1) * (-pkin(5) * t125 - pkin(10) * t146 + qJ(6) * t126 - t166 * t232 + t192) - g(2) * (pkin(5) * t129 + pkin(10) * t150 + qJ(6) * t128 + t168 * t232 + t199)) * MDP(32) + t212 * (g(1) * t125 - g(2) * t129) + t233 * (g(1) * t126 + g(2) * t128) + t234 * (-g(1) * t145 + g(2) * t149) + t236 * (-g(1) * t146 + g(2) * t150) - t235 * (g(1) * t166 - g(2) * t168); (-g(1) * (pkin(9) * t169 + t209) - g(2) * (pkin(9) * t167 + t210) - g(3) * t206) * MDP(21) + (-g(1) * (pkin(5) * t135 - pkin(10) * t222 + qJ(6) * t134 + t169 * t232 + t209) - g(2) * (pkin(5) * t133 - pkin(10) * t223 + qJ(6) * t132 + t167 * t232 + t210) - g(3) * (t154 * pkin(5) + t153 * qJ(6) + (pkin(4) * t187 + pkin(10) * t214) * t184 + t206)) * MDP(32) + t212 * (-g(1) * t135 - g(2) * t133 - g(3) * t154) + t233 * (g(1) * t134 + g(2) * t132 + g(3) * t153) + t235 * (g(1) * t169 + g(2) * t167 + g(3) * t221) + (t234 * t186 + t236 * t190 - MDP(9)) * (-g(1) * t168 - g(2) * t166 + g(3) * t218); (-g(1) * (qJ(4) * t150 - t141) - g(2) * (qJ(4) * t146 - t139) - g(3) * (qJ(4) * t165 - t159)) * MDP(21) + (-g(1) * (-pkin(10) * t149 - t141) - g(2) * (-pkin(10) * t145 - t139) - g(3) * (-pkin(10) * t164 - t159)) * MDP(32) - t236 * t122 + ((-pkin(5) * t185 + qJ(6) * t189 - qJ(4)) * MDP(32) + t234 - t233 * t189 - t212 * t185) * (g(1) * t150 + g(2) * t146 + g(3) * t165); -(MDP(21) + MDP(32)) * t122; (-g(1) * (-pkin(5) * t128 + qJ(6) * t129) - g(2) * (pkin(5) * t126 + qJ(6) * t125) - g(3) * (pkin(5) * t143 - qJ(6) * t144)) * MDP(32) + t212 * t198 + t233 * (g(1) * t129 + g(2) * t125 - g(3) * t144); -t198 * MDP(32);];
taug  = t1;
