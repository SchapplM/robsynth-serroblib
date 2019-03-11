% Calculate Gravitation load on the joints for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:11:09
% EndTime: 2019-03-09 17:11:12
% DurationCPUTime: 1.14s
% Computational Cost: add. (632->145), mult. (1156->224), div. (0->0), fcn. (1375->12), ass. (0->74)
t237 = MDP(10) - MDP(18);
t216 = MDP(25) + MDP(27);
t236 = MDP(26) - MDP(29);
t187 = sin(qJ(2));
t188 = sin(qJ(1));
t191 = cos(qJ(2));
t192 = cos(qJ(1));
t232 = cos(pkin(6));
t209 = t192 * t232;
t161 = t187 * t209 + t188 * t191;
t182 = qJ(3) + pkin(11);
t179 = sin(t182);
t180 = cos(t182);
t183 = sin(pkin(6));
t221 = t183 * t192;
t136 = t161 * t180 - t179 * t221;
t160 = t187 * t188 - t191 * t209;
t185 = sin(qJ(5));
t189 = cos(qJ(5));
t122 = t136 * t185 - t160 * t189;
t123 = t136 * t189 + t160 * t185;
t186 = sin(qJ(3));
t190 = cos(qJ(3));
t212 = t190 * t221;
t200 = t161 * t186 + t212;
t210 = t188 * t232;
t163 = -t187 * t210 + t191 * t192;
t228 = t163 * t186;
t227 = t180 * t185;
t226 = t180 * t189;
t225 = t183 * t187;
t224 = t183 * t188;
t223 = t183 * t190;
t222 = t183 * t191;
t184 = -qJ(4) - pkin(9);
t220 = t184 * t187;
t219 = t189 * t191;
t178 = pkin(3) * t190 + pkin(2);
t218 = -t160 * t178 - t161 * t184;
t162 = t192 * t187 + t191 * t210;
t217 = -t162 * t178 - t163 * t184;
t215 = t186 * t224;
t214 = t186 * t225;
t213 = t188 * t223;
t211 = t185 * t222;
t173 = t186 * t221;
t208 = t232 * t190;
t207 = t161 * t190 - t173;
t206 = pkin(4) * t180 + pkin(10) * t179;
t203 = t192 * pkin(1) + pkin(3) * t215 + pkin(8) * t224 - t162 * t184 + t163 * t178;
t201 = -t161 * t179 - t180 * t221;
t199 = -pkin(1) * t188 + pkin(3) * t173 + pkin(8) * t221 + t160 * t184 - t161 * t178;
t140 = t163 * t180 + t179 * t224;
t126 = t140 * t185 - t162 * t189;
t150 = t232 * t179 + t180 * t225;
t133 = t150 * t185 + t183 * t219;
t198 = g(1) * t126 + g(2) * t122 + g(3) * t133;
t194 = -g(1) * t162 - g(2) * t160 + g(3) * t222;
t193 = g(1) * t163 + g(2) * t161 + g(3) * t225;
t177 = pkin(3) * t208;
t168 = pkin(3) * t213;
t164 = t178 * t222;
t144 = (t180 * t219 + t185 * t187) * t183;
t143 = t180 * t211 - t189 * t225;
t142 = t163 * t190 + t215;
t141 = t213 - t228;
t139 = t163 * t179 - t180 * t224;
t134 = t150 * t189 - t211;
t131 = -t162 * t226 + t163 * t185;
t130 = -t162 * t227 - t163 * t189;
t129 = -t160 * t226 + t161 * t185;
t128 = -t160 * t227 - t161 * t189;
t127 = t140 * t189 + t162 * t185;
t1 = [(g(1) * t188 - g(2) * t192) * MDP(2) + (g(1) * t192 + g(2) * t188) * MDP(3) + (g(1) * t161 - g(2) * t163) * MDP(9) + (g(1) * t207 - g(2) * t142) * MDP(16) + (-g(1) * t200 - g(2) * t141) * MDP(17) + (-g(1) * t199 - g(2) * t203) * MDP(19) + (-g(1) * t201 - g(2) * t139) * MDP(28) + (-g(1) * (-pkin(4) * t136 - pkin(5) * t123 + pkin(10) * t201 - qJ(6) * t122 + t199) - g(2) * (pkin(4) * t140 + pkin(5) * t127 + pkin(10) * t139 + qJ(6) * t126 + t203)) * MDP(30) + t216 * (g(1) * t123 - g(2) * t127) + t236 * (-g(1) * t122 + g(2) * t126) - t237 * (g(1) * t160 - g(2) * t162); (-g(1) * t217 - g(2) * t218 - g(3) * (-t183 * t220 + t164)) * MDP(19) + (-g(1) * (pkin(5) * t131 + qJ(6) * t130 - t206 * t162 + t217) - g(2) * (pkin(5) * t129 + qJ(6) * t128 - t206 * t160 + t218) + (-pkin(5) * t144 - qJ(6) * t143 - t164 - (t206 * t191 - t220) * t183) * g(3)) * MDP(30) + t236 * (g(1) * t130 + g(2) * t128 + g(3) * t143) + t237 * t193 + t216 * (-g(1) * t131 - g(2) * t129 - g(3) * t144) + (-MDP(16) * t190 + MDP(17) * t186 - t179 * MDP(28) - MDP(9)) * t194; (-g(1) * t141 + g(2) * t200 - g(3) * (t208 - t214)) * MDP(16) + (g(1) * t142 + g(2) * t207 - g(3) * (-t232 * t186 - t187 * t223)) * MDP(17) + (-g(1) * t168 - g(3) * t177 + (g(2) * t212 + t193 * t186) * pkin(3)) * MDP(19) + (-g(1) * t140 - g(2) * t136 - g(3) * t150) * MDP(28) + (-g(1) * (-pkin(3) * t228 + pkin(10) * t140 + t168) - g(2) * (-t200 * pkin(3) + t136 * pkin(10)) - g(3) * (-pkin(3) * t214 + pkin(10) * t150 + t177)) * MDP(30) + ((pkin(5) * t189 + qJ(6) * t185 + pkin(4)) * MDP(30) + t216 * t189 - t236 * t185) * (-g(3) * (-t179 * t225 + t232 * t180) - g(2) * t201 + g(1) * t139); (MDP(19) + MDP(30)) * t194; (-g(1) * (-pkin(5) * t126 + qJ(6) * t127) - g(2) * (-pkin(5) * t122 + qJ(6) * t123) - g(3) * (-pkin(5) * t133 + qJ(6) * t134)) * MDP(30) + t216 * t198 + t236 * (g(1) * t127 + g(2) * t123 + g(3) * t134); -t198 * MDP(30);];
taug  = t1;
