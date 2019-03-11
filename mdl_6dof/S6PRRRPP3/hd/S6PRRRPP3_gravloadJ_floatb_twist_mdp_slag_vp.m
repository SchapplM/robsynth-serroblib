% Calculate Gravitation load on the joints for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:58:31
% EndTime: 2019-03-08 22:58:32
% DurationCPUTime: 0.69s
% Computational Cost: add. (524->112), mult. (1348->173), div. (0->0), fcn. (1637->10), ass. (0->67)
t226 = MDP(18) - MDP(21) - MDP(24);
t225 = MDP(20) - MDP(17) - MDP(25);
t182 = cos(qJ(3));
t229 = -pkin(3) * t182 - pkin(2);
t178 = sin(qJ(4));
t228 = qJ(5) * t178 + pkin(3);
t227 = MDP(11) - MDP(19) - MDP(23);
t224 = pkin(5) + pkin(9);
t222 = cos(pkin(6));
t221 = cos(pkin(10));
t176 = sin(pkin(10));
t180 = sin(qJ(2));
t183 = cos(qJ(2));
t194 = t222 * t221;
t162 = t176 * t183 + t180 * t194;
t179 = sin(qJ(3));
t177 = sin(pkin(6));
t198 = t177 * t221;
t139 = -t162 * t179 - t182 * t198;
t181 = cos(qJ(4));
t219 = t139 * t181;
t199 = t176 * t222;
t164 = -t180 * t199 + t221 * t183;
t213 = t177 * t182;
t141 = -t164 * t179 + t176 * t213;
t218 = t141 * t181;
t161 = t176 * t180 - t183 * t194;
t217 = t161 * t179;
t163 = t221 * t180 + t183 * t199;
t216 = t163 * t179;
t214 = t177 * t180;
t165 = -t179 * t214 + t222 * t182;
t215 = t165 * t181;
t212 = t177 * t183;
t211 = t178 * t182;
t210 = t181 * t182;
t209 = t182 * t183;
t205 = MDP(22) + MDP(26);
t204 = t179 * t212;
t203 = t177 * t209;
t202 = pkin(4) * t219 + t228 * t139;
t201 = pkin(4) * t218 + t228 * t141;
t200 = pkin(4) * t215 + t228 * t165;
t140 = t162 * t182 - t179 * t198;
t121 = t140 * t178 - t161 * t181;
t122 = t140 * t181 + t161 * t178;
t197 = -t121 * pkin(4) + qJ(5) * t122;
t142 = t176 * t177 * t179 + t164 * t182;
t123 = t142 * t178 - t163 * t181;
t124 = t142 * t181 + t163 * t178;
t196 = -t123 * pkin(4) + qJ(5) * t124;
t166 = t222 * t179 + t180 * t213;
t143 = t166 * t178 + t181 * t212;
t144 = t166 * t181 - t178 * t212;
t195 = -t143 * pkin(4) + qJ(5) * t144;
t193 = g(1) * t123 + g(2) * t121 + g(3) * t143;
t192 = g(1) * t124 + g(2) * t122 + g(3) * t144;
t146 = t178 * t203 - t181 * t214;
t147 = (t178 * t180 + t181 * t209) * t177;
t187 = pkin(2) * t212 + pkin(3) * t203 + t147 * pkin(4) + pkin(8) * t214 + pkin(9) * t204 + t146 * qJ(5);
t128 = -t161 * t211 - t162 * t181;
t129 = -t161 * t210 + t162 * t178;
t185 = t129 * pkin(4) + pkin(8) * t162 - pkin(9) * t217 + qJ(5) * t128 + t229 * t161;
t130 = -t163 * t211 - t164 * t181;
t131 = -t163 * t210 + t164 * t178;
t184 = t131 * pkin(4) + pkin(8) * t164 - pkin(9) * t216 + qJ(5) * t130 + t229 * t163;
t1 = [(-MDP(1) - t205) * g(3); (g(1) * t164 + g(2) * t162 + g(3) * t214) * MDP(4) + (-g(1) * t184 - g(2) * t185 - g(3) * t187) * MDP(22) + (-g(1) * (-pkin(5) * t216 + qJ(6) * t131 + t184) - g(2) * (-pkin(5) * t217 + qJ(6) * t129 + t185) - g(3) * (pkin(5) * t204 + t147 * qJ(6) + t187)) * MDP(26) + t226 * (g(1) * t130 + g(2) * t128 + g(3) * t146) + t225 * (g(1) * t131 + g(2) * t129 + g(3) * t147) + (-MDP(10) * t182 + t179 * t227 - MDP(3)) * (-g(1) * t163 - g(2) * t161 + g(3) * t212); (-g(1) * (pkin(9) * t142 + t201) - g(2) * (pkin(9) * t140 + t202) - g(3) * (pkin(9) * t166 + t200)) * MDP(22) + (-g(1) * (qJ(6) * t218 + t224 * t142 + t201) - g(2) * (qJ(6) * t219 + t224 * t140 + t202) - g(3) * (qJ(6) * t215 + t224 * t166 + t200)) * MDP(26) + t227 * (g(1) * t142 + g(2) * t140 + g(3) * t166) + (t178 * t226 + t181 * t225 - MDP(10)) * (g(1) * t141 + g(2) * t139 + g(3) * t165); (-g(1) * t196 - g(2) * t197 - g(3) * t195) * MDP(22) + (-g(1) * (-qJ(6) * t123 + t196) - g(2) * (-qJ(6) * t121 + t197) - g(3) * (-qJ(6) * t143 + t195)) * MDP(26) - t225 * t193 + t226 * t192; -t205 * t193; -t192 * MDP(26);];
taug  = t1;
