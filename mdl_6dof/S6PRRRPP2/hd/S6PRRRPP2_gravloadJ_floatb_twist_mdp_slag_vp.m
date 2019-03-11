% Calculate Gravitation load on the joints for
% S6PRRRPP2
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
%   see S6PRRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:52:41
% EndTime: 2019-03-08 22:52:43
% DurationCPUTime: 0.61s
% Computational Cost: add. (519->112), mult. (1335->173), div. (0->0), fcn. (1619->10), ass. (0->67)
t183 = cos(qJ(3));
t228 = -pkin(3) * t183 - pkin(2);
t179 = sin(qJ(4));
t227 = qJ(5) * t179 + pkin(3);
t226 = MDP(18) - MDP(21) - MDP(24);
t225 = MDP(11) - MDP(20) + MDP(25);
t223 = pkin(9) - qJ(6);
t222 = cos(pkin(6));
t221 = cos(pkin(10));
t177 = sin(pkin(10));
t181 = sin(qJ(2));
t184 = cos(qJ(2));
t194 = t222 * t221;
t163 = t177 * t184 + t181 * t194;
t180 = sin(qJ(3));
t178 = sin(pkin(6));
t199 = t178 * t221;
t140 = -t163 * t180 - t183 * t199;
t182 = cos(qJ(4));
t219 = t140 * t182;
t200 = t177 * t222;
t165 = -t181 * t200 + t221 * t184;
t213 = t178 * t183;
t142 = -t165 * t180 + t177 * t213;
t218 = t142 * t182;
t162 = t177 * t181 - t184 * t194;
t217 = t162 * t180;
t164 = t221 * t181 + t184 * t200;
t216 = t164 * t180;
t214 = t178 * t181;
t166 = -t180 * t214 + t222 * t183;
t215 = t166 * t182;
t212 = t178 * t184;
t211 = t179 * t183;
t210 = t182 * t183;
t209 = t183 * t184;
t206 = MDP(22) + MDP(26);
t205 = t180 * t212;
t204 = t179 * t212;
t203 = pkin(4) * t219 + t227 * t140;
t202 = pkin(4) * t218 + t227 * t142;
t201 = pkin(4) * t215 + t227 * t166;
t198 = MDP(17) + MDP(19) + MDP(23);
t141 = t163 * t183 - t180 * t199;
t122 = t141 * t179 - t162 * t182;
t123 = t141 * t182 + t162 * t179;
t197 = -t122 * pkin(4) + qJ(5) * t123;
t143 = t177 * t178 * t180 + t165 * t183;
t124 = t143 * t179 - t164 * t182;
t125 = t143 * t182 + t164 * t179;
t196 = -t124 * pkin(4) + qJ(5) * t125;
t167 = t222 * t180 + t181 * t213;
t144 = t167 * t179 + t182 * t212;
t145 = t167 * t182 - t204;
t195 = -t144 * pkin(4) + qJ(5) * t145;
t193 = g(1) * t124 + g(2) * t122 + g(3) * t144;
t190 = g(1) * t142 + g(2) * t140 + g(3) * t166;
t147 = -t182 * t214 + t183 * t204;
t148 = (t179 * t181 + t182 * t209) * t178;
t188 = t178 * pkin(3) * t209 + pkin(2) * t212 + t148 * pkin(4) + pkin(8) * t214 + pkin(9) * t205 + t147 * qJ(5);
t129 = -t162 * t211 - t163 * t182;
t130 = -t162 * t210 + t163 * t179;
t186 = t130 * pkin(4) + t163 * pkin(8) - pkin(9) * t217 + t129 * qJ(5) + t228 * t162;
t131 = -t164 * t211 - t165 * t182;
t132 = -t164 * t210 + t165 * t179;
t185 = t132 * pkin(4) + t165 * pkin(8) - pkin(9) * t216 + t131 * qJ(5) + t228 * t164;
t1 = [(-MDP(1) - t206) * g(3); (g(1) * t165 + g(2) * t163 + g(3) * t214) * MDP(4) + (-g(1) * t185 - g(2) * t186 - g(3) * t188) * MDP(22) + (-g(1) * (t132 * pkin(5) + qJ(6) * t216 + t185) - g(2) * (t130 * pkin(5) + qJ(6) * t217 + t186) - g(3) * (t148 * pkin(5) - qJ(6) * t205 + t188)) * MDP(26) + t198 * (-g(1) * t132 - g(2) * t130 - g(3) * t148) + t226 * (g(1) * t131 + g(2) * t129 + g(3) * t147) + (-MDP(10) * t183 + t225 * t180 - MDP(3)) * (-g(1) * t164 - g(2) * t162 + g(3) * t212); (-g(1) * (t143 * pkin(9) + t202) - g(2) * (t141 * pkin(9) + t203) - g(3) * (pkin(9) * t167 + t201)) * MDP(22) + (-g(1) * (pkin(5) * t218 + t223 * t143 + t202) - g(2) * (pkin(5) * t219 + t223 * t141 + t203) - g(3) * (pkin(5) * t215 + t223 * t167 + t201)) * MDP(26) + t225 * (g(1) * t143 + g(2) * t141 + g(3) * t167) + (t226 * t179 - t198 * t182 - MDP(10)) * t190; (-g(1) * t196 - g(2) * t197 - g(3) * t195) * MDP(22) + (-g(1) * (-pkin(5) * t124 + t196) - g(2) * (-pkin(5) * t122 + t197) - g(3) * (-t144 * pkin(5) + t195)) * MDP(26) + t198 * t193 + t226 * (g(1) * t125 + g(2) * t123 + g(3) * t145); -t206 * t193; -t190 * MDP(26);];
taug  = t1;
