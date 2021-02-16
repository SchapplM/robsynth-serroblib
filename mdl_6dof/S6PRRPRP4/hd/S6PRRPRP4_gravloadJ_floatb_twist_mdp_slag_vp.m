% Calculate Gravitation load on the joints for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:12:51
% EndTime: 2021-01-16 03:12:57
% DurationCPUTime: 1.39s
% Computational Cost: add. (359->129), mult. (858->213), div. (0->0), fcn. (951->10), ass. (0->75)
t167 = cos(qJ(2));
t162 = sin(qJ(5));
t154 = t162 * pkin(5) + qJ(4);
t157 = qJ(6) + pkin(3) + pkin(9);
t163 = sin(qJ(3));
t166 = cos(qJ(3));
t185 = t154 * t163 + t157 * t166;
t164 = sin(qJ(2));
t165 = cos(qJ(5));
t190 = t165 * pkin(5) + pkin(4) + pkin(8);
t232 = t190 * t164;
t238 = t232 + (pkin(2) + t185) * t167;
t237 = MDP(10) - MDP(13) + MDP(25);
t159 = sin(pkin(6));
t204 = t163 * t164;
t161 = cos(pkin(6));
t207 = t161 * t166;
t142 = t159 * t204 - t207;
t212 = t159 * t167;
t236 = g(3) * (t142 * t165 + t162 * t212);
t160 = cos(pkin(10));
t234 = t160 * t164;
t199 = t166 * t164;
t191 = t161 * t199;
t205 = t163 * t159;
t144 = t191 - t205;
t158 = sin(pkin(10));
t198 = t166 * t167;
t209 = t161 * t163;
t218 = g(3) * (t159 * t199 + t209);
t230 = -t218 - g(2) * (t160 * t144 + t158 * t198);
t229 = -MDP(11) + MDP(14);
t227 = MDP(21) + MDP(23);
t226 = MDP(22) + MDP(24);
t200 = t166 * t159;
t225 = (t154 * t207 - t157 * t209) * t164 - t154 * t205 - t157 * t200;
t224 = (pkin(3) * t209 - qJ(4) * t207) * t164 + pkin(3) * t200 + qJ(4) * t205;
t188 = t166 * pkin(3) + qJ(4) * t163;
t223 = (pkin(2) + t188) * t167 + pkin(8) * t164;
t219 = g(3) * t142;
t217 = g(3) * t159;
t215 = pkin(2) * t234;
t214 = t158 * t161;
t213 = t159 * t164;
t211 = t160 * t161;
t210 = t161 * t162;
t208 = t161 * t164;
t206 = t161 * t167;
t203 = t163 * t165;
t202 = t163 * t167;
t201 = t165 * t167;
t197 = MDP(15) + MDP(26);
t196 = t158 * t200;
t193 = t160 * t200;
t192 = t161 * t203;
t141 = t161 * t204 + t200;
t131 = t160 * t141 + t158 * t202;
t133 = -t158 * t141 + t160 * t202;
t187 = -t163 * pkin(3) + qJ(4) * t166;
t186 = t154 * t166 - t157 * t163;
t182 = -t158 * t164 + t160 * t206;
t140 = t158 * t206 + t234;
t181 = t167 * t187;
t180 = t188 * t161;
t179 = t167 * t186;
t169 = -t185 * t164 + t190 * t167;
t155 = t158 * pkin(2);
t151 = pkin(2) * t211;
t150 = t160 * t198;
t145 = t158 * t205;
t139 = t158 * t167 + t160 * t208;
t137 = t158 * t208 - t160 * t167;
t136 = t140 * t163;
t135 = t182 * t163;
t1 = [(-MDP(1) - t197) * g(3); (-g(1) * (-t215 + (pkin(8) * t167 - t188 * t164) * t160 - t223 * t214) - g(2) * (-(-pkin(8) * t211 + t188 * t158 + t155) * t164 + (t158 * pkin(8) + t160 * t180 + t151) * t167) - t223 * t217) * MDP(15) + (-g(1) * (t169 * t160 - t238 * t214 - t215) - g(2) * (t151 * t167 - t155 * t164 + (t185 * t167 + t232) * t211 + t169 * t158) - t238 * t217) * MDP(26) + (MDP(4) - MDP(12)) * (-g(1) * t137 + g(2) * t139 + g(3) * t213) + t226 * (-g(1) * (-t136 * t165 + t137 * t162) - g(2) * (t135 * t165 - t139 * t162) - (-t162 * t164 + t163 * t201) * t217) + t227 * (-g(1) * (-t136 * t162 - t137 * t165) - g(2) * (t135 * t162 + t139 * t165) - (t162 * t202 + t164 * t165) * t217) + (t229 * t163 + t237 * t166 + MDP(3)) * (g(1) * t140 - g(2) * t182 - g(3) * t212); (-g(3) * (t187 * t213 + t180) + (-g(1) * t181 + g(2) * t224) * t160 + (-g(1) * t224 - g(2) * t181) * t158) * MDP(15) + (-g(3) * (t185 * t161 + t186 * t213) + (-g(1) * t179 - g(2) * t225) * t160 + (g(1) * t225 - g(2) * t179) * t158) * MDP(26) + t229 * (g(1) * (t137 * t166 - t145) - g(2) * (t139 * t166 - t160 * t205) - t218) + t237 * (g(2) * (t139 * t163 + t193) + t219 - g(1) * (t137 * t163 + t196)) + (t227 * t162 + t226 * t165) * (-g(1) * (-t158 * t144 + t150) + t230); t197 * (-g(1) * t133 - g(2) * t131 - t219); (-g(1) * ((-t158 * t210 + t160 * t203) * t167 + (-t158 * t192 - t160 * t162) * t164 - t165 * t196) - g(2) * ((t158 * t203 + t160 * t210) * t167 + (-t158 * t162 + t160 * t192) * t164 + t165 * t193) - t236) * pkin(5) * MDP(26) + t227 * (-g(1) * (t133 * t165 - t140 * t162) - g(2) * (t131 * t165 + t162 * t182) - t236) + t226 * (-g(1) * (-t133 * t162 - t140 * t165) - g(2) * (-t131 * t162 + t165 * t182) - g(3) * (-t142 * t162 + t159 * t201)); (-g(1) * (-t158 * t191 + t145 + t150) + t230) * MDP(26);];
taug = t1;
