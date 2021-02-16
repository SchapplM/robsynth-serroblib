% Calculate Gravitation load on the joints for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:06:22
% EndTime: 2021-01-16 01:06:28
% DurationCPUTime: 1.17s
% Computational Cost: add. (446->130), mult. (488->214), div. (0->0), fcn. (512->18), ass. (0->90)
t143 = qJ(4) + pkin(12);
t141 = cos(t143);
t148 = cos(pkin(6));
t134 = t148 * t141;
t139 = sin(t143);
t144 = qJ(2) + pkin(11);
t140 = sin(t144);
t146 = sin(pkin(6));
t142 = cos(t144);
t147 = cos(pkin(10));
t133 = t147 * t142;
t145 = sin(pkin(10));
t205 = t145 * t148;
t214 = t140 * t205 - t133;
t132 = t145 * t142;
t194 = t147 * t148;
t215 = t140 * t194 + t132;
t216 = -(g(1) * t214 - g(2) * t215) * t139 - g(3) * (-t139 * t140 * t146 + t134);
t138 = pkin(6) - t144;
t212 = sin(t138);
t210 = g(3) * t146;
t209 = t139 * t148;
t208 = t140 * t145;
t207 = t141 * t146;
t206 = t145 * t146;
t150 = sin(qJ(6));
t204 = t145 * t150;
t152 = sin(qJ(2));
t203 = t145 * t152;
t153 = cos(qJ(6));
t202 = t145 * t153;
t201 = t146 * t147;
t200 = t146 * t150;
t151 = sin(qJ(4));
t199 = t146 * t151;
t198 = t146 * t153;
t154 = cos(qJ(4));
t197 = t146 * t154;
t155 = cos(qJ(2));
t196 = t146 * t155;
t195 = t147 * t140;
t193 = t147 * t150;
t192 = t147 * t153;
t191 = t148 * t150;
t190 = t148 * t152;
t189 = t148 * t153;
t188 = t148 * t154;
t187 = t148 * t155;
t186 = MDP(16) + MDP(5);
t185 = g(3) * t200;
t184 = g(3) * t198;
t183 = cos(t138) / 0.2e1;
t180 = t141 * t204;
t179 = t141 * t202;
t178 = t141 * t193;
t177 = t141 * t192;
t176 = t145 * t197;
t175 = t145 * t191;
t174 = t145 * t189;
t173 = t147 * t197;
t172 = t147 * t191;
t171 = t147 * t189;
t170 = t147 * t187;
t169 = pkin(6) + t144;
t164 = sin(t169) / 0.2e1;
t128 = t164 - t212 / 0.2e1;
t168 = t128 * t147 + t132;
t167 = -t128 * t145 + t133;
t166 = g(1) * t145 - g(2) * t147;
t165 = cos(t169);
t163 = t166 * t146;
t160 = -t145 * t187 - t147 * t152;
t157 = t165 / 0.2e1 + t183;
t108 = -t147 * t157 + t208;
t111 = t145 * t157 + t195;
t129 = t212 / 0.2e1 + t164;
t159 = g(1) * t111 + g(2) * t108 - g(3) * t129;
t156 = -g(1) * t160 - g(3) * t196;
t149 = -qJ(5) - pkin(8);
t137 = pkin(4) * t154 + pkin(3);
t131 = pkin(2) * t170;
t130 = t183 - t165 / 0.2e1;
t121 = t141 * t174 - t193;
t120 = t141 * t171 + t204;
t119 = t141 * t175 + t192;
t118 = -t141 * t172 + t202;
t117 = t140 * t207 + t209;
t115 = t139 * t198 + t142 * t191;
t114 = -t139 * t200 + t142 * t189;
t1 = [(-MDP(1) - t186) * g(3); (-g(2) * (t170 - t203) + t156) * MDP(3) + (-g(1) * (t145 * t190 - t147 * t155) - g(2) * (-t145 * t155 - t147 * t190) + t152 * t210) * MDP(4) + (-g(2) * t131 + (g(2) * t203 + t156) * pkin(2)) * MDP(5) + (-g(1) * t167 - g(2) * t168 - g(3) * t130) * MDP(15) + (-g(1) * (pkin(2) * t160 - t111 * t137 - t149 * t167) - g(2) * (-pkin(2) * t203 - t108 * t137 - t149 * t168 + t131) - g(3) * (pkin(2) * t196 + t129 * t137 - t130 * t149)) * MDP(16) + ((g(1) * t121 - g(2) * t120 - t141 * t184) * t142 + (-g(1) * (-t175 - t177) - g(2) * (t172 - t179) - t185) * t140) * MDP(22) + ((-g(1) * t119 - g(2) * t118 + t141 * t185) * t142 + (-g(1) * (-t174 + t178) - g(2) * (t171 + t180) - t184) * t140) * MDP(23) + (MDP(13) * t141 - MDP(14) * t139) * t159 + (-MDP(11) * t154 + MDP(12) * t151) * (-g(1) * (t142 * t205 + t195) + g(2) * (t142 * t194 - t208) + t142 * t210); t186 * (-g(3) * t148 - t163); (-g(1) * (t151 * t214 + t176) - g(2) * (-t151 * t215 - t173) - g(3) * (-t140 * t199 + t188)) * MDP(11) + (-g(1) * (-t145 * t199 + t154 * t214) - g(2) * (t147 * t199 - t154 * t215) - g(3) * (-t140 * t197 - t148 * t151)) * MDP(12) + (-g(1) * (-t139 * t167 + t141 * t206) - g(2) * (-t139 * t168 - t141 * t201) - g(3) * (-t130 * t139 + t134)) * MDP(13) + (-g(1) * (-t139 * t206 - t141 * t167) - g(2) * (t139 * t201 - t141 * t168) - g(3) * (-t130 * t141 - t209)) * MDP(14) + (-g(1) * (-t151 * t167 + t176) - g(2) * (-t151 * t168 - t173) - g(3) * (-t130 * t151 + t188)) * pkin(4) * MDP(16) + (-t141 * t163 + t216) * t153 * MDP(22) + (t166 * t207 - t216) * t150 * MDP(23); -t159 * MDP(16); (-g(1) * (t114 * t145 + t119 * t140 - t142 * t178) - g(2) * (-t114 * t147 + t118 * t140 - t142 * t180) - g(3) * (-t117 * t150 - t142 * t198)) * MDP(22) + (-g(1) * (-t115 * t145 + t121 * t140 - t142 * t177) - g(2) * (t115 * t147 - t120 * t140 - t142 * t179) - g(3) * (-t117 * t153 + t142 * t200)) * MDP(23);];
taug = t1;
