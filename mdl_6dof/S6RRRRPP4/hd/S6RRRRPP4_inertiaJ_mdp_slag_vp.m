% Calculate joint inertia matrix for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRRPP4_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:02:38
% EndTime: 2019-03-09 21:02:41
% DurationCPUTime: 1.05s
% Computational Cost: add. (1911->233), mult. (3635->326), div. (0->0), fcn. (3829->8), ass. (0->87)
t173 = sin(qJ(3));
t176 = cos(qJ(3));
t172 = sin(qJ(4));
t175 = cos(qJ(4));
t207 = pkin(8) + pkin(9);
t189 = t207 * t176;
t190 = t207 * t173;
t133 = -t172 * t190 + t175 * t189;
t150 = t172 * t173 - t175 * t176;
t123 = -t150 * qJ(5) + t133;
t170 = sin(pkin(10));
t171 = cos(pkin(10));
t132 = -t172 * t189 - t175 * t190;
t151 = t172 * t176 + t175 * t173;
t180 = -t151 * qJ(5) + t132;
t114 = t170 * t123 - t171 * t180;
t116 = t171 * t123 + t170 * t180;
t181 = t151 * MDP(20) - t150 * MDP(21) + t132 * MDP(23) - t133 * MDP(24) - t114 * MDP(27) + t116 * MDP(29);
t219 = t181 - (t173 * MDP(16) + t176 * MDP(17)) * pkin(8) + t173 * MDP(13) + t176 * MDP(14);
t194 = t175 * MDP(23);
t218 = (-t172 * MDP(24) + t194) * pkin(3);
t174 = sin(qJ(2));
t142 = t151 * t174;
t143 = t150 * t174;
t217 = t143 * MDP(20) + t142 * MDP(21);
t177 = cos(qJ(2));
t153 = -t177 * pkin(2) - t174 * pkin(8) - pkin(1);
t149 = t176 * t153;
t202 = pkin(9) * t174;
t205 = pkin(7) * t173;
t127 = -t176 * t202 + t149 + (-pkin(3) - t205) * t177;
t203 = pkin(7) * t177;
t191 = t176 * t203;
t129 = t191 + (t153 - t202) * t173;
t118 = t175 * t127 - t172 * t129;
t108 = -t177 * pkin(4) + t143 * qJ(5) + t118;
t119 = t172 * t127 + t175 * t129;
t117 = -t142 * qJ(5) + t119;
t101 = t171 * t108 - t170 * t117;
t102 = t170 * t108 + t171 * t117;
t216 = t118 * MDP(23) - t119 * MDP(24) + t101 * MDP(27) + t102 * MDP(29) - t217;
t214 = -2 * MDP(19);
t213 = 0.2e1 * MDP(23);
t212 = 0.2e1 * MDP(24);
t211 = 2 * MDP(25);
t210 = 0.2e1 * MDP(27);
t209 = 2 * MDP(28);
t208 = 0.2e1 * MDP(29);
t206 = pkin(3) * t172;
t204 = pkin(7) * t176;
t201 = t173 * t176;
t162 = t175 * pkin(3) + pkin(4);
t145 = t170 * t162 + t171 * t206;
t152 = (pkin(3) * t173 + pkin(7)) * t174;
t125 = t171 * t150 + t170 * t151;
t126 = -t170 * t150 + t171 * t151;
t163 = -t176 * pkin(3) - pkin(2);
t137 = t150 * pkin(4) + t163;
t112 = t125 * pkin(5) - t126 * qJ(6) + t137;
t200 = t112 * MDP(30);
t197 = t125 * MDP(27);
t196 = t126 * MDP(29);
t195 = t143 * MDP(18);
t193 = MDP(15) + MDP(22);
t192 = t114 ^ 2 + t116 ^ 2;
t188 = MDP(12) * t201;
t120 = t171 * t142 - t170 * t143;
t121 = -t170 * t142 - t171 * t143;
t187 = t114 * t121 - t116 * t120;
t155 = t171 * t162;
t144 = -t170 * t206 + t155;
t128 = t142 * pkin(4) + t152;
t185 = t176 * MDP(13) - t173 * MDP(14);
t168 = t176 ^ 2;
t167 = t174 ^ 2;
t166 = t173 ^ 2;
t164 = t170 * pkin(4);
t159 = t171 * pkin(4) + pkin(5);
t158 = t164 + qJ(6);
t141 = -pkin(5) - t144;
t140 = qJ(6) + t145;
t136 = t173 * t153 + t191;
t135 = -t173 * t203 + t149;
t104 = t120 * pkin(5) - t121 * qJ(6) + t128;
t100 = t177 * pkin(5) - t101;
t99 = -t177 * qJ(6) + t102;
t1 = [-0.2e1 * pkin(1) * t174 * MDP(10) + MDP(1) + (t101 ^ 2 + t102 ^ 2 + t128 ^ 2) * MDP(26) + (t100 ^ 2 + t104 ^ 2 + t99 ^ 2) * MDP(30) + t193 * t177 ^ 2 - (t142 * t214 - t195) * t143 + (t168 * MDP(11) + MDP(4) - 0.2e1 * t188) * t167 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t185) * t174 + t217) * t177 + 0.2e1 * (t136 * t177 + t167 * t204) * MDP(17) + 0.2e1 * (-t135 * t177 + t167 * t205) * MDP(16) + (-t118 * t177 + t152 * t142) * t213 + (t119 * t177 - t152 * t143) * t212 + (-t104 * t121 - t99 * t177) * t208 + (t100 * t177 + t104 * t120) * t210 + (t100 * t121 - t99 * t120) * t209 + (-t101 * t121 - t102 * t120) * t211; -t151 * t195 + (-t151 * t142 + t143 * t150) * MDP(19) + (t163 * t142 + t152 * t150) * MDP(23) + (-t163 * t143 + t152 * t151) * MDP(24) + (-t101 * t126 - t102 * t125 + t187) * MDP(25) + (-t101 * t114 + t102 * t116 + t128 * t137) * MDP(26) + (t104 * t125 + t112 * t120) * MDP(27) + (t100 * t126 - t99 * t125 + t187) * MDP(28) + (-t104 * t126 - t112 * t121) * MDP(29) + (t100 * t114 + t104 * t112 + t99 * t116) * MDP(30) + (-pkin(7) * MDP(10) + MDP(7) - t219) * t177 + (MDP(6) - pkin(7) * MDP(9) + MDP(11) * t201 + (-t166 + t168) * MDP(12) + (-pkin(2) * t173 - t204) * MDP(16) + (-pkin(2) * t176 + t205) * MDP(17)) * t174; MDP(8) + t166 * MDP(11) + 0.2e1 * t188 + t163 * t150 * t213 + (t137 ^ 2 + t192) * MDP(26) + t192 * MDP(30) + (-0.2e1 * t196 + 0.2e1 * t197 + t200) * t112 + 0.2e1 * (t176 * MDP(16) - t173 * MDP(17)) * pkin(2) + (MDP(18) * t151 + t150 * t214 + t163 * t212) * t151 + (t211 + t209) * (t114 * t126 - t116 * t125); t135 * MDP(16) - t136 * MDP(17) + (-t145 * t120 - t144 * t121) * MDP(25) + (t101 * t144 + t102 * t145) * MDP(26) + (-t140 * t120 + t141 * t121) * MDP(28) + (t100 * t141 + t99 * t140) * MDP(30) + t185 * t174 + ((-pkin(5) + t141) * MDP(27) + (-qJ(6) - t140) * MDP(29) - t218 - t193) * t177 + t216; (-t145 * t125 - t144 * t126) * MDP(25) + (-t114 * t144 + t116 * t145) * MDP(26) + (-t140 * t125 + t141 * t126) * MDP(28) + (t114 * t141 + t116 * t140) * MDP(30) + t219; (t144 ^ 2 + t145 ^ 2) * MDP(26) + (t140 ^ 2 + t141 ^ 2) * MDP(30) + 0.2e1 * t218 - 0.2e1 * t141 * MDP(27) + t140 * t208 + t193; (-t158 * t120 - t159 * t121) * MDP(28) + (-t100 * t159 + t99 * t158) * MDP(30) + (-MDP(22) + (-pkin(5) - t159) * MDP(27) + (-qJ(6) - t158) * MDP(29)) * t177 + ((-t120 * t170 - t121 * t171) * MDP(25) + (t101 * t171 + t102 * t170) * MDP(26)) * pkin(4) + t216; (-t158 * t125 - t159 * t126) * MDP(28) + (-t114 * t159 + t116 * t158) * MDP(30) + ((-t125 * t170 - t126 * t171) * MDP(25) + (-t114 * t171 + t116 * t170) * MDP(26)) * pkin(4) + t181; MDP(22) + (0.2e1 * pkin(5) + t155) * MDP(27) + (t164 + 0.2e1 * qJ(6) + t145) * MDP(29) + (t140 * t158 - t141 * t159) * MDP(30) + ((t144 * t171 + t145 * t170) * MDP(26) + t171 * MDP(27)) * pkin(4) + (t194 + (-MDP(27) * t170 - MDP(24)) * t172) * pkin(3); MDP(22) + (t158 ^ 2 + t159 ^ 2) * MDP(30) + (t170 ^ 2 + t171 ^ 2) * MDP(26) * pkin(4) ^ 2 + t159 * t210 + t158 * t208; t128 * MDP(26) + t120 * MDP(27) - t121 * MDP(29) + t104 * MDP(30); t137 * MDP(26) - t196 + t197 + t200; 0; 0; MDP(26) + MDP(30); t177 * MDP(27) + t121 * MDP(28) + t100 * MDP(30); t126 * MDP(28) + t114 * MDP(30); t141 * MDP(30) - MDP(27); -t159 * MDP(30) - MDP(27); 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
