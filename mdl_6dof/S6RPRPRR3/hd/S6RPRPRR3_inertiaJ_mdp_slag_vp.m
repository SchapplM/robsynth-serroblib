% Calculate joint inertia matrix for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR3_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:42:28
% EndTime: 2019-03-09 03:42:30
% DurationCPUTime: 0.58s
% Computational Cost: add. (772->142), mult. (1513->209), div. (0->0), fcn. (1602->10), ass. (0->77)
t192 = MDP(15) * qJ(4);
t146 = sin(qJ(6));
t149 = cos(qJ(6));
t157 = (MDP(28) * t149 - MDP(29) * t146) * pkin(5);
t142 = sin(pkin(11));
t144 = cos(pkin(11));
t173 = t142 ^ 2 + t144 ^ 2;
t191 = t173 * t192;
t182 = pkin(8) + qJ(4);
t130 = t182 * t142;
t131 = t182 * t144;
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t115 = -t150 * t130 - t131 * t147;
t116 = -t130 * t147 + t131 * t150;
t127 = t142 * t147 - t150 * t144;
t128 = t142 * t150 + t144 * t147;
t105 = -pkin(9) * t128 + t115;
t106 = -pkin(9) * t127 + t116;
t111 = t149 * t127 + t128 * t146;
t112 = -t127 * t146 + t128 * t149;
t163 = t112 * MDP(25) - t111 * MDP(26) + (t105 * t149 - t106 * t146) * MDP(28) - (t105 * t146 + t106 * t149) * MDP(29);
t190 = t128 * MDP(18) - t127 * MDP(19) + t115 * MDP(21) - t116 * MDP(22) + t163;
t148 = sin(qJ(3));
t119 = t128 * t148;
t120 = t127 * t148;
t102 = t149 * t119 - t120 * t146;
t103 = -t119 * t146 - t120 * t149;
t179 = -t102 * MDP(28) - t103 * MDP(29);
t189 = -t119 * MDP(21) + t120 * MDP(22) + t179;
t181 = t103 * MDP(25) - t102 * MDP(26);
t151 = cos(qJ(3));
t145 = cos(pkin(10));
t136 = -pkin(1) * t145 - pkin(2);
t125 = -pkin(3) * t151 - qJ(4) * t148 + t136;
t122 = t144 * t125;
t143 = sin(pkin(10));
t135 = pkin(1) * t143 + pkin(7);
t178 = t135 * t142;
t104 = -pkin(8) * t144 * t148 + t122 + (-pkin(4) - t178) * t151;
t177 = t135 * t151;
t114 = t142 * t125 + t144 * t177;
t176 = t142 * t148;
t107 = -pkin(8) * t176 + t114;
t95 = t150 * t104 - t107 * t147;
t89 = -pkin(5) * t151 + pkin(9) * t120 + t95;
t96 = t104 * t147 + t107 * t150;
t90 = -pkin(9) * t119 + t96;
t86 = -t146 * t90 + t149 * t89;
t87 = t146 * t89 + t149 * t90;
t188 = t86 * MDP(28) - t87 * MDP(29) + t181;
t186 = -2 * MDP(17);
t185 = 0.2e1 * MDP(22);
t184 = -2 * MDP(24);
t183 = 0.2e1 * MDP(29);
t180 = pkin(3) * MDP(15);
t123 = pkin(4) * t176 + t148 * t135;
t172 = MDP(12) * t144;
t171 = MDP(13) * t142;
t170 = MDP(16) * t128;
t169 = MDP(21) * t127;
t168 = MDP(23) * t112;
t167 = MDP(28) * t111;
t166 = MDP(20) + MDP(27);
t137 = -pkin(4) * t144 - pkin(3);
t164 = t173 * MDP(14);
t113 = -t142 * t177 + t122;
t162 = -t113 * t142 + t114 * t144;
t161 = MDP(12) * t142 + MDP(13) * t144;
t160 = -t120 * MDP(18) - t119 * MDP(19);
t156 = t171 - t172 - t180;
t154 = MDP(22) * t128 + MDP(29) * t112 + t156 + t167 + t169;
t141 = t151 ^ 2;
t140 = t148 ^ 2;
t118 = pkin(5) * t127 + t137;
t108 = pkin(5) * t119 + t123;
t1 = [t140 * MDP(5) + (t140 * t135 ^ 2 + t113 ^ 2 + t114 ^ 2) * MDP(15) + MDP(1) + (t143 ^ 2 + t145 ^ 2) * MDP(4) * pkin(1) ^ 2 + t166 * t141 - (-MDP(16) * t120 + t119 * t186) * t120 + (MDP(23) * t103 + t102 * t184) * t103 + 0.2e1 * (-t136 * MDP(10) + t148 * MDP(6) - t160 - t181) * t151 + 0.2e1 * (-t113 * t151 + t140 * t178) * MDP(12) + 0.2e1 * (t135 * t140 * t144 + t114 * t151) * MDP(13) + 0.2e1 * (t119 * t123 - t151 * t95) * MDP(21) + (-t120 * t123 + t151 * t96) * t185 + 0.2e1 * (t102 * t108 - t151 * t86) * MDP(28) + (t103 * t108 + t151 * t87) * t183 + 0.2e1 * (t136 * MDP(11) + (-t113 * t144 - t114 * t142) * MDP(14)) * t148; (t162 - t177) * t148 * MDP(15); MDP(4) + (t173 * t140 + t141) * MDP(15); -t120 * t170 + (-t119 * t128 + t120 * t127) * MDP(17) + (t119 * t137 + t123 * t127) * MDP(21) + (-t120 * t137 + t123 * t128) * MDP(22) + t103 * t168 + (-t102 * t112 - t103 * t111) * MDP(24) + (t102 * t118 + t108 * t111) * MDP(28) + (t103 * t118 + t108 * t112) * MDP(29) + (MDP(7) - t161 * pkin(3) + (-MDP(10) + t156) * t135) * t148 + (-t135 * MDP(11) + t161 * qJ(4) + MDP(8) - t190) * t151 + (MDP(14) + t192) * t162; (-MDP(11) + t164 + t191) * t148 + (MDP(10) - t154) * t151; 0.2e1 * t137 * t169 + 0.2e1 * t118 * t167 + MDP(9) + (-0.2e1 * t171 + 0.2e1 * t172 + t180) * pkin(3) + (t127 * t186 + t137 * t185 + t170) * t128 + (t111 * t184 + t118 * t183 + t168) * t112 + (0.2e1 * t164 + t191) * qJ(4); (MDP(15) * t135 + t161) * t148 - t189; -t151 * MDP(15); t154; MDP(15); t95 * MDP(21) - t96 * MDP(22) + (-t166 - t157) * t151 + t160 + t188; t189; t190; 0; 0.2e1 * t157 + t166; -t151 * MDP(27) + t188; t179; t163; 0; MDP(27) + t157; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
