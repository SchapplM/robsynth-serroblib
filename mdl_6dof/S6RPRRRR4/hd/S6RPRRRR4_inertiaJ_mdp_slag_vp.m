% Calculate joint inertia matrix for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR4_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:05:56
% EndTime: 2019-03-09 07:05:58
% DurationCPUTime: 0.49s
% Computational Cost: add. (1376->142), mult. (2575->186), div. (0->0), fcn. (3151->10), ass. (0->84)
t144 = sin(qJ(6));
t148 = cos(qJ(6));
t179 = t144 * MDP(31) + t148 * MDP(32);
t170 = t148 * MDP(34);
t158 = -t144 * MDP(35) + t170;
t153 = 0.2e1 * t158;
t142 = sin(pkin(11));
t143 = cos(pkin(11));
t147 = sin(qJ(3));
t189 = cos(qJ(3));
t122 = t147 * t142 - t189 * t143;
t123 = t189 * t142 + t147 * t143;
t146 = sin(qJ(4));
t150 = cos(qJ(4));
t113 = t150 * t122 + t146 * t123;
t130 = -t143 * pkin(2) - pkin(1);
t115 = t122 * pkin(3) + t130;
t106 = t113 * pkin(4) + t115;
t192 = 0.2e1 * t106;
t191 = 0.2e1 * t115;
t190 = 0.2e1 * t130;
t188 = pkin(1) * MDP(7);
t187 = pkin(3) * t146;
t186 = pkin(5) * t144;
t149 = cos(qJ(5));
t185 = t149 * pkin(4);
t184 = pkin(7) + qJ(2);
t145 = sin(qJ(5));
t114 = -t146 * t122 + t150 * t123;
t124 = t184 * t142;
t125 = t184 * t143;
t165 = -t189 * t124 - t147 * t125;
t108 = -t123 * pkin(8) + t165;
t156 = t147 * t124 - t189 * t125;
t109 = -t122 * pkin(8) - t156;
t166 = t150 * t108 - t146 * t109;
t94 = -t114 * pkin(9) + t166;
t160 = -t146 * t108 - t150 * t109;
t95 = -t113 * pkin(9) - t160;
t90 = t145 * t95 - t149 * t94;
t87 = t90 * t144;
t183 = t90 * t148;
t182 = t142 * MDP(5);
t181 = t143 * MDP(4);
t180 = t144 * t148;
t104 = t149 * t113 + t145 * t114;
t177 = t104 * MDP(33);
t105 = -t145 * t113 + t149 * t114;
t176 = t105 * MDP(28);
t175 = t113 * MDP(20);
t133 = t150 * pkin(3) + pkin(4);
t126 = t149 * t133;
t119 = -t145 * t187 + t126;
t174 = t119 * MDP(27);
t120 = -t145 * t133 - t149 * t187;
t173 = t120 * MDP(28);
t172 = t122 * MDP(13);
t171 = t145 * MDP(28);
t169 = t150 * MDP(20);
t168 = MDP(30) * t180;
t140 = t144 ^ 2;
t167 = t140 * MDP(29) + MDP(26) + 0.2e1 * t168;
t164 = MDP(19) + t167;
t163 = -pkin(5) * t105 - pkin(10) * t104;
t117 = -pkin(5) - t119;
t118 = pkin(10) - t120;
t162 = -t104 * t118 + t105 * t117;
t131 = t145 * pkin(4) + pkin(10);
t132 = -pkin(5) - t185;
t161 = -t104 * t131 + t105 * t132;
t159 = MDP(31) * t148 - MDP(32) * t144;
t157 = -MDP(34) * t144 - MDP(35) * t148;
t155 = (t149 * MDP(27) - t171) * pkin(4);
t141 = t148 ^ 2;
t91 = t145 * t94 + t149 * t95;
t154 = -t90 * MDP(27) - t91 * MDP(28) + (MDP(24) + (-t140 + t141) * MDP(30) + MDP(29) * t180) * t105 + (-MDP(25) + t179) * t104;
t152 = t114 * MDP(17) - t113 * MDP(18) + t166 * MDP(20) + t160 * MDP(21) + t154;
t137 = pkin(5) * t148;
t127 = t132 * t144;
t116 = t117 * t144;
t92 = t104 * pkin(5) - t105 * pkin(10) + t106;
t86 = t144 * t92 + t148 * t91;
t85 = -t144 * t91 + t148 * t92;
t1 = [t172 * t190 + t175 * t191 + t176 * t192 + MDP(1) + (0.2e1 * t181 - 0.2e1 * t182 + t188) * pkin(1) + (MDP(14) * t190 + MDP(8) * t123 - 0.2e1 * t122 * MDP(9)) * t123 + (MDP(15) * t114 - 0.2e1 * t113 * MDP(16) + MDP(21) * t191) * t114 + (t141 * MDP(29) + MDP(22) - 0.2e1 * t168) * t105 ^ 2 + (MDP(27) * t192 + t177 + 0.2e1 * (-MDP(23) + t159) * t105) * t104 + 0.2e1 * (t85 * t104 + t105 * t87) * MDP(34) + 0.2e1 * (-t86 * t104 + t105 * t183) * MDP(35) + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t142 ^ 2 + t143 ^ 2) * qJ(2); t172 + t123 * MDP(14) + t175 + t114 * MDP(21) + t176 - t181 + t182 - t188 + (MDP(27) + t158) * t104; MDP(7); (t162 * t144 - t183) * MDP(34) + (t162 * t148 + t87) * MDP(35) + t152 + t165 * MDP(13) + t156 * MDP(14) - t122 * MDP(11) + t123 * MDP(10); 0; MDP(12) - t117 * t153 + 0.2e1 * (-t146 * MDP(21) + t169) * pkin(3) + 0.2e1 * t174 + 0.2e1 * t173 + t164; (t161 * t144 - t183) * MDP(34) + (t161 * t148 + t87) * MDP(35) + t152; 0; (t126 + t185) * MDP(27) + (t127 + t116) * MDP(35) + (-t117 - t132) * t170 + (-pkin(4) - t133) * t171 + (t169 + (-MDP(27) * t145 - MDP(28) * t149 - MDP(21)) * t146) * pkin(3) + t164; -t132 * t153 + 0.2e1 * t155 + t164; (t163 * t144 - t183) * MDP(34) + (t163 * t148 + t87) * MDP(35) + t154; 0; t174 + t173 + (-t117 * t148 + t137) * MDP(34) + (t116 - t186) * MDP(35) + t167; (-t132 * t148 + t137) * MDP(34) + (t127 - t186) * MDP(35) + t155 + t167; pkin(5) * t153 + t167; t85 * MDP(34) - t86 * MDP(35) + t159 * t105 + t177; t158; t157 * t118 + t179; t157 * t131 + t179; t157 * pkin(10) + t179; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
