% Calculate joint inertia matrix for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% MDP [34x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR7_inertiaJ_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:17:57
% EndTime: 2019-03-09 07:17:58
% DurationCPUTime: 0.46s
% Computational Cost: add. (806->135), mult. (1359->175), div. (0->0), fcn. (1542->8), ass. (0->73)
t140 = sin(qJ(4));
t141 = sin(qJ(3));
t144 = cos(qJ(4));
t145 = cos(qJ(3));
t120 = t140 * t145 + t141 * t144;
t121 = -t140 * t141 + t144 * t145;
t139 = sin(qJ(5));
t143 = cos(qJ(5));
t101 = t120 * t139 - t121 * t143;
t187 = t101 ^ 2;
t138 = sin(qJ(6));
t142 = cos(qJ(6));
t175 = MDP(30) * t138 + MDP(31) * t142;
t173 = MDP(33) * t142;
t186 = -MDP(34) * t138 + t173;
t149 = 0.2e1 * t186;
t185 = MDP(33) * t138 + MDP(34) * t142;
t127 = pkin(3) * t141 + qJ(2);
t108 = pkin(4) * t120 + t127;
t184 = 0.2e1 * t108;
t183 = 0.2e1 * t127;
t182 = (pkin(1) * MDP(6));
t181 = pkin(3) * t140;
t180 = pkin(5) * t138;
t179 = t143 * pkin(4);
t146 = -pkin(1) - pkin(7);
t178 = -pkin(8) + t146;
t122 = t178 * t141;
t123 = t178 * t145;
t163 = -t122 * t140 + t123 * t144;
t90 = -pkin(9) * t121 + t163;
t156 = -t122 * t144 - t123 * t140;
t91 = -pkin(9) * t120 - t156;
t86 = t139 * t91 - t143 * t90;
t83 = t86 * t138;
t177 = t86 * t142;
t176 = t138 * t142;
t157 = t120 * t143 + t121 * t139;
t170 = t157 * MDP(32);
t130 = pkin(3) * t144 + pkin(4);
t124 = t143 * t130;
t111 = -t139 * t181 + t124;
t169 = t111 * MDP(26);
t112 = -t130 * t139 - t143 * t181;
t168 = t112 * MDP(27);
t167 = t139 * MDP(27);
t166 = t144 * MDP(19);
t165 = MDP(29) * t176;
t136 = t138 ^ 2;
t164 = MDP(28) * t136 + MDP(25) + 0.2e1 * t165;
t162 = MDP(18) + t164;
t161 = pkin(5) * t101 - pkin(10) * t157;
t109 = -pkin(5) - t111;
t110 = pkin(10) - t112;
t159 = -t101 * t109 - t110 * t157;
t128 = pkin(4) * t139 + pkin(10);
t129 = -pkin(5) - t179;
t158 = -t101 * t129 - t128 * t157;
t155 = MDP(30) * t142 - MDP(31) * t138;
t152 = -MDP(27) * t157 + (-MDP(26) - t186) * t101;
t151 = (MDP(26) * t143 - t167) * pkin(4);
t137 = t142 ^ 2;
t87 = t139 * t90 + t143 * t91;
t150 = -t86 * MDP(26) - t87 * MDP(27) + (-MDP(24) + t175) * t157 + (-(-t136 + t137) * MDP(29) - MDP(28) * t176 - MDP(23)) * t101;
t148 = MDP(19) * t121 - MDP(20) * t120 + t152;
t147 = MDP(16) * t121 - MDP(17) * t120 + MDP(19) * t163 + MDP(20) * t156 + t150;
t135 = pkin(5) * t142;
t125 = t129 * t138;
t107 = t109 * t138;
t88 = pkin(5) * t157 + pkin(10) * t101 + t108;
t82 = t138 * t88 + t142 * t87;
t81 = -t138 * t87 + t142 * t88;
t1 = [t120 * MDP(19) * t183 - t101 * MDP(27) * t184 + MDP(1) + (MDP(7) * t145 - 0.2e1 * MDP(8) * t141) * t145 + ((-2 * MDP(4) + t182) * pkin(1)) + (MDP(26) * t184 + t170 - 0.2e1 * (-MDP(22) + t155) * t101) * t157 + 0.2e1 * (-t101 * t83 + t157 * t81) * MDP(33) + 0.2e1 * (-t101 * t177 - t157 * t82) * MDP(34) + (MDP(14) * t121 - 0.2e1 * MDP(15) * t120 + MDP(20) * t183) * t121 + (MDP(28) * t137 + MDP(21) - 0.2e1 * t165) * t187 + (0.2e1 * MDP(12) * t141 + 0.2e1 * MDP(13) * t145 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2); MDP(4) - t182 + t185 * (-t157 ^ 2 - t187); MDP(6); (t142 * t159 + t83) * MDP(34) + (t138 * t159 - t177) * MDP(33) + (-MDP(13) * t146 - MDP(10)) * t141 + (MDP(12) * t146 + MDP(9)) * t145 + t147; MDP(12) * t145 - MDP(13) * t141 + t148; MDP(11) - t109 * t149 + 0.2e1 * (-MDP(20) * t140 + t166) * pkin(3) + 0.2e1 * t169 + 0.2e1 * t168 + t162; (t142 * t158 + t83) * MDP(34) + (t138 * t158 - t177) * MDP(33) + t147; t148; (t124 + t179) * MDP(26) + (t125 + t107) * MDP(34) + (-t109 - t129) * t173 + (-pkin(4) - t130) * t167 + (t166 + (-MDP(26) * t139 - MDP(27) * t143 - MDP(20)) * t140) * pkin(3) + t162; -t129 * t149 + 0.2e1 * t151 + t162; (t138 * t161 - t177) * MDP(33) + (t142 * t161 + t83) * MDP(34) + t150; t152; t169 + t168 + (-t109 * t142 + t135) * MDP(33) + (t107 - t180) * MDP(34) + t164; (-t129 * t142 + t135) * MDP(33) + (t125 - t180) * MDP(34) + t151 + t164; pkin(5) * t149 + t164; t81 * MDP(33) - t82 * MDP(34) - t101 * t155 + t170; -t185 * t157; -t110 * t185 + t175; -t128 * t185 + t175; -pkin(10) * t185 + t175; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
