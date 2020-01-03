% Calculate joint inertia matrix for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR7_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:18
% EndTime: 2019-12-31 21:17:20
% DurationCPUTime: 0.50s
% Computational Cost: add. (659->135), mult. (1244->190), div. (0->0), fcn. (1315->8), ass. (0->75)
t150 = sin(qJ(3));
t136 = t150 * pkin(2) + qJ(4);
t147 = sin(pkin(9));
t148 = cos(pkin(9));
t171 = t147 ^ 2 + t148 ^ 2;
t173 = t171 * t136;
t149 = sin(qJ(5));
t152 = cos(qJ(5));
t123 = t149 * t147 - t152 * t148;
t124 = t152 * t147 + t149 * t148;
t175 = t124 * MDP(24) - t123 * MDP(25);
t191 = t123 * MDP(27) + t124 * MDP(28);
t190 = t148 * MDP(18) - t147 * MDP(19);
t153 = cos(qJ(2));
t140 = -t153 * pkin(2) - pkin(1);
t189 = 0.2e1 * t140;
t188 = 0.2e1 * t153;
t187 = 2 * MDP(20);
t186 = -2 * MDP(23);
t185 = -pkin(7) - pkin(6);
t184 = cos(qJ(3));
t183 = t148 * pkin(4);
t144 = t148 * pkin(8);
t181 = pkin(3) * MDP(21);
t151 = sin(qJ(2));
t125 = t150 * t151 - t184 * t153;
t126 = t150 * t153 + t184 * t151;
t102 = t125 * pkin(3) - t126 * qJ(4) + t140;
t132 = t185 * t151;
t133 = t185 * t153;
t113 = t150 * t132 - t184 * t133;
t91 = t147 * t102 + t148 * t113;
t97 = t123 * t126;
t180 = MDP(22) * t97;
t112 = -t184 * t132 - t150 * t133;
t179 = t112 * t148;
t178 = t126 * t147;
t96 = t124 * t126;
t177 = t96 * MDP(25);
t176 = t97 * MDP(24);
t172 = t171 * qJ(4);
t170 = t125 * MDP(26);
t139 = -t184 * pkin(2) - pkin(3);
t169 = t139 * MDP(21);
t167 = MDP(15) + (MDP(22) * t124 + t123 * t186) * t124;
t90 = t148 * t102 - t147 * t113;
t166 = t171 * MDP(21);
t165 = -pkin(3) * t126 - qJ(4) * t125;
t164 = -t90 * t147 + t91 * t148;
t88 = t125 * pkin(4) - t126 * t144 + t90;
t89 = -pkin(8) * t178 + t91;
t163 = (-t149 * t89 + t152 * t88) * MDP(27) - (t149 * t88 + t152 * t89) * MDP(28);
t162 = t96 * MDP(27) - t97 * MDP(28);
t161 = -t125 * t136 + t126 * t139;
t160 = 0.2e1 * t190;
t159 = t147 * MDP(18) + t148 * MDP(19);
t158 = -t190 + t191;
t157 = 0.2e1 * t191;
t156 = (t184 * MDP(16) - t150 * MDP(17)) * pkin(2);
t155 = t164 * MDP(20) - t112 * MDP(16) - t113 * MDP(17) + t126 * MDP(13) + (t97 * t123 - t124 * t96) * MDP(23) - t124 * t180 + (-MDP(14) + t175) * t125;
t137 = -pkin(3) - t183;
t129 = t148 * qJ(4) + t144;
t128 = (-pkin(8) - qJ(4)) * t147;
t127 = t139 - t183;
t120 = t148 * t136 + t144;
t119 = (-pkin(8) - t136) * t147;
t109 = t149 * t128 + t152 * t129;
t108 = t152 * t128 - t149 * t129;
t105 = t112 * t147;
t101 = t149 * t119 + t152 * t120;
t100 = t152 * t119 - t149 * t120;
t95 = pkin(4) * t178 + t112;
t93 = t95 * t124;
t92 = t95 * t123;
t1 = [MDP(1) + pkin(1) * MDP(9) * t188 + (t112 ^ 2 + t90 ^ 2 + t91 ^ 2) * MDP(21) - (t96 * t186 - t180) * t97 + (MDP(11) * t126 + MDP(17) * t189) * t126 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t151 + MDP(5) * t188) * t151 + (-0.2e1 * t126 * MDP(12) + MDP(16) * t189 + t170 - 0.2e1 * t176 - 0.2e1 * t177) * t125 + 0.2e1 * t162 * t95 + 0.2e1 * (t90 * MDP(18) - t91 * MDP(19) + t163) * t125 + 0.2e1 * ((-t147 * t91 - t148 * t90) * MDP(20) + t159 * t112) * t126; (-t153 * MDP(10) - t151 * MDP(9)) * pkin(6) + t151 * MDP(6) + t153 * MDP(7) + (t100 * t125 + t127 * t96 + t92) * MDP(27) + (-t101 * t125 - t127 * t97 + t93) * MDP(28) + (t161 * t147 - t179) * MDP(18) + (t112 * t139 + t164 * t136) * MDP(21) + (t161 * t148 + t105) * MDP(19) + t155; MDP(8) + t173 * t187 + t136 ^ 2 * t166 + (-t160 + t169) * t139 + t127 * t157 + 0.2e1 * t156 + t167; (t165 * t147 - t179) * MDP(18) + (t165 * t148 + t105) * MDP(19) + (-t112 * pkin(3) + t164 * qJ(4)) * MDP(21) + (t108 * t125 + t137 * t96 + t92) * MDP(27) + (-t109 * t125 - t137 * t97 + t93) * MDP(28) + t155; (t172 + t173) * MDP(20) + (-t139 * pkin(3) + qJ(4) * t173) * MDP(21) + t156 + t167 + t190 * (pkin(3) - t139) + t191 * (t127 + t137); t172 * t187 + qJ(4) ^ 2 * t166 + t137 * t157 + (t160 + t181) * pkin(3) + t167; t112 * MDP(21) + t159 * t126 + t162; t158 + t169; t158 - t181; MDP(21); t163 + t170 - t176 - t177; t100 * MDP(27) - t101 * MDP(28) + t175; t108 * MDP(27) - t109 * MDP(28) + t175; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
