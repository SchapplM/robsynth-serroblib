% Calculate joint inertia matrix for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR9_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:54
% EndTime: 2019-12-31 22:29:56
% DurationCPUTime: 0.57s
% Computational Cost: add. (727->146), mult. (1463->206), div. (0->0), fcn. (1529->8), ass. (0->78)
t143 = sin(qJ(3));
t147 = cos(qJ(3));
t180 = pkin(7) + pkin(8);
t129 = t180 * t143;
t130 = t180 * t147;
t142 = sin(qJ(4));
t146 = cos(qJ(4));
t108 = -t146 * t129 - t142 * t130;
t109 = -t142 * t129 + t146 * t130;
t124 = t142 * t143 - t146 * t147;
t125 = t142 * t147 + t146 * t143;
t141 = sin(qJ(5));
t145 = cos(qJ(5));
t101 = t145 * t124 + t141 * t125;
t102 = -t141 * t124 + t145 * t125;
t96 = -t125 * pkin(9) + t108;
t97 = -t124 * pkin(9) + t109;
t158 = t102 * MDP(27) - t101 * MDP(28) + (-t141 * t97 + t145 * t96) * MDP(30) - (t141 * t96 + t145 * t97) * MDP(31);
t150 = t125 * MDP(20) - t124 * MDP(21) + t108 * MDP(23) - t109 * MDP(24) + t158;
t193 = t150 - (t143 * MDP(16) + t147 * MDP(17)) * pkin(7) + t143 * MDP(13) + t147 * MDP(14);
t163 = t146 * MDP(23);
t192 = pkin(3) * (-t142 * MDP(24) + t163);
t164 = t141 * MDP(31);
t151 = (t145 * MDP(30) - t164) * pkin(4);
t144 = sin(qJ(2));
t116 = t125 * t144;
t117 = t124 * t144;
t94 = t145 * t116 - t141 * t117;
t95 = -t141 * t116 - t145 * t117;
t174 = t95 * MDP(27) - t94 * MDP(28);
t191 = -t117 * MDP(20) - t116 * MDP(21);
t148 = cos(qJ(2));
t128 = -t148 * pkin(2) - t144 * pkin(7) - pkin(1);
t123 = t147 * t128;
t175 = pkin(8) * t144;
t178 = pkin(6) * t143;
t103 = -t147 * t175 + t123 + (-pkin(3) - t178) * t148;
t176 = pkin(6) * t148;
t161 = t147 * t176;
t105 = t161 + (t128 - t175) * t143;
t90 = t146 * t103 - t142 * t105;
t84 = -t148 * pkin(4) + t117 * pkin(9) + t90;
t91 = t142 * t103 + t146 * t105;
t89 = -t116 * pkin(9) + t91;
t81 = -t141 * t89 + t145 * t84;
t82 = t141 * t84 + t145 * t89;
t189 = t81 * MDP(30) - t82 * MDP(31) + t174;
t190 = t90 * MDP(23) - t91 * MDP(24) + t189 + t191;
t186 = -2 * MDP(19);
t185 = 0.2e1 * MDP(23);
t184 = 0.2e1 * MDP(24);
t183 = -2 * MDP(26);
t182 = 0.2e1 * MDP(30);
t181 = 0.2e1 * MDP(31);
t179 = pkin(3) * t142;
t177 = pkin(6) * t147;
t173 = t143 * t147;
t127 = (pkin(3) * t143 + pkin(6)) * t144;
t168 = t102 * MDP(25);
t134 = t146 * pkin(3) + pkin(4);
t131 = t145 * t134;
t167 = (-t141 * t179 + t131) * MDP(30);
t166 = (t141 * t134 + t145 * t179) * MDP(31);
t165 = t125 * MDP(18);
t162 = MDP(22) + MDP(29);
t135 = -t147 * pkin(3) - pkin(2);
t160 = MDP(12) * t173;
t159 = MDP(15) + t162;
t157 = t147 * MDP(13) - t143 * MDP(14);
t152 = MDP(29) - t166 + t167;
t139 = t147 ^ 2;
t138 = t144 ^ 2;
t137 = t143 ^ 2;
t113 = t124 * pkin(4) + t135;
t112 = t143 * t128 + t161;
t111 = -t143 * t176 + t123;
t104 = t116 * pkin(4) + t127;
t1 = [-0.2e1 * pkin(1) * t144 * MDP(10) + MDP(1) + (t95 * MDP(25) + t94 * t183) * t95 - (-t117 * MDP(18) + t116 * t186) * t117 + t159 * t148 ^ 2 + (t139 * MDP(11) + MDP(4) - 0.2e1 * t160) * t138 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t157) * t144 - t191 - t174) * t148 + 0.2e1 * (t112 * t148 + t138 * t177) * MDP(17) + 0.2e1 * (-t111 * t148 + t138 * t178) * MDP(16) + (t127 * t116 - t90 * t148) * t185 + (-t127 * t117 + t91 * t148) * t184 + (t104 * t95 + t82 * t148) * t181 + (t104 * t94 - t81 * t148) * t182; -t117 * t165 + (-t125 * t116 + t117 * t124) * MDP(19) + (t135 * t116 + t127 * t124) * MDP(23) + (-t135 * t117 + t127 * t125) * MDP(24) + t95 * t168 + (-t95 * t101 - t102 * t94) * MDP(26) + (t104 * t101 + t113 * t94) * MDP(30) + (t104 * t102 + t113 * t95) * MDP(31) + (MDP(6) - pkin(6) * MDP(9) + MDP(11) * t173 + (-t137 + t139) * MDP(12) + (-pkin(2) * t143 - t177) * MDP(16) + (-pkin(2) * t147 + t178) * MDP(17)) * t144 + (-pkin(6) * MDP(10) + MDP(7) - t193) * t148; 0.2e1 * t160 + t135 * t124 * t185 + t113 * t101 * t182 + t137 * MDP(11) + MDP(8) + 0.2e1 * (t147 * MDP(16) - t143 * MDP(17)) * pkin(2) + (t124 * t186 + t135 * t184 + t165) * t125 + (t101 * t183 + t113 * t181 + t168) * t102; t111 * MDP(16) - t112 * MDP(17) + t157 * t144 + (-MDP(15) - MDP(22) - t152 - t192) * t148 + t190; t193; t159 - 0.2e1 * t166 + 0.2e1 * t167 + 0.2e1 * t192; (-t162 - t151) * t148 + t190; t150; (t145 * pkin(4) + t131) * MDP(30) + (-pkin(4) - t134) * t164 + (t163 + (-MDP(30) * t141 - MDP(31) * t145 - MDP(24)) * t142) * pkin(3) + t162; 0.2e1 * t151 + t162; -t148 * MDP(29) + t189; t158; t152; MDP(29) + t151; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
