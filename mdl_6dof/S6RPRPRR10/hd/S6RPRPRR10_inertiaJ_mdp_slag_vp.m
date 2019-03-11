% Calculate joint inertia matrix for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRPRR10_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:09:22
% EndTime: 2019-03-09 04:09:24
% DurationCPUTime: 0.65s
% Computational Cost: add. (746->166), mult. (1450->240), div. (0->0), fcn. (1543->8), ass. (0->83)
t192 = MDP(17) * qJ(4);
t143 = sin(pkin(10));
t144 = cos(pkin(10));
t146 = sin(qJ(5));
t149 = cos(qJ(5));
t130 = t143 * t149 + t144 * t146;
t150 = cos(qJ(3));
t121 = t130 * t150;
t129 = t143 * t146 - t149 * t144;
t123 = t129 * t150;
t145 = sin(qJ(6));
t148 = cos(qJ(6));
t101 = t148 * t121 - t123 * t145;
t103 = -t121 * t145 - t123 * t148;
t191 = t103 * MDP(27) - t101 * MDP(28);
t174 = t143 ^ 2 + t144 ^ 2;
t190 = t174 * t192;
t181 = pkin(8) + qJ(4);
t132 = t181 * t143;
t133 = t181 * t144;
t114 = -t149 * t132 - t133 * t146;
t115 = -t132 * t146 + t133 * t149;
t104 = -pkin(9) * t130 + t114;
t105 = -pkin(9) * t129 + t115;
t110 = t148 * t129 + t130 * t145;
t111 = -t129 * t145 + t130 * t148;
t160 = t111 * MDP(27) - t110 * MDP(28) + (t104 * t148 - t105 * t145) * MDP(30) - (t104 * t145 + t105 * t148) * MDP(31);
t189 = MDP(20) * t130 - MDP(21) * t129 + MDP(23) * t114 - MDP(24) * t115 + t160;
t188 = 0.2e1 * t150;
t187 = -2 * MDP(19);
t186 = 0.2e1 * MDP(24);
t185 = -2 * MDP(26);
t184 = 0.2e1 * MDP(31);
t183 = (pkin(1) * MDP(6));
t147 = sin(qJ(3));
t182 = pkin(5) * t147;
t120 = t130 * t147;
t122 = t129 * t147;
t100 = -t120 * t148 + t122 * t145;
t102 = -t120 * t145 - t122 * t148;
t180 = t100 * MDP(30) - t102 * MDP(31);
t179 = pkin(3) * MDP(17);
t131 = pkin(3) * t147 - qJ(4) * t150 + qJ(2);
t126 = t144 * t131;
t151 = -pkin(1) - pkin(7);
t176 = t143 * t151;
t109 = -pkin(8) * t144 * t150 + t126 + (pkin(4) - t176) * t147;
t175 = t147 * t151;
t117 = t143 * t131 + t144 * t175;
t177 = t143 * t150;
t113 = -pkin(8) * t177 + t117;
t95 = t109 * t146 + t113 * t149;
t91 = -pkin(9) * t121 + t95;
t178 = t148 * t91;
t141 = t147 ^ 2;
t142 = t150 ^ 2;
t173 = -t141 - t142;
t172 = MDP(17) * t147;
t171 = t110 * MDP(30);
t170 = t111 * MDP(25);
t169 = t129 * MDP(23);
t168 = t130 * MDP(18);
t167 = t143 * MDP(15);
t166 = t144 * MDP(14);
t165 = t151 * MDP(17);
t164 = MDP(22) + MDP(29);
t163 = t147 * MDP(29) + t191;
t137 = -pkin(4) * t144 - pkin(3);
t94 = t149 * t109 - t113 * t146;
t88 = pkin(9) * t123 + t182 + t94;
t85 = -t145 * t91 + t148 * t88;
t127 = pkin(4) * t177 - t150 * t151;
t161 = t174 * MDP(16);
t158 = -t143 * MDP(14) - t144 * MDP(15);
t157 = -t123 * MDP(20) - t121 * MDP(21);
t156 = (MDP(30) * t148 - MDP(31) * t145) * pkin(5);
t155 = t166 - t167 + t179;
t153 = t130 * MDP(24) + t111 * MDP(31) - t155 + t169 + t171;
t119 = pkin(5) * t129 + t137;
t116 = -t143 * t175 + t126;
t112 = pkin(5) * t121 + t127;
t86 = t145 * t88 + t178;
t1 = [(t142 * t151 ^ 2 + t116 ^ 2 + t117 ^ 2) * MDP(17) + t142 * MDP(7) + MDP(1) + t164 * t141 - (-t123 * MDP(18) + t121 * t187) * t123 + (t103 * MDP(25) + t101 * t185) * t103 + ((-2 * MDP(4) + t183) * pkin(1)) + (MDP(13) * t188 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (qJ(2) * MDP(12) - t150 * MDP(8) + t157 + t191) * t147 + 0.2e1 * (t101 * t112 + t147 * t85) * MDP(30) + (t103 * t112 - t147 * t86) * t184 + 0.2e1 * (t121 * t127 + t147 * t94) * MDP(23) + (-t123 * t127 - t147 * t95) * t186 + 0.2e1 * (t116 * t147 - t142 * t176) * MDP(14) + 0.2e1 * (-t142 * t144 * t151 - t117 * t147) * MDP(15) + (-t116 * t144 - t117 * t143) * MDP(16) * t188; MDP(4) - t183 + t142 * t165 + (-t120 * t147 - t121 * t150) * MDP(23) + (t122 * t147 + t123 * t150) * MDP(24) + (t100 * t147 - t101 * t150) * MDP(30) + (-t102 * t147 - t103 * t150) * MDP(31) + (t173 * MDP(15) + t117 * t172) * t144 + (t173 * MDP(14) - t116 * t172) * t143; MDP(6) + (t174 * t141 + t142) * MDP(17); -t123 * t168 + (-t121 * t130 + t123 * t129) * MDP(19) + (t121 * t137 + t127 * t129) * MDP(23) + (-t123 * t137 + t127 * t130) * MDP(24) + t103 * t170 + (-t101 * t111 - t103 * t110) * MDP(26) + (t101 * t119 + t110 * t112) * MDP(30) + (t103 * t119 + t111 * t112) * MDP(31) + (MDP(9) + t158 * pkin(3) + (MDP(12) + t155) * t151) * t150 + (-MDP(13) * t151 + t158 * qJ(4) - MDP(10) + t189) * t147 + (MDP(16) + t192) * (-t116 * t143 + t117 * t144); (-MDP(13) + t161 + t190) * t147 + (MDP(12) - t153) * t150; 0.2e1 * t137 * t169 + 0.2e1 * t119 * t171 + MDP(11) + (0.2e1 * t166 - 0.2e1 * t167 + t179) * pkin(3) + (t129 * t187 + t137 * t186 + t168) * t130 + (t110 * t185 + t119 * t184 + t170) * t111 + (0.2e1 * t161 + t190) * qJ(4); t121 * MDP(23) - t123 * MDP(24) + t101 * MDP(30) + t103 * MDP(31) + (-t158 - t165) * t150; -t150 * MDP(17); t153; MDP(17); t147 * MDP(22) + t94 * MDP(23) - t95 * MDP(24) + (t148 * t182 + t85) * MDP(30) + (-t178 + (-t88 - t182) * t145) * MDP(31) + t157 + t163; -t120 * MDP(23) + t122 * MDP(24) + t180; t189; 0; 0.2e1 * t156 + t164; t85 * MDP(30) - t86 * MDP(31) + t163; t180; t160; 0; MDP(29) + t156; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
