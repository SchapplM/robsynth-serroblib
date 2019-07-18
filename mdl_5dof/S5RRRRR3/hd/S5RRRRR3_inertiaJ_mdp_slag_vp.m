% Calculate joint inertia matrix for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR3_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:18:44
% EndTime: 2019-07-18 17:18:46
% DurationCPUTime: 0.41s
% Computational Cost: add. (437->107), mult. (913->144), div. (0->0), fcn. (1015->8), ass. (0->64)
t132 = sin(qJ(4));
t136 = cos(qJ(4));
t135 = cos(qJ(5));
t131 = sin(qJ(5));
t168 = t131 * t132;
t115 = -t135 * t136 + t168;
t166 = t135 * t132;
t117 = t131 * t136 + t166;
t163 = t117 * MDP(27) - t115 * MDP(28);
t148 = t132 * MDP(20) + t136 * MDP(21) + t163;
t133 = sin(qJ(3));
t123 = t133 * pkin(1) + pkin(5);
t98 = t117 * t123;
t99 = t115 * t123;
t180 = -t98 * MDP(30) + t99 * MDP(31);
t111 = t117 * pkin(5);
t112 = t115 * pkin(5);
t179 = -t111 * MDP(30) + t112 * MDP(31);
t155 = t136 * MDP(23);
t143 = -t132 * MDP(24) + t155;
t178 = t115 * MDP(30) + t117 * MDP(31);
t177 = MDP(23) * t132 + MDP(24) * t136;
t176 = -2 * MDP(26);
t175 = cos(qJ(2));
t174 = cos(qJ(3));
t173 = pkin(3) * t132;
t134 = sin(qJ(2));
t116 = t133 * t134 - t174 * t175;
t172 = t116 * pkin(3);
t171 = t136 * pkin(3);
t118 = t133 * t175 + t174 * t134;
t88 = t115 * t118;
t169 = MDP(25) * t88;
t167 = t132 * t136;
t150 = t175 * pkin(1);
t90 = t116 * pkin(2) - t118 * pkin(5) - t150;
t84 = t136 * t90 + t172;
t82 = t135 * t84;
t165 = (-t90 * t168 + t82) * MDP(30);
t164 = (t131 * t84 + t90 * t166) * MDP(31);
t87 = t117 * t118;
t85 = t87 * MDP(28);
t86 = t88 * MDP(27);
t159 = MDP(31) * t131;
t154 = MDP(22) + MDP(29);
t153 = 0.2e1 * t118;
t152 = t118 * t173;
t151 = t116 * MDP(29) - t85 - t86;
t149 = MDP(19) * t167;
t124 = -t174 * pkin(1) - pkin(2);
t129 = t132 ^ 2;
t146 = t129 * MDP(18) + MDP(15) + 0.2e1 * t149 + (MDP(25) * t117 + t115 * t176) * t117;
t144 = t136 * MDP(20) - t132 * MDP(21);
t141 = (MDP(30) * t135 - t159) * pkin(3);
t140 = 0.2e1 * t143;
t139 = 0.2e1 * t178;
t138 = (t174 * MDP(16) - t133 * MDP(17)) * pkin(1);
t130 = t136 ^ 2;
t137 = (t88 * t115 - t117 * t87) * MDP(26) - t117 * t169 + (MDP(18) * t167 + MDP(13) + (-t129 + t130) * MDP(19)) * t118 + (-MDP(14) + t148) * t116;
t125 = -pkin(2) - t171;
t119 = t124 - t171;
t92 = t117 * t152;
t91 = t115 * t152;
t1 = [MDP(1) - (t87 * t176 - t169) * t88 + (MDP(4) * t134 + 0.2e1 * t175 * MDP(5)) * t134 + (-MDP(17) * t150 + (t87 * MDP(30) - t88 * MDP(31)) * t173) * t153 + (t130 * MDP(18) + MDP(11) - 0.2e1 * t149) * t118 ^ 2 + (-0.2e1 * MDP(16) * t150 - 0.2e1 * t86 - 0.2e1 * t85 + t90 * t140 + 0.2e1 * t165 - 0.2e1 * t164 + (-MDP(12) + t144) * t153 + t154 * t116) * t116; t137 + t134 * MDP(6) + (-t98 * t116 + t119 * t87 + t91) * MDP(30) + (t99 * t116 - t119 * t88 + t92) * MDP(31) + t175 * MDP(7) + t177 * (-t116 * t123 + t118 * t124); t119 * t139 - 0.2e1 * t124 * t143 + MDP(8) + 0.2e1 * t138 + t146; t137 + (-t111 * t116 + t125 * t87 + t91) * MDP(30) + (t112 * t116 - t125 * t88 + t92) * MDP(31) + t177 * (-pkin(2) * t118 - pkin(5) * t116); t138 + t146 + t143 * (pkin(2) - t124) + t178 * (t119 + t125); pkin(2) * t140 + t125 * t139 + t146; t116 * MDP(22) + (t135 * t172 + t82) * MDP(30) + (-t84 - t172) * t159 + (t155 + (-MDP(30) * t131 - MDP(31) * t135 - MDP(24)) * t132) * t90 + t144 * t118 + t151; -t177 * t123 + t148 + t180; -t177 * pkin(5) + t148 + t179; 0.2e1 * t141 + t154; t151 - t164 + t165; t163 + t180; t163 + t179; MDP(29) + t141; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq  = res;
