% Calculate joint inertia matrix for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR8_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:37
% EndTime: 2019-03-09 03:59:38
% DurationCPUTime: 0.51s
% Computational Cost: add. (657->131), mult. (1154->193), div. (0->0), fcn. (1285->8), ass. (0->69)
t129 = sin(pkin(10));
t122 = pkin(3) * t129 + pkin(8);
t132 = sin(qJ(5));
t134 = cos(qJ(5));
t150 = t134 * MDP(22);
t151 = t132 * MDP(21);
t140 = t150 + t151;
t162 = pkin(9) + t122;
t108 = t162 * t132;
t109 = t162 * t134;
t131 = sin(qJ(6));
t133 = cos(qJ(6));
t115 = t131 * t132 - t133 * t134;
t116 = t131 * t134 + t132 * t133;
t144 = t116 * MDP(25) - t115 * MDP(26) + (-t108 * t133 - t109 * t131) * MDP(28) - (-t108 * t131 + t109 * t133) * MDP(29);
t173 = MDP(18) * t132 + MDP(19) * t134 - t140 * t122 + t144;
t130 = cos(pkin(10));
t135 = cos(qJ(3));
t165 = sin(qJ(3));
t112 = -t129 * t165 + t130 * t135;
t89 = t116 * t112;
t91 = t115 * t112;
t172 = -t91 * MDP(25) - t89 * MDP(26);
t171 = (MDP(18) * t134 - MDP(19) * t132) * t112;
t136 = -pkin(1) - pkin(7);
t170 = -qJ(4) + t136;
t141 = t134 * MDP(21) - t132 * MDP(22);
t104 = t115 * MDP(28);
t153 = -t116 * MDP(29) - t104;
t169 = t141 + t153;
t167 = -2 * MDP(24);
t166 = 0.2e1 * MDP(29);
t164 = (pkin(1) * MDP(6));
t113 = -t129 * t135 - t130 * t165;
t163 = pkin(5) * t113;
t161 = (MDP(28) * t116 - MDP(29) * t115) * t113;
t160 = MDP(15) * pkin(3);
t117 = t170 * t165;
t145 = t170 * t135;
t101 = t130 * t117 + t129 * t145;
t158 = t101 * t134;
t124 = t165 * pkin(3) + qJ(2);
t98 = -pkin(4) * t113 - pkin(8) * t112 + t124;
t81 = t158 + (-pkin(9) * t112 + t98) * t132;
t159 = t133 * t81;
t157 = t112 * t132;
t156 = t112 * t134;
t155 = t113 * t129;
t154 = t132 * t134;
t152 = t116 * MDP(23);
t149 = MDP(20) + MDP(27);
t148 = -t113 * MDP(27) + t172;
t123 = -pkin(3) * t130 - pkin(4);
t147 = MDP(17) * t154;
t146 = t165 * MDP(13);
t82 = -t101 * t132 + t134 * t98;
t80 = -pkin(9) * t156 - t163 + t82;
t77 = -t131 * t81 + t133 * t80;
t99 = t117 * t129 - t130 * t145;
t139 = (MDP(28) * t133 - MDP(29) * t131) * pkin(5);
t128 = t134 ^ 2;
t127 = t132 ^ 2;
t118 = -pkin(5) * t134 + t123;
t111 = t113 ^ 2;
t110 = t112 ^ 2;
t84 = pkin(5) * t157 + t99;
t83 = t132 * t98 + t158;
t78 = t131 * t80 + t159;
t1 = [(t101 ^ 2 + t124 ^ 2 + t99 ^ 2) * MDP(15) + MDP(1) - (-t91 * MDP(23) + t89 * t167) * t91 + (MDP(7) * t135 - 0.2e1 * t165 * MDP(8)) * t135 + t149 * t111 + (t128 * MDP(16) - 0.2e1 * t147) * t110 + ((-2 * MDP(4) + t164) * pkin(1)) + (0.2e1 * t165 * MDP(12) + 0.2e1 * t135 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (-t171 - t172) * t113 + 0.2e1 * (t101 * t113 + t112 * t99) * MDP(14) + 0.2e1 * (-t113 * t82 + t99 * t157) * MDP(21) + 0.2e1 * (t113 * t83 + t99 * t156) * MDP(22) + 0.2e1 * (-t113 * t77 + t84 * t89) * MDP(28) + (t113 * t78 - t84 * t91) * t166; -t111 * t151 - t164 + MDP(4) + (-t101 * MDP(15) + (-MDP(14) - t150) * t113 - t161) * t113 + (-t99 * MDP(15) - t89 * MDP(28) + t91 * MDP(29) + (-MDP(14) - t140) * t112) * t112; MDP(6) + (t111 + t110) * MDP(15); -t165 * MDP(10) - t136 * t146 - t91 * t152 + (t115 * t91 - t116 * t89) * MDP(24) + (t115 * t84 + t118 * t89) * MDP(28) + (t116 * t84 - t118 * t91) * MDP(29) - t141 * t99 + (MDP(12) * t136 + MDP(9)) * t135 + (MDP(16) * t154 + (-t127 + t128) * MDP(17) + t140 * t123) * t112 - t173 * t113 + ((-t112 * t130 + t155) * MDP(14) + (t101 * t129 - t130 * t99) * MDP(15)) * pkin(3); -t155 * t160 + t135 * MDP(12) - t146 + (t130 * t160 + t169) * t112; 0.2e1 * t147 + 0.2e1 * t118 * t104 + t127 * MDP(16) + MDP(11) + (t129 ^ 2 + t130 ^ 2) * MDP(15) * pkin(3) ^ 2 - 0.2e1 * t141 * t123 + (t115 * t167 + t118 * t166 + t152) * t116; MDP(15) * t124 - t113 * t169; 0; 0; MDP(15); -t113 * MDP(20) + t82 * MDP(21) - t83 * MDP(22) + (-t133 * t163 + t77) * MDP(28) + (-t159 + (-t80 + t163) * t131) * MDP(29) + t171 + t148; t140 * t113 + t161; t173; t169; 0.2e1 * t139 + t149; t77 * MDP(28) - t78 * MDP(29) + t148; t161; t144; t153; MDP(27) + t139; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
