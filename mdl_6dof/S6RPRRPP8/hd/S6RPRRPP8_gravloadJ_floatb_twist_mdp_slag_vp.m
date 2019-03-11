% Calculate Gravitation load on the joints for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:55:53
% EndTime: 2019-03-09 04:55:55
% DurationCPUTime: 0.61s
% Computational Cost: add. (216->94), mult. (488->129), div. (0->0), fcn. (477->6), ass. (0->47)
t158 = MDP(20) - MDP(23) - MDP(26);
t157 = -MDP(19) + MDP(22) - MDP(27);
t120 = cos(qJ(3));
t121 = cos(qJ(1));
t139 = t120 * t121;
t117 = sin(qJ(3));
t148 = g(3) * t117;
t156 = -g(2) * t139 - t148;
t153 = MDP(13) - MDP(21) - MDP(25);
t152 = -pkin(1) - pkin(7);
t151 = -pkin(5) - pkin(8);
t118 = sin(qJ(1));
t150 = g(1) * t118;
t149 = g(2) * t121;
t147 = -pkin(4) - qJ(6);
t116 = sin(qJ(4));
t146 = t116 * t120;
t145 = t117 * t118;
t144 = t117 * t121;
t143 = t118 * t116;
t119 = cos(qJ(4));
t142 = t118 * t119;
t141 = t118 * t120;
t140 = t119 * t120;
t138 = t121 * t119;
t137 = t121 * pkin(1) + t118 * qJ(2);
t136 = MDP(24) + MDP(28);
t134 = t116 * t141;
t133 = t118 * t140;
t132 = -qJ(5) * t116 - pkin(3);
t89 = t117 * t143 - t138;
t90 = t116 * t121 + t117 * t142;
t131 = -t89 * pkin(4) + qJ(5) * t90;
t91 = t116 * t144 + t142;
t92 = t117 * t138 - t143;
t130 = t91 * pkin(4) - qJ(5) * t92;
t129 = pkin(3) * t141 + pkin(4) * t133 + pkin(8) * t145 + qJ(5) * t134;
t98 = -t149 + t150;
t126 = -pkin(4) * t119 + t132;
t125 = pkin(3) * t145 + t90 * pkin(4) + t121 * pkin(7) + qJ(5) * t89 + t137;
t124 = g(1) * t89 - g(2) * t91 + g(3) * t146;
t123 = g(1) * t90 - g(2) * t92 + g(3) * t140;
t111 = t121 * qJ(2);
t122 = pkin(3) * t144 + t92 * pkin(4) - pkin(8) * t139 + t91 * qJ(5) + t111;
t112 = t120 * pkin(8);
t100 = qJ(5) * t140;
t1 = [(-g(1) * (-t118 * pkin(1) + t111) - g(2) * t137) * MDP(6) + (-g(1) * (t152 * t118 + t122) - g(2) * (-pkin(8) * t141 + t125)) * MDP(24) + (-g(1) * (-pkin(5) * t139 + t92 * qJ(6) + t122) - g(2) * (qJ(6) * t90 + t125) + (-g(2) * t151 * t120 - g(1) * t152) * t118) * MDP(28) + (MDP(2) - MDP(4)) * t98 + (-t117 * MDP(12) - t153 * t120 + MDP(3) - MDP(5)) * (g(1) * t121 + g(2) * t118) + t158 * (g(1) * t91 + g(2) * t89) + t157 * (g(1) * t92 + g(2) * t90); (-MDP(6) - t136) * t98; (-t98 * t120 + t148) * MDP(12) + (-g(1) * t129 - g(3) * t112 - t126 * t148 - (-pkin(8) * t117 + t126 * t120) * t149) * MDP(24) + (-g(1) * (qJ(6) * t133 + t129) - g(3) * (pkin(5) * t120 + t112) + (-pkin(5) * t150 - g(3) * (-qJ(6) * t119 + t126)) * t117 - (t151 * t117 + (t147 * t119 + t132) * t120) * t149) * MDP(28) + t153 * (g(1) * t145 - g(2) * t144 + g(3) * t120) + t157 * (g(1) * t133 + t156 * t119) + t158 * (g(1) * t134 + t156 * t116); (-g(1) * t131 - g(2) * t130 - g(3) * (-pkin(4) * t146 + t100)) * MDP(24) + (-g(1) * (-qJ(6) * t89 + t131) - g(2) * (qJ(6) * t91 + t130) - g(3) * (t147 * t146 + t100)) * MDP(28) - t157 * t124 + t158 * t123; -t136 * t124; -t123 * MDP(28);];
taug  = t1;
