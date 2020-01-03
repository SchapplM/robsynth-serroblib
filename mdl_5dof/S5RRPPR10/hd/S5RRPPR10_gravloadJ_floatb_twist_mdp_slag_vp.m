% Calculate Gravitation load on the joints for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:41
% EndTime: 2019-12-31 19:44:43
% DurationCPUTime: 0.43s
% Computational Cost: add. (148->66), mult. (361->99), div. (0->0), fcn. (366->8), ass. (0->40)
t85 = sin(qJ(2));
t88 = cos(qJ(2));
t102 = t88 * pkin(2) + t85 * qJ(3);
t117 = -pkin(1) - t102;
t86 = sin(qJ(1));
t89 = cos(qJ(1));
t98 = g(1) * t89 + g(2) * t86;
t113 = t85 * t98;
t63 = -g(3) * t88 + t113;
t116 = MDP(11) + MDP(15);
t115 = -MDP(12) + MDP(17);
t114 = MDP(10) - MDP(13) - MDP(16);
t111 = pkin(2) * t85;
t110 = g(1) * t86;
t107 = g(3) * t85;
t105 = t86 * t88;
t82 = sin(pkin(8));
t104 = t89 * t82;
t83 = cos(pkin(8));
t103 = t89 * t83;
t101 = qJ(3) * t88;
t100 = t86 * pkin(6) - t117 * t89;
t96 = pkin(3) * t83 + qJ(4) * t82;
t65 = t105 * t82 + t103;
t66 = t105 * t83 - t104;
t84 = sin(qJ(5));
t87 = cos(qJ(5));
t95 = t65 * t87 - t66 * t84;
t94 = t65 * t84 + t66 * t87;
t93 = t82 * t87 - t83 * t84;
t92 = t82 * t84 + t83 * t87;
t90 = t117 * t110;
t79 = t89 * pkin(6);
t72 = t89 * t101;
t70 = t86 * t101;
t68 = t103 * t88 + t86 * t82;
t67 = t104 * t88 - t86 * t83;
t59 = t67 * t84 + t68 * t87;
t58 = t67 * t87 - t68 * t84;
t1 = [t98 * MDP(3) + (-g(1) * t79 - g(2) * t100 - t90) * MDP(14) + (-g(1) * (-t66 * pkin(3) - t65 * qJ(4) + t79) - g(2) * (t68 * pkin(3) + t67 * qJ(4) + t100) - t90) * MDP(18) + (g(1) * t94 - g(2) * t59) * MDP(24) + (g(1) * t95 - g(2) * t58) * MDP(25) + t115 * (g(1) * t65 - g(2) * t67) + t116 * (g(1) * t66 - g(2) * t68) + (MDP(9) * t88 - t114 * t85 + MDP(2)) * (-g(2) * t89 + t110); (-g(1) * (-t111 * t89 + t72) - g(2) * (-t111 * t86 + t70) - g(3) * t102) * MDP(14) + (-g(1) * t72 - g(2) * t70 - g(3) * (t88 * t96 + t102) + (pkin(2) + t96) * t113) * MDP(18) + t114 * (t88 * t98 + t107) + (t92 * MDP(24) + t93 * MDP(25) + t115 * t82 + t116 * t83 + MDP(9)) * t63; (-MDP(14) - MDP(18)) * t63; (-g(1) * t67 - g(2) * t65 - t107 * t82) * MDP(18); (-g(1) * t58 - g(2) * t95) * MDP(24) + (g(1) * t59 + g(2) * t94) * MDP(25) + (-MDP(24) * t93 + MDP(25) * t92) * t107;];
taug = t1;
