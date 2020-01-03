% Calculate Gravitation load on the joints for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:00
% EndTime: 2019-12-31 21:02:02
% DurationCPUTime: 0.37s
% Computational Cost: add. (217->82), mult. (349->118), div. (0->0), fcn. (328->8), ass. (0->41)
t89 = cos(qJ(3));
t77 = t89 * pkin(3) + pkin(2);
t90 = cos(qJ(2));
t70 = t90 * t77;
t85 = -qJ(4) - pkin(7);
t87 = sin(qJ(2));
t115 = -t87 * t85 + t70;
t114 = -pkin(1) - t115;
t113 = MDP(10) - MDP(18) - MDP(21);
t88 = sin(qJ(1));
t91 = cos(qJ(1));
t99 = g(1) * t91 + g(2) * t88;
t62 = -g(3) * t90 + t99 * t87;
t110 = g(3) * t87;
t86 = sin(qJ(3));
t107 = t88 * t86;
t106 = t88 * t89;
t105 = t88 * t90;
t84 = qJ(3) + pkin(8);
t78 = sin(t84);
t104 = t91 * t78;
t79 = cos(t84);
t103 = t91 * t79;
t102 = t91 * t86;
t101 = t91 * t89;
t100 = t90 * t102;
t97 = pkin(4) * t79 + qJ(5) * t78;
t64 = t105 * t86 + t101;
t94 = pkin(3) * t107 + t88 * pkin(6) - t114 * t91;
t58 = t105 * t78 + t103;
t60 = t104 * t90 - t88 * t79;
t93 = g(1) * t60 + g(2) * t58 + t110 * t78;
t92 = pkin(3) * t102 + t91 * pkin(6) + t114 * t88;
t63 = t99 * t90 + t110;
t75 = pkin(3) * t106;
t67 = t101 * t90 + t107;
t66 = -t100 + t106;
t65 = -t105 * t89 + t102;
t61 = t103 * t90 + t88 * t78;
t59 = t105 * t79 - t104;
t1 = [t99 * MDP(3) + (-g(1) * t65 - g(2) * t67) * MDP(16) + (-g(1) * t64 - g(2) * t66) * MDP(17) + (-g(1) * t92 - g(2) * t94) * MDP(19) + (g(1) * t59 - g(2) * t61) * MDP(20) + (g(1) * t58 - g(2) * t60) * MDP(22) + (-g(1) * (-t59 * pkin(4) - t58 * qJ(5) + t92) - g(2) * (t61 * pkin(4) + t60 * qJ(5) + t94)) * MDP(23) + (t90 * MDP(9) - t113 * t87 + MDP(2)) * (g(1) * t88 - g(2) * t91); (-g(3) * t115 + t99 * (t77 * t87 + t85 * t90)) * MDP(19) + (-g(3) * t70 + (-g(3) * t97 + t85 * t99) * t90 + (g(3) * t85 + t99 * (t77 + t97)) * t87) * MDP(23) + t113 * t63 + (t89 * MDP(16) - t86 * MDP(17) + t79 * MDP(20) + t78 * MDP(22) + MDP(9)) * t62; (-g(1) * t66 + g(2) * t64 + t110 * t86) * MDP(16) + (g(1) * t67 - g(2) * t65 + t110 * t89) * MDP(17) + (-g(1) * t75 + (g(2) * t101 + t63 * t86) * pkin(3)) * MDP(19) + t93 * MDP(20) + (-g(1) * t61 - g(2) * t59 - t110 * t79) * MDP(22) + (-g(1) * (-pkin(3) * t100 - t60 * pkin(4) + t61 * qJ(5) + t75) - g(2) * (-pkin(3) * t64 - t58 * pkin(4) + t59 * qJ(5)) - (-pkin(3) * t86 - pkin(4) * t78 + qJ(5) * t79) * t110) * MDP(23); (-MDP(19) - MDP(23)) * t62; -t93 * MDP(23);];
taug = t1;
