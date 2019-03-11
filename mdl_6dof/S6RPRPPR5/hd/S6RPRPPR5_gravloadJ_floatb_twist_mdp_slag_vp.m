% Calculate Gravitation load on the joints for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:38
% EndTime: 2019-03-09 02:51:40
% DurationCPUTime: 0.60s
% Computational Cost: add. (240->74), mult. (278->99), div. (0->0), fcn. (241->10), ass. (0->44)
t87 = cos(pkin(9));
t83 = pkin(9) + qJ(3);
t78 = sin(t83);
t80 = cos(t83);
t97 = t80 * pkin(3) + t78 * qJ(4);
t116 = pkin(2) * t87 + pkin(1) + t97;
t90 = cos(qJ(1));
t115 = g(2) * t90;
t114 = MDP(14) - MDP(17);
t113 = MDP(13) - MDP(16) + MDP(21);
t111 = g(1) * t90;
t89 = sin(qJ(1));
t71 = g(2) * t89 + t111;
t112 = t71 * t78;
t59 = g(3) * t78 + t71 * t80;
t108 = g(3) * t80;
t88 = -pkin(7) - qJ(2);
t107 = pkin(4) - t88;
t106 = t78 * t90;
t82 = pkin(10) + qJ(6);
t79 = cos(t82);
t105 = t79 * t90;
t84 = sin(pkin(10));
t104 = t84 * t90;
t86 = cos(pkin(10));
t103 = t86 * t90;
t77 = sin(t82);
t102 = t89 * t77;
t101 = t89 * t79;
t100 = t89 * t84;
t99 = t89 * t86;
t96 = qJ(4) * t80;
t95 = qJ(5) * t80;
t94 = -MDP(18) - MDP(22);
t93 = t116 * t115;
t70 = g(1) * t89 - t115;
t68 = t90 * t96;
t66 = t89 * t96;
t63 = -t78 * t102 + t105;
t62 = t78 * t101 + t77 * t90;
t61 = t77 * t106 + t101;
t60 = t78 * t105 - t102;
t58 = -t108 + t112;
t1 = [(-g(1) * (-t89 * pkin(1) + qJ(2) * t90) - g(2) * (pkin(1) * t90 + t89 * qJ(2))) * MDP(7) + (t88 * t111 - t93 + (g(1) * t116 + g(2) * t88) * t89) * MDP(18) + (-g(1) * (-t78 * t100 + t103) - g(2) * (t78 * t104 + t99)) * MDP(19) + (-g(1) * (-t78 * t99 - t104) - g(2) * (t78 * t103 - t100)) * MDP(20) + (-t93 + (-g(1) * t107 - g(2) * t95) * t90 + (-g(1) * (-t116 - t95) - g(2) * t107) * t89) * MDP(22) + (-g(1) * t63 - g(2) * t61) * MDP(28) + (g(1) * t62 - g(2) * t60) * MDP(29) + (MDP(3) - MDP(6) - MDP(15)) * t71 + (-t114 * t78 + MDP(4) * t87 - MDP(5) * sin(pkin(9)) + MDP(2) + t113 * t80) * t70; (-MDP(7) + t94) * t70; (-g(1) * (-pkin(3) * t106 + t68) - g(2) * (-pkin(3) * t78 * t89 + t66) - g(3) * t97) * MDP(18) + (-g(1) * t68 - g(2) * t66 - g(3) * (t95 + t97) + (pkin(3) + qJ(5)) * t112) * MDP(22) + t113 * t58 + (-MDP(19) * t84 - MDP(20) * t86 - MDP(28) * t77 - MDP(29) * t79 + t114) * t59; t94 * t58; -t59 * MDP(22); (-g(1) * t60 - g(2) * t62 + t79 * t108) * MDP(28) + (g(1) * t61 - g(2) * t63 - t77 * t108) * MDP(29);];
taug  = t1;
