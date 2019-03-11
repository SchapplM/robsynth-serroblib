% Calculate Gravitation load on the joints for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:58:14
% EndTime: 2019-03-10 00:58:15
% DurationCPUTime: 0.26s
% Computational Cost: add. (370->62), mult. (336->82), div. (0->0), fcn. (290->10), ass. (0->36)
t118 = MDP(24) - MDP(32);
t89 = cos(qJ(2));
t83 = qJ(2) + qJ(3);
t80 = qJ(4) + t83;
t75 = sin(t80);
t76 = cos(t80);
t88 = cos(qJ(5));
t77 = t88 * pkin(5) + pkin(4);
t84 = -qJ(6) - pkin(10);
t100 = -t75 * t84 + t76 * t77;
t79 = cos(t83);
t99 = pkin(3) * t79 + t100;
t117 = t89 * pkin(2) + t99;
t78 = sin(t83);
t95 = t75 * t77 + t76 * t84;
t116 = pkin(3) * t78 + t95;
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t97 = g(1) * t90 + g(2) * t87;
t62 = -g(3) * t76 + t97 * t75;
t109 = g(3) * t75;
t85 = sin(qJ(5));
t106 = t87 * t85;
t105 = t87 * t88;
t104 = t90 * t85;
t103 = t90 * t88;
t101 = pkin(5) * t85 + pkin(7) + pkin(8) + pkin(9);
t98 = t118 * (t97 * t76 + t109) + (MDP(30) * t88 - MDP(31) * t85 + MDP(23)) * t62;
t69 = -t76 * t104 + t105;
t67 = t76 * t106 + t103;
t94 = -pkin(1) - t117;
t92 = (-g(3) * t79 + t97 * t78) * MDP(16) + (g(3) * t78 + t97 * t79) * MDP(17) + t98;
t86 = sin(qJ(2));
t70 = t76 * t103 + t106;
t68 = -t76 * t105 + t104;
t1 = [t97 * MDP(3) + (-g(1) * t68 - g(2) * t70) * MDP(30) + (-g(1) * t67 - g(2) * t69) * MDP(31) + ((-g(1) * t101 + g(2) * t94) * t90 + (-g(1) * t94 - g(2) * t101) * t87) * MDP(33) + (-t86 * MDP(10) + MDP(16) * t79 - MDP(17) * t78 + t76 * MDP(23) + t89 * MDP(9) - t118 * t75 + MDP(2)) * (g(1) * t87 - g(2) * t90); (-g(3) * t89 + t97 * t86) * MDP(9) + (g(3) * t86 + t97 * t89) * MDP(10) + (-g(3) * t117 + t97 * (t86 * pkin(2) + t116)) * MDP(33) + t92; (-g(3) * t99 + t97 * t116) * MDP(33) + t92; (-g(3) * t100 + t97 * t95) * MDP(33) + t98; (g(1) * t70 - g(2) * t68 + t88 * t109) * MDP(31) + (pkin(5) * MDP(33) + MDP(30)) * (-g(1) * t69 + g(2) * t67 + t85 * t109); -t62 * MDP(33);];
taug  = t1;
