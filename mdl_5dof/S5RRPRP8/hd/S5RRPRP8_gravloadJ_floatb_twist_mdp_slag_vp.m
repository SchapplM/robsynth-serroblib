% Calculate Gravitation load on the joints for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:27
% EndTime: 2019-12-31 20:04:28
% DurationCPUTime: 0.26s
% Computational Cost: add. (116->61), mult. (257->87), div. (0->0), fcn. (234->6), ass. (0->38)
t77 = sin(qJ(2));
t69 = t77 * qJ(3);
t80 = cos(qJ(2));
t90 = t80 * pkin(2) + t69;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t62 = g(1) * t81 + g(2) * t78;
t103 = MDP(9) + MDP(11);
t102 = MDP(10) - MDP(13);
t79 = cos(qJ(4));
t76 = sin(qJ(4));
t92 = t80 * t76;
t85 = -t77 * t79 + t92;
t53 = t85 * t78;
t95 = t77 * t76;
t59 = t80 * t79 + t95;
t54 = t59 * t78;
t91 = t80 * t81;
t87 = t76 * t91;
t94 = t77 * t81;
t55 = -t79 * t94 + t87;
t56 = t59 * t81;
t97 = g(3) * t59;
t101 = (g(1) * t55 + g(2) * t53 + t97) * MDP(20) + (g(1) * t56 + g(2) * t54 - g(3) * t85) * MDP(21);
t100 = g(1) * t78;
t68 = t79 * pkin(4) + pkin(3);
t93 = t80 * t68;
t89 = qJ(3) * t80;
t88 = pkin(4) * t95;
t86 = pkin(2) * t91 + t78 * pkin(6) + (pkin(1) + t69) * t81;
t61 = -g(2) * t81 + t100;
t84 = -pkin(1) - t90;
t75 = -qJ(5) - pkin(7);
t72 = t81 * pkin(6);
t66 = t81 * t89;
t64 = t78 * t89;
t51 = -g(3) * t80 + t62 * t77;
t1 = [(-g(1) * t72 - g(2) * t86 - t84 * t100) * MDP(14) + (g(1) * t54 - g(2) * t56) * MDP(20) + (-g(1) * t53 + g(2) * t55) * MDP(21) + (-g(1) * (t81 * t75 + t72) - g(2) * (t68 * t91 + t81 * t88 + t86) + (-g(1) * (t84 - t88 - t93) - g(2) * t75) * t78) * MDP(23) + (MDP(3) - MDP(12) + MDP(22)) * t62 + (-t102 * t77 + t103 * t80 + MDP(2)) * t61; (-g(1) * (-pkin(2) * t94 + t66) - g(2) * (-t78 * t77 * pkin(2) + t64) - g(3) * t90) * MDP(14) + (-g(1) * (pkin(4) * t87 + t66) - g(2) * (t78 * pkin(4) * t92 + t64) - g(3) * (t90 + t93) + (-g(3) * pkin(4) * t76 + t62 * (pkin(2) + t68)) * t77) * MDP(23) + t102 * (g(3) * t77 + t62 * t80) + t103 * t51 - t101; (-MDP(14) - MDP(23)) * t51; (t62 * t85 + t97) * pkin(4) * MDP(23) + t101; t61 * MDP(23);];
taug = t1;
