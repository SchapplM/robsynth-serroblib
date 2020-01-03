% Calculate Gravitation load on the joints for
% S5RRPRP11
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:05
% EndTime: 2019-12-31 20:14:06
% DurationCPUTime: 0.40s
% Computational Cost: add. (144->69), mult. (328->95), div. (0->0), fcn. (308->6), ass. (0->36)
t77 = sin(qJ(2));
t70 = t77 * qJ(3);
t80 = cos(qJ(2));
t106 = t80 * pkin(2) + t70;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t64 = g(1) * t81 + g(2) * t78;
t105 = MDP(10) - MDP(13);
t104 = MDP(20) + MDP(22);
t103 = MDP(21) - MDP(24);
t102 = MDP(9) - MDP(12) + MDP(23);
t100 = pkin(2) * t77;
t99 = g(1) * t78;
t95 = g(3) * t80;
t76 = sin(qJ(4));
t94 = t78 * t76;
t79 = cos(qJ(4));
t93 = t78 * t79;
t92 = t80 * t81;
t91 = t81 * t76;
t90 = t81 * t79;
t89 = qJ(3) * t80;
t88 = g(3) * t106;
t87 = pkin(2) * t92 + t78 * pkin(6) + (pkin(1) + t70) * t81;
t84 = pkin(4) * t76 - qJ(5) * t79;
t83 = -pkin(1) - t106;
t58 = -t77 * t90 + t94;
t60 = t77 * t93 + t91;
t52 = g(1) * t58 - g(2) * t60 + t79 * t95;
t73 = t81 * pkin(6);
t68 = t81 * t89;
t66 = t78 * t89;
t61 = -t77 * t94 + t90;
t59 = t77 * t91 + t93;
t56 = t64 * t77 - t95;
t1 = [(-g(1) * t73 - g(2) * t87 - t83 * t99) * MDP(14) + (-g(1) * (t81 * pkin(3) + t61 * pkin(4) + t60 * qJ(5) + t73) - g(2) * (t59 * pkin(4) + pkin(7) * t92 + t58 * qJ(5) + t87) + (-g(1) * (-t80 * pkin(7) + t83) - g(2) * pkin(3)) * t78) * MDP(25) + t103 * (g(1) * t60 + g(2) * t58) + (MDP(3) - MDP(11)) * t64 + t104 * (-g(1) * t61 - g(2) * t59) + (t102 * t80 - t105 * t77 + MDP(2)) * (-g(2) * t81 + t99); (-g(1) * (-t100 * t81 + t68) - g(2) * (-t100 * t78 + t66) - t88) * MDP(14) + (-g(1) * t68 - g(2) * t66 - t88 + (-g(3) * pkin(7) - t64 * t84) * t80 + (-g(3) * t84 + t64 * (pkin(2) + pkin(7))) * t77) * MDP(25) + t102 * t56 + (-t103 * t79 - t104 * t76 + t105) * (g(3) * t77 + t64 * t80); (-MDP(14) - MDP(25)) * t56; (-g(1) * (-t58 * pkin(4) + t59 * qJ(5)) - g(2) * (t60 * pkin(4) - t61 * qJ(5)) - (-pkin(4) * t79 - qJ(5) * t76) * t95) * MDP(25) - t103 * (-g(1) * t59 + g(2) * t61 + t76 * t95) + t104 * t52; -t52 * MDP(25);];
taug = t1;
