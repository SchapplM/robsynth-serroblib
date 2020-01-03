% Calculate Gravitation load on the joints for
% S5RRPRP10
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
%   see S5RRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:06
% EndTime: 2019-12-31 20:11:07
% DurationCPUTime: 0.32s
% Computational Cost: add. (111->62), mult. (241->87), div. (0->0), fcn. (209->6), ass. (0->38)
t77 = sin(qJ(2));
t69 = t77 * qJ(3);
t80 = cos(qJ(2));
t90 = t80 * pkin(2) + t69;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t63 = g(1) * t81 + g(2) * t78;
t109 = MDP(10) - MDP(13);
t107 = MDP(9) - MDP(12) + MDP(22);
t56 = g(3) * t77 + t63 * t80;
t76 = sin(qJ(4));
t106 = pkin(4) * t76;
t104 = g(1) * t78;
t99 = g(3) * t80;
t75 = -qJ(5) - pkin(7);
t97 = t75 * t80;
t96 = t76 * t80;
t95 = t77 * t81;
t94 = t78 * t76;
t79 = cos(qJ(4));
t93 = t78 * t79;
t92 = t79 * t81;
t91 = t80 * t81;
t89 = qJ(3) * t80;
t88 = pkin(4) * t96;
t86 = t76 * t95;
t85 = pkin(2) * t91 + t78 * pkin(6) + (pkin(1) + t69) * t81;
t83 = -pkin(1) - t90;
t57 = t77 * t92 - t94;
t59 = t76 * t81 + t77 * t93;
t72 = t81 * pkin(6);
t68 = pkin(4) * t79 + pkin(3);
t66 = t81 * t89;
t64 = t78 * t89;
t60 = -t77 * t94 + t92;
t58 = t86 + t93;
t55 = t63 * t77 - t99;
t1 = [(-g(1) * t72 - g(2) * t85 - t83 * t104) * MDP(14) + (-g(1) * t60 - g(2) * t58) * MDP(20) + (g(1) * t59 - g(2) * t57) * MDP(21) + (-g(1) * (t68 * t81 + t72) - g(2) * (pkin(4) * t86 - t75 * t91 + t85) + (-g(1) * (-t77 * t106 + t83 + t97) - g(2) * t68) * t78) * MDP(23) + (MDP(3) - MDP(11)) * t63 + (t107 * t80 - t109 * t77 + MDP(2)) * (-g(2) * t81 + t104); (-g(1) * (-pkin(2) * t95 + t66) - g(2) * (-pkin(2) * t77 * t78 + t64) - g(3) * t90) * MDP(14) + (-g(1) * (t81 * t88 + t66) - g(2) * (t78 * t88 + t64) - g(3) * (t90 - t97) + (-g(3) * t106 + t63 * (pkin(2) - t75)) * t77) * MDP(23) + t107 * t55 + (-t76 * MDP(20) - MDP(21) * t79 + t109) * t56; (-MDP(14) - MDP(23)) * t55; (g(1) * t58 - g(2) * t60 - g(3) * t96) * MDP(21) + (pkin(4) * MDP(23) + MDP(20)) * (-g(1) * t57 - g(2) * t59 + t79 * t99); -t56 * MDP(23);];
taug = t1;
