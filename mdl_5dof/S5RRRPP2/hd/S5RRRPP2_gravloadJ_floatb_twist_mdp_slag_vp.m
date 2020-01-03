% Calculate Gravitation load on the joints for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:59
% EndTime: 2019-12-31 20:52:00
% DurationCPUTime: 0.19s
% Computational Cost: add. (232->54), mult. (250->65), div. (0->0), fcn. (199->6), ass. (0->33)
t99 = MDP(12) + MDP(14) + MDP(18);
t98 = MDP(13) - MDP(16) - MDP(19);
t75 = qJ(1) + qJ(2);
t69 = sin(t75);
t70 = cos(t75);
t60 = g(1) * t70 + g(2) * t69;
t76 = sin(qJ(3));
t97 = t60 * t76;
t71 = t76 * qJ(4);
t78 = cos(qJ(3));
t88 = t78 * pkin(3) + t71;
t95 = pkin(3) * t76;
t94 = g(1) * t69;
t77 = sin(qJ(1));
t91 = t77 * pkin(1);
t90 = t78 * pkin(4);
t89 = t70 * t78;
t87 = qJ(4) * t78;
t67 = t70 * pkin(7);
t86 = -t70 * qJ(5) + t67;
t85 = pkin(3) * t89 + t69 * pkin(7) + (pkin(2) + t71) * t70;
t79 = cos(qJ(1));
t84 = t79 * pkin(1) + t85;
t59 = -g(2) * t70 + t94;
t83 = -pkin(2) - t88;
t82 = t83 * t94;
t81 = (-MDP(15) + MDP(20) + MDP(6)) * t60 + (-t98 * t76 + t99 * t78 + MDP(5)) * t59;
t80 = (-g(1) * (t83 - t90) + g(2) * qJ(5)) * t69;
t64 = pkin(4) * t89;
t63 = t70 * t87;
t61 = t69 * t87;
t45 = -g(3) * t78 + t97;
t1 = [t81 + (-g(1) * (t86 - t91) - g(2) * (t64 + t84) + t80) * MDP(21) + (-g(1) * (t67 - t91) - g(2) * t84 - t82) * MDP(17) + (g(1) * t77 - g(2) * t79) * MDP(2) + (g(1) * t79 + g(2) * t77) * MDP(3); (-g(1) * t67 - g(2) * t85 - t82) * MDP(17) + (-g(1) * t86 - g(2) * (t64 + t85) + t80) * MDP(21) + t81; (-g(1) * (-t70 * t95 + t63) - g(2) * (-t69 * t95 + t61) - g(3) * t88) * MDP(17) + (-g(1) * t63 - g(2) * t61 - g(3) * (t88 + t90) + (pkin(3) + pkin(4)) * t97) * MDP(21) + t98 * (g(3) * t76 + t60 * t78) + t99 * t45; (-MDP(17) - MDP(21)) * t45; t59 * MDP(21);];
taug = t1;
