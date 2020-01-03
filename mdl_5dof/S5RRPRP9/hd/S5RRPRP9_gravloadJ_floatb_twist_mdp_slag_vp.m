% Calculate Gravitation load on the joints for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:28
% EndTime: 2019-12-31 20:07:30
% DurationCPUTime: 0.46s
% Computational Cost: add. (233->75), mult. (348->107), div. (0->0), fcn. (331->8), ass. (0->34)
t72 = cos(pkin(8));
t63 = t72 * pkin(3) + pkin(2);
t73 = -pkin(7) - qJ(3);
t74 = sin(qJ(2));
t76 = cos(qJ(2));
t101 = t76 * t63 - t74 * t73;
t75 = sin(qJ(1));
t77 = cos(qJ(1));
t83 = g(1) * t77 + g(2) * t75;
t100 = MDP(20) + MDP(22);
t99 = MDP(21) - MDP(24);
t58 = -g(3) * t76 + t74 * t83;
t98 = MDP(10) - MDP(13) - MDP(23);
t97 = g(1) * t75;
t94 = g(3) * t74;
t91 = t75 * t76;
t70 = pkin(8) + qJ(4);
t64 = sin(t70);
t89 = t77 * t64;
t65 = cos(t70);
t88 = t77 * t65;
t71 = sin(pkin(8));
t87 = t77 * t71;
t86 = t77 * t72;
t85 = t77 * pkin(1) + t75 * pkin(6);
t81 = t76 * pkin(2) + t74 * qJ(3);
t79 = pkin(4) * t65 + qJ(5) * t64 + t63;
t54 = t64 * t91 + t88;
t56 = -t75 * t65 + t76 * t89;
t50 = g(1) * t56 + g(2) * t54 + t64 * t94;
t67 = t77 * pkin(6);
t57 = t75 * t64 + t76 * t88;
t55 = t65 * t91 - t89;
t1 = [t83 * MDP(3) + (-g(1) * (-t72 * t91 + t87) - g(2) * (t75 * t71 + t76 * t86)) * MDP(11) + (-g(1) * (t71 * t91 + t86) - g(2) * (t75 * t72 - t76 * t87)) * MDP(12) + (-g(1) * t67 - g(2) * (t77 * t81 + t85) - (-pkin(1) - t81) * t97) * MDP(14) + (-g(1) * (pkin(3) * t87 - t55 * pkin(4) - t54 * qJ(5) + t67) - g(2) * (t57 * pkin(4) + t56 * qJ(5) + t101 * t77 + t85) + (-g(1) * (-pkin(1) - t101) - g(2) * pkin(3) * t71) * t75) * MDP(25) - t99 * (g(1) * t54 - g(2) * t56) + t100 * (g(1) * t55 - g(2) * t57) + (t76 * MDP(9) - t98 * t74 + MDP(2)) * (-g(2) * t77 + t97); (-g(3) * t81 + t83 * (pkin(2) * t74 - qJ(3) * t76)) * MDP(14) + ((-g(3) * t79 + t73 * t83) * t76 + (g(3) * t73 + t83 * t79) * t74) * MDP(25) + t98 * (t76 * t83 + t94) + (MDP(11) * t72 - t71 * MDP(12) + t100 * t65 - t99 * t64 + MDP(9)) * t58; (-MDP(14) - MDP(25)) * t58; (-g(1) * (-t56 * pkin(4) + t57 * qJ(5)) - g(2) * (-t54 * pkin(4) + t55 * qJ(5)) - (-pkin(4) * t64 + qJ(5) * t65) * t94) * MDP(25) + t99 * (g(1) * t57 + g(2) * t55 + t65 * t94) + t100 * t50; -t50 * MDP(25);];
taug = t1;
