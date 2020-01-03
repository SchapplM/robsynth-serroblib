% Calculate Gravitation load on the joints for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:48
% EndTime: 2019-12-31 19:41:49
% DurationCPUTime: 0.32s
% Computational Cost: add. (100->54), mult. (220->75), div. (0->0), fcn. (189->6), ass. (0->34)
t90 = MDP(9) + MDP(11) - MDP(16);
t89 = -MDP(10) + MDP(13) + MDP(15);
t67 = sin(qJ(1));
t70 = cos(qJ(1));
t54 = g(1) * t70 + g(2) * t67;
t66 = sin(qJ(2));
t88 = t54 * t66;
t59 = t66 * qJ(3);
t69 = cos(qJ(2));
t74 = t69 * pkin(2) + t59;
t86 = pkin(2) * t66;
t85 = g(1) * t67;
t81 = g(3) * t69;
t80 = t69 * pkin(3);
t65 = sin(qJ(5));
t79 = t67 * t65;
t68 = cos(qJ(5));
t78 = t67 * t68;
t77 = t69 * t70;
t76 = t70 * t65;
t75 = t70 * t68;
t73 = qJ(3) * t69;
t72 = pkin(2) * t77 + t67 * pkin(6) + (pkin(1) + t59) * t70;
t53 = -g(2) * t70 + t85;
t71 = -pkin(1) - t74;
t62 = t70 * pkin(6);
t57 = t70 * t73;
t55 = t67 * t73;
t50 = t66 * t75 - t79;
t49 = -t66 * t76 - t78;
t48 = -t66 * t78 - t76;
t47 = t66 * t79 - t75;
t45 = -t81 + t88;
t1 = [(-g(1) * t62 - g(2) * t72 - t71 * t85) * MDP(14) + (-g(1) * (-t70 * qJ(4) + t62) - g(2) * (pkin(3) * t77 + t72) + (-g(1) * (t71 - t80) + g(2) * qJ(4)) * t67) * MDP(18) + (-g(1) * t48 - g(2) * t50) * MDP(24) + (-g(1) * t47 - g(2) * t49) * MDP(25) + (MDP(3) - MDP(12) + MDP(17)) * t54 + (t89 * t66 + t90 * t69 + MDP(2)) * t53; (-g(1) * (-t70 * t86 + t57) - g(2) * (-t67 * t86 + t55) - g(3) * t74) * MDP(14) + (-g(1) * t57 - g(2) * t55 - g(3) * (t74 + t80) + (pkin(2) + pkin(3)) * t88) * MDP(18) + t90 * t45 + (-MDP(24) * t68 + MDP(25) * t65 - t89) * (g(3) * t66 + t54 * t69); (-MDP(14) - MDP(18)) * t45; t53 * MDP(18); (-g(1) * t49 + g(2) * t47 - t65 * t81) * MDP(24) + (g(1) * t50 - g(2) * t48 - t68 * t81) * MDP(25);];
taug = t1;
