% Calculate Gravitation load on the joints for
% S5RRRPP4
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
%   see S5RRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:46
% EndTime: 2019-12-31 20:55:47
% DurationCPUTime: 0.20s
% Computational Cost: add. (212->59), mult. (206->75), div. (0->0), fcn. (161->8), ass. (0->32)
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t47 = g(1) * t63 + g(2) * t61;
t59 = qJ(2) + qJ(3);
t54 = sin(t59);
t73 = t47 * t54;
t53 = pkin(8) + t59;
t50 = sin(t53);
t51 = cos(t53);
t64 = t51 * pkin(4) + t50 * qJ(5);
t72 = MDP(19) + MDP(23);
t71 = pkin(4) * t50;
t55 = cos(t59);
t70 = g(3) * t55;
t52 = pkin(3) * t55;
t62 = cos(qJ(2));
t56 = t62 * pkin(2);
t69 = t52 + t56;
t68 = qJ(5) * t51;
t67 = t52 + t64;
t60 = sin(qJ(2));
t43 = -t60 * pkin(2) - pkin(3) * t54;
t66 = t43 - t71;
t38 = -g(3) * t51 + t47 * t50;
t65 = t38 * MDP(20) + (-g(3) * t50 - t47 * t51) * MDP(22) + (-t70 + t73) * MDP(16) + (g(3) * t54 + t47 * t55) * MDP(17);
t46 = g(1) * t61 - g(2) * t63;
t58 = -qJ(4) - pkin(7) - pkin(6);
t45 = t63 * t68;
t44 = t61 * t68;
t42 = pkin(1) + t69;
t41 = t63 * t42;
t1 = [(-g(1) * (-t61 * t42 - t63 * t58) - g(2) * (-t61 * t58 + t41)) * MDP(19) + (-g(2) * t41 + (g(1) * t58 - g(2) * t64) * t63 + (-g(1) * (-t42 - t64) + g(2) * t58) * t61) * MDP(23) + (MDP(3) - MDP(18) - MDP(21)) * t47 + (-t60 * MDP(10) + MDP(16) * t55 - MDP(17) * t54 + t51 * MDP(20) + t50 * MDP(22) + t62 * MDP(9) + MDP(2)) * t46; (-g(3) * t62 + t47 * t60) * MDP(9) + (g(3) * t60 + t47 * t62) * MDP(10) + (-g(3) * t69 - t47 * t43) * MDP(19) + (-g(1) * (t66 * t63 + t45) - g(2) * (t66 * t61 + t44) - g(3) * (t56 + t67)) * MDP(23) + t65; (-g(1) * (-t63 * t71 + t45) - g(2) * (-t61 * t71 + t44) - g(3) * t67) * MDP(23) + (-MDP(19) * t70 + t72 * t73) * pkin(3) + t65; -t72 * t46; -t38 * MDP(23);];
taug = t1;
