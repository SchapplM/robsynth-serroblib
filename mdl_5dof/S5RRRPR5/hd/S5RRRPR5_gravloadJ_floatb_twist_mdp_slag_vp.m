% Calculate Gravitation load on the joints for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:12
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:10:44
% EndTime: 2021-01-15 23:10:45
% DurationCPUTime: 0.20s
% Computational Cost: add. (196->50), mult. (210->71), div. (0->0), fcn. (179->10), ass. (0->29)
t59 = qJ(2) + qJ(3);
t54 = pkin(9) + t59;
t51 = sin(t54);
t76 = g(3) * t51;
t60 = sin(qJ(5));
t65 = cos(qJ(1));
t74 = t60 * t65;
t62 = sin(qJ(1));
t73 = t62 * t60;
t63 = cos(qJ(5));
t72 = t62 * t63;
t71 = t63 * t65;
t56 = cos(t59);
t64 = cos(qJ(2));
t70 = t64 * pkin(2) + pkin(3) * t56;
t69 = g(1) * t65 + g(2) * t62;
t49 = g(1) * t62 - g(2) * t65;
t52 = cos(t54);
t55 = sin(t59);
t66 = -g(3) * t56 + t69 * t55;
t68 = (t69 * t52 + t76) * MDP(19) + t66 * MDP(16) + (g(3) * t55 + t69 * t56) * MDP(17) + (t63 * MDP(27) - t60 * MDP(28) + MDP(18)) * (-g(3) * t52 + t69 * t51);
t61 = sin(qJ(2));
t58 = -qJ(4) - pkin(7) - pkin(6);
t47 = pkin(1) + t70;
t46 = t52 * t71 + t73;
t45 = -t52 * t74 + t72;
t44 = -t52 * t72 + t74;
t43 = t52 * t73 + t71;
t1 = [(-g(1) * (-t62 * t47 - t58 * t65) - g(2) * (t47 * t65 - t62 * t58)) * MDP(21) + (-g(1) * t44 - g(2) * t46) * MDP(27) + (-g(1) * t43 - g(2) * t45) * MDP(28) + (MDP(3) - MDP(20)) * t69 + (-t61 * MDP(10) + MDP(16) * t56 - MDP(17) * t55 + MDP(18) * t52 - MDP(19) * t51 + t64 * MDP(9) + MDP(2)) * t49; (-g(3) * t64 + t69 * t61) * MDP(9) + (g(3) * t61 + t69 * t64) * MDP(10) + (-g(3) * t70 - t69 * (-pkin(2) * t61 - pkin(3) * t55)) * MDP(21) + t68; t66 * MDP(21) * pkin(3) + t68; -t49 * MDP(21); (-g(1) * t45 + g(2) * t43 + t60 * t76) * MDP(27) + (g(1) * t46 - g(2) * t44 + t63 * t76) * MDP(28);];
taug = t1;
