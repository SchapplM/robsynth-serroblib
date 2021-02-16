% Calculate Gravitation load on the joints for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:13:19
% EndTime: 2021-01-15 21:13:20
% DurationCPUTime: 0.13s
% Computational Cost: add. (114->39), mult. (184->56), div. (0->0), fcn. (159->8), ass. (0->25)
t52 = cos(qJ(2));
t66 = pkin(1) * t52;
t47 = qJ(2) + qJ(4);
t45 = sin(t47);
t65 = g(3) * t45;
t48 = sin(qJ(5));
t50 = sin(qJ(1));
t63 = t50 * t48;
t51 = cos(qJ(5));
t62 = t50 * t51;
t53 = cos(qJ(1));
t61 = t53 * t48;
t60 = t53 * t51;
t59 = MDP(9) + MDP(11);
t58 = MDP(10) + MDP(12);
t46 = cos(t47);
t56 = g(1) * t53 + g(2) * t50;
t57 = (t46 * t56 + t65) * MDP(21) + (t51 * MDP(27) - t48 * MDP(28) + MDP(20)) * (-g(3) * t46 + t45 * t56);
t43 = g(1) * t50 - g(2) * t53;
t49 = sin(qJ(2));
t40 = t46 * t60 + t63;
t39 = -t46 * t61 + t62;
t38 = -t46 * t62 + t61;
t37 = t46 * t63 + t60;
t1 = [(-g(1) * (t53 * qJ(3) - t50 * t66) - g(2) * (t50 * qJ(3) + t53 * t66)) * MDP(14) + (-g(1) * t38 - g(2) * t40) * MDP(27) + (-g(1) * t37 - g(2) * t39) * MDP(28) + (MDP(3) - MDP(13)) * t56 + (MDP(20) * t46 - MDP(21) * t45 - t58 * t49 + t59 * t52 + MDP(2)) * t43; t58 * (g(3) * t49 + t52 * t56) + t57 + (MDP(14) * pkin(1) + t59) * (-g(3) * t52 + t49 * t56); -t43 * MDP(14); t57; (-g(1) * t39 + g(2) * t37 + t48 * t65) * MDP(27) + (g(1) * t40 - g(2) * t38 + t51 * t65) * MDP(28);];
taug = t1;
