% Calculate Gravitation load on the joints for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:37
% EndTime: 2019-12-31 18:45:38
% DurationCPUTime: 0.26s
% Computational Cost: add. (199->56), mult. (260->83), div. (0->0), fcn. (249->8), ass. (0->25)
t77 = MDP(11) - MDP(20);
t75 = MDP(17) + MDP(19);
t74 = MDP(18) - MDP(21);
t52 = qJ(1) + pkin(8);
t50 = sin(t52);
t51 = cos(t52);
t76 = -g(1) * t51 - g(2) * t50;
t54 = sin(qJ(3));
t71 = g(3) * t54;
t70 = t54 * pkin(7);
t53 = sin(qJ(4));
t57 = cos(qJ(3));
t69 = t53 * t57;
t56 = cos(qJ(4));
t68 = t56 * t57;
t63 = pkin(3) * t57 + pkin(2) + t70;
t62 = pkin(4) * t56 + qJ(5) * t53 + pkin(3);
t43 = t50 * t69 + t51 * t56;
t45 = -t50 * t56 + t51 * t69;
t38 = g(1) * t45 + g(2) * t43 + t53 * t71;
t58 = cos(qJ(1));
t55 = sin(qJ(1));
t46 = t50 * t53 + t51 * t68;
t44 = t50 * t68 - t51 * t53;
t1 = [(g(1) * t58 + g(2) * t55) * MDP(3) + (-g(1) * (-pkin(1) * t55 - pkin(4) * t44 - t43 * qJ(5)) - g(2) * (pkin(1) * t58 + t46 * pkin(4) + t45 * qJ(5)) + (-g(1) * pkin(6) - g(2) * t63) * t51 + (-g(2) * pkin(6) + g(1) * t63) * t50) * MDP(22) - t74 * (g(1) * t43 - g(2) * t45) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t55 - g(2) * t58) + t75 * (g(1) * t44 - g(2) * t46) + (t57 * MDP(10) - t77 * t54) * (g(1) * t50 - g(2) * t51); (-MDP(22) - MDP(4)) * g(3); (-g(3) * (t62 * t57 + t70) + t76 * (pkin(7) * t57 - t62 * t54)) * MDP(22) + t77 * (-t57 * t76 + t71) + (-t74 * t53 + t75 * t56 + MDP(10)) * (-g(3) * t57 - t76 * t54); (-g(1) * (-pkin(4) * t45 + qJ(5) * t46) - g(2) * (-pkin(4) * t43 + qJ(5) * t44) - (-pkin(4) * t53 + qJ(5) * t56) * t71) * MDP(22) + t74 * (g(1) * t46 + g(2) * t44 + t56 * t71) + t75 * t38; -t38 * MDP(22);];
taug = t1;
