% Calculate Gravitation load on the joints for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:15
% EndTime: 2019-12-31 18:30:16
% DurationCPUTime: 0.30s
% Computational Cost: add. (165->56), mult. (196->80), div. (0->0), fcn. (174->10), ass. (0->32)
t63 = sin(qJ(1));
t64 = cos(qJ(1));
t49 = g(1) * t64 + g(2) * t63;
t81 = MDP(14) - MDP(17);
t57 = pkin(8) + qJ(3);
t52 = sin(t57);
t54 = cos(t57);
t41 = -g(3) * t54 + t49 * t52;
t78 = g(3) * t52;
t56 = pkin(9) + qJ(5);
t51 = sin(t56);
t76 = t63 * t51;
t53 = cos(t56);
t75 = t63 * t53;
t58 = sin(pkin(9));
t74 = t63 * t58;
t60 = cos(pkin(9));
t73 = t63 * t60;
t72 = t64 * t51;
t71 = t64 * t53;
t70 = t64 * t58;
t69 = t64 * t60;
t48 = g(1) * t63 - g(2) * t64;
t68 = t54 * pkin(3) + t52 * qJ(4);
t61 = cos(pkin(8));
t66 = t61 * pkin(2) + pkin(1) + t68;
t62 = -pkin(6) - qJ(2);
t46 = t54 * t71 + t76;
t45 = -t54 * t72 + t75;
t44 = -t54 * t75 + t72;
t43 = t54 * t76 + t71;
t1 = [(-g(1) * (-t63 * pkin(1) + t64 * qJ(2)) - g(2) * (t64 * pkin(1) + t63 * qJ(2))) * MDP(7) + (-g(1) * (-t54 * t73 + t70) - g(2) * (t54 * t69 + t74)) * MDP(15) + (-g(1) * (t54 * t74 + t69) - g(2) * (-t54 * t70 + t73)) * MDP(16) + ((g(1) * t62 - g(2) * t66) * t64 + (g(1) * t66 + g(2) * t62) * t63) * MDP(18) + (-g(1) * t44 - g(2) * t46) * MDP(24) + (-g(1) * t43 - g(2) * t45) * MDP(25) + (MDP(3) - MDP(6)) * t49 + (t54 * MDP(13) + MDP(4) * t61 - MDP(5) * sin(pkin(8)) - t81 * t52 + MDP(2)) * t48; (-MDP(18) - MDP(7)) * t48; (-g(3) * t68 + t49 * (pkin(3) * t52 - qJ(4) * t54)) * MDP(18) + t81 * (t49 * t54 + t78) + (MDP(15) * t60 - MDP(16) * t58 + MDP(24) * t53 - MDP(25) * t51 + MDP(13)) * t41; -t41 * MDP(18); (-g(1) * t45 + g(2) * t43 + t51 * t78) * MDP(24) + (g(1) * t46 - g(2) * t44 + t53 * t78) * MDP(25);];
taug = t1;
