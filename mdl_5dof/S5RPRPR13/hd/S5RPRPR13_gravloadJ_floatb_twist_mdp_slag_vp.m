% Calculate Gravitation load on the joints for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:53
% EndTime: 2019-12-31 18:32:54
% DurationCPUTime: 0.20s
% Computational Cost: add. (129->47), mult. (174->64), div. (0->0), fcn. (149->8), ass. (0->26)
t56 = sin(qJ(1));
t58 = cos(qJ(1));
t45 = g(1) * t58 + g(2) * t56;
t73 = MDP(13) - MDP(16);
t72 = MDP(14) - MDP(17);
t51 = pkin(8) + qJ(3);
t49 = cos(t51);
t67 = g(3) * t49;
t55 = sin(qJ(5));
t66 = t56 * t55;
t57 = cos(qJ(5));
t65 = t56 * t57;
t64 = t58 * t55;
t63 = t58 * t57;
t44 = g(1) * t56 - g(2) * t58;
t48 = sin(t51);
t62 = t49 * pkin(3) + t48 * qJ(4);
t53 = cos(pkin(8));
t60 = t53 * pkin(2) + pkin(1) + t62;
t54 = -pkin(6) - qJ(2);
t43 = -t48 * t66 + t63;
t42 = t48 * t65 + t64;
t41 = t48 * t64 + t65;
t40 = t48 * t63 - t66;
t36 = t45 * t48 - t67;
t1 = [(-g(1) * (-t56 * pkin(1) + t58 * qJ(2)) - g(2) * (t58 * pkin(1) + t56 * qJ(2))) * MDP(7) + ((g(1) * t54 - g(2) * t60) * t58 + (g(1) * t60 + g(2) * t54) * t56) * MDP(18) + (-g(1) * t43 - g(2) * t41) * MDP(24) + (g(1) * t42 - g(2) * t40) * MDP(25) + (MDP(3) - MDP(6) - MDP(15)) * t45 + (MDP(4) * t53 - MDP(5) * sin(pkin(8)) - t72 * t48 + t73 * t49 + MDP(2)) * t44; (-MDP(18) - MDP(7)) * t44; (-g(3) * t62 + t45 * (pkin(3) * t48 - qJ(4) * t49)) * MDP(18) + t73 * t36 + (-MDP(24) * t55 - MDP(25) * t57 + t72) * (g(3) * t48 + t45 * t49); -t36 * MDP(18); (-g(1) * t40 - g(2) * t42 + t57 * t67) * MDP(24) + (g(1) * t41 - g(2) * t43 - t55 * t67) * MDP(25);];
taug = t1;
