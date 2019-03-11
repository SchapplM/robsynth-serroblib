% Calculate Gravitation load on the joints for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:26
% EndTime: 2019-03-08 18:26:26
% DurationCPUTime: 0.14s
% Computational Cost: add. (88->41), mult. (197->56), div. (0->0), fcn. (199->6), ass. (0->24)
t51 = sin(pkin(4));
t71 = pkin(3) * t51;
t70 = g(3) * t51;
t50 = sin(pkin(6));
t54 = sin(qJ(1));
t69 = t54 * t50;
t52 = cos(pkin(6));
t68 = t54 * t52;
t55 = cos(qJ(1));
t67 = t55 * t50;
t66 = t55 * t52;
t64 = qJ(2) * t51;
t65 = t55 * pkin(1) + t54 * t64;
t63 = MDP(11) + MDP(15);
t62 = -t54 * pkin(1) + t55 * t64;
t58 = g(1) * t54 - g(2) * t55;
t53 = cos(pkin(4));
t43 = t53 * t68 + t67;
t44 = -t53 * t69 + t66;
t57 = t44 * pkin(2) + t43 * qJ(3) + t65;
t41 = -t53 * t66 + t69;
t42 = t53 * t67 + t68;
t56 = -t42 * pkin(2) - t41 * qJ(3) + t62;
t1 = [t58 * MDP(2) + (-g(1) * t62 - g(2) * t65) * MDP(7) + (-g(1) * t56 - g(2) * t57) * MDP(11) + (-g(1) * (-t42 * qJ(4) + t55 * t71 + t56) - g(2) * (t44 * qJ(4) + t54 * t71 + t57)) * MDP(15) + (-MDP(5) + MDP(10) + MDP(13)) * (g(1) * t41 - g(2) * t43) + (-MDP(9) + MDP(4) + MDP(14)) * (g(1) * t42 - g(2) * t44) + (MDP(3) - (MDP(6) + MDP(8) + MDP(12)) * t51) * (g(1) * t55 + g(2) * t54); (MDP(7) + t63) * (-g(3) * t53 - t58 * t51); t63 * (-g(1) * t43 - g(2) * t41 + t52 * t70); (-g(1) * t44 - g(2) * t42 - t50 * t70) * MDP(15);];
taug  = t1;
