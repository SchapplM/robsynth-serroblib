% Calculate Gravitation load on the joints for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:37
% EndTime: 2019-12-31 17:00:38
% DurationCPUTime: 0.20s
% Computational Cost: add. (78->38), mult. (165->49), div. (0->0), fcn. (132->4), ass. (0->21)
t47 = sin(qJ(2));
t49 = cos(qJ(2));
t56 = t49 * pkin(2) + t47 * qJ(3);
t51 = -pkin(1) - t56;
t64 = MDP(9) - MDP(12) + MDP(17);
t63 = MDP(10) - MDP(13) - MDP(16);
t48 = sin(qJ(1));
t50 = cos(qJ(1));
t36 = g(1) * t50 + g(2) * t48;
t62 = t36 * t47;
t61 = pkin(2) * t47;
t60 = g(1) * t48;
t55 = qJ(3) * t49;
t54 = t49 * qJ(4);
t53 = t48 * pkin(5) - t51 * t50;
t44 = t50 * pkin(5);
t39 = t50 * t55;
t37 = t48 * t55;
t33 = g(3) * t47 + t36 * t49;
t32 = -g(3) * t49 + t62;
t1 = [(-g(1) * t44 - g(2) * t53 - t51 * t60) * MDP(14) + (-g(1) * (t50 * pkin(3) + t44) - g(2) * (t50 * t54 + t53) + (-g(1) * (t51 - t54) - g(2) * pkin(3)) * t48) * MDP(18) + (MDP(3) - MDP(11) - MDP(15)) * t36 + (-t63 * t47 + t64 * t49 + MDP(2)) * (-g(2) * t50 + t60); (-g(1) * (-t50 * t61 + t39) - g(2) * (-t48 * t61 + t37) - g(3) * t56) * MDP(14) + (-g(1) * t39 - g(2) * t37 - g(3) * (t54 + t56) + (pkin(2) + qJ(4)) * t62) * MDP(18) + t63 * t33 + t64 * t32; (-MDP(14) - MDP(18)) * t32; -t33 * MDP(18);];
taug = t1;
