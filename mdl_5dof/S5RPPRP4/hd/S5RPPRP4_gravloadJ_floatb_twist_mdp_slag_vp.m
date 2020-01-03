% Calculate Gravitation load on the joints for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:20
% EndTime: 2019-12-31 17:52:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (78->37), mult. (141->48), div. (0->0), fcn. (133->6), ass. (0->19)
t57 = cos(qJ(1));
t56 = sin(qJ(1));
t55 = t57 * pkin(1) + t56 * qJ(2);
t54 = cos(pkin(7));
t53 = sin(pkin(7));
t52 = MDP(18) + MDP(9);
t51 = t57 * pkin(2) + t55;
t50 = -pkin(1) * t56 + t57 * qJ(2);
t30 = -t53 * t56 - t54 * t57;
t31 = t53 * t57 - t54 * t56;
t49 = g(1) * t31 - g(2) * t30;
t48 = g(1) * t30 + g(2) * t31;
t47 = -pkin(2) * t56 + t50;
t45 = cos(qJ(4));
t44 = sin(qJ(4));
t43 = -qJ(5) - pkin(6);
t37 = pkin(4) * t45 + pkin(3);
t32 = g(1) * t56 - g(2) * t57;
t1 = [(-g(1) * t50 - g(2) * t55) * MDP(6) + (-g(1) * t47 - g(2) * t51) * MDP(9) + (-g(1) * (-t30 * t43 + t31 * t37 + t47) - g(2) * (-t30 * t37 - t31 * t43 + t51)) * MDP(18) + (MDP(3) - MDP(5)) * (g(1) * t57 + g(2) * t56) + (MDP(2) + MDP(4)) * t32 + (MDP(8) - MDP(17)) * t48 + (-MDP(15) * t45 + t44 * MDP(16) - MDP(7)) * t49; (-MDP(6) - t52) * t32; t52 * g(3); (-g(3) * t44 - t45 * t48) * MDP(16) + (MDP(18) * pkin(4) + MDP(15)) * (g(3) * t45 - t44 * t48); -t49 * MDP(18);];
taug = t1;
