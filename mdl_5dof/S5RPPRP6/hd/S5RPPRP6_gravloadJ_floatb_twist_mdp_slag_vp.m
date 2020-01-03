% Calculate Gravitation load on the joints for
% S5RPPRP6
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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:17
% EndTime: 2019-12-31 17:55:18
% DurationCPUTime: 0.22s
% Computational Cost: add. (107->42), mult. (144->47), div. (0->0), fcn. (109->6), ass. (0->19)
t47 = sin(qJ(1));
t48 = cos(qJ(1));
t59 = -g(1) * t47 + g(2) * t48;
t58 = MDP(16) + MDP(18);
t57 = MDP(17) - MDP(20);
t54 = t48 * pkin(1) + t47 * qJ(2);
t53 = -MDP(10) - MDP(21);
t52 = g(2) * t54;
t34 = g(1) * t48 + g(2) * t47;
t43 = pkin(7) + qJ(4);
t37 = sin(t43);
t38 = cos(t43);
t50 = -t37 * pkin(4) + t38 * qJ(5);
t44 = sin(pkin(7));
t49 = pkin(3) * t44 - t50;
t46 = -pkin(6) - qJ(3);
t40 = t48 * qJ(2);
t30 = -g(3) * t37 - t59 * t38;
t1 = [(-g(1) * (-t47 * pkin(1) + t40) - t52) * MDP(6) + (-g(1) * (t40 + (-pkin(1) - qJ(3)) * t47) - g(2) * (t48 * qJ(3) + t54)) * MDP(10) + (-g(1) * t40 - t52 + (-g(1) * t49 + g(2) * t46) * t48 + (-g(1) * (-pkin(1) + t46) - g(2) * t49) * t47) * MDP(21) - (MDP(2) - MDP(4) + MDP(9) + MDP(19)) * t59 + (-t44 * MDP(7) - MDP(8) * cos(pkin(7)) - t58 * t37 - t57 * t38 + MDP(3) - MDP(5)) * t34; -(-MDP(6) + t53) * t59; t53 * t34; (-g(3) * t50 + t59 * (pkin(4) * t38 + qJ(5) * t37)) * MDP(21) - t58 * t30 + t57 * (g(3) * t38 - t37 * t59); t30 * MDP(21);];
taug = t1;
