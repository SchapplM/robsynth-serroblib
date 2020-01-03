% Calculate Gravitation load on the joints for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPPPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:30
% EndTime: 2019-12-31 17:46:31
% DurationCPUTime: 0.18s
% Computational Cost: add. (87->37), mult. (144->48), div. (0->0), fcn. (140->8), ass. (0->18)
t57 = cos(qJ(1));
t56 = sin(qJ(1));
t55 = t57 * pkin(1) + t56 * qJ(2);
t54 = cos(pkin(7));
t53 = sin(pkin(7));
t52 = MDP(13) + MDP(9);
t51 = t57 * pkin(2) + t55;
t50 = -t56 * pkin(1) + t57 * qJ(2);
t30 = -t56 * t53 - t57 * t54;
t31 = t57 * t53 - t56 * t54;
t49 = g(1) * t31 - g(2) * t30;
t48 = g(1) * t30 + g(2) * t31;
t47 = -t56 * pkin(2) + t50;
t44 = pkin(8) + qJ(5);
t38 = cos(t44);
t37 = sin(t44);
t32 = g(1) * t56 - g(2) * t57;
t1 = [(-g(1) * t50 - g(2) * t55) * MDP(6) + (-g(1) * t47 - g(2) * t51) * MDP(9) + (-g(1) * (t31 * pkin(3) + t30 * qJ(4) + t47) - g(2) * (-t30 * pkin(3) + t31 * qJ(4) + t51)) * MDP(13) + (MDP(3) - MDP(5)) * (g(1) * t57 + g(2) * t56) + (MDP(2) + MDP(4)) * t32 + (MDP(8) - MDP(12)) * t48 + (-MDP(10) * cos(pkin(8)) + MDP(11) * sin(pkin(8)) - t38 * MDP(19) + t37 * MDP(20) - MDP(7)) * t49; (-MDP(6) - t52) * t32; t52 * g(3); -t49 * MDP(13); (g(3) * t38 - t48 * t37) * MDP(19) + (-g(3) * t37 - t48 * t38) * MDP(20);];
taug = t1;
