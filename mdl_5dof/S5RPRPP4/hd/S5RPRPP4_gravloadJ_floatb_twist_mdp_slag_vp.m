% Calculate Gravitation load on the joints for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:50
% EndTime: 2019-12-31 18:14:51
% DurationCPUTime: 0.18s
% Computational Cost: add. (103->44), mult. (152->54), div. (0->0), fcn. (115->6), ass. (0->19)
t45 = sin(qJ(1));
t47 = cos(qJ(1));
t59 = -g(1) * t45 + g(2) * t47;
t44 = sin(qJ(3));
t56 = t44 * pkin(3);
t55 = t47 * pkin(1) + t45 * qJ(2);
t54 = -MDP(15) - MDP(19);
t53 = -t45 * pkin(1) + t47 * qJ(2);
t32 = g(1) * t47 + g(2) * t45;
t42 = qJ(3) + pkin(7);
t36 = sin(t42);
t37 = cos(t42);
t52 = t36 * pkin(4) - t37 * qJ(5);
t43 = -qJ(4) - pkin(6);
t51 = t45 * t43 + t47 * t56 + t53;
t50 = -t47 * t43 + t45 * t56 + t55;
t46 = cos(qJ(3));
t30 = -g(3) * t36 - t37 * t59;
t1 = [(-g(1) * t53 - g(2) * t55) * MDP(6) + (-g(1) * t51 - g(2) * t50) * MDP(15) + (-g(1) * (t52 * t47 + t51) - g(2) * (t52 * t45 + t50)) * MDP(19) - (MDP(2) - MDP(4) + MDP(14) + MDP(17)) * t59 + (-t44 * MDP(12) - t46 * MDP(13) - t36 * MDP(16) + t37 * MDP(18) + MDP(3) - MDP(5)) * t32; -(-MDP(6) + t54) * t59; (g(3) * t46 - t44 * t59) * MDP(13) - t30 * MDP(16) + (-g(3) * t37 + t36 * t59) * MDP(18) + (-g(3) * (-t52 - t56) + t59 * (pkin(3) * t46 + pkin(4) * t37 + qJ(5) * t36)) * MDP(19) + (pkin(3) * MDP(15) + MDP(12)) * (g(3) * t44 + t46 * t59); t54 * t32; t30 * MDP(19);];
taug = t1;
