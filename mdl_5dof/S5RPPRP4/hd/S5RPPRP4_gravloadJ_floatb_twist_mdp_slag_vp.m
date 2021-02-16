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
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 17:13:16
% EndTime: 2021-01-15 17:13:17
% DurationCPUTime: 0.17s
% Computational Cost: add. (92->40), mult. (167->52), div. (0->0), fcn. (155->6), ass. (0->22)
t62 = MDP(13) + MDP(15);
t61 = MDP(14) + MDP(16);
t45 = sin(pkin(7));
t46 = cos(pkin(7));
t47 = qJ(5) + pkin(6);
t50 = cos(qJ(4));
t56 = pkin(4) * t50 + pkin(3);
t60 = t45 * t56 - t47 * t46 + qJ(2);
t52 = pkin(1) + pkin(2);
t57 = MDP(18) + MDP(7);
t49 = sin(qJ(1));
t51 = cos(qJ(1));
t35 = -t49 * t45 - t51 * t46;
t36 = t45 * t51 - t49 * t46;
t54 = g(1) * t36 - g(2) * t35;
t53 = -g(1) * t35 - g(2) * t36;
t48 = sin(qJ(4));
t43 = t51 * qJ(2);
t42 = t49 * qJ(2);
t37 = g(1) * t49 - g(2) * t51;
t34 = t47 * t45 + t46 * t56 + t52;
t1 = [(-g(1) * (-pkin(1) * t49 + t43) - g(2) * (pkin(1) * t51 + t42)) * MDP(6) + (-g(1) * (-t49 * t52 + t43) - g(2) * (t51 * t52 + t42)) * MDP(7) + t53 * MDP(17) + (-g(1) * (-t34 * t49 + t51 * t60) - g(2) * (t34 * t51 + t49 * t60)) * MDP(18) + (MDP(3) - MDP(5)) * (g(1) * t51 + g(2) * t49) + (MDP(2) + MDP(4)) * t37 + (t61 * t48 - t50 * t62) * t54; (-MDP(6) - t57) * t37; t57 * g(3); t61 * (-g(3) * t48 + t50 * t53) + (MDP(18) * pkin(4) + t62) * (g(3) * t50 + t48 * t53); -t54 * MDP(18);];
taug = t1;
