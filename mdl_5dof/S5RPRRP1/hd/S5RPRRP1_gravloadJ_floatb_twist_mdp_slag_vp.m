% Calculate Gravitation load on the joints for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:20
% EndTime: 2021-01-15 12:27:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (120->38), mult. (157->47), div. (0->0), fcn. (119->6), ass. (0->17)
t55 = MDP(19) + MDP(21);
t54 = MDP(20) + MDP(22);
t40 = qJ(3) + qJ(4);
t36 = sin(t40);
t41 = sin(qJ(3));
t53 = -t41 * pkin(3) - pkin(4) * t36;
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t52 = -g(1) * t42 + g(2) * t44;
t47 = qJ(2) - t53;
t37 = cos(t40);
t30 = g(3) * t36 + t37 * t52;
t46 = t55 * t30 + t54 * (g(3) * t37 - t36 * t52);
t34 = g(1) * t44 + g(2) * t42;
t43 = cos(qJ(3));
t38 = qJ(5) + pkin(1) + pkin(6) + pkin(7);
t1 = [(-g(1) * (-t42 * pkin(1) + t44 * qJ(2)) - g(2) * (t44 * pkin(1) + t42 * qJ(2))) * MDP(6) + (-g(1) * (-t38 * t42 + t47 * t44) - g(2) * (t38 * t44 + t47 * t42)) * MDP(24) - (MDP(2) - MDP(4) + MDP(23)) * t52 + (-t41 * MDP(12) - t43 * MDP(13) - t55 * t36 - t54 * t37 + MDP(3) - MDP(5)) * t34; -(-MDP(24) - MDP(6)) * t52; (g(3) * t41 + t43 * t52) * MDP(12) + (g(3) * t43 - t41 * t52) * MDP(13) + (-g(3) * t53 + t52 * (pkin(3) * t43 + pkin(4) * t37)) * MDP(24) + t46; t30 * MDP(24) * pkin(4) + t46; -t34 * MDP(24);];
taug = t1;
