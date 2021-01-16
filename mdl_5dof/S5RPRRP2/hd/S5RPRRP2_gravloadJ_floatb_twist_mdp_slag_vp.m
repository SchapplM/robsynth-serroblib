% Calculate Gravitation load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:36:28
% EndTime: 2021-01-15 12:36:28
% DurationCPUTime: 0.08s
% Computational Cost: add. (166->32), mult. (129->39), div. (0->0), fcn. (93->8), ass. (0->18)
t56 = -MDP(13) - MDP(15);
t55 = MDP(14) + MDP(16);
t44 = qJ(1) + pkin(8);
t43 = qJ(3) + t44;
t40 = sin(t43);
t41 = cos(t43);
t48 = cos(qJ(4));
t42 = t48 * pkin(4) + pkin(3);
t45 = -qJ(5) - pkin(7);
t54 = t40 * t42 + t41 * t45;
t53 = -t40 * t45 + t41 * t42;
t36 = g(2) * t41 + g(3) * t40;
t35 = g(2) * t40 - g(3) * t41;
t46 = sin(qJ(4));
t51 = (-MDP(17) + MDP(7)) * t35 + (t55 * t46 + t56 * t48 - MDP(6)) * t36;
t49 = cos(qJ(1));
t47 = sin(qJ(1));
t1 = [(g(2) * t47 - g(3) * t49) * MDP(3) + (-g(2) * (pkin(2) * cos(t44) + t49 * pkin(1) + t53) - g(3) * (pkin(2) * sin(t44) + t47 * pkin(1) + t54)) * MDP(18) + t51 + (pkin(1) * MDP(4) + MDP(2)) * (-g(2) * t49 - g(3) * t47); (-MDP(18) - MDP(4)) * g(1); (-g(2) * t53 - g(3) * t54) * MDP(18) + t51; t55 * (g(1) * t46 + t35 * t48) + (MDP(18) * pkin(4) - t56) * (-g(1) * t48 + t35 * t46); t36 * MDP(18);];
taug = t1;
