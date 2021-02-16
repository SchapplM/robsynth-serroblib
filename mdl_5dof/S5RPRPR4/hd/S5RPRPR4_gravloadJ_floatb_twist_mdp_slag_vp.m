% Calculate Gravitation load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:45
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:44:15
% EndTime: 2021-01-15 11:44:16
% DurationCPUTime: 0.14s
% Computational Cost: add. (136->37), mult. (117->50), div. (0->0), fcn. (87->10), ass. (0->19)
t32 = qJ(3) + pkin(9);
t31 = qJ(5) + t32;
t24 = sin(t31);
t25 = cos(t31);
t33 = qJ(1) + pkin(8);
t28 = sin(t33);
t30 = cos(t33);
t41 = g(2) * t28 - g(3) * t30;
t43 = (-g(1) * t25 + t24 * t41) * MDP(21) + (g(1) * t24 + t25 * t41) * MDP(22);
t42 = g(2) * t30 + g(3) * t28;
t38 = cos(qJ(1));
t37 = cos(qJ(3));
t36 = sin(qJ(1));
t35 = sin(qJ(3));
t34 = -qJ(4) - pkin(6);
t29 = cos(t32);
t27 = sin(t32);
t26 = t37 * pkin(3) + pkin(2);
t1 = [(g(2) * t36 - g(3) * t38) * MDP(3) - t41 * MDP(14) + (-g(2) * (t38 * pkin(1) + t30 * t26 - t28 * t34) - g(3) * (t36 * pkin(1) + t28 * t26 + t30 * t34)) * MDP(15) + (MDP(4) * pkin(1) + MDP(2)) * (-g(2) * t38 - g(3) * t36) + (-MDP(10) * t37 + MDP(11) * t35 - MDP(12) * t29 + MDP(13) * t27 - MDP(21) * t25 + MDP(22) * t24) * t42; (-MDP(15) - MDP(4)) * g(1); (g(1) * t35 + t37 * t41) * MDP(11) + (-g(1) * t29 + t27 * t41) * MDP(12) + (g(1) * t27 + t29 * t41) * MDP(13) + t43 + (MDP(15) * pkin(3) + MDP(10)) * (-g(1) * t37 + t35 * t41); t42 * MDP(15); t43;];
taug = t1;
