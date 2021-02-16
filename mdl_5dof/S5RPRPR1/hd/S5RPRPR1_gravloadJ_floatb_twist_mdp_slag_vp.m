% Calculate Gravitation load on the joints for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:33:52
% EndTime: 2021-01-15 11:33:53
% DurationCPUTime: 0.14s
% Computational Cost: add. (106->38), mult. (127->48), div. (0->0), fcn. (95->8), ass. (0->16)
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t21 = g(1) * t33 - g(2) * t35;
t31 = qJ(3) + pkin(8);
t28 = qJ(5) + t31;
t23 = sin(t28);
t24 = cos(t28);
t37 = (g(3) * t23 - t21 * t24) * MDP(23) + (g(3) * t24 + t21 * t23) * MDP(24);
t22 = g(1) * t35 + g(2) * t33;
t34 = cos(qJ(3));
t32 = sin(qJ(3));
t30 = pkin(1) + pkin(6) + qJ(4);
t27 = cos(t31);
t26 = sin(t31);
t25 = pkin(3) * t32 + qJ(2);
t1 = [(-g(1) * (-pkin(1) * t33 + qJ(2) * t35) - g(2) * (pkin(1) * t35 + qJ(2) * t33)) * MDP(6) + (-g(1) * (t25 * t35 - t30 * t33) - g(2) * (t25 * t33 + t30 * t35)) * MDP(17) + (MDP(2) - MDP(4) + MDP(16)) * t21 + (-MDP(12) * t32 - t34 * MDP(13) - t26 * MDP(14) - t27 * MDP(15) - MDP(23) * t23 - MDP(24) * t24 + MDP(3) - MDP(5)) * t22; (-MDP(17) - MDP(6)) * t21; (g(3) * t34 + t21 * t32) * MDP(13) + (g(3) * t26 - t21 * t27) * MDP(14) + (g(3) * t27 + t21 * t26) * MDP(15) + t37 + (MDP(17) * pkin(3) + MDP(12)) * (g(3) * t32 - t21 * t34); -t22 * MDP(17); t37;];
taug = t1;
