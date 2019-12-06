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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:43
% EndTime: 2019-12-05 17:47:43
% DurationCPUTime: 0.09s
% Computational Cost: add. (86->33), mult. (109->40), div. (0->0), fcn. (81->6), ass. (0->15)
t32 = sin(qJ(3));
t39 = pkin(3) * t32;
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t22 = g(1) * t33 - g(2) * t35;
t26 = qJ(3) + pkin(8) + qJ(5);
t24 = sin(t26);
t25 = cos(t26);
t38 = (g(3) * t24 - t22 * t25) * MDP(21) + (g(3) * t25 + t22 * t24) * MDP(22);
t37 = t35 * pkin(1) + t33 * qJ(2);
t23 = g(1) * t35 + g(2) * t33;
t34 = cos(qJ(3));
t31 = -qJ(4) - pkin(6);
t28 = t35 * qJ(2);
t1 = [(-g(1) * (-t33 * pkin(1) + t28) - g(2) * t37) * MDP(6) + (-g(1) * (t35 * t39 + t28 + (-pkin(1) + t31) * t33) - g(2) * (-t35 * t31 + t33 * t39 + t37)) * MDP(15) + (MDP(2) - MDP(4) + MDP(14)) * t22 + (-t32 * MDP(12) - t34 * MDP(13) - MDP(21) * t24 - MDP(22) * t25 + MDP(3) - MDP(5)) * t23; (-MDP(15) - MDP(6)) * t22; (g(3) * t34 + t22 * t32) * MDP(13) + t38 + (MDP(15) * pkin(3) + MDP(12)) * (g(3) * t32 - t22 * t34); -t23 * MDP(15); t38;];
taug = t1;
