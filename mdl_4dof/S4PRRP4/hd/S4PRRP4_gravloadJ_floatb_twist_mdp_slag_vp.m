% Calculate Gravitation load on the joints for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:58
% EndTime: 2019-12-31 16:27:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (77->23), mult. (90->30), div. (0->0), fcn. (69->4), ass. (0->12)
t24 = pkin(6) + qJ(2);
t22 = sin(t24);
t23 = cos(t24);
t21 = g(1) * t23 + g(2) * t22;
t34 = MDP(10) + MDP(12);
t33 = MDP(11) - MDP(14);
t25 = sin(qJ(3));
t26 = cos(qJ(3));
t29 = t26 * pkin(3) + t25 * qJ(4);
t27 = pkin(2) + t29;
t17 = -g(3) * t26 + t21 * t25;
t1 = [(-MDP(1) - MDP(15)) * g(3); ((-g(1) * pkin(5) - g(2) * t27) * t23 + (-g(2) * pkin(5) + g(1) * t27) * t22) * MDP(15) + (MDP(4) - MDP(13)) * t21 + (-t33 * t25 + t34 * t26 + MDP(3)) * (g(1) * t22 - g(2) * t23); (-g(3) * t29 + t21 * (pkin(3) * t25 - qJ(4) * t26)) * MDP(15) + t33 * (g(3) * t25 + t21 * t26) + t34 * t17; -t17 * MDP(15);];
taug = t1;
