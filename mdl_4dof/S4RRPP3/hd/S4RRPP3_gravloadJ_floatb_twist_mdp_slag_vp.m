% Calculate Gravitation load on the joints for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRPP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:51
% EndTime: 2019-12-31 16:57:51
% DurationCPUTime: 0.12s
% Computational Cost: add. (86->38), mult. (122->49), div. (0->0), fcn. (93->6), ass. (0->16)
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t23 = g(1) * t35 + g(2) * t33;
t22 = g(1) * t33 - g(2) * t35;
t30 = qJ(2) + pkin(6);
t26 = sin(t30);
t27 = cos(t30);
t38 = t27 * pkin(3) + t26 * qJ(4);
t34 = cos(qJ(2));
t32 = sin(qJ(2));
t31 = -qJ(3) - pkin(5);
t28 = t34 * pkin(2);
t25 = t28 + pkin(1);
t24 = t35 * t25;
t21 = -g(3) * t27 + t23 * t26;
t1 = [(-g(1) * (-t33 * t25 - t35 * t31) - g(2) * (-t33 * t31 + t24)) * MDP(12) + (-g(2) * t24 + (g(1) * t31 - g(2) * t38) * t35 + (-g(1) * (-t25 - t38) + g(2) * t31) * t33) * MDP(16) + (MDP(3) - MDP(11) - MDP(14)) * t23 + (-t32 * MDP(10) + t27 * MDP(13) + t26 * MDP(15) + t34 * MDP(9) + MDP(2)) * t22; (g(3) * t32 + t23 * t34) * MDP(10) + t21 * MDP(13) + (-g(3) * t26 - t23 * t27) * MDP(15) + (-g(3) * (t28 + t38) + t23 * (pkin(2) * t32 + pkin(3) * t26 - qJ(4) * t27)) * MDP(16) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t34 + t23 * t32); (-MDP(12) - MDP(16)) * t22; -t21 * MDP(16);];
taug = t1;
