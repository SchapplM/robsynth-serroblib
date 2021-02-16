% Calculate Gravitation load on the joints for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:29
% EndTime: 2021-01-15 14:56:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (71->24), mult. (143->32), div. (0->0), fcn. (149->6), ass. (0->16)
t44 = cos(qJ(3));
t43 = sin(qJ(3));
t42 = cos(pkin(7));
t41 = sin(pkin(7));
t40 = MDP(16) + MDP(2);
t39 = MDP(11) + MDP(13);
t38 = MDP(12) + MDP(14);
t25 = -t41 * t43 - t42 * t44;
t26 = -t41 * t44 + t42 * t43;
t37 = g(1) * t26 - g(2) * t25;
t36 = g(1) * t25 + g(2) * t26;
t34 = cos(qJ(4));
t33 = sin(qJ(4));
t32 = qJ(5) + pkin(6);
t31 = t34 * pkin(4) + pkin(3);
t1 = [(-MDP(1) - t40) * g(3); t40 * (-g(1) * t41 + g(2) * t42); (-g(1) * (-t25 * t32 - t26 * t31) - g(2) * (t25 * t31 - t26 * t32)) * MDP(16) + (-MDP(5) + MDP(15)) * t36 + (-t38 * t33 + t39 * t34 + MDP(4)) * t37; t38 * (-g(3) * t33 - t36 * t34) + (MDP(16) * pkin(4) + t39) * (g(3) * t34 - t36 * t33); -t37 * MDP(16);];
taug = t1;
