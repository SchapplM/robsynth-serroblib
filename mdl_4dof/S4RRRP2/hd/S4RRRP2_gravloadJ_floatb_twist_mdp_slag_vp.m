% Calculate Gravitation load on the joints for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:09
% EndTime: 2019-12-31 17:13:09
% DurationCPUTime: 0.06s
% Computational Cost: add. (85->26), mult. (88->35), div. (0->0), fcn. (63->6), ass. (0->15)
t33 = cos(qJ(3));
t26 = t33 * pkin(3) + pkin(2);
t29 = qJ(1) + qJ(2);
t27 = sin(t29);
t28 = cos(t29);
t30 = -qJ(4) - pkin(6);
t38 = t28 * t26 - t27 * t30;
t23 = g(1) * t27 - g(2) * t28;
t24 = g(1) * t28 + g(2) * t27;
t31 = sin(qJ(3));
t37 = (-MDP(14) + MDP(6)) * t24 + (t33 * MDP(12) - t31 * MDP(13) + MDP(5)) * t23;
t36 = -t27 * t26 - t28 * t30;
t34 = cos(qJ(1));
t32 = sin(qJ(1));
t1 = [(g(1) * t32 - g(2) * t34) * MDP(2) + (g(1) * t34 + g(2) * t32) * MDP(3) + (-g(1) * (-t32 * pkin(1) + t36) - g(2) * (t34 * pkin(1) + t38)) * MDP(15) + t37; (-g(1) * t36 - g(2) * t38) * MDP(15) + t37; (g(3) * t31 + t24 * t33) * MDP(13) + (MDP(15) * pkin(3) + MDP(12)) * (-g(3) * t33 + t24 * t31); -t23 * MDP(15);];
taug = t1;
