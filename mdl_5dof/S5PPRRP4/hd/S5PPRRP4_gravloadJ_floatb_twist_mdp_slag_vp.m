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
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:41
% EndTime: 2019-12-31 17:34:42
% DurationCPUTime: 0.08s
% Computational Cost: add. (53->22), mult. (105->32), div. (0->0), fcn. (107->6), ass. (0->14)
t35 = cos(qJ(3));
t34 = sin(qJ(3));
t33 = cos(pkin(7));
t32 = sin(pkin(7));
t31 = MDP(14) + MDP(2);
t18 = -t32 * t34 - t33 * t35;
t19 = -t32 * t35 + t33 * t34;
t30 = g(1) * t19 - g(2) * t18;
t29 = g(1) * t18 + g(2) * t19;
t27 = cos(qJ(4));
t26 = sin(qJ(4));
t25 = -qJ(5) - pkin(6);
t24 = t27 * pkin(4) + pkin(3);
t1 = [(-MDP(1) - t31) * g(3); t31 * (-g(1) * t32 + g(2) * t33); (-g(1) * (t18 * t25 - t19 * t24) - g(2) * (t18 * t24 + t19 * t25)) * MDP(14) + (-MDP(5) + MDP(13)) * t29 + (t27 * MDP(11) - t26 * MDP(12) + MDP(4)) * t30; (-g(3) * t26 - t29 * t27) * MDP(12) + (MDP(14) * pkin(4) + MDP(11)) * (g(3) * t27 - t29 * t26); -t30 * MDP(14);];
taug = t1;
