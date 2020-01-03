% Calculate Gravitation load on the joints for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:24
% EndTime: 2019-12-31 17:32:24
% DurationCPUTime: 0.06s
% Computational Cost: add. (62->22), mult. (108->32), div. (0->0), fcn. (114->8), ass. (0->13)
t35 = cos(qJ(3));
t34 = sin(qJ(3));
t33 = cos(pkin(7));
t32 = sin(pkin(7));
t31 = MDP(2) + MDP(9);
t18 = -t32 * t34 - t33 * t35;
t19 = -t32 * t35 + t33 * t34;
t30 = g(1) * t19 - g(2) * t18;
t29 = g(1) * t18 + g(2) * t19;
t26 = pkin(8) + qJ(5);
t25 = cos(t26);
t24 = sin(t26);
t1 = [(-MDP(1) - t31) * g(3); t31 * (-g(1) * t32 + g(2) * t33); (-g(1) * (-t19 * pkin(3) - t18 * qJ(4)) - g(2) * (t18 * pkin(3) - t19 * qJ(4))) * MDP(9) + (-MDP(5) + MDP(8)) * t29 + (t25 * MDP(15) - t24 * MDP(16) + MDP(6) * cos(pkin(8)) - MDP(7) * sin(pkin(8)) + MDP(4)) * t30; -t30 * MDP(9); (g(3) * t25 - t29 * t24) * MDP(15) + (-g(3) * t24 - t29 * t25) * MDP(16);];
taug = t1;
