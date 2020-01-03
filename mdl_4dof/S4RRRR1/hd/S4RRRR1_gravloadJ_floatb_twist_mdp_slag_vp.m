% Calculate Gravitation load on the joints for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:14
% EndTime: 2019-12-31 17:22:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (102->18), mult. (78->26), div. (0->0), fcn. (58->8), ass. (0->14)
t28 = qJ(1) + qJ(2);
t27 = qJ(3) + t28;
t23 = sin(t27);
t24 = cos(t27);
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t35 = g(1) * t24 + g(2) * t23;
t36 = t35 * MDP(9) + (t31 * MDP(15) - t29 * MDP(16) + MDP(8)) * (g(1) * t23 - g(2) * t24);
t25 = sin(t28);
t26 = cos(t28);
t33 = (g(1) * t25 - g(2) * t26) * MDP(5) + (g(1) * t26 + g(2) * t25) * MDP(6) + t36;
t32 = cos(qJ(1));
t30 = sin(qJ(1));
t1 = [(g(1) * t30 - g(2) * t32) * MDP(2) + (g(1) * t32 + g(2) * t30) * MDP(3) + t33; t33; t36; (-g(3) * t31 + t35 * t29) * MDP(15) + (g(3) * t29 + t35 * t31) * MDP(16);];
taug = t1;
