% Calculate Gravitation load on the joints for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:16
% EndTime: 2019-12-31 17:45:16
% DurationCPUTime: 0.11s
% Computational Cost: add. (88->31), mult. (84->40), div. (0->0), fcn. (58->8), ass. (0->14)
t41 = -MDP(11) - MDP(7);
t33 = qJ(1) + pkin(7);
t28 = sin(t33);
t30 = cos(t33);
t37 = cos(qJ(1));
t40 = t37 * pkin(1) + t30 * pkin(2) + t28 * qJ(3);
t36 = sin(qJ(1));
t39 = -t36 * pkin(1) + t30 * qJ(3);
t22 = -g(1) * t30 - g(2) * t28;
t21 = g(1) * t28 - g(2) * t30;
t32 = pkin(8) + qJ(5);
t29 = cos(t32);
t27 = sin(t32);
t1 = [(g(1) * t37 + g(2) * t36) * MDP(3) + (-g(1) * (-t28 * pkin(2) + t39) - g(2) * t40) * MDP(7) + (-g(1) * ((-pkin(2) - qJ(4)) * t28 + t39) - g(2) * (t30 * qJ(4) + t40)) * MDP(11) + (-MDP(5) + MDP(10)) * t21 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t36 - g(2) * t37) + (t27 * MDP(17) + t29 * MDP(18) + MDP(8) * sin(pkin(8)) + MDP(9) * cos(pkin(8)) + MDP(6)) * t22; (-MDP(4) + t41) * g(3); t41 * t21; t22 * MDP(11); (g(3) * t27 - t21 * t29) * MDP(17) + (g(3) * t29 + t21 * t27) * MDP(18);];
taug = t1;
