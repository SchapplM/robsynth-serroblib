% Calculate Gravitation load on the joints for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:04:41
% EndTime: 2021-01-15 17:04:42
% DurationCPUTime: 0.12s
% Computational Cost: add. (99->34), mult. (109->42), div. (0->0), fcn. (77->6), ass. (0->17)
t39 = sin(qJ(4));
t50 = pkin(4) * t39;
t49 = -MDP(18) - MDP(7);
t48 = MDP(13) + MDP(15);
t47 = MDP(14) + MDP(16);
t37 = qJ(1) + pkin(7);
t34 = sin(t37);
t35 = cos(t37);
t42 = cos(qJ(1));
t46 = t42 * pkin(1) + t35 * pkin(2) + t34 * qJ(3);
t40 = sin(qJ(1));
t45 = -t40 * pkin(1) + t35 * qJ(3);
t29 = -g(1) * t35 - g(2) * t34;
t28 = g(1) * t34 - g(2) * t35;
t41 = cos(qJ(4));
t38 = -qJ(5) - pkin(6);
t1 = [(g(1) * t42 + g(2) * t40) * MDP(3) + (-g(1) * (-t34 * pkin(2) + t45) - g(2) * t46) * MDP(7) + (-g(1) * (t35 * t50 + (-pkin(2) + t38) * t34 + t45) - g(2) * (t34 * t50 - t35 * t38 + t46)) * MDP(18) + (-MDP(5) + MDP(17)) * t28 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t40 - g(2) * t42) + (t48 * t39 + t47 * t41 + MDP(6)) * t29; (-MDP(4) + t49) * g(3); t49 * t28; t47 * (g(3) * t41 + t28 * t39) + (MDP(18) * pkin(4) + t48) * (g(3) * t39 - t28 * t41); t29 * MDP(18);];
taug = t1;
