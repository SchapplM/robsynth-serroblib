% Calculate Gravitation load on the joints for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:40
% EndTime: 2019-12-31 20:15:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (135->31), mult. (125->41), div. (0->0), fcn. (94->8), ass. (0->16)
t46 = qJ(1) + qJ(2);
t42 = sin(t46);
t44 = cos(t46);
t35 = g(1) * t42 - g(2) * t44;
t45 = qJ(4) + qJ(5);
t41 = sin(t45);
t43 = cos(t45);
t54 = (g(3) * t41 - t35 * t43) * MDP(22) + (g(3) * t43 + t35 * t41) * MDP(23);
t53 = t44 * pkin(2) + t42 * qJ(3);
t52 = -t42 * pkin(2) + t44 * qJ(3);
t47 = sin(qJ(4));
t49 = cos(qJ(4));
t51 = (MDP(5) - MDP(7)) * t35 + (-t47 * MDP(15) - t49 * MDP(16) - t41 * MDP(22) - t43 * MDP(23) + MDP(6) - MDP(8)) * (g(1) * t44 + g(2) * t42);
t50 = cos(qJ(1));
t48 = sin(qJ(1));
t1 = [(g(1) * t48 - g(2) * t50) * MDP(2) + (g(1) * t50 + g(2) * t48) * MDP(3) + (-g(1) * (-t48 * pkin(1) + t52) - g(2) * (t50 * pkin(1) + t53)) * MDP(9) + t51; (-g(1) * t52 - g(2) * t53) * MDP(9) + t51; -t35 * MDP(9); (g(3) * t47 - t35 * t49) * MDP(15) + (g(3) * t49 + t35 * t47) * MDP(16) + t54; t54;];
taug = t1;
