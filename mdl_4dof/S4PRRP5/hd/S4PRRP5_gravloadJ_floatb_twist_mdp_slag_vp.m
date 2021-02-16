% Calculate Gravitation load on the joints for
% S4PRRP5
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
%   see S4PRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:36
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:36:09
% EndTime: 2021-01-14 22:36:10
% DurationCPUTime: 0.12s
% Computational Cost: add. (59->27), mult. (139->39), div. (0->0), fcn. (125->6), ass. (0->16)
t47 = MDP(10) + MDP(12);
t46 = MDP(11) + MDP(13);
t31 = sin(pkin(6));
t32 = cos(pkin(6));
t39 = g(1) * t32 + g(2) * t31;
t35 = sin(qJ(2));
t37 = cos(qJ(2));
t28 = -g(3) * t37 + t39 * t35;
t43 = g(3) * t35;
t34 = sin(qJ(3));
t41 = t34 * t37;
t36 = cos(qJ(3));
t40 = t36 * t37;
t33 = qJ(4) + pkin(5);
t30 = t36 * pkin(3) + pkin(2);
t1 = [(-MDP(1) - MDP(15)) * g(3); (-g(3) * (t37 * t30 + t35 * t33) + t39 * (t30 * t35 - t33 * t37)) * MDP(15) + (MDP(4) - MDP(14)) * (t39 * t37 + t43) + (-t46 * t34 + t47 * t36 + MDP(3)) * t28; t46 * (-g(1) * (-t31 * t34 - t32 * t40) - g(2) * (-t31 * t40 + t32 * t34) + t36 * t43) + (MDP(15) * pkin(3) + t47) * (-g(1) * (t31 * t36 - t32 * t41) - g(2) * (-t31 * t41 - t32 * t36) + t34 * t43); -t28 * MDP(15);];
taug = t1;
