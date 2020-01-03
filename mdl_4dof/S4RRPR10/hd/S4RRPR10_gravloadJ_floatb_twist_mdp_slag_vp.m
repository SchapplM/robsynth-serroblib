% Calculate Gravitation load on the joints for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:55
% EndTime: 2019-12-31 17:11:56
% DurationCPUTime: 0.14s
% Computational Cost: add. (63->36), mult. (146->53), div. (0->0), fcn. (129->6), ass. (0->22)
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t38 = g(1) * t45 + g(2) * t42;
t61 = MDP(9) - MDP(12);
t60 = MDP(10) - MDP(13);
t44 = cos(qJ(2));
t55 = g(3) * t44;
t40 = sin(qJ(4));
t54 = t42 * t40;
t43 = cos(qJ(4));
t53 = t42 * t43;
t52 = t45 * t40;
t51 = t45 * t43;
t41 = sin(qJ(2));
t49 = pkin(2) * t44 + qJ(3) * t41;
t47 = pkin(1) + t49;
t35 = -t41 * t54 + t51;
t34 = t41 * t53 + t52;
t33 = t41 * t52 + t53;
t32 = t41 * t51 - t54;
t30 = t38 * t41 - t55;
t1 = [((-g(1) * pkin(5) - g(2) * t47) * t45 + (-g(2) * pkin(5) + g(1) * t47) * t42) * MDP(14) + (-g(1) * t35 - g(2) * t33) * MDP(20) + (g(1) * t34 - g(2) * t32) * MDP(21) + (MDP(3) - MDP(11)) * t38 + (-t41 * t60 + t44 * t61 + MDP(2)) * (g(1) * t42 - g(2) * t45); (-g(3) * t49 + t38 * (pkin(2) * t41 - qJ(3) * t44)) * MDP(14) + t61 * t30 + (-MDP(20) * t40 - MDP(21) * t43 + t60) * (g(3) * t41 + t38 * t44); -t30 * MDP(14); (-g(1) * t32 - g(2) * t34 + t43 * t55) * MDP(20) + (g(1) * t33 - g(2) * t35 - t40 * t55) * MDP(21);];
taug = t1;
