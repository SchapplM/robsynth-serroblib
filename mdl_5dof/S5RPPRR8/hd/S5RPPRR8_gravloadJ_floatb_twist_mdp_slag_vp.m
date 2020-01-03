% Calculate Gravitation load on the joints for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:14
% EndTime: 2019-12-31 18:01:14
% DurationCPUTime: 0.13s
% Computational Cost: add. (108->32), mult. (133->48), div. (0->0), fcn. (132->8), ass. (0->19)
t46 = cos(qJ(1));
t52 = pkin(8) + qJ(4);
t50 = sin(t52);
t51 = cos(t52);
t56 = sin(qJ(1));
t29 = -t46 * t51 - t56 * t50;
t30 = t46 * t50 - t56 * t51;
t44 = sin(qJ(5));
t45 = cos(qJ(5));
t47 = g(1) * t29 + g(2) * t30;
t59 = (t45 * MDP(18) - t44 * MDP(19) + MDP(11)) * (g(1) * t30 - g(2) * t29) - t47 * MDP(12);
t55 = t46 * pkin(1) + t56 * qJ(2);
t49 = -t56 * pkin(1) + t46 * qJ(2);
t43 = cos(pkin(8));
t42 = sin(pkin(8));
t33 = g(1) * t56 - g(2) * t46;
t32 = t56 * t42 + t46 * t43;
t31 = t46 * t42 - t56 * t43;
t1 = [(-g(1) * t49 - g(2) * t55) * MDP(6) + (-g(1) * t31 - g(2) * t32) * MDP(7) + (-g(1) * t32 + g(2) * t31) * MDP(8) + (-g(1) * (-t56 * pkin(2) + t49) - g(2) * (t46 * pkin(2) + t55)) * MDP(9) + (MDP(3) - MDP(5)) * (g(1) * t46 + g(2) * t56) + (MDP(2) + MDP(4)) * t33 - t59; (-MDP(6) - MDP(9)) * t33; g(3) * MDP(9); t59; (g(3) * t45 - t47 * t44) * MDP(18) + (-g(3) * t44 - t47 * t45) * MDP(19);];
taug = t1;
