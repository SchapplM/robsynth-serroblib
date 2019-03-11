% Calculate minimal parameter regressor of potential energy for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:20
% EndTime: 2019-03-08 18:25:20
% DurationCPUTime: 0.04s
% Computational Cost: add. (31->10), mult. (13->13), div. (0->0), fcn. (12->6), ass. (0->10)
t55 = pkin(7) + qJ(2);
t54 = qJ(3) + t55;
t53 = cos(t55);
t52 = sin(t55);
t51 = qJ(4) + t54;
t50 = cos(t54);
t49 = sin(t54);
t48 = cos(t51);
t47 = sin(t51);
t1 = [-g(3) * qJ(1), 0, -g(1) * t53 - g(2) * t52, g(1) * t52 - g(2) * t53, 0, -g(1) * t50 - g(2) * t49, g(1) * t49 - g(2) * t50, 0, -g(1) * t48 - g(2) * t47, g(1) * t47 - g(2) * t48;];
U_reg  = t1;
