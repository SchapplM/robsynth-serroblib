% Calculate minimal parameter regressor of potential energy for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% U_reg [1x8]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:48
% EndTime: 2019-03-08 18:15:48
% DurationCPUTime: 0.03s
% Computational Cost: add. (25->14), mult. (32->23), div. (0->0), fcn. (36->6), ass. (0->13)
t55 = g(3) * qJ(1);
t54 = cos(qJ(3));
t53 = sin(qJ(3));
t52 = cos(pkin(6));
t51 = sin(pkin(6));
t50 = qJ(3) + qJ(4);
t49 = cos(t50);
t48 = sin(t50);
t47 = t51 * t54 - t52 * t53;
t46 = -t51 * t53 - t52 * t54;
t45 = -t52 * t48 + t51 * t49;
t44 = -t51 * t48 - t52 * t49;
t1 = [-t55, -g(1) * (t52 * pkin(1) + t51 * qJ(2)) - g(2) * (t51 * pkin(1) - t52 * qJ(2)) - t55, 0, g(1) * t46 - g(2) * t47, -g(1) * t47 - g(2) * t46, 0, g(1) * t44 - g(2) * t45, -g(1) * t45 - g(2) * t44;];
U_reg  = t1;
