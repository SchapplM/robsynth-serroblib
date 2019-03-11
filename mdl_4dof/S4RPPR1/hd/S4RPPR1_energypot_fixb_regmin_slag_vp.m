% Calculate minimal parameter regressor of potential energy for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:31
% EndTime: 2019-03-08 18:27:31
% DurationCPUTime: 0.05s
% Computational Cost: add. (36->17), mult. (34->26), div. (0->0), fcn. (32->6), ass. (0->12)
t64 = g(3) * (qJ(2) + pkin(4));
t60 = sin(qJ(1));
t62 = cos(qJ(1));
t63 = -g(1) * t62 - g(2) * t60;
t61 = cos(qJ(4));
t59 = sin(qJ(4));
t57 = qJ(1) + pkin(6);
t56 = cos(t57);
t55 = sin(t57);
t54 = t55 * t61 - t56 * t59;
t53 = -t55 * t59 - t56 * t61;
t1 = [0, t63, g(1) * t60 - g(2) * t62, pkin(1) * t63 - t64, -g(1) * t56 - g(2) * t55, -g(1) * t55 + g(2) * t56, -g(1) * (pkin(1) * t62 + pkin(2) * t56 + qJ(3) * t55) - g(2) * (pkin(1) * t60 + pkin(2) * t55 - qJ(3) * t56) - t64, 0, g(1) * t53 - g(2) * t54, -g(1) * t54 - g(2) * t53;];
U_reg  = t1;
