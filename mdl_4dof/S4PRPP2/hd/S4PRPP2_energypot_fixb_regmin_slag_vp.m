% Calculate minimal parameter regressor of potential energy for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% 
% Output:
% U_reg [1x8]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:18:59
% EndTime: 2019-03-08 18:18:59
% DurationCPUTime: 0.02s
% Computational Cost: add. (27->17), mult. (23->20), div. (0->0), fcn. (16->4), ass. (0->9)
t58 = g(3) * (qJ(3) + pkin(4));
t55 = cos(qJ(2));
t57 = t55 * pkin(2) + pkin(1);
t54 = sin(qJ(2));
t56 = t54 * pkin(2) + qJ(1);
t52 = qJ(2) + pkin(5);
t49 = cos(t52);
t48 = sin(t52);
t1 = [-g(2) * qJ(1), 0, -g(1) * t55 - g(2) * t54, g(1) * t54 - g(2) * t55, -g(1) * t57 - g(2) * t56 - t58, -g(1) * t49 - g(2) * t48, -g(1) * t48 + g(2) * t49, -g(1) * (t49 * pkin(3) + t48 * qJ(4) + t57) - g(2) * (t48 * pkin(3) - t49 * qJ(4) + t56) - t58;];
U_reg  = t1;
