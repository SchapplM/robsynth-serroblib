% Calculate minimal parameter regressor of potential energy for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:20
% EndTime: 2021-01-15 10:46:20
% DurationCPUTime: 0.05s
% Computational Cost: add. (42->19), mult. (46->25), div. (0->0), fcn. (43->8), ass. (0->15)
t53 = qJ(2) + pkin(7);
t56 = sin(qJ(1));
t58 = cos(qJ(1));
t59 = g(1) * t58 + g(2) * t56;
t57 = cos(qJ(2));
t55 = sin(qJ(2));
t54 = pkin(5) + qJ(3);
t52 = qJ(4) + t53;
t51 = cos(t53);
t50 = sin(t53);
t49 = t57 * pkin(2) + pkin(1);
t48 = cos(t52);
t47 = sin(t52);
t46 = g(1) * t56 - g(2) * t58;
t1 = [0, -t59, t46, 0, 0, 0, 0, 0, -g(3) * t55 - t59 * t57, -g(3) * t57 + t59 * t55, -g(3) * t50 - t59 * t51, -g(3) * t51 + t59 * t50, -t46, -g(1) * (t58 * t49 + t54 * t56) - g(2) * (t56 * t49 - t58 * t54) - g(3) * (t55 * pkin(2) + pkin(4)), 0, 0, 0, 0, 0, -g(3) * t47 - t59 * t48, -g(3) * t48 + t59 * t47;];
U_reg = t1;
