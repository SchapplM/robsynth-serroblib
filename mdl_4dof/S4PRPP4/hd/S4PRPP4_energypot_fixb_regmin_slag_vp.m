% Calculate minimal parameter regressor of potential energy for
% S4PRPP4
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:10
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PRPP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:29
% EndTime: 2018-11-14 14:09:29
% DurationCPUTime: 0.03s
% Computational Cost: add. (27->17), mult. (23->20), div. (0->0), fcn. (16->4), ass. (0->9)
t58 = g(3) * (-qJ(3) - pkin(4));
t55 = cos(qJ(2));
t57 = t55 * pkin(2) + pkin(1);
t54 = sin(qJ(2));
t56 = t54 * pkin(2) + qJ(1);
t52 = qJ(2) + pkin(5);
t49 = cos(t52);
t48 = sin(t52);
t1 = [-g(1) * qJ(1), 0, -g(1) * t54 - g(2) * t55, -g(1) * t55 + g(2) * t54, -g(1) * t56 - g(2) * t57 - t58, -g(1) * t48 - g(2) * t49, g(1) * t49 - g(2) * t48, -g(1) * (t48 * pkin(3) - t49 * qJ(4) + t56) - g(2) * (t49 * pkin(3) + t48 * qJ(4) + t57) - t58;];
U_reg  = t1;
