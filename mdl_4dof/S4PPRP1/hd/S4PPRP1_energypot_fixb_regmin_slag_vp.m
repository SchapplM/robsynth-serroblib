% Calculate minimal parameter regressor of potential energy for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% 
% Output:
% U_reg [1x8]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:39
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:39:14
% EndTime: 2018-11-14 13:39:14
% DurationCPUTime: 0.03s
% Computational Cost: add. (32->19), mult. (53->24), div. (0->0), fcn. (58->4), ass. (0->12)
t56 = cos(qJ(3));
t55 = sin(qJ(3));
t54 = g(3) * qJ(1);
t49 = cos(pkin(5));
t52 = sin(pkin(5));
t53 = t49 * pkin(1) + t52 * qJ(2);
t51 = t52 * pkin(1) - t49 * qJ(2);
t41 = -t49 * t56 - t52 * t55;
t42 = t49 * t55 - t52 * t56;
t50 = g(1) * t42 - g(2) * t41;
t40 = g(1) * t41 + g(2) * t42;
t1 = [-t54, -g(1) * t53 - g(2) * t51 - t54, 0, t40, t50, t40, -t50, -g(1) * (t49 * pkin(2) - t41 * pkin(3) + t42 * qJ(4) + t53) - g(2) * (t52 * pkin(2) - t42 * pkin(3) - t41 * qJ(4) + t51) - g(3) * (-pkin(4) + qJ(1));];
U_reg  = t1;
