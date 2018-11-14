% Calculate minimal parameter regressor of potential energy for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:50
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:49:34
% EndTime: 2018-11-14 13:49:34
% DurationCPUTime: 0.05s
% Computational Cost: add. (26->21), mult. (42->28), div. (0->0), fcn. (40->4), ass. (0->13)
t55 = sin(qJ(3));
t56 = sin(qJ(1));
t60 = t56 * t55;
t58 = cos(qJ(1));
t59 = pkin(1) * t58 + qJ(2) * t56;
t57 = cos(qJ(3));
t53 = t56 * pkin(1);
t51 = pkin(3) * t57 + pkin(2);
t50 = -g(1) * t58 - g(2) * t56;
t49 = g(1) * t56 - g(2) * t58;
t48 = -t55 * t58 + t56 * t57;
t47 = -t57 * t58 - t60;
t1 = [0, t50, t49, t50, -t49, -g(1) * t59 - g(2) * (-qJ(2) * t58 + t53) - g(3) * pkin(4), 0, g(1) * t47 - g(2) * t48, -g(1) * t48 - g(2) * t47, -g(1) * (pkin(3) * t60 + t51 * t58 + t59) - g(2) * (t56 * t51 + t53 + (-pkin(3) * t55 - qJ(2)) * t58) - g(3) * (-qJ(4) - pkin(5) + pkin(4));];
U_reg  = t1;
