% Calculate minimal parameter regressor of potential energy for
% S4PRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PRPP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP5_energypot_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:10:22
% EndTime: 2018-11-14 14:10:22
% DurationCPUTime: 0.02s
% Computational Cost: add. (22->15), mult. (29->17), div. (0->0), fcn. (22->2), ass. (0->7)
t37 = sin(qJ(2));
t38 = cos(qJ(2));
t40 = t38 * pkin(2) + t37 * qJ(3) + pkin(1);
t39 = t37 * pkin(2) - t38 * qJ(3) + qJ(1);
t33 = g(1) * t38 - g(2) * t37;
t32 = -g(1) * t37 - g(2) * t38;
t1 = [-g(1) * qJ(1), 0, t32, -t33, t32, t33, g(3) * pkin(4) - g(1) * t39 - g(2) * t40, t32, t33, -g(1) * (t37 * pkin(3) + t39) - g(2) * (t38 * pkin(3) + t40) - g(3) * (qJ(4) - pkin(4));];
U_reg  = t1;
