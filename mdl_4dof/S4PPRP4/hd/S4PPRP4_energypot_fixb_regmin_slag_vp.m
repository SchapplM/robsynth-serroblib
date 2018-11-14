% Calculate minimal parameter regressor of potential energy for
% S4PPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% 
% Output:
% U_reg [1x6]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:06
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PPRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP4_energypot_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:06:13
% EndTime: 2018-11-14 14:06:13
% DurationCPUTime: 0.02s
% Computational Cost: add. (12->12), mult. (13->12), div. (0->0), fcn. (6->2), ass. (0->4)
t28 = g(1) * qJ(1);
t27 = cos(qJ(3));
t26 = sin(qJ(3));
t1 = [-t28, g(3) * pkin(1) + g(2) * qJ(2) - t28, 0, -g(1) * t27 + g(2) * t26, g(1) * t26 + g(2) * t27, -g(1) * (t27 * pkin(3) + pkin(2) + qJ(1)) - g(2) * (-t26 * pkin(3) - qJ(2)) - g(3) * (-qJ(4) - pkin(4) - pkin(1));];
U_reg  = t1;
