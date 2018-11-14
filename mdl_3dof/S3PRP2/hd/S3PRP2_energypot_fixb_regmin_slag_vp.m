% Calculate minimal parameter regressor of potential energy for
% S3PRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% 
% Output:
% U_reg [1x7]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:07
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S3PRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP2_energypot_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP2_energypot_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:07:17
% EndTime: 2018-11-14 10:07:17
% DurationCPUTime: 0.02s
% Computational Cost: add. (11->10), mult. (16->12), div. (0->0), fcn. (12->2), ass. (0->5)
t22 = cos(qJ(2));
t21 = sin(qJ(2));
t20 = g(1) * t22 - g(2) * t21;
t19 = -g(1) * t21 - g(2) * t22;
t1 = [-g(1) * qJ(1), 0, t19, -t20, t19, t20, -g(1) * (t21 * pkin(2) - t22 * qJ(3) + qJ(1)) - g(2) * (t22 * pkin(2) + t21 * qJ(3) + pkin(1)) + g(3) * pkin(3);];
U_reg  = t1;
