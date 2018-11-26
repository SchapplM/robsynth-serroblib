% Calculate minimal parameter regressor of potential energy for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S2RR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_energypot_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_energypot_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:48:47
% EndTime: 2018-11-16 16:48:47
% DurationCPUTime: 0.02s
% Computational Cost: add. (6->5), mult. (14->8), div. (0->0), fcn. (14->4), ass. (0->6)
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t29 = g(1) * t28 - g(3) * t26;
t27 = cos(qJ(2));
t25 = sin(qJ(2));
t1 = [0, -t29, g(1) * t26 + g(3) * t28, 0, 0, 0, 0, 0, -g(2) * t25 - t29 * t27, -g(2) * t27 + t29 * t25;];
U_reg  = t1;
