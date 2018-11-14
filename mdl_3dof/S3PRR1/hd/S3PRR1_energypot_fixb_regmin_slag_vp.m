% Calculate minimal parameter regressor of potential energy for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% 
% Output:
% U_reg [1x7]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S3PRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_energypot_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_energypot_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:57
% EndTime: 2018-11-14 10:11:57
% DurationCPUTime: 0.02s
% Computational Cost: add. (9->6), mult. (9->9), div. (0->0), fcn. (8->4), ass. (0->6)
t30 = cos(qJ(2));
t29 = sin(qJ(2));
t28 = qJ(2) + qJ(3);
t27 = cos(t28);
t26 = sin(t28);
t1 = [-g(2) * qJ(1), 0, -g(1) * t30 - g(2) * t29, g(1) * t29 - g(2) * t30, 0, -g(1) * t27 - g(2) * t26, g(1) * t26 - g(2) * t27;];
U_reg  = t1;
