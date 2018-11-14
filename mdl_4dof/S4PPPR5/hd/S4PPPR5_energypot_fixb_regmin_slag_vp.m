% Calculate minimal parameter regressor of potential energy for
% S4PPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
% 
% Output:
% U_reg [1x6]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:05
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PPPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR5_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:05:11
% EndTime: 2018-11-14 14:05:11
% DurationCPUTime: 0.02s
% Computational Cost: add. (15->13), mult. (23->17), div. (0->0), fcn. (20->4), ass. (0->9)
t47 = g(1) * qJ(1);
t46 = cos(qJ(4));
t45 = sin(qJ(4));
t44 = g(3) * qJ(2);
t43 = cos(pkin(5));
t42 = sin(pkin(5));
t41 = t42 * t46 - t43 * t45;
t40 = -t42 * t45 - t43 * t46;
t1 = [-t47, -g(2) * pkin(1) + t44 - t47, -g(1) * (t42 * pkin(2) - t43 * qJ(3) + qJ(1)) - g(2) * (t43 * pkin(2) + t42 * qJ(3) + pkin(1)) + t44, 0, -g(1) * t41 + g(2) * t40, -g(1) * t40 - g(2) * t41;];
U_reg  = t1;
