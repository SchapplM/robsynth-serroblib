% Calculate minimal parameter regressor of potential energy for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% U_reg [1x12]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:34
% EndTime: 2019-07-18 18:16:34
% DurationCPUTime: 0.05s
% Computational Cost: add. (39->17), mult. (33->25), div. (0->0), fcn. (34->6), ass. (0->12)
t63 = cos(qJ(1));
t62 = cos(qJ(4));
t61 = sin(qJ(1));
t60 = sin(qJ(4));
t59 = qJ(1) + qJ(2);
t58 = cos(t59);
t57 = sin(t59);
t56 = -g(1) * t58 - g(2) * t57;
t55 = g(1) * t57 - g(2) * t58;
t54 = t57 * t62 - t58 * t60;
t53 = -t57 * t60 - t58 * t62;
t1 = [0, -g(1) * t63 - g(2) * t61, g(1) * t61 - g(2) * t63, 0, t56, t55, t56, -t55, -g(1) * (t63 * pkin(1) + t58 * pkin(2) + t57 * qJ(3)) - g(2) * (t61 * pkin(1) + t57 * pkin(2) - t58 * qJ(3)) - g(3) * (pkin(5) + pkin(4)), 0, g(1) * t53 - g(2) * t54, -g(1) * t54 - g(2) * t53;];
U_reg  = t1;
