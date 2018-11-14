% Calculate minimal parameter regressor of potential energy for
% S4RRPP2
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RRPP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:30
% EndTime: 2018-11-14 13:52:30
% DurationCPUTime: 0.03s
% Computational Cost: add. (47->19), mult. (36->22), div. (0->0), fcn. (30->4), ass. (0->11)
t62 = pkin(5) + pkin(4);
t57 = qJ(1) + qJ(2);
t53 = sin(t57);
t54 = cos(t57);
t59 = cos(qJ(1));
t61 = t59 * pkin(1) + t54 * pkin(2) + t53 * qJ(3);
t58 = sin(qJ(1));
t60 = t58 * pkin(1) + t53 * pkin(2) - t54 * qJ(3);
t48 = -g(1) * t54 - g(2) * t53;
t47 = g(1) * t53 - g(2) * t54;
t1 = [0, -g(1) * t59 - g(2) * t58, g(1) * t58 - g(2) * t59, 0, t48, t47, t48, -t47, -g(1) * t61 - g(2) * t60 - g(3) * t62, t48, -t47, -g(1) * (t54 * pkin(3) + t61) - g(2) * (t53 * pkin(3) + t60) - g(3) * (-qJ(4) + t62);];
U_reg  = t1;
