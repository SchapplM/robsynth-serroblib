% Calculate potential energy for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S3RRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_energypot_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_energypot_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_energypot_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRP1_energypot_fixb_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:06
% EndTime: 2018-11-14 10:15:06
% DurationCPUTime: 0.05s
% Computational Cost: add. (42->28), mult. (36->32), div. (0->0), fcn. (20->4), ass. (0->11)
t23 = pkin(3) + pkin(4);
t22 = rSges(4,1) + pkin(2);
t21 = rSges(4,3) + qJ(3);
t20 = cos(qJ(1));
t19 = sin(qJ(1));
t18 = qJ(1) + qJ(2);
t17 = t20 * pkin(1);
t16 = t19 * pkin(1);
t15 = cos(t18);
t14 = sin(t18);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t20 * rSges(2,1) - t19 * rSges(2,2)) + g(2) * (t19 * rSges(2,1) + t20 * rSges(2,2)) + g(3) * (pkin(3) + rSges(2,3))) - m(3) * (g(1) * (t15 * rSges(3,1) - t14 * rSges(3,2) + t17) + g(2) * (t14 * rSges(3,1) + t15 * rSges(3,2) + t16) + g(3) * (rSges(3,3) + t23)) - m(4) * (g(1) * t17 + g(2) * t16 + g(3) * (rSges(4,2) + t23) + (g(1) * t22 - g(2) * t21) * t15 + (g(1) * t21 + g(2) * t22) * t14);
U  = t1;
