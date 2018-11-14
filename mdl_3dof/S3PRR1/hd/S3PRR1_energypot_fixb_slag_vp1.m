% Calculate potential energy for
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
% Datum: 2018-11-14 10:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S3PRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_energypot_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_energypot_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_energypot_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRR1_energypot_fixb_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:56
% EndTime: 2018-11-14 10:11:56
% DurationCPUTime: 0.03s
% Computational Cost: add. (29->26), mult. (26->26), div. (0->0), fcn. (10->4), ass. (0->6)
t13 = cos(qJ(2));
t12 = sin(qJ(2));
t11 = qJ(2) + qJ(3);
t10 = cos(t11);
t9 = sin(t11);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * rSges(2,1) + g(2) * (qJ(1) + rSges(2,3)) - g(3) * rSges(2,2)) - m(3) * (g(1) * (t13 * rSges(3,1) - t12 * rSges(3,2) + pkin(1)) + g(2) * (t12 * rSges(3,1) + t13 * rSges(3,2) + qJ(1)) + g(3) * (pkin(3) + rSges(3,3))) - m(4) * (g(1) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t13 * pkin(2) + pkin(1)) + g(2) * (t9 * rSges(4,1) + t10 * rSges(4,2) + t12 * pkin(2) + qJ(1)) + g(3) * (pkin(4) + pkin(3) + rSges(4,3)));
U  = t1;
