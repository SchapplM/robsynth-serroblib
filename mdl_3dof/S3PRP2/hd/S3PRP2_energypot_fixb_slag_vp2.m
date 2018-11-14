% Calculate potential energy for
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
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:07
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S3PRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP2_energypot_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP2_energypot_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP2_energypot_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRP2_energypot_fixb_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:07:06
% EndTime: 2018-11-14 10:07:06
% DurationCPUTime: 0.03s
% Computational Cost: add. (26->20), mult. (31->12), div. (0->0), fcn. (12->2), ass. (0->6)
t14 = m(3) + m(4);
t13 = -m(4) * pkin(2) - mrSges(3,1) - mrSges(4,1);
t12 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t11 = cos(qJ(2));
t10 = sin(qJ(2));
t1 = (t14 * pkin(3) - mrSges(2,2) + mrSges(4,2) - mrSges(1,3) + mrSges(3,3)) * g(3) + (-t14 * pkin(1) - t12 * t10 + t13 * t11 - mrSges(2,1) - mrSges(1,2)) * g(2) + (-mrSges(1,1) - mrSges(2,3) + t12 * t11 + t13 * t10 + (-m(2) - t14) * qJ(1)) * g(1);
U  = t1;
