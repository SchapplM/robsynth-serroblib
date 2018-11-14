% Calculate potential energy for
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
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:06
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4PPRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP4_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP4_energypot_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP4_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP4_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:06:02
% EndTime: 2018-11-14 14:06:02
% DurationCPUTime: 0.04s
% Computational Cost: add. (35->28), mult. (35->14), div. (0->0), fcn. (10->2), ass. (0->8)
t18 = -m(4) - m(5);
t17 = -pkin(4) - pkin(1);
t16 = mrSges(4,2) + mrSges(5,2);
t15 = m(3) - t18;
t14 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t13 = cos(qJ(3));
t12 = sin(qJ(3));
t1 = (-mrSges(1,3) + mrSges(2,1) + m(3) * pkin(1) - mrSges(3,2) - m(4) * t17 + mrSges(4,3) - m(5) * (-qJ(4) + t17) + mrSges(5,3)) * g(3) + (t15 * qJ(2) + t14 * t12 + t16 * t13 - mrSges(1,2) - mrSges(2,2) + mrSges(3,3)) * g(2) + (-mrSges(1,1) - mrSges(3,1) - mrSges(2,3) + t18 * pkin(2) - t14 * t13 + t16 * t12 + (-m(2) - t15) * qJ(1)) * g(1);
U  = t1;
