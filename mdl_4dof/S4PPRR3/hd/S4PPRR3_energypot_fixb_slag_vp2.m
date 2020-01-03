% Calculate potential energy for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:19
% EndTime: 2019-12-31 16:17:20
% DurationCPUTime: 0.15s
% Computational Cost: add. (65->37), mult. (97->33), div. (0->0), fcn. (86->6), ass. (0->16)
t56 = -m(4) - m(5);
t55 = -mrSges(2,1) - mrSges(3,1);
t42 = sin(qJ(4));
t43 = cos(qJ(4));
t54 = m(5) * pkin(3) + t43 * mrSges(5,1) - t42 * mrSges(5,2) + mrSges(4,1);
t53 = mrSges(2,2) - mrSges(3,3);
t52 = m(5) * pkin(5) - mrSges(4,2) + mrSges(5,3);
t51 = cos(qJ(3));
t50 = sin(qJ(3));
t40 = cos(pkin(6));
t48 = sin(pkin(6));
t49 = t40 * pkin(1) + t48 * qJ(2);
t46 = t48 * pkin(1) - qJ(2) * t40;
t31 = t40 * t50 - t48 * t51;
t30 = -t40 * t51 - t48 * t50;
t1 = (t42 * mrSges(5,1) + t43 * mrSges(5,2) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + t56 * (-pkin(4) + qJ(1)) + (-m(2) - m(3)) * qJ(1)) * g(3) + (-m(3) * t46 - mrSges(1,2) + t55 * t48 + t56 * (t48 * pkin(2) + t46) - t53 * t40 + t54 * t31 + t52 * t30) * g(2) + (-m(3) * t49 - mrSges(1,1) + t53 * t48 + t56 * (t40 * pkin(2) + t49) + t55 * t40 - t52 * t31 + t54 * t30) * g(1);
U = t1;
