% Calculate potential energy for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:18
% EndTime: 2019-12-31 17:46:19
% DurationCPUTime: 0.24s
% Computational Cost: add. (106->48), mult. (148->41), div. (0->0), fcn. (141->8), ass. (0->19)
t68 = -mrSges(2,1) - mrSges(3,1);
t67 = mrSges(2,2) - mrSges(3,3);
t66 = -m(4) - m(5) - m(6);
t49 = pkin(8) + qJ(5);
t42 = sin(t49);
t43 = cos(t49);
t50 = sin(pkin(8));
t51 = cos(pkin(8));
t65 = mrSges(4,1) + m(5) * pkin(3) + t51 * mrSges(5,1) - t50 * mrSges(5,2) + m(6) * (t51 * pkin(4) + pkin(3)) + t43 * mrSges(6,1) - t42 * mrSges(6,2);
t64 = m(5) * qJ(4) - m(6) * (-pkin(6) - qJ(4)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t63 = sin(qJ(1));
t54 = cos(qJ(1));
t62 = t54 * pkin(1) + t63 * qJ(2);
t61 = cos(pkin(7));
t60 = sin(pkin(7));
t58 = t63 * pkin(1) - t54 * qJ(2);
t37 = t54 * t60 - t63 * t61;
t36 = -t54 * t61 - t63 * t60;
t1 = (t42 * mrSges(6,1) + t51 * mrSges(5,2) + t43 * mrSges(6,2) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t50 + t66 * (-qJ(3) + pkin(5)) + (-m(2) - m(3)) * pkin(5)) * g(3) + (-m(3) * t58 - mrSges(1,2) + t68 * t63 - t67 * t54 + t66 * (t63 * pkin(2) + t58) + t65 * t37 + t64 * t36) * g(2) + (-m(3) * t62 - mrSges(1,1) + t67 * t63 + t68 * t54 + t66 * (t54 * pkin(2) + t62) - t64 * t37 + t65 * t36) * g(1);
U = t1;
