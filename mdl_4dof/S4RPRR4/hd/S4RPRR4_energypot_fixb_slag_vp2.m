% Calculate potential energy for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:12
% EndTime: 2019-12-31 16:50:12
% DurationCPUTime: 0.18s
% Computational Cost: add. (92->38), mult. (94->33), div. (0->0), fcn. (73->8), ass. (0->17)
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t58 = -m(5) * pkin(3) - t44 * mrSges(5,1) + t41 * mrSges(5,2) - mrSges(4,1);
t57 = -m(5) * pkin(6) + mrSges(4,2) - mrSges(5,3);
t56 = m(4) + m(5);
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t55 = t57 * t42 + t58 * t45 - mrSges(3,1);
t54 = -t41 * mrSges(5,1) - t44 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t43 = sin(qJ(1));
t37 = t43 * pkin(1);
t46 = cos(qJ(1));
t38 = t46 * pkin(1);
t39 = qJ(1) + pkin(7);
t36 = cos(t39);
t35 = sin(t39);
t1 = (-m(2) * pkin(4) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t58 * t42 - t57 * t45 + (-m(3) - t56) * (qJ(2) + pkin(4))) * g(3) + (-m(3) * t37 - t43 * mrSges(2,1) - t46 * mrSges(2,2) - mrSges(1,2) - t56 * (t35 * pkin(2) + t37) + (t56 * pkin(5) - t54) * t36 + t55 * t35) * g(2) + (-m(3) * t38 - t46 * mrSges(2,1) + t43 * mrSges(2,2) - mrSges(1,1) - t56 * (t36 * pkin(2) + t35 * pkin(5) + t38) + t55 * t36 + t54 * t35) * g(1);
U = t1;
