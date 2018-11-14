% Calculate Gravitation load on the joints for
% S4PRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4PRPP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP5_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:10:12
% EndTime: 2018-11-14 14:10:12
% DurationCPUTime: 0.11s
% Computational Cost: add. (28->17), mult. (45->18), div. (0->0), fcn. (26->2), ass. (0->6)
t10 = m(4) + m(5);
t9 = mrSges(3,1) + mrSges(4,1) + mrSges(5,1);
t8 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t5 = cos(qJ(2));
t4 = sin(qJ(2));
t1 = [(-m(2) - m(3) - t10) * g(1) ((m(4) * pkin(2) - m(5) * (-pkin(2) - pkin(3)) + t9) * t4 + (-t10 * qJ(3) + t8) * t5) * g(2) + (-t10 * (t5 * pkin(2) + t4 * qJ(3)) + (-m(5) * pkin(3) - t9) * t5 + t8 * t4) * g(1), t10 * (g(1) * t5 - g(2) * t4) -g(3) * m(5)];
taug  = t1(:);
