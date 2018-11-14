% Calculate Gravitation load on the joints for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2018-11-14 14:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4PRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:12:28
% EndTime: 2018-11-14 14:12:28
% DurationCPUTime: 0.08s
% Computational Cost: add. (49->16), mult. (49->16), div. (0->0), fcn. (28->4), ass. (0->12)
t14 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t11 = m(4) + m(5);
t13 = t11 * pkin(2) + mrSges(3,1);
t10 = mrSges(4,2) + mrSges(5,2);
t5 = qJ(2) + qJ(3);
t3 = sin(t5);
t4 = cos(t5);
t9 = t10 * t3 - t14 * t4;
t8 = t10 * t4 + t14 * t3;
t7 = cos(qJ(2));
t6 = sin(qJ(2));
t1 = [(-m(2) - m(3) - t11) * g(1) (t7 * mrSges(3,2) + t13 * t6 + t8) * g(2) + (t6 * mrSges(3,2) - t13 * t7 + t9) * g(1), t9 * g(1) + t8 * g(2), g(3) * m(5)];
taug  = t1(:);
