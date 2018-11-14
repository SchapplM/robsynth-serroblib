% Calculate joint inertia matrix for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:03
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP2_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP2_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:03:29
% EndTime: 2018-11-14 14:03:29
% DurationCPUTime: 0.09s
% Computational Cost: add. (61->33), mult. (127->41), div. (0->0), fcn. (89->4), ass. (0->18)
t11 = cos(qJ(3));
t17 = -mrSges(4,2) - mrSges(5,2);
t9 = sin(qJ(3));
t22 = (t11 * mrSges(4,1) + t17 * t9) * pkin(2);
t21 = 2 * mrSges(5,1);
t20 = m(5) * pkin(3);
t19 = pkin(2) * t11;
t16 = Ifges(4,3) + Ifges(5,3);
t10 = sin(qJ(2));
t12 = cos(qJ(2));
t5 = t10 * t11 + t12 * t9;
t6 = -t10 * t9 + t11 * t12;
t14 = t17 * t5 + (mrSges(4,1) + mrSges(5,1)) * t6;
t13 = pkin(2) ^ 2;
t8 = t9 ^ 2 * t13;
t7 = pkin(3) + t19;
t2 = t9 * pkin(2) * t5;
t1 = [m(2) + m(3) * (t10 ^ 2 + t12 ^ 2) + 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t5 ^ 2 + t6 ^ 2); t12 * mrSges(3,1) - t10 * mrSges(3,2) + m(4) * (t6 * t19 + t2) + m(5) * (t6 * t7 + t2) + t14; t7 * t21 + Ifges(3,3) + m(4) * (t11 ^ 2 * t13 + t8) + m(5) * (t7 ^ 2 + t8) + 0.2e1 * t22 + t16; t6 * t20 + t14; t7 * t20 + (pkin(3) + t7) * mrSges(5,1) + t22 + t16; (t21 + t20) * pkin(3) + t16; 0; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
