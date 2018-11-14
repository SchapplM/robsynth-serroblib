% Calculate joint inertia matrix for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RRPP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:29
% EndTime: 2018-11-14 13:51:30
% DurationCPUTime: 0.11s
% Computational Cost: add. (96->44), mult. (176->55), div. (0->0), fcn. (91->4), ass. (0->20)
t23 = -2 * mrSges(5,1);
t22 = 2 * mrSges(5,3);
t10 = sin(pkin(6));
t11 = cos(pkin(6));
t12 = sin(qJ(2));
t21 = pkin(1) * t12;
t13 = cos(qJ(2));
t9 = t13 * pkin(1) + pkin(2);
t4 = t10 * t9 + t11 * t21;
t3 = -t10 * t21 + t11 * t9;
t20 = t3 * mrSges(4,1);
t19 = t4 * mrSges(4,2);
t18 = Ifges(5,2) + Ifges(3,3) + Ifges(4,3);
t17 = t11 * mrSges(4,1) - t10 * mrSges(4,2);
t16 = (t13 * mrSges(3,1) - t12 * mrSges(3,2)) * pkin(1);
t8 = -t11 * pkin(2) - pkin(3);
t7 = t10 * pkin(2) + qJ(4);
t2 = -pkin(3) - t3;
t1 = qJ(4) + t4;
t5 = [0.2e1 * t20 + t2 * t23 - 0.2e1 * t19 + t1 * t22 + Ifges(2,3) + 0.2e1 * t16 + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + m(3) * (t12 ^ 2 + t13 ^ 2) * pkin(1) ^ 2 + t18; m(5) * (t7 * t1 + t8 * t2) - t19 + t20 + t16 + (t7 + t1) * mrSges(5,3) + (-t8 - t2) * mrSges(5,1) + (m(4) * (t10 * t4 + t11 * t3) + t17) * pkin(2) + t18; t8 * t23 + t7 * t22 + m(5) * (t7 ^ 2 + t8 ^ 2) + t18 + (0.2e1 * t17 + m(4) * (t10 ^ 2 + t11 ^ 2) * pkin(2)) * pkin(2); 0; 0; m(4) + m(5); m(5) * t2 - mrSges(5,1); m(5) * t8 - mrSges(5,1); 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1) t5(2) t5(4) t5(7); t5(2) t5(3) t5(5) t5(8); t5(4) t5(5) t5(6) t5(9); t5(7) t5(8) t5(9) t5(10);];
Mq  = res;
