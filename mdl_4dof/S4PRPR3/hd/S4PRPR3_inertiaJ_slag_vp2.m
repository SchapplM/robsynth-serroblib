% Calculate joint inertia matrix for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2018-11-14 14:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR3_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:11:15
% EndTime: 2018-11-14 14:11:15
% DurationCPUTime: 0.12s
% Computational Cost: add. (75->35), mult. (160->51), div. (0->0), fcn. (146->6), ass. (0->19)
t20 = m(4) * pkin(2);
t9 = sin(pkin(6));
t19 = pkin(2) * t9;
t11 = sin(qJ(4));
t13 = cos(qJ(4));
t10 = cos(pkin(6));
t8 = pkin(2) * t10 + pkin(3);
t4 = -t11 * t19 + t13 * t8;
t18 = t4 * mrSges(5,1);
t5 = t11 * t8 + t13 * t19;
t17 = t5 * mrSges(5,2);
t12 = sin(qJ(2));
t14 = cos(qJ(2));
t6 = t10 * t14 - t12 * t9;
t7 = t10 * t12 + t14 * t9;
t2 = -t11 * t7 + t13 * t6;
t3 = t11 * t6 + t13 * t7;
t16 = t2 * mrSges(5,1) - t3 * mrSges(5,2);
t1 = [m(2) + m(3) * (t12 ^ 2 + t14 ^ 2) + m(4) * (t6 ^ 2 + t7 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2); t14 * mrSges(3,1) + t6 * mrSges(4,1) - t12 * mrSges(3,2) - t7 * mrSges(4,2) + m(5) * (t2 * t4 + t3 * t5) + (t10 * t6 + t7 * t9) * t20 + t16; 0.2e1 * t18 - 0.2e1 * t17 + Ifges(3,3) + Ifges(4,3) + Ifges(5,3) + m(5) * (t4 ^ 2 + t5 ^ 2) + (0.2e1 * t10 * mrSges(4,1) - 0.2e1 * t9 * mrSges(4,2) + (t10 ^ 2 + t9 ^ 2) * t20) * pkin(2); 0; 0; m(4) + m(5); t16; Ifges(5,3) - t17 + t18; 0; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
