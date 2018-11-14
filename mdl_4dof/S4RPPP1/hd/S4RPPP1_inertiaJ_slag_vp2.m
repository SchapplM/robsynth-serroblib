% Calculate joint inertia matrix for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RPPP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPP1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:25
% EndTime: 2018-11-14 13:45:25
% DurationCPUTime: 0.23s
% Computational Cost: add. (122->73), mult. (309->90), div. (0->0), fcn. (236->4), ass. (0->31)
t36 = 2 * mrSges(3,1);
t35 = -2 * mrSges(4,3);
t34 = -2 * mrSges(5,3);
t33 = m(3) * pkin(1);
t32 = m(4) + m(5);
t22 = cos(pkin(4));
t31 = pkin(1) * t22;
t19 = sin(pkin(6));
t20 = sin(pkin(4));
t30 = t19 * t20;
t21 = cos(pkin(6));
t29 = t20 * t21;
t26 = qJ(2) * t20;
t8 = t19 * t31 + t21 * t26;
t28 = mrSges(5,1) * t29 + t22 * mrSges(5,2);
t27 = mrSges(4,1) * t30 + t22 * mrSges(4,2);
t25 = 0.2e1 * t22;
t24 = -pkin(1) * t21 - pkin(2);
t23 = -qJ(3) * t19 - pkin(1);
t4 = -t22 * qJ(3) - t8;
t14 = mrSges(4,2) * t29;
t12 = mrSges(5,1) * t30;
t11 = mrSges(3,2) * t30;
t9 = t19 * t26;
t7 = t21 * t31 - t9;
t6 = (-pkin(2) * t21 + t23) * t20;
t5 = t24 * t22 + t9;
t3 = ((-pkin(2) - qJ(4)) * t21 + t23) * t20;
t2 = pkin(3) * t29 - t4;
t1 = pkin(3) * t30 + t9 + (-qJ(4) + t24) * t22;
t10 = [Ifges(2,3) + 0.2e1 * t2 * t28 + 0.2e1 * t6 * t14 + 0.2e1 * t1 * t12 + 0.2e1 * t5 * t27 + m(4) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(3) * (t7 ^ 2 + t8 ^ 2) + (t7 * t36 - 0.2e1 * t8 * mrSges(3,2) + t4 * t35 + t1 * t34 + (Ifges(4,1) + Ifges(3,3) + Ifges(5,1)) * t22) * t22 + ((t20 * t33 - 0.2e1 * t11) * pkin(1) + (-0.2e1 * t3 * mrSges(5,2) - 0.2e1 * t7 * mrSges(3,3) + t6 * t35 + (Ifges(5,3) + Ifges(4,2) + Ifges(3,1)) * t30 + (-Ifges(4,4) + Ifges(3,5) + Ifges(5,5)) * t25) * t19 + (-0.2e1 * t4 * mrSges(4,1) + 0.2e1 * t8 * mrSges(3,3) + t3 * t34 + (pkin(1) * t36 + (Ifges(5,2) + Ifges(4,3) + Ifges(3,2)) * t21) * t20 + 0.2e1 * (Ifges(3,4) + Ifges(4,6) - Ifges(5,6)) * t30 + (-Ifges(5,4) - Ifges(4,5) + Ifges(3,6)) * t25) * t21) * t20; t14 + t11 + m(4) * t6 + m(5) * t3 + (-t33 + (-mrSges(3,1) - mrSges(5,3)) * t21 + (-mrSges(5,2) - mrSges(4,3)) * t19) * t20; m(3) + t32; m(4) * t5 + m(5) * t1 - t22 * mrSges(5,3) + t12 + t27; 0; t32; m(5) * t2 + t28; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1) t10(2) t10(4) t10(7); t10(2) t10(3) t10(5) t10(8); t10(4) t10(5) t10(6) t10(9); t10(7) t10(8) t10(9) t10(10);];
Mq  = res;
