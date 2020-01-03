% Calculate joint inertia matrix for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:18
% EndTime: 2019-12-31 18:14:19
% DurationCPUTime: 0.31s
% Computational Cost: add. (301->97), mult. (523->124), div. (0->0), fcn. (425->4), ass. (0->33)
t43 = m(5) + m(6);
t47 = (-mrSges(6,2) - mrSges(5,3));
t27 = cos(qJ(3));
t46 = t27 ^ 2;
t45 = 2 * mrSges(5,1);
t44 = m(5) * pkin(3);
t24 = sin(pkin(7));
t25 = cos(pkin(7));
t26 = sin(qJ(3));
t12 = t24 * t26 - t25 * t27;
t13 = t24 * t27 + t25 * t26;
t37 = t12 ^ 2 + t13 ^ 2;
t42 = 2 * mrSges(6,3);
t40 = t26 ^ 2 + t46;
t19 = t26 * pkin(3) + qJ(2);
t28 = -pkin(1) - pkin(6);
t39 = -qJ(4) + t28;
t15 = t39 * t26;
t33 = t39 * t27;
t4 = t24 * t15 - t25 * t33;
t6 = t25 * t15 + t24 * t33;
t38 = t4 ^ 2 + t6 ^ 2;
t35 = m(4) * t40;
t34 = t40 * mrSges(4,3);
t32 = 2 * t47;
t31 = t12 * t25 - t13 * t24;
t29 = qJ(2) ^ 2;
t18 = -t25 * pkin(3) - pkin(4);
t16 = t24 * pkin(3) + qJ(5);
t8 = t12 * mrSges(5,2);
t7 = t13 * mrSges(6,1);
t2 = t13 * pkin(4) + t12 * qJ(5) + t19;
t1 = [Ifges(3,1) + Ifges(2,3) + 0.2e1 * t2 * t7 - 0.2e1 * t19 * t8 + Ifges(4,1) * t46 - (2 * pkin(1) * mrSges(3,2)) + (t19 * t45 + (Ifges(6,3) + Ifges(5,2)) * t13 + t6 * t32) * t13 + m(5) * (t19 ^ 2 + t38) + m(6) * (t2 ^ 2 + t38) + m(4) * (t40 * t28 ^ 2 + t29) + m(3) * ((pkin(1) ^ 2) + t29) - 0.2e1 * t28 * t34 + (-0.2e1 * Ifges(4,4) * t27 + Ifges(4,2) * t26) * t26 + 0.2e1 * (t26 * mrSges(4,1) + t27 * mrSges(4,2) + mrSges(3,3)) * qJ(2) + (t2 * t42 + t4 * t32 + 0.2e1 * (Ifges(5,4) - Ifges(6,5)) * t13 + (Ifges(5,1) + Ifges(6,1)) * t12) * t12; -m(3) * pkin(1) + mrSges(3,2) - t34 + t28 * t35 + t43 * (t12 * t4 + t13 * t6) + t47 * t37; t43 * t37 + m(3) + t35; -t4 * mrSges(5,1) - t4 * mrSges(6,1) + t6 * mrSges(6,3) - t6 * mrSges(5,2) + m(6) * (t16 * t6 + t18 * t4) + (t28 * mrSges(4,1) + Ifges(4,5)) * t27 + (-t28 * mrSges(4,2) - Ifges(4,6)) * t26 + (-t16 * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t13 + (-t18 * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t12 + (m(5) * (t24 * t6 - t25 * t4) + t31 * mrSges(5,3)) * pkin(3); t27 * mrSges(4,1) - t26 * mrSges(4,2) + (-mrSges(5,2) + mrSges(6,3)) * t13 + (-mrSges(5,1) - mrSges(6,1)) * t12 + m(6) * (t18 * t12 + t16 * t13) - t31 * t44; -0.2e1 * t18 * mrSges(6,1) + t16 * t42 + Ifges(6,2) + Ifges(4,3) + Ifges(5,3) + m(6) * (t16 ^ 2 + t18 ^ 2) + (t25 * t45 - 0.2e1 * t24 * mrSges(5,2) + (t24 ^ 2 + t25 ^ 2) * t44) * pkin(3); m(5) * t19 + m(6) * t2 + t13 * mrSges(5,1) + t12 * mrSges(6,3) + t7 - t8; 0; 0; t43; m(6) * t4 - t12 * mrSges(6,2); m(6) * t12; m(6) * t18 - mrSges(6,1); 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
