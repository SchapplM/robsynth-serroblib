% Calculate joint inertia matrix for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:21
% EndTime: 2019-12-31 17:53:22
% DurationCPUTime: 0.29s
% Computational Cost: add. (223->86), mult. (442->103), div. (0->0), fcn. (351->4), ass. (0->29)
t19 = sin(pkin(7));
t20 = cos(pkin(7));
t41 = t19 ^ 2 + t20 ^ 2;
t40 = m(6) + m(5);
t39 = -mrSges(5,2) + mrSges(6,3);
t32 = mrSges(6,2) + mrSges(5,3);
t38 = -m(6) * pkin(4) - mrSges(6,1);
t37 = mrSges(5,1) - t38;
t36 = m(6) * qJ(5) + t39;
t26 = t19 * qJ(3) + pkin(1);
t7 = (pkin(2) + pkin(3)) * t20 + t26;
t35 = 0.2e1 * t7;
t34 = 2 * mrSges(6,1);
t11 = -t20 * pkin(2) - t26;
t33 = -0.2e1 * t11;
t31 = -pkin(6) + qJ(2);
t30 = t41 * qJ(2) ^ 2;
t12 = t31 * t20;
t21 = sin(qJ(4));
t22 = cos(qJ(4));
t27 = t31 * t19;
t4 = t21 * t12 - t22 * t27;
t6 = t22 * t12 + t21 * t27;
t29 = t4 ^ 2 + t6 ^ 2;
t16 = t19 * mrSges(3,2);
t10 = t19 * t22 - t20 * t21;
t9 = t19 * t21 + t20 * t22;
t1 = t9 * pkin(4) - t10 * qJ(5) + t7;
t2 = [-0.2e1 * pkin(1) * t16 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t33 + (Ifges(4,3) + Ifges(3,2)) * t20) * t20 + (mrSges(4,3) * t33 + (Ifges(4,1) + Ifges(3,1)) * t19 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t20) * t19 + (mrSges(5,1) * t35 + t1 * t34 + (Ifges(6,3) + Ifges(5,2)) * t9 - 0.2e1 * t32 * t6) * t9 + (mrSges(5,2) * t35 - 0.2e1 * t1 * mrSges(6,3) + (Ifges(6,1) + Ifges(5,1)) * t10 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t9 + 0.2e1 * t32 * t4) * t10 + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * qJ(2) * t41 + m(6) * (t1 ^ 2 + t29) + m(5) * (t7 ^ 2 + t29) + m(3) * (pkin(1) ^ 2 + t30) + m(4) * (t11 ^ 2 + t30); -m(3) * pkin(1) - t19 * mrSges(4,3) + t16 + (-mrSges(5,1) - mrSges(6,1)) * t9 + (-mrSges(4,1) - mrSges(3,1)) * t20 + t39 * t10 + m(4) * t11 - m(5) * t7 - m(6) * t1; m(3) + m(4) + t40; (m(4) * qJ(2) + mrSges(4,2)) * t19 + t40 * (t21 * t6 - t22 * t4) + t32 * (-t22 * t10 - t21 * t9); 0; m(4) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t21 ^ 2 + t22 ^ 2); (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t9 + (-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t10 + t36 * t6 - t37 * t4; 0; t36 * t21 + t37 * t22; Ifges(6,2) + Ifges(5,3) + pkin(4) * t34 + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); m(6) * t4 + t10 * mrSges(6,2); 0; -m(6) * t22; t38; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
