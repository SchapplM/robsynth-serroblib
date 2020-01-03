% Calculate joint inertia matrix for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:06
% EndTime: 2019-12-31 17:49:06
% DurationCPUTime: 0.21s
% Computational Cost: add. (252->73), mult. (475->97), div. (0->0), fcn. (399->6), ass. (0->31)
t43 = m(5) + m(6);
t42 = -m(6) * pkin(4) - mrSges(6,1);
t25 = cos(pkin(8));
t22 = t25 ^ 2;
t26 = cos(pkin(7));
t19 = -t26 * pkin(1) - pkin(2);
t15 = -t25 * pkin(3) + t19;
t41 = 0.2e1 * t15;
t24 = sin(pkin(7));
t17 = t24 * pkin(1) + qJ(3);
t40 = pkin(6) + t17;
t39 = cos(qJ(4));
t23 = sin(pkin(8));
t27 = sin(qJ(4));
t38 = t27 * t23;
t37 = mrSges(6,2) + mrSges(5,3);
t36 = t23 ^ 2 + t22;
t32 = t39 * t23;
t8 = t40 * t25;
t4 = t27 * t8 + t40 * t32;
t6 = -t40 * t38 + t39 * t8;
t35 = t4 ^ 2 + t6 ^ 2;
t31 = -t25 * mrSges(4,1) + t23 * mrSges(4,2);
t12 = -t39 * t25 + t38;
t14 = t27 * t25 + t32;
t30 = -pkin(4) * t12 + t14 * qJ(5);
t10 = t14 * mrSges(5,2);
t9 = t12 * mrSges(6,1);
t29 = -t12 * mrSges(5,1) + t14 * mrSges(6,3) - t10 - t9;
t2 = t15 - t30;
t1 = [Ifges(2,3) + Ifges(3,3) + t10 * t41 + 0.2e1 * t2 * t9 + 0.2e1 * t19 * t31 + Ifges(4,2) * t22 + (Ifges(4,1) * t23 + 0.2e1 * Ifges(4,4) * t25) * t23 + 0.2e1 * (t26 * mrSges(3,1) - t24 * mrSges(3,2)) * pkin(1) + 0.2e1 * t36 * t17 * mrSges(4,3) + (-0.2e1 * t2 * mrSges(6,3) + (Ifges(5,1) + Ifges(6,1)) * t14 + 0.2e1 * t37 * t4) * t14 + (mrSges(5,1) * t41 + (Ifges(5,2) + Ifges(6,3)) * t12 - 0.2e1 * t37 * t6 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t14) * t12 + m(5) * (t15 ^ 2 + t35) + m(6) * (t2 ^ 2 + t35) + m(4) * (t36 * t17 ^ 2 + t19 ^ 2) + m(3) * (t24 ^ 2 + t26 ^ 2) * pkin(1) ^ 2; t43 * (t4 * t12 + t6 * t14); m(4) * t36 + m(3) + t43 * (t12 ^ 2 + t14 ^ 2); m(4) * t19 + m(5) * t15 + m(6) * t2 - t29 + t31; 0; m(4) + t43; (-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t14 + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t12 + (m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t6 + (-mrSges(5,1) + t42) * t4; m(6) * t30 + t29; 0; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); m(6) * t4 + t14 * mrSges(6,2); m(6) * t12; 0; t42; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
