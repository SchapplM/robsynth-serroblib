% Calculate joint inertia matrix for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:21
% EndTime: 2019-12-31 17:46:22
% DurationCPUTime: 0.22s
% Computational Cost: add. (267->75), mult. (442->110), div. (0->0), fcn. (368->6), ass. (0->32)
t28 = sin(pkin(8));
t30 = cos(pkin(8));
t18 = t30 * mrSges(5,1) - t28 * mrSges(5,2);
t32 = sin(qJ(5));
t33 = cos(qJ(5));
t11 = t32 * t28 - t33 * t30;
t13 = t33 * t28 + t32 * t30;
t3 = -t11 * mrSges(6,1) - t13 * mrSges(6,2);
t41 = -t18 - t3;
t26 = t30 ^ 2;
t40 = 2 * mrSges(6,3);
t39 = m(5) + m(6);
t29 = sin(pkin(7));
t31 = cos(pkin(7));
t34 = -pkin(1) - pkin(2);
t17 = t31 * qJ(2) + t29 * t34;
t14 = -qJ(4) + t17;
t38 = pkin(6) - t14;
t36 = t28 ^ 2 + t26;
t35 = t36 * mrSges(5,3);
t16 = -t29 * qJ(2) + t31 * t34;
t15 = pkin(3) - t16;
t27 = t31 ^ 2;
t25 = t29 ^ 2;
t8 = t30 * pkin(4) + t15;
t7 = t11 * t29;
t6 = t13 * t29;
t5 = t38 * t30;
t4 = t38 * t28;
t2 = t32 * t4 - t33 * t5;
t1 = t32 * t5 + t33 * t4;
t9 = [t26 * Ifges(5,2) + (2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t16 * mrSges(4,1) + 0.2e1 * t17 * mrSges(4,2) + 0.2e1 * qJ(2) * mrSges(3,3) + 0.2e1 * t15 * t18 + 0.2e1 * t8 * t3 + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) + (Ifges(5,1) * t28 + 0.2e1 * Ifges(5,4) * t30) * t28 + (Ifges(6,1) * t13 + t1 * t40) * t13 - 0.2e1 * t14 * t35 + (-0.2e1 * Ifges(6,4) * t13 + Ifges(6,2) * t11 + t2 * t40) * t11 + m(6) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) + m(5) * (t14 ^ 2 * t36 + t15 ^ 2) + m(4) * (t16 ^ 2 + t17 ^ 2) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2); -m(3) * pkin(1) - mrSges(3,1) + (-t7 * t11 - t6 * t13) * mrSges(6,3) + (-mrSges(4,1) + t41) * t31 + (mrSges(4,2) - t35) * t29 + m(6) * (-t6 * t1 - t7 * t2 - t31 * t8) + m(5) * (t14 * t29 * t36 - t31 * t15) + m(4) * (t31 * t16 + t29 * t17); m(3) + m(4) * (t25 + t27) + m(5) * (t25 * t36 + t27) + m(6) * (t6 ^ 2 + t7 ^ 2 + t27); m(6) * (-t11 * t1 + t13 * t2); m(6) * (t11 * t6 - t13 * t7); m(4) + m(5) * t36 + m(6) * (t11 ^ 2 + t13 ^ 2); m(5) * t15 + m(6) * t8 - t41; -t39 * t31; 0; t39; t1 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,5) * t13 + Ifges(6,6) * t11; -t6 * mrSges(6,1) + t7 * mrSges(6,2); t3; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
