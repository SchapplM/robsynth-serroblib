% Calculate joint inertia matrix for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:29
% EndTime: 2019-12-31 17:56:30
% DurationCPUTime: 0.20s
% Computational Cost: add. (217->69), mult. (347->92), div. (0->0), fcn. (208->6), ass. (0->28)
t23 = cos(qJ(5));
t37 = t23 ^ 2;
t21 = sin(qJ(5));
t8 = -t23 * mrSges(6,1) + t21 * mrSges(6,2);
t36 = -mrSges(5,1) + t8;
t32 = t21 ^ 2 + t37;
t30 = t32 * mrSges(6,3);
t35 = -0.2e1 * t8;
t20 = cos(pkin(8));
t13 = -t20 * pkin(1) - pkin(2);
t11 = -pkin(3) + t13;
t19 = sin(pkin(8));
t12 = t19 * pkin(1) + qJ(3);
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t5 = t22 * t11 + t24 * t12;
t4 = t24 * t11 - t22 * t12;
t34 = t4 * mrSges(5,1);
t33 = t5 * mrSges(5,2);
t3 = -pkin(7) + t5;
t31 = t32 * t3;
t29 = t32 * t22;
t28 = -mrSges(6,1) * t21 - mrSges(6,2) * t23;
t27 = Ifges(6,2) * t37 + Ifges(5,3) + (Ifges(6,1) * t21 + 0.2e1 * Ifges(6,4) * t23) * t21;
t18 = t24 ^ 2;
t16 = t22 ^ 2;
t2 = pkin(4) - t4;
t1 = [-0.2e1 * t13 * mrSges(4,1) - 0.2e1 * t34 + 0.2e1 * t33 + 0.2e1 * t12 * mrSges(4,3) + t2 * t35 + Ifges(4,2) + Ifges(2,3) + Ifges(3,3) + m(6) * (t32 * t3 ^ 2 + t2 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(3) * (t19 ^ 2 + t20 ^ 2) * pkin(1) ^ 2 + t27 + 0.2e1 * (t20 * mrSges(3,1) - t19 * mrSges(3,2)) * pkin(1) - 0.2e1 * t3 * t30; 0; m(6) * t32 + m(3) + m(4) + m(5); -mrSges(4,1) + t36 * t24 + (mrSges(5,2) - t30) * t22 + m(6) * (-t24 * t2 + t3 * t29) + m(5) * (t22 * t5 + t24 * t4) + m(4) * t13; 0; m(4) + m(5) * (t16 + t18) + m(6) * (t32 * t16 + t18); m(6) * (-pkin(4) * t2 + pkin(7) * t31) - t33 + t34 + (t2 + pkin(4)) * t8 + (-t32 * pkin(7) + t31) * mrSges(6,3) - t27; 0; -t22 * mrSges(5,2) + (m(6) * pkin(7) + mrSges(6,3)) * t29 + (m(6) * pkin(4) - t36) * t24; pkin(4) * t35 + m(6) * (t32 * pkin(7) ^ 2 + pkin(4) ^ 2) + 0.2e1 * pkin(7) * t30 + t27; (-mrSges(6,2) * t3 - Ifges(6,6)) * t23 + (-mrSges(6,1) * t3 - Ifges(6,5)) * t21; t8; t28 * t22; Ifges(6,5) * t21 + Ifges(6,6) * t23 + t28 * pkin(7); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
