% Calculate joint inertia matrix for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:48
% EndTime: 2019-12-31 16:20:48
% DurationCPUTime: 0.11s
% Computational Cost: add. (85->38), mult. (178->60), div. (0->0), fcn. (147->4), ass. (0->17)
t13 = cos(pkin(7));
t11 = t13 ^ 2;
t19 = pkin(5) + qJ(3);
t12 = sin(pkin(7));
t18 = t12 ^ 2 + t11;
t17 = -t13 * mrSges(4,1) + t12 * mrSges(4,2);
t14 = sin(qJ(4));
t15 = cos(qJ(4));
t4 = -t14 * t12 + t15 * t13;
t5 = t15 * t12 + t14 * t13;
t1 = -t4 * mrSges(5,1) + t5 * mrSges(5,2);
t8 = -t13 * pkin(3) - pkin(2);
t7 = t19 * t13;
t6 = t19 * t12;
t3 = -t14 * t6 + t15 * t7;
t2 = -t14 * t7 - t15 * t6;
t9 = [m(2) + m(3) + m(4) * t18 + m(5) * (t4 ^ 2 + t5 ^ 2); m(5) * (t2 * t4 + t3 * t5); 0.2e1 * t8 * t1 - 0.2e1 * pkin(2) * t17 + Ifges(4,2) * t11 + Ifges(3,3) + (Ifges(4,1) * t12 + 0.2e1 * Ifges(4,4) * t13) * t12 + (-0.2e1 * t2 * mrSges(5,3) + Ifges(5,1) * t5) * t5 + 0.2e1 * t18 * qJ(3) * mrSges(4,3) + m(5) * (t2 ^ 2 + t3 ^ 2 + t8 ^ 2) + m(4) * (t18 * qJ(3) ^ 2 + pkin(2) ^ 2) + (0.2e1 * t3 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t5 + Ifges(5,2) * t4) * t4; 0; -m(4) * pkin(2) + m(5) * t8 + t1 + t17; m(4) + m(5); -t1; t2 * mrSges(5,1) - t3 * mrSges(5,2) + Ifges(5,5) * t5 + Ifges(5,6) * t4; 0; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1), t9(2), t9(4), t9(7); t9(2), t9(3), t9(5), t9(8); t9(4), t9(5), t9(6), t9(9); t9(7), t9(8), t9(9), t9(10);];
Mq = res;
