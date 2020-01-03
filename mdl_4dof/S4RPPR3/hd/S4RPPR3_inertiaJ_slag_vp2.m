% Calculate joint inertia matrix for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:50
% EndTime: 2019-12-31 16:37:50
% DurationCPUTime: 0.13s
% Computational Cost: add. (115->44), mult. (224->71), div. (0->0), fcn. (177->6), ass. (0->21)
t17 = cos(pkin(7));
t14 = t17 ^ 2;
t16 = sin(pkin(6));
t10 = t16 * pkin(1) + qJ(3);
t24 = pkin(5) + t10;
t15 = sin(pkin(7));
t23 = t15 ^ 2 + t14;
t18 = cos(pkin(6));
t11 = -t18 * pkin(1) - pkin(2);
t22 = -t17 * mrSges(4,1) + t15 * mrSges(4,2);
t19 = sin(qJ(4));
t20 = cos(qJ(4));
t6 = -t19 * t15 + t20 * t17;
t7 = t20 * t15 + t19 * t17;
t3 = -t6 * mrSges(5,1) + t7 * mrSges(5,2);
t8 = -t17 * pkin(3) + t11;
t5 = t24 * t17;
t4 = t24 * t15;
t2 = -t19 * t4 + t20 * t5;
t1 = -t19 * t5 - t20 * t4;
t9 = [Ifges(2,3) + Ifges(3,3) + 0.2e1 * t8 * t3 + 0.2e1 * t11 * t22 + Ifges(4,2) * t14 + (Ifges(4,1) * t15 + 0.2e1 * Ifges(4,4) * t17) * t15 + (-0.2e1 * t1 * mrSges(5,3) + Ifges(5,1) * t7) * t7 + (0.2e1 * t2 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t7 + Ifges(5,2) * t6) * t6 + m(5) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) + m(4) * (t23 * t10 ^ 2 + t11 ^ 2) + m(3) * (t16 ^ 2 + t18 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t18 * mrSges(3,1) - t16 * mrSges(3,2)) * pkin(1) + 0.2e1 * t23 * t10 * mrSges(4,3); m(5) * (t1 * t6 + t2 * t7); m(3) + m(4) * t23 + m(5) * (t6 ^ 2 + t7 ^ 2); m(4) * t11 + m(5) * t8 + t22 + t3; 0; m(4) + m(5); t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t7 + Ifges(5,6) * t6; -t3; 0; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1), t9(2), t9(4), t9(7); t9(2), t9(3), t9(5), t9(8); t9(4), t9(5), t9(6), t9(9); t9(7), t9(8), t9(9), t9(10);];
Mq = res;
