% Calculate joint inertia matrix for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR2_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR2_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:08
% EndTime: 2019-12-31 16:48:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (113->43), mult. (224->60), div. (0->0), fcn. (136->6), ass. (0->23)
t20 = cos(qJ(4));
t33 = t20 ^ 2;
t16 = sin(pkin(7));
t32 = pkin(1) * t16;
t17 = cos(pkin(7));
t11 = t17 * pkin(1) + pkin(2);
t19 = sin(qJ(3));
t21 = cos(qJ(3));
t4 = t21 * t11 - t19 * t32;
t31 = t4 * mrSges(4,1);
t5 = t19 * t11 + t21 * t32;
t30 = t5 * mrSges(4,2);
t18 = sin(qJ(4));
t29 = Ifges(5,5) * t18 + Ifges(5,6) * t20;
t28 = t18 ^ 2 + t33;
t27 = Ifges(5,2) * t33 + Ifges(4,3) + (Ifges(5,1) * t18 + 0.2e1 * Ifges(5,4) * t20) * t18;
t3 = pkin(6) + t5;
t26 = t28 * t3;
t25 = -mrSges(5,1) * t18 - mrSges(5,2) * t20;
t24 = 0.2e1 * t28 * mrSges(5,3);
t8 = -t20 * mrSges(5,1) + t18 * mrSges(5,2);
t2 = -pkin(3) - t4;
t1 = [0.2e1 * t31 - 0.2e1 * t30 + 0.2e1 * t2 * t8 + Ifges(2,3) + Ifges(3,3) + t3 * t24 + m(5) * (t28 * t3 ^ 2 + t2 ^ 2) + m(4) * (t4 ^ 2 + t5 ^ 2) + t27 + (0.2e1 * t17 * mrSges(3,1) - 0.2e1 * t16 * mrSges(3,2) + m(3) * (t16 ^ 2 + t17 ^ 2) * pkin(1)) * pkin(1); 0; m(5) * t28 + m(3) + m(4); m(5) * (-pkin(3) * t2 + pkin(6) * t26) - t30 + t31 + (-pkin(3) + t2) * t8 + (t28 * pkin(6) + t26) * mrSges(5,3) + t27; 0; -0.2e1 * pkin(3) * t8 + m(5) * (t28 * pkin(6) ^ 2 + pkin(3) ^ 2) + pkin(6) * t24 + t27; t25 * t3 + t29; -t8; t25 * pkin(6) + t29; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
