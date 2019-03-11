% Calculate joint inertia matrix for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR2_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR2_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:16:49
% EndTime: 2019-03-08 18:16:50
% DurationCPUTime: 0.07s
% Computational Cost: add. (49->23), mult. (108->34), div. (0->0), fcn. (106->6), ass. (0->13)
t10 = cos(qJ(4));
t11 = cos(qJ(3));
t6 = sin(pkin(6));
t7 = cos(pkin(6));
t9 = sin(qJ(3));
t4 = t11 * t6 + t7 * t9;
t5 = t11 * t7 - t6 * t9;
t8 = sin(qJ(4));
t2 = t10 * t5 - t4 * t8;
t3 = t10 * t4 + t5 * t8;
t14 = t2 * mrSges(5,1) - t3 * mrSges(5,2);
t13 = (mrSges(5,1) * t10 - mrSges(5,2) * t8) * pkin(3);
t1 = [m(2) + m(3) * (t6 ^ 2 + t7 ^ 2) + m(4) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2); 0; m(3) + m(4) + m(5); t5 * mrSges(4,1) - t4 * mrSges(4,2) + m(5) * (t10 * t2 + t3 * t8) * pkin(3) + t14; 0; Ifges(4,3) + Ifges(5,3) + m(5) * (t10 ^ 2 + t8 ^ 2) * pkin(3) ^ 2 + 0.2e1 * t13; t14; 0; Ifges(5,3) + t13; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
