% Calculate joint inertia matrix for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3PRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_inertiaJ_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_inertiaJ_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_inertiaJ_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRR1_inertiaJ_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRR1_inertiaJ_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:50
% EndTime: 2019-03-08 18:03:50
% DurationCPUTime: 0.08s
% Computational Cost: add. (25->17), mult. (62->27), div. (0->0), fcn. (44->4), ass. (0->9)
t4 = sin(qJ(3));
t5 = sin(qJ(2));
t6 = cos(qJ(3));
t7 = cos(qJ(2));
t2 = t4 * t7 + t6 * t5;
t3 = -t4 * t5 + t6 * t7;
t10 = t3 * mrSges(4,1) - t2 * mrSges(4,2);
t9 = (t6 * mrSges(4,1) - t4 * mrSges(4,2)) * pkin(2);
t1 = [m(2) + m(3) * (t5 ^ 2 + t7 ^ 2) + m(4) * (t2 ^ 2 + t3 ^ 2); t7 * mrSges(3,1) - t5 * mrSges(3,2) + m(4) * (t2 * t4 + t3 * t6) * pkin(2) + t10; Ifges(3,3) + Ifges(4,3) + m(4) * (t4 ^ 2 + t6 ^ 2) * pkin(2) ^ 2 + 0.2e1 * t9; t10; Ifges(4,3) + t9; Ifges(4,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1) t1(2) t1(4); t1(2) t1(3) t1(5); t1(4) t1(5) t1(6);];
Mq  = res;
