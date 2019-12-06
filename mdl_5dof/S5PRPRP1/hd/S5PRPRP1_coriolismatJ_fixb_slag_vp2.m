% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:21
% EndTime: 2019-12-05 15:28:22
% DurationCPUTime: 0.33s
% Computational Cost: add. (1012->61), mult. (2194->82), div. (0->0), fcn. (2136->4), ass. (0->38)
t77 = -Ifges(6,5) + Ifges(5,4);
t45 = sin(pkin(8));
t46 = cos(pkin(8));
t47 = sin(qJ(4));
t71 = cos(qJ(4));
t37 = t47 * t45 - t71 * t46;
t39 = t71 * t45 + t47 * t46;
t60 = -t46 * pkin(3) - pkin(2);
t17 = t37 * pkin(4) - t39 * qJ(5) + t60;
t76 = m(6) * t17 + t37 * mrSges(6,1) - t39 * mrSges(6,3);
t42 = m(6) * qJ(5) + mrSges(6,3);
t74 = 0.2e1 * m(6);
t73 = t37 ^ 2;
t72 = pkin(4) * t39;
t69 = pkin(6) + qJ(3);
t67 = t37 * qJ(5);
t23 = t67 + t72;
t32 = t37 * mrSges(6,3);
t58 = t39 * mrSges(5,1) - t37 * mrSges(5,2);
t48 = t39 * mrSges(6,1) + t32 + t58;
t7 = (t23 / 0.4e1 + t72 / 0.4e1 + t67 / 0.4e1) * t74 + t48;
t66 = t7 * qJD(2);
t64 = qJD(4) * t39;
t10 = t76 * t39;
t63 = t10 * qJD(2);
t18 = m(6) * t39;
t62 = t18 * qJD(2);
t61 = t42 * qJD(4);
t59 = t69 * t45;
t1 = t60 * t58 + t17 * t32 + (t17 * mrSges(6,1) + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,3)) * t37 - t77 * t39) * t39 + t77 * t73 + t76 * t23;
t51 = t1 * qJD(2);
t40 = t69 * t46;
t25 = t47 * t40 + t71 * t59;
t26 = t71 * t40 - t47 * t59;
t4 = (t39 ^ 2 + t73) * (mrSges(5,3) + mrSges(6,2)) + (m(4) * qJ(3) + mrSges(4,3)) * (t45 ^ 2 + t46 ^ 2) + (m(6) + m(5)) * (t25 * t39 - t26 * t37);
t49 = t4 * qJD(2);
t12 = m(6) * t26 - t37 * mrSges(6,2);
t2 = [0, 0, 0, -t48 * qJD(4) + (-t23 * qJD(4) / 0.2e1 + t39 * qJD(5) / 0.2e1) * t74, m(6) * t64; 0, t4 * qJD(3) + t1 * qJD(4) - t10 * qJD(5), t49, t12 * qJD(5) + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t64 + t51 + ((-m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1)) * t26 + (mrSges(5,2) - t42) * t25 + (pkin(4) * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t37) * qJD(4), t12 * qJD(4) - t63; 0, t7 * qJD(4) - t18 * qJD(5) - t49, 0, t66, -t62; 0, -t7 * qJD(3) - t51, -t66, t42 * qJD(5), t61; 0, t18 * qJD(3) + t63, t62, -t61, 0;];
Cq = t2;
