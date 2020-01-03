% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:07
% EndTime: 2019-12-31 17:49:08
% DurationCPUTime: 0.34s
% Computational Cost: add. (1238->63), mult. (2420->84), div. (0->0), fcn. (2362->6), ass. (0->39)
t79 = Ifges(5,4) - Ifges(6,5);
t45 = sin(pkin(8));
t46 = cos(pkin(8));
t48 = sin(qJ(4));
t72 = cos(qJ(4));
t38 = t48 * t45 - t46 * t72;
t40 = t45 * t72 + t48 * t46;
t52 = -cos(pkin(7)) * pkin(1) - t46 * pkin(3) - pkin(2);
t17 = t38 * pkin(4) - t40 * qJ(5) + t52;
t78 = m(6) * t17 + t38 * mrSges(6,1) - t40 * mrSges(6,3);
t42 = m(6) * qJ(5) + mrSges(6,3);
t75 = 0.2e1 * m(6);
t74 = t38 ^ 2;
t73 = pkin(4) * t40;
t70 = t38 * qJ(5);
t25 = t70 + t73;
t34 = t38 * mrSges(6,3);
t63 = t40 * mrSges(5,1) - t38 * mrSges(5,2);
t50 = t40 * mrSges(6,1) + t34 + t63;
t7 = (t25 / 0.4e1 + t73 / 0.4e1 + t70 / 0.4e1) * t75 + t50;
t69 = t7 * qJD(1);
t67 = qJD(4) * t40;
t10 = t78 * t40;
t66 = t10 * qJD(1);
t20 = m(6) * t40;
t65 = t20 * qJD(1);
t64 = t42 * qJD(4);
t56 = sin(pkin(7)) * pkin(1) + qJ(3);
t1 = t52 * t63 + t17 * t34 + (t17 * mrSges(6,1) + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,3)) * t38 - t79 * t40) * t40 + t79 * t74 + t78 * t25;
t55 = t1 * qJD(1);
t51 = pkin(6) + t56;
t33 = t51 * t46;
t49 = t51 * t45;
t18 = t48 * t33 + t49 * t72;
t19 = t33 * t72 - t48 * t49;
t4 = (t40 ^ 2 + t74) * (mrSges(6,2) + mrSges(5,3)) + (m(5) + m(6)) * (t18 * t40 - t19 * t38) + (m(4) * t56 + mrSges(4,3)) * (t45 ^ 2 + t46 ^ 2);
t53 = t4 * qJD(1);
t12 = m(6) * t19 - t38 * mrSges(6,2);
t2 = [t4 * qJD(3) + t1 * qJD(4) - t10 * qJD(5), 0, t53, t12 * qJD(5) + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t67 + t55 + ((-m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1)) * t19 + (mrSges(5,2) - t42) * t18 + (pkin(4) * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t38) * qJD(4), t12 * qJD(4) - t66; 0, 0, 0, -t50 * qJD(4) + (-t25 * qJD(4) / 0.2e1 + t40 * qJD(5) / 0.2e1) * t75, m(6) * t67; t7 * qJD(4) - t20 * qJD(5) - t53, 0, 0, t69, -t65; -t7 * qJD(3) - t55, 0, -t69, t42 * qJD(5), t64; t20 * qJD(3) + t66, 0, t65, -t64, 0;];
Cq = t2;
