% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:28
% EndTime: 2019-12-31 17:34:28
% DurationCPUTime: 0.27s
% Computational Cost: add. (303->65), mult. (810->84), div. (0->0), fcn. (509->4), ass. (0->32)
t39 = sin(qJ(4));
t37 = t39 ^ 2;
t41 = cos(qJ(4));
t38 = t41 ^ 2;
t54 = t37 + t38;
t56 = -qJ(5) - pkin(6);
t27 = t56 * t39;
t28 = t56 * t41;
t64 = m(6) * (-t27 * t39 - t28 * t41);
t63 = m(6) * pkin(4);
t49 = mrSges(6,1) + t63;
t40 = sin(qJ(3));
t61 = m(6) * t40;
t60 = t41 * mrSges(6,2);
t57 = Ifges(5,4) + Ifges(6,4);
t55 = t39 * mrSges(5,1) + t41 * mrSges(5,2);
t42 = cos(qJ(3));
t4 = 0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * (-0.1e1 + t54) * t42 * t40;
t53 = t4 * qJD(2);
t15 = t39 * t49 + t60;
t52 = t15 * qJD(3);
t48 = t41 * pkin(4) + pkin(3);
t47 = -t41 * mrSges(6,1) + t39 * mrSges(6,2);
t46 = t39 * mrSges(6,1) + t60;
t1 = pkin(3) * t55 - t57 * t38 + t48 * t46 + (-pkin(4) * t47 + t48 * t63 + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2)) * t41 + t57 * t39) * t39;
t44 = t1 * qJD(3);
t12 = (t37 / 0.2e1 + t38 / 0.2e1 - 0.1e1 / 0.2e1) * t61;
t7 = mrSges(6,3) * t54 + t64;
t43 = t12 * qJD(2) + t7 * qJD(3);
t11 = (t54 + 0.1e1) * t61 / 0.2e1;
t3 = (-t46 / 0.2e1 - t55 / 0.2e1 + (-mrSges(5,2) / 0.2e1 - mrSges(6,2) / 0.2e1) * t41 + (-t63 - mrSges(5,1) / 0.2e1 - mrSges(6,1) / 0.2e1) * t39) * t42;
t2 = [0, 0, 0, (t15 + t55) * qJD(4), 0; 0, t4 * qJD(3), t3 * qJD(4) + t11 * qJD(5) + t53 + ((-m(5) * pkin(3) - m(6) * t48 - t41 * mrSges(5,1) + t39 * mrSges(5,2) - mrSges(4,1) + t47) * t40 + (-mrSges(4,2) + t64 + (m(5) * pkin(6) + mrSges(5,3) + mrSges(6,3)) * t54) * t42) * qJD(3), t3 * qJD(3) + ((mrSges(5,2) + mrSges(6,2)) * t39 + (-mrSges(5,1) - t49) * t41) * qJD(4) * t40, t11 * qJD(3); 0, t12 * qJD(5) - t53, -t1 * qJD(4) + t7 * qJD(5), -t44 + (-t27 * mrSges(6,2) + (mrSges(5,2) * pkin(6) - Ifges(5,6) - Ifges(6,6)) * t39 + (-mrSges(5,1) * pkin(6) - mrSges(6,3) * pkin(4) + Ifges(5,5) + Ifges(6,5)) * t41 + t49 * t28) * qJD(4), t43; 0, 0, -t15 * qJD(5) + t44, 0, -t52; 0, -t12 * qJD(3), t15 * qJD(4) - t43, t52, 0;];
Cq = t2;
