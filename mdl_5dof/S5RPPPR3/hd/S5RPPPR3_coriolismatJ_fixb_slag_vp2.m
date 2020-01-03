% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:53
% EndTime: 2019-12-31 17:43:53
% DurationCPUTime: 0.24s
% Computational Cost: add. (754->50), mult. (1527->74), div. (0->0), fcn. (1415->6), ass. (0->29)
t44 = sin(pkin(8));
t45 = cos(pkin(8));
t47 = cos(qJ(5));
t63 = sin(qJ(5));
t35 = -t44 * t63 - t45 * t47;
t36 = t44 * t47 - t45 * t63;
t21 = t36 * mrSges(6,1) + t35 * mrSges(6,2);
t69 = t21 * qJD(5);
t67 = -t63 * t35 / 0.2e1 + t47 * t36 / 0.2e1;
t66 = t36 ^ 2;
t37 = sin(pkin(7)) * pkin(1) + qJ(3);
t64 = -pkin(6) + t37;
t55 = cos(pkin(7)) * pkin(1) + pkin(2);
t24 = t44 * qJ(4) + (pkin(3) + pkin(4)) * t45 + t55;
t7 = (-t35 * mrSges(6,1) + t36 * mrSges(6,2) + m(6) * t24 - m(5) * (-t45 * pkin(3) - t55) + t45 * mrSges(5,1) + (m(5) * qJ(4) + mrSges(5,3)) * t44) * t44;
t59 = t7 * qJD(1);
t58 = t21 * qJD(1);
t14 = m(5) * t44 + (t44 / 0.2e1 + t67) * m(6);
t57 = t14 * qJD(1);
t1 = -t66 * Ifges(6,4) + t24 * t21 + (Ifges(6,4) * t35 + (Ifges(6,1) - Ifges(6,2)) * t36) * t35;
t50 = t1 * qJD(1);
t32 = t64 * t44;
t33 = t64 * t45;
t18 = t47 * t32 - t33 * t63;
t19 = t32 * t63 + t47 * t33;
t3 = (-t35 ^ 2 - t66) * mrSges(6,3) + m(6) * (t18 * t36 - t19 * t35) + (mrSges(4,3) + mrSges(5,2) + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t37) * (t44 ^ 2 + t45 ^ 2);
t49 = t3 * qJD(1);
t15 = (-t44 / 0.2e1 + t67) * m(6);
t2 = [t3 * qJD(3) + t7 * qJD(4) + t1 * qJD(5), 0, t15 * qJD(4) + t49, t15 * qJD(3) + t59, (-t19 * mrSges(6,1) - t18 * mrSges(6,2) + Ifges(6,5) * t35 - Ifges(6,6) * t36) * qJD(5) + t50; 0, 0, 0, 0, -t69; -t14 * qJD(4) - t49 - t69, 0, 0, -t57, -t58; t14 * qJD(3) - t59, 0, t57, 0, (-mrSges(6,1) * t63 - t47 * mrSges(6,2)) * qJD(5); qJD(3) * t21 - t50, 0, t58, 0, 0;];
Cq = t2;
