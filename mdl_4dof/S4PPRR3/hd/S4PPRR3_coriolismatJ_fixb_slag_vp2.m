% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:08
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:08:16
% EndTime: 2018-11-14 14:08:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (274->14), mult. (622->24), div. (0->0), fcn. (738->6), ass. (0->15)
t25 = sin(pkin(6));
t26 = cos(pkin(6));
t32 = sin(qJ(3));
t33 = cos(qJ(3));
t15 = -t32 * t25 + t33 * t26;
t16 = -t33 * t25 - t32 * t26;
t21 = sin(qJ(4));
t22 = cos(qJ(4));
t13 = t21 * t15 - t22 * t16;
t24 = t22 * t15 + t21 * t16;
t2 = -t13 * mrSges(5,1) - t24 * mrSges(5,2);
t18 = (mrSges(5,1) * t21 + mrSges(5,2) * t22) * pkin(3);
t23 = t18 * qJD(3);
t17 = t18 * qJD(4);
t1 = [0, 0 (t16 * mrSges(4,1) - t15 * mrSges(4,2) + m(5) * (-t13 * t22 + t21 * t24) * pkin(3) + t2) * qJD(3) + t2 * qJD(4) (qJD(3) + qJD(4)) * t2; 0, 0, 0, 0; 0, 0, -t17, -t17 - t23; 0, 0, t23, 0;];
Cq  = t1;
