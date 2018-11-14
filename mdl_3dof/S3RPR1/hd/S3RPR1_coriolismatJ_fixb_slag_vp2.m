% Calculate matrix of centrifugal and coriolis load on the joints for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% Cq [3x3]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S3RPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:23
% EndTime: 2018-11-14 10:14:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (123->15), mult. (189->19), div. (0->0), fcn. (110->2), ass. (0->14)
t10 = sin(qJ(3));
t11 = cos(qJ(3));
t13 = t10 * mrSges(4,1) + t11 * mrSges(4,2);
t19 = t13 * qJD(3);
t12 = -pkin(1) - pkin(2);
t8 = -t10 * qJ(2) + t11 * t12;
t9 = t11 * qJ(2) + t10 * t12;
t1 = t9 * mrSges(4,1) + t8 * mrSges(4,2);
t18 = t1 * qJD(1);
t17 = t1 * qJD(3);
t2 = mrSges(3,3) + m(3) * qJ(2) + m(4) * (-t8 * t10 + t9 * t11) + t13;
t16 = t2 * qJD(1);
t15 = t13 * qJD(1);
t3 = [t2 * qJD(2) + t17, t16, -t17 + t18; -t16 + t19, 0, t15 - t19; -qJD(2) * t13 - t18, -t15, 0;];
Cq  = t3;
