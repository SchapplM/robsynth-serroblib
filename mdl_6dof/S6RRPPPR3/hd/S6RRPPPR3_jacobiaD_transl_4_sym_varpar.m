% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:03
% EndTime: 2019-02-26 21:23:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (59->24), mult. (168->36), div. (0->0), fcn. (115->4), ass. (0->15)
t19 = sin(qJ(2));
t21 = cos(qJ(2));
t28 = pkin(2) + pkin(3) - r_i_i_C(2);
t33 = r_i_i_C(1) + qJ(3);
t34 = t28 * t19 - t33 * t21;
t35 = t34 * qJD(2) - t19 * qJD(3);
t20 = sin(qJ(1));
t32 = qJD(1) * t20;
t22 = cos(qJ(1));
t31 = qJD(1) * t22;
t30 = qJD(2) * t22;
t27 = pkin(7) - r_i_i_C(3) - qJ(4);
t26 = -t33 * t19 - t28 * t21;
t24 = -pkin(1) + t26;
t1 = [-t22 * qJD(4) + t35 * t20 + (-t27 * t20 + t24 * t22) * qJD(1) (t28 * t32 - t33 * t30) * t19 + (-t33 * t32 + (-t28 * qJD(2) + qJD(3)) * t22) * t21, -t19 * t32 + t21 * t30, -t31, 0, 0; -t20 * qJD(4) - t35 * t22 + (t24 * t20 + t27 * t22) * qJD(1), -t34 * t31 + (t26 * qJD(2) + qJD(3) * t21) * t20, t20 * qJD(2) * t21 + t19 * t31, -t32, 0, 0; 0, -t35, qJD(2) * t19, 0, 0, 0;];
JaD_transl  = t1;
