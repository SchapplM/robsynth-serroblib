% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR1_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_3_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR1_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiaD_transl_3_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:25
% EndTime: 2019-07-18 13:26:25
% DurationCPUTime: 0.07s
% Computational Cost: add. (19->16), mult. (64->31), div. (0->0), fcn. (42->4), ass. (0->13)
t27 = r_i_i_C(3) + qJ(2);
t17 = sin(qJ(1));
t26 = qJD(1) * t17;
t19 = cos(qJ(1));
t25 = qJD(1) * t19;
t24 = qJD(3) * t17;
t23 = qJD(3) * t19;
t16 = sin(qJ(3));
t18 = cos(qJ(3));
t22 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16;
t21 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
t20 = t21 * qJD(3);
t1 = [t19 * qJD(2) + t21 * t24 + (-t17 * t27 + t19 * t22) * qJD(1), t25, (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0; t17 * qJD(2) - t19 * t20 + (t17 * t22 + t19 * t27) * qJD(1), t26, (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0; 0, 0, -t20, 0, 0;];
JaD_transl  = t1;
