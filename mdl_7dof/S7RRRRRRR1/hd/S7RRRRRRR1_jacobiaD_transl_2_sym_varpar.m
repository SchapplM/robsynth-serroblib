% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JaD_transl [3x7]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S7RRRRRRR1_jacobiaD_transl_2_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_2_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_2_sym_varpar: qJD has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobiaD_transl_2_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_2_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:20
% EndTime: 2019-02-26 22:54:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (15->13), mult. (56->29), div. (0->0), fcn. (36->4), ass. (0->12)
t16 = sin(qJ(1));
t25 = qJD(1) * t16;
t18 = cos(qJ(1));
t24 = qJD(1) * t18;
t23 = qJD(2) * t16;
t22 = qJD(2) * t18;
t15 = sin(qJ(2));
t17 = cos(qJ(2));
t21 = -r_i_i_C(1) * t17 + r_i_i_C(2) * t15;
t20 = r_i_i_C(1) * t15 + r_i_i_C(2) * t17;
t19 = t20 * qJD(2);
t1 = [t20 * t23 + (-r_i_i_C(3) * t16 + t21 * t18) * qJD(1) (t15 * t22 + t17 * t25) * r_i_i_C(2) + (t15 * t25 - t17 * t22) * r_i_i_C(1), 0, 0, 0, 0, 0; -t18 * t19 + (r_i_i_C(3) * t18 + t21 * t16) * qJD(1) (t15 * t23 - t17 * t24) * r_i_i_C(2) + (-t15 * t24 - t17 * t23) * r_i_i_C(1), 0, 0, 0, 0, 0; 0, -t19, 0, 0, 0, 0, 0;];
JaD_transl  = t1;
