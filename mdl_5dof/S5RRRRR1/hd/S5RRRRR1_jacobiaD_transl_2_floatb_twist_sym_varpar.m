% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S5RRRRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR1_jacobiaD_transl_2_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_jacobiaD_transl_2_floatb_twist_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_jacobiaD_transl_2_floatb_twist_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR1_jacobiaD_transl_2_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_jacobiaD_transl_2_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:52:29
% EndTime: 2018-11-16 14:52:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (17->13), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
t17 = sin(qJ(1));
t26 = qJD(1) * t17;
t19 = cos(qJ(1));
t25 = qJD(1) * t19;
t24 = qJD(2) * t17;
t23 = qJD(2) * t19;
t16 = sin(qJ(2));
t18 = cos(qJ(2));
t22 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
t21 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16 - pkin(1);
t20 = t22 * qJD(2);
t1 = [t17 * t20 + (r_i_i_C(3) * t17 + t19 * t21) * qJD(1) (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0, 0; -t22 * t23 + (-r_i_i_C(3) * t19 + t17 * t21) * qJD(1) (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0, 0; 0, t20, 0, 0, 0;];
JaD_transl  = t1;
