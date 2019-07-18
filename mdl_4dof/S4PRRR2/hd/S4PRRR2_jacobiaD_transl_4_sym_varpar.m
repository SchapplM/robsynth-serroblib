% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4PRRR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PRRR2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_jacobiaD_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_jacobiaD_transl_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRR2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_jacobiaD_transl_4_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:29
% EndTime: 2019-07-18 13:27:29
% DurationCPUTime: 0.02s
% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
t24 = qJD(2) + qJD(3);
t33 = pkin(2) * t24;
t21 = qJD(4) + t24;
t25 = qJ(2) + qJ(3);
t23 = qJ(4) + t25;
t32 = sin(t23) * t21;
t31 = cos(t23) * t21;
t30 = r_i_i_C(1) * t32 + r_i_i_C(2) * t31;
t29 = pkin(1) * qJD(2);
t28 = sin(t25) * t33 + t30;
t27 = -r_i_i_C(1) * t31 + r_i_i_C(2) * t32;
t26 = -cos(t25) * t33 + t27;
t1 = [0, sin(qJ(2)) * t29 + t28, t28, t30; 0, 0, 0, 0; 0, -cos(qJ(2)) * t29 + t26, t26, t27;];
JaD_transl  = t1;
