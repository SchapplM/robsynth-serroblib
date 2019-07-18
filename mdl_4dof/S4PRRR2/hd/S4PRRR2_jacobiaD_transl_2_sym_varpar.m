% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
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

function JaD_transl = S4PRRR2_jacobiaD_transl_2_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_jacobiaD_transl_2_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_jacobiaD_transl_2_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRR2_jacobiaD_transl_2_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_jacobiaD_transl_2_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:29
% EndTime: 2019-07-18 13:27:29
% DurationCPUTime: 0.01s
% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
t6 = cos(qJ(2));
t5 = sin(qJ(2));
t1 = [0, (r_i_i_C(1) * t5 + r_i_i_C(2) * t6) * qJD(2), 0, 0; 0, 0, 0, 0; 0, (-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(2), 0, 0;];
JaD_transl  = t1;
