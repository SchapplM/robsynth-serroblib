% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR1_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiaD_transl_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_jacobiaD_transl_3_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR1_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiaD_transl_3_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:14
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.06s
% Computational Cost: add. (15->12), mult. (56->29), div. (0->0), fcn. (36->4), ass. (0->12)
t66 = sin(qJ(2));
t75 = qJD(2) * t66;
t68 = cos(qJ(2));
t74 = qJD(2) * t68;
t73 = qJD(3) * t66;
t72 = qJD(3) * t68;
t65 = sin(qJ(3));
t67 = cos(qJ(3));
t71 = -r_i_i_C(1) * t67 + r_i_i_C(2) * t65;
t70 = r_i_i_C(1) * t65 + r_i_i_C(2) * t67;
t69 = t70 * qJD(3);
t1 = [0, t66 * t69 + (-r_i_i_C(3) * t66 + t71 * t68) * qJD(2), (t65 * t72 + t67 * t75) * r_i_i_C(2) + (t65 * t75 - t67 * t72) * r_i_i_C(1), 0, 0; 0, 0, t69, 0, 0; 0, -t70 * t72 + (r_i_i_C(3) * t68 + t71 * t66) * qJD(2), (t65 * t73 - t67 * t74) * r_i_i_C(2) + (-t65 * t74 - t67 * t73) * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
