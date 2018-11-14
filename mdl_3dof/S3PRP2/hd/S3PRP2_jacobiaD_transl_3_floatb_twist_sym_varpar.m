% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S3PRP2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
%
% Output:
% JaD_transl [3x3]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:07
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S3PRP2_jacobiaD_transl_3_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP2_jacobiaD_transl_3_floatb_twist_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP2_jacobiaD_transl_3_floatb_twist_sym_varpar: qJD has to be [3x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S3PRP2_jacobiaD_transl_3_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP2_jacobiaD_transl_3_floatb_twist_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:07:28
% EndTime: 2018-11-14 10:07:28
% DurationCPUTime: 0.02s
% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
t10 = -pkin(2) - r_i_i_C(1);
t9 = r_i_i_C(3) + qJ(3);
t8 = cos(qJ(2));
t7 = sin(qJ(2));
t1 = [0, t7 * qJD(3) + (t10 * t7 + t9 * t8) * qJD(2), qJD(2) * t7; 0, t8 * qJD(3) + (t10 * t8 - t9 * t7) * qJD(2), qJD(2) * t8; 0, 0, 0;];
JaD_transl  = t1;
