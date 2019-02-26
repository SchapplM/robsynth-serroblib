% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S3RPR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
%
% Output:
% JaD_transl [3x3]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S3RPR1_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_jacobiaD_transl_3_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_jacobiaD_transl_3_sym_varpar: qJD has to be [3x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S3RPR1_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_jacobiaD_transl_3_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:15:17
% EndTime: 2019-02-26 19:15:17
% DurationCPUTime: 0.03s
% Computational Cost: add. (36->15), mult. (88->23), div. (0->0), fcn. (74->4), ass. (0->14)
t87 = -pkin(1) - pkin(2);
t78 = sin(qJ(1));
t86 = qJD(1) * t78;
t80 = cos(qJ(1));
t85 = qJD(1) * t80;
t84 = qJD(3) * t78;
t83 = qJD(3) * t80;
t77 = sin(qJ(3));
t79 = cos(qJ(3));
t71 = -t77 * t84 - t79 * t83 + (t77 * t78 + t79 * t80) * qJD(1);
t72 = (t84 - t86) * t79 + (-t83 + t85) * t77;
t82 = -t72 * r_i_i_C(1) - t71 * r_i_i_C(2);
t81 = t71 * r_i_i_C(1) - t72 * r_i_i_C(2);
t1 = [t80 * qJD(2) + (-qJ(2) * t78 + t87 * t80) * qJD(1) - t81, t85, t81; t78 * qJD(2) + (qJ(2) * t80 + t87 * t78) * qJD(1) - t82, t86, t82; 0, 0, 0;];
JaD_transl  = t1;
