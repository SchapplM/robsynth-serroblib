% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RPRRRR12
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR12_jacobiaD_transl_2_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_transl_2_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_transl_2_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR12_jacobiaD_transl_2_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiaD_transl_2_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:13
% EndTime: 2019-02-26 21:21:13
% DurationCPUTime: 0.05s
% Computational Cost: add. (14->13), mult. (48->27), div. (0->0), fcn. (38->6), ass. (0->12)
t115 = sin(pkin(6));
t125 = t115 * (r_i_i_C(3) + qJ(2));
t117 = cos(pkin(6));
t118 = sin(qJ(1));
t123 = t117 * t118;
t119 = cos(qJ(1));
t122 = t117 * t119;
t121 = qJD(1) * t115;
t120 = t115 * qJD(2);
t116 = cos(pkin(14));
t114 = sin(pkin(14));
t1 = [t119 * t120 + ((t114 * t123 - t116 * t119) * r_i_i_C(1) + (t114 * t119 + t116 * t123) * r_i_i_C(2) - t119 * pkin(1) - t118 * t125) * qJD(1), t119 * t121, 0, 0, 0, 0; t118 * t120 + ((-t114 * t122 - t116 * t118) * r_i_i_C(1) + (t114 * t118 - t116 * t122) * r_i_i_C(2) - t118 * pkin(1) + t119 * t125) * qJD(1), t118 * t121, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JaD_transl  = t1;
