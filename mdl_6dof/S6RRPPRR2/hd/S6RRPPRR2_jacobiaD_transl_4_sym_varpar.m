% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:45
% EndTime: 2019-02-26 21:28:45
% DurationCPUTime: 0.14s
% Computational Cost: add. (123->22), mult. (204->32), div. (0->0), fcn. (152->8), ass. (0->19)
t176 = sin(pkin(11));
t177 = cos(pkin(11));
t175 = qJ(2) + pkin(10);
t173 = sin(t175);
t174 = cos(t175);
t190 = r_i_i_C(1) * t177 - r_i_i_C(2) * t176 + pkin(3);
t195 = r_i_i_C(3) + qJ(4);
t186 = t190 * t173 - t195 * t174 + sin(qJ(2)) * pkin(2);
t183 = -t186 * qJD(2) + t173 * qJD(4);
t200 = t183 + (r_i_i_C(1) * t176 + r_i_i_C(2) * t177 + pkin(7) + qJ(3)) * qJD(1);
t199 = -t195 * t173 - t190 * t174 - cos(qJ(2)) * pkin(2);
t180 = sin(qJ(1));
t194 = qJD(1) * t180;
t182 = cos(qJ(1));
t193 = qJD(1) * t182;
t192 = qJD(2) * t174;
t185 = qJD(3) + (-pkin(1) + t199) * qJD(1);
t184 = t199 * qJD(2) + qJD(4) * t174;
t1 = [-t200 * t180 + t185 * t182, t184 * t182 + t186 * t194, t193, -t173 * t194 + t182 * t192, 0, 0; t185 * t180 + t200 * t182, t184 * t180 - t186 * t193, t194, t173 * t193 + t180 * t192, 0, 0; 0, t183, 0, qJD(2) * t173, 0, 0;];
JaD_transl  = t1;
